#pragma once
// Host side of the kind-3 triangle-pair evaluator (emsurf_blocks.hpp): a set
// of blocks of row edges by column edges, each written to its own output or
// to the set's scratch, launched in groups whose scratch (the pair sums, and
// the outputs kept there) stays within a budget.
//
//   set.clear();
//   set.add(rows, cols, out, rs, cs);   // out == nullptr: the output in the scratch
//   set.plan(edge_table, budget);       // triangles, groups, scratch
//   set.out(i);                         // a scratch output lasts until the next group runs
//   set.launch(g, spec, stream);        // the image is uploaded at the first launch
//   set.release();                      // later launches are ordered after these
//
// A set is kept and reused (efie_block_sets()): its pinned metadata image
// only grows.

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "emsurf_blocks.hpp"
#include "gpu_runtime.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace fmm {
namespace gpu {

class EfieBlockSet {
public:
    void clear() {
        blocks_.clear();
        items_.clear();
        group_.assign(1, 0);
        max_pairs_.clear();
        max_entries_.clear();
        md_ = nullptr;
    }
    size_t size() const { return blocks_.size(); }

    // out(r, c) = out[r * rs + c * cs] for the row edges r and column edges c
    // (0-based edge ids of the mesh tables)
    void add(std::vector<int> rows, std::vector<int> cols, dcomplex* out, int64_t rs, int64_t cs) {
        Block b;
        b.rows = std::move(rows);
        b.cols = std::move(cols);
        b.out = out;
        b.rs = rs;
        b.cs = cs;
        blocks_.push_back(std::move(b));
    }

    // Each block's distinct triangles (its edges' two triangles: entries 2
    // and 3 of the edge table's six per edge, -1 for none), then groups within
    // `budget` bytes of scratch and the scratch itself.
    void plan(const std::vector<int>& edge_table, size_t budget) {
        if (scratch_ != nullptr) throw std::runtime_error("EfieBlockSet: plan before the last release");
        #pragma omp parallel for schedule(dynamic)
        for (int64_t i = 0; i < static_cast<int64_t>(blocks_.size()); ++i) {
            Block& b = blocks_[static_cast<size_t>(i)];
            side(edge_table, b.rows, b.row_tri, b.tr);
            side(edge_table, b.cols, b.col_tri, b.tc);
        }
        const size_t D = sizeof(dcomplex);
        std::vector<size_t> sums_at(blocks_.size()), out_at(blocks_.size(), SIZE_MAX);
        size_t group_bytes = 0, scratch_bytes = 0;
        int64_t max_pairs = 0, max_entries = 0;
        for (size_t i = 0; i < blocks_.size(); ++i) {
            const Block& b = blocks_[i];
            const size_t sums = align_up(b.tr.size() * b.tc.size() * kEfieSums * D);
            const size_t outs = b.out == nullptr ? align_up(b.rows.size() * b.cols.size() * D) : 0;
            if (group_bytes > 0 && group_bytes + sums + outs > budget) {
                group_.push_back(i);
                max_pairs_.push_back(max_pairs);
                max_entries_.push_back(max_entries);
                group_bytes = 0;
                max_pairs = max_entries = 0;
            }
            sums_at[i] = group_bytes;
            if (b.out == nullptr) out_at[i] = group_bytes + sums;
            group_bytes += sums + outs;
            scratch_bytes = std::max(scratch_bytes, group_bytes);
            max_pairs = std::max<int64_t>(max_pairs, static_cast<int64_t>(b.tr.size() * b.tc.size()));
            max_entries = std::max<int64_t>(max_entries, static_cast<int64_t>(b.rows.size() * b.cols.size()));
        }
        group_.push_back(blocks_.size());
        max_pairs_.push_back(max_pairs);
        max_entries_.push_back(max_entries);
        if (blocks_.empty()) return;
        scratch_ = DeviceHeap::instance().alloc(std::max<size_t>(scratch_bytes, 1));
        meta_.clear();
        items_.resize(blocks_.size());
        for (size_t i = 0; i < blocks_.size(); ++i) {
            const Block& b = blocks_[i];
            EfieBlockItem& e = items_[i];
            e.nrow = static_cast<int>(b.rows.size());
            e.ncol = static_cast<int>(b.cols.size());
            e.ntr = static_cast<int>(b.tr.size());
            e.ntc = static_cast<int>(b.tc.size());
            e.row_edges = static_cast<int64_t>(meta_.append(b.rows));
            e.col_edges = static_cast<int64_t>(meta_.append(b.cols));
            e.row_tri = static_cast<int64_t>(meta_.append(b.row_tri));
            e.col_tri = static_cast<int64_t>(meta_.append(b.col_tri));
            e.tr = static_cast<int64_t>(meta_.append(b.tr));
            e.tc = static_cast<int64_t>(meta_.append(b.tc));
            e.M = reinterpret_cast<dcomplex*>(scratch_ + sums_at[i]);
            e.out = b.out != nullptr ? b.out : reinterpret_cast<dcomplex*>(scratch_ + out_at[i]);
            e.rs = b.rs;
            e.cs = b.cs;
        }
        off_items_ = meta_.append(items_);
    }

    size_t groups() const { return group_.size() - 1; }
    size_t group_begin(size_t g) const { return group_[g]; }
    size_t group_end(size_t g) const { return group_[g + 1]; }
    dcomplex* out(size_t i) const { return items_[i].out; }

    void launch(size_t g, const KernelSpec& spec, cudaStream_t stream) {
        if (group_end(g) == group_begin(g)) return;
        if (md_ == nullptr) md_ = meta_.upload(meta_device_, stream);
        launch_efie_blocks(reinterpret_cast<const EfieBlockItem*>(md_ + off_items_) + group_begin(g),
                           static_cast<int>(group_end(g) - group_begin(g)), max_pairs_[g], max_entries_[g], md_, spec,
                           stream);
    }
    void launch_all(const KernelSpec& spec, cudaStream_t stream) {
        for (size_t g = 0; g < groups(); ++g) launch(g, spec, stream);
    }
    void release() {
        DeviceHeap::instance().free(scratch_);  // later launches are ordered after these
        scratch_ = nullptr;
    }

private:
    struct Block {
        std::vector<int> rows, cols, row_tri, col_tri, tr, tc;
        dcomplex* out = nullptr;
        int64_t rs = 0, cs = 0;
    };

    // the distinct triangles of `edges` (sorted), and each edge's two as
    // indices into them
    static void side(const std::vector<int>& edge_table, const std::vector<int>& edges, std::vector<int>& tri_of,
                     std::vector<int>& tris) {
        tris.clear();
        for (int e : edges) {
            for (int a = 0; a < 2; ++a) {
                const int t = edge_table[6 * static_cast<size_t>(e) + 2 + static_cast<size_t>(a)];
                if (t >= 0) tris.push_back(t);
            }
        }
        std::sort(tris.begin(), tris.end());
        tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
        tri_of.resize(2 * edges.size());
        for (size_t i = 0; i < edges.size(); ++i) {
            for (int a = 0; a < 2; ++a) {
                const int t = edge_table[6 * static_cast<size_t>(edges[i]) + 2 + static_cast<size_t>(a)];
                tri_of[2 * i + static_cast<size_t>(a)] =
                    t < 0 ? -1 : static_cast<int>(std::lower_bound(tris.begin(), tris.end(), t) - tris.begin());
            }
        }
    }

    std::vector<Block> blocks_;
    std::vector<EfieBlockItem> items_;
    std::vector<size_t> group_{0};
    std::vector<int64_t> max_pairs_, max_entries_;
    char* scratch_ = nullptr;
    MetaBuilder meta_;
    DeviceBuffer meta_device_;
    size_t off_items_ = 0;
    char* md_ = nullptr;
};

// The backend's sets, kept for the run like pinned_pool(): the sketch rows
// and the other blocks of the factorization's levels, and the compression's.
struct EfieBlockSets {
    EfieBlockSet rows, eval, compression;
};
inline EfieBlockSets& efie_block_sets() {
    static EfieBlockSets* sets = new EfieBlockSets;
    return *sets;
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
