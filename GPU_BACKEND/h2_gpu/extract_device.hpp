#pragma once
// =============================================================================
// The entries of a compression-only H2 on the device (H2_extract;
// unstructured/extraction.hpp extract_compressed_blocks).
//
// The host keeps every decision: it resolves the request to its leaf boxes,
// groups the points, builds the task list and carries the interpolation
// vectors up the tree (plan_extraction), exactly as for the host path.  This
// takes that plan and runs the values on the device:
//   - a near task gathers its entries out of the leaf block K(B_t, B_s)
//     straight into their place in the output (one kernel, no intermediate);
//   - a far task gathers u_t and u_s out of the plan's vectors, takes
//     w = K(S_t, S_s) u_s and g = u_t^T w as two batched GEMMs, and scatters
//     g into the output.
// The blocks are the ones the compression already kept on the device for the
// matvec (h2_matvec_store.hpp): the extraction reads that store and uploads
// nothing but the plan's vectors, its index lists and its items.  So the
// device path needs the store (store.usable, every rank kept every level) and
// adds no residency of its own.
//
// One rank of the caller's communicator only (every intersection is then this
// rank's, so there is no routing and the output image is written in place and
// downloaded once).  That is the case the cost sits in: a STRUMPACK H2 front
// builds its children's contribution blocks on self communicators, and those
// calls are 32.1 s of the 65.2 s a rank the F22 compression takes at 200^3.
// A call over several ranks falls back to the host path, which still routes
// the blocks to their owners.
// =============================================================================

#ifdef H2_HAVE_GPU

#include "device_kernels.hpp"
#include "h2_matvec.hpp"
#include "h2_matvec_store.hpp"
// PartialDeviceWork is the backend's transient-work helper (transient device
// allocations, transfers through the two pinned staging buffers, the
// operator's metadata image).  It is shared with H2_partial because an
// operator is either a partial factorization or a compression-only H2, never
// both.
#include "partial_device.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// H2_extract auto (0): the host, for now, whatever the call's size.
//
// MEASURED at 64^3 (262144 points, leaf 64, 698092 entries in one call, one
// A100, 16 threads), with the 8-step profile of core/extraction_profile.hpp:
//   the plan, which BOTH paths run:  resolve 0.069, tasks 0.008,
//                                    vectors up 0.097, column side 0.002
//   what the device replaces:        compute 0.205, pack 0.007, route 0.005,
//                                    unpack 0.002  (0.219 s on the host)
//   the device instead:              0.33 s (the request's index lists built
//                                    per point, the uploads, the kernels,
//                                    the output down)
// So the device loses 1.5x on the 47% of a 0.436 s call it replaces, and the
// whole call would go 0.436 -> 0.55 s.  (`vectors up` was 0.297 s of a 0.62 s
// call before the window-restricted GEMM form of plan_extraction; closing
// that raised the device's share of the call from 35% to 47% without making
// it win.)
//
// The device path therefore stays available (H2_extract device) for the work
// that would make it pay -- the request built inside plan_extraction instead
// of walked again per point -- but auto does not choose it.  Raise this to
// 1 << 16 (the crossing H2_partial measured for its own entries) once the
// device path wins in situ.
inline constexpr int64_t kExtractDeviceEntries = std::numeric_limits<int64_t>::max();

// One task of the plan: an intersection's row group by column group.  The
// index lists are offsets into the request's `index` pool (positions, not
// bytes), each as long as the count it indexes.
struct ExtractDeviceTask {
    int level = 0;
    int near = 0;        // 1: a near block of the leaf level
    int64_t t = -1;      // the row side's box of `level`
    int64_t s = -1;      // the column side's
    int nr = 0, nc = 0;
    int64_t out_at = 0;  // the intersection's first entry of the output
    int64_t out_ld = 0;  // its rows
    int rvec = -1;       // far: the row and column sides' vector blocks
    int cvec = -1;
    int64_t orow = 0, ocol = 0;  // the output's rows and columns
    int64_t rpos = 0, cpos = 0;  // near: the points' positions in the two boxes
                                 // far: the points' columns of rvec and cvec
};

// The interpolation vectors of one box of a level, as plan_extraction built
// them: column j of a k x n column-major block, on the host.
template<typename DataType>
struct ExtractDeviceVectors {
    int k = 0, n = 0;
    const DataType* u = nullptr;
};

template<typename DataType>
struct ExtractDeviceRequest {
    int64_t entries = 0;  // the output this rank owns, in the host's layout
    std::vector<ExtractDeviceVectors<DataType>> vectors;
    std::vector<ExtractDeviceTask> tasks;
    std::vector<int> index;
};

// Whether a compression-only H2's entries can be computed on the device: the
// compression kept its blocks there for the matvec.  (device_matvec_usable
// also refuses while a matvec check holds the device path off, which has
// nothing to do with the entries, so the store is read directly.)
// `tree` must be the active operator: device_matvec_store() serves whichever
// operator is active and inserts an empty store for it, so asking before
// activate_operator answers for the wrong matrix -- measured, as a device
// path that declined every call while its blocks sat on the card.
inline bool compressed_extract_ready(const void* tree) {
    const DeviceMatvecStore& store = device_matvec_store();
    return store.usable && store.tree == tree;
}

namespace extract_device_detail {

// The device block of a task, with the strides that read it: a near block of
// two local boxes of a symmetric kernel is kept once, for the pair's larger
// Morton index, so the other direction reads it transposed (rs and cs
// swapped).  Returns false when the task's block is not in the store.
struct ResolvedBlock {
    const void* K = nullptr;
    int64_t rs = 1, cs = 1;
    int rows = 0;        // its leading dimension (the GEMM reads K this way)
    bool empty = false;  // an empty skeleton or no coupling: the block is zero
};

template<typename CoordType, typename DataType>
class BlockFinder {
public:
    BlockFinder(ParallelTree<CoordType, DataType>* tree, const DeviceMatvecStore& store)
        : tree_(tree), store_(store), index_(static_cast<size_t>(tree->num_levels)),
          pair_of_(static_cast<size_t>(tree->num_levels)), ready_(static_cast<size_t>(tree->num_levels), false) {}

    bool find(const ExtractDeviceTask& task, ResolvedBlock& out) {
        const size_t l = static_cast<size_t>(task.level);
        if (l >= store_.levels.size() || !store_.levels[l].active) return false;
        prepare(l);
        const auto bit = index_[l].find(task.t);
        if (bit == index_[l].end()) return false;
        const H2MatvecLevel& lv = store_.levels[l];
        const H2Box& tb = lv.boxes[bit->second];
        const int i0 = task.near ? tb.near0 : tb.block0, ni = task.near ? tb.nnear : tb.nblocks;
        const int rows = task.near ? tb.n : tb.k;
        if (rows == 0) {  // an empty skeleton: no coupling, zeros
            out = ResolvedBlock{nullptr, 1, 1, 0, true};
            return true;
        }
        for (int i = i0; i < i0 + ni; ++i) {
            if (lv.block_source[static_cast<size_t>(i)] != task.s) continue;
            const void* K = lv.block_K[static_cast<size_t>(i)];
            if (K != nullptr) {
                out = ResolvedBlock{K, 1, rows, rows, false};
                return true;
            }
            // the pair's kept block, found by the partial vector it was given
            const auto pit = pair_of_[l].find(lv.block_partial[static_cast<size_t>(i)]);
            if (pit == pair_of_[l].end()) return false;
            const H2Pair& pr = lv.pairs[static_cast<size_t>(pit->second.first)];
            if (pr.K == nullptr) return false;
            out = pit->second.second ? ResolvedBlock{pr.K, pr.rows, 1, pr.rows, false}   // transposed
                                     : ResolvedBlock{pr.K, 1, pr.rows, pr.rows, false};
            return true;
        }
        // The two paths part on a missing block: the row side of a far task
        // need not couple to the column side (the host leaves those entries
        // zero), while a leaf near block it asks for must be there.
        if (!task.near) {
            out = ResolvedBlock{nullptr, 1, 1, 0, true};
            return true;
        }
        return false;
    }

private:
    void prepare(size_t l) {
        if (ready_[l]) return;
        const auto& boxes = tree_->levels[l].local_boxes;
        for (size_t b = 0; b < boxes.size(); ++b) index_[l].emplace(boxes[b].morton_index, b);
        const H2MatvecLevel& lv = store_.levels[l];
        for (size_t q = 0; q < lv.pairs.size(); ++q) {
            pair_of_[l].emplace(lv.pairs[q].pa, std::make_pair(q, false));
            pair_of_[l].emplace(lv.pairs[q].pb, std::make_pair(q, true));
        }
        ready_[l] = true;
    }

    ParallelTree<CoordType, DataType>* tree_;
    const DeviceMatvecStore& store_;
    std::vector<std::unordered_map<int64_t, size_t>> index_;
    std::vector<std::unordered_map<int64_t, std::pair<size_t, bool>>> pair_of_;
    std::vector<char> ready_;
};

}  // namespace extract_device_detail

// The request's values into `out` (the host's layout, request.entries of
// them).  Returns false, having written nothing, when the store cannot serve
// a task: the caller then runs the host path, and `why` says what stopped it
// (an explicit H2_extract device must not fall back silently).
template<typename CoordType, typename DataType>
bool extract_compressed_blocks_device(ParallelTree<CoordType, DataType>* tree,
                                      const ExtractDeviceRequest<DataType>& request, DataType* out,
                                      std::string* why = nullptr) {
    using S = typename DeviceScalar<DataType>::type;
    using extract_device_detail::BlockFinder;
    using extract_device_detail::ResolvedBlock;
    constexpr size_t D = sizeof(DataType);
    auto no = [&](const std::string& reason) {
        if (why != nullptr) *why = reason;
        return false;
    };
    OperatorCall call(tree);
    if (!compressed_extract_ready(tree)) {
        const DeviceMatvecStore& s = device_matvec_store();
        return no(s.levels.empty() ? "the compression kept no blocks on the device"
                                   : s.tree != tree ? "the device blocks belong to another operator"
                                                    : "a rank could not keep all of its levels");
    }
    const DeviceMatvecStore& store = device_matvec_store();

    // ---- the tasks' blocks, before anything is uploaded
    BlockFinder<CoordType, DataType> finder(tree, store);
    std::vector<ResolvedBlock> blocks(request.tasks.size());
    for (size_t q = 0; q < request.tasks.size(); ++q) {
        if (!finder.find(request.tasks[q], blocks[q])) {
            const ExtractDeviceTask& t = request.tasks[q];
            return no("level " + std::to_string(t.level) + (t.near ? " near" : " coupling") + " block of box " +
                      std::to_string(t.t) + " from " + std::to_string(t.s) + " is not in the device store");
        }
    }

    PartialDeviceWork work;
    const cudaStream_t stream = work.stream();

    // ---- the index lists and the plan's vectors, each in one buffer (the
    // lists are read through device pointers, so they are not in the
    // metadata image and every group of tasks reuses them)
    int* d_index = nullptr;
    if (!request.index.empty()) {
        d_index = reinterpret_cast<int*>(work.alloc(request.index.size() * sizeof(int)));
        work.upload(d_index, request.index.data(), request.index.size() * sizeof(int));
    }
    std::vector<int64_t> vec_at(request.vectors.size() + 1, 0);
    for (size_t v = 0; v < request.vectors.size(); ++v) {
        vec_at[v + 1] = vec_at[v] + static_cast<int64_t>(request.vectors[v].k) * request.vectors[v].n;
    }
    S* d_vec = nullptr;
    if (vec_at.back() > 0) {
        // The plan holds a block per box of a level, so a call at 64^3 has
        // thousands of them; each upload ends in a stream synchronize, so they
        // are packed into one buffer and sent once.  Measured: 4672 blocks a
        // call one at a time made the device path 10% slower than the host's.
        std::vector<DataType> packed(static_cast<size_t>(vec_at.back()));
        const int64_t nv = static_cast<int64_t>(request.vectors.size());
#pragma omp parallel for schedule(dynamic)
        for (int64_t v = 0; v < nv; ++v) {
            const auto& bv = request.vectors[static_cast<size_t>(v)];
            const size_t n = static_cast<size_t>(bv.k) * bv.n;
            if (n > 0) std::copy_n(bv.u, n, packed.begin() + vec_at[static_cast<size_t>(v)]);
        }
        d_vec = reinterpret_cast<S*>(work.alloc(static_cast<size_t>(vec_at.back()) * D));
        work.upload(d_vec, packed.data(), static_cast<size_t>(vec_at.back()) * D);
    }
    // ---- the output: the tasks write their own rectangles, and a task whose
    // block is empty writes nothing, so it starts at zero (as the host's does)
    S* d_out = reinterpret_cast<S*>(work.alloc(static_cast<size_t>(std::max<int64_t>(request.entries, 1)) * D));
    work.zero(d_out, static_cast<size_t>(std::max<int64_t>(request.entries, 1)) * D);

    auto list = [&](int64_t at) {
        IndexList l;
        l.ptr = d_index + at;
        return l;
    };

    // ---- the near tasks: one gather each, into the output
    {
        std::vector<BlockCopyItemT<S>> items;
        int max_m = 0, max_n = 0;
        for (size_t q = 0; q < request.tasks.size(); ++q) {
            const ExtractDeviceTask& task = request.tasks[q];
            if (!task.near || blocks[q].empty || task.nr == 0 || task.nc == 0) continue;
            BlockCopyItemT<S> it{};
            it.out = d_out + task.out_at;
            it.ors = 1;
            it.ocs = task.out_ld;
            it.out_rows = list(task.orow);
            it.out_cols = list(task.ocol);
            it.src = reinterpret_cast<const S*>(blocks[q].K);
            it.rs = blocks[q].rs;
            it.cs = blocks[q].cs;
            it.src_rows = list(task.rpos);
            it.src_cols = list(task.cpos);
            it.m = task.nr;
            it.n = task.nc;
            items.push_back(it);
            max_m = std::max(max_m, task.nr);
            max_n = std::max(max_n, task.nc);
        }
        if (!items.empty()) {
            MetaBuilder& meta = work.meta();
            const size_t off = meta.append(items);
            char* md = work.upload_meta();
            launch_block_copy(reinterpret_cast<const BlockCopyItemT<S>*>(md + off), static_cast<int>(items.size()),
                              max_m, max_n, md, stream);
            work.sync();  // (the far tasks' first image overwrites this one)
        }
    }

    // ---- the far tasks: u_t and u_s gathered, w = K u_s, g = u_t^T w, g
    // scattered into the output.  Their intermediates are transient, so the
    // tasks run in groups that fit the work's budget.
    std::vector<size_t> far;
    for (size_t q = 0; q < request.tasks.size(); ++q) {
        const ExtractDeviceTask& task = request.tasks[q];
        if (task.near || blocks[q].empty || task.nr == 0 || task.nc == 0) continue;
        if (task.rvec < 0 || task.cvec < 0) return no("a far task without interpolation vectors");
        if (blocks[q].rows != request.vectors[static_cast<size_t>(task.rvec)].k) {
            return no("level " + std::to_string(task.level) + " box " + std::to_string(task.t) + ": the device block has " +
                      std::to_string(blocks[q].rows) + " rows, the plan's skeleton " +
                      std::to_string(request.vectors[static_cast<size_t>(task.rvec)].k));
        }
        far.push_back(q);
    }
    const size_t budget = work.budget();
    for (size_t at = 0; at < far.size();) {
        size_t take = 0, bytes = 0;
        while (at + take < far.size()) {
            const ExtractDeviceTask& task = request.tasks[far[at + take]];
            const int64_t kt = request.vectors[static_cast<size_t>(task.rvec)].k;
            const int64_t ks = request.vectors[static_cast<size_t>(task.cvec)].k;
            const size_t need = static_cast<size_t>(kt * task.nr + ks * task.nc + kt * task.nc +
                                                    static_cast<int64_t>(task.nr) * task.nc) * D;
            if (take > 0 && bytes + need > budget) break;
            bytes += need;
            ++take;
        }
        char* d = work.alloc(bytes);
        std::vector<BlockCopyItemT<S>> gathers, scatters;
        VBatch<DataType> couple, apply;  // w = K u_s (N, N); g = u_t^T w (T, N)
        int gm = 0, gn = 0, sm = 0, sn = 0;
        size_t off = 0;
        for (size_t i = 0; i < take; ++i) {
            const ExtractDeviceTask& task = request.tasks[far[at + i]];
            const auto& rv = request.vectors[static_cast<size_t>(task.rvec)];
            const auto& cv = request.vectors[static_cast<size_t>(task.cvec)];
            S* ut = reinterpret_cast<S*>(d + off);
            off += static_cast<size_t>(rv.k) * task.nr * D;
            S* us = reinterpret_cast<S*>(d + off);
            off += static_cast<size_t>(cv.k) * task.nc * D;
            S* w = reinterpret_cast<S*>(d + off);
            off += static_cast<size_t>(rv.k) * task.nc * D;
            S* g = reinterpret_cast<S*>(d + off);
            off += static_cast<size_t>(task.nr) * task.nc * D;
            auto gather = [&](S* dst, int k, int n, int64_t vat, int64_t cols) {
                BlockCopyItemT<S> it{};
                it.out = dst;
                it.ors = 1;
                it.ocs = k;
                it.src = d_vec + vat;
                it.rs = 1;
                it.cs = k;
                it.src_cols = list(cols);
                it.m = k;
                it.n = n;
                gathers.push_back(it);
                gm = std::max(gm, k);
                gn = std::max(gn, n);
            };
            gather(ut, rv.k, task.nr, vec_at[static_cast<size_t>(task.rvec)], task.rpos);
            gather(us, cv.k, task.nc, vec_at[static_cast<size_t>(task.cvec)], task.cpos);
            couple.entries.push_back({reinterpret_cast<const DataType*>(blocks[far[at + i]].K),
                                      reinterpret_cast<const DataType*>(us), reinterpret_cast<DataType*>(w), rv.k,
                                      task.nc, cv.k, rv.k, cv.k, rv.k});
            apply.entries.push_back({reinterpret_cast<const DataType*>(ut), reinterpret_cast<const DataType*>(w),
                                     reinterpret_cast<DataType*>(g), task.nr, task.nc, rv.k, rv.k, rv.k, task.nr});
            BlockCopyItemT<S> it{};
            it.out = d_out + task.out_at;
            it.ors = 1;
            it.ocs = task.out_ld;
            it.out_rows = list(task.orow);
            it.out_cols = list(task.ocol);
            it.src = g;
            it.rs = 1;
            it.cs = task.nr;
            it.m = task.nr;
            it.n = task.nc;
            scatters.push_back(it);
            sm = std::max(sm, task.nr);
            sn = std::max(sn, task.nc);
        }
        MetaBuilder& meta = work.meta();
        const size_t gather_off = meta.append(gathers);
        const size_t scatter_off = meta.append(scatters);
        couple.stage(meta);
        apply.stage(meta);
        char* md = work.upload_meta();
        launch_block_copy(reinterpret_cast<const BlockCopyItemT<S>*>(md + gather_off),
                          static_cast<int>(gathers.size()), gm, gn, md, stream);
        couple.gemm(md, MagmaNoTrans, MagmaNoTrans, DataType(1), DataType(0), work.queue());
        apply.gemm(md, MagmaTrans, MagmaNoTrans, DataType(1), DataType(0), work.queue());
        launch_block_copy(reinterpret_cast<const BlockCopyItemT<S>*>(md + scatter_off),
                          static_cast<int>(scatters.size()), sm, sn, md, stream);
        work.sync();  // (the next group reuses these bytes and the image)
        work.release(d);
        at += take;
    }

    if (request.entries > 0) work.download(out, d_out, static_cast<size_t>(request.entries) * D);
    return true;
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
