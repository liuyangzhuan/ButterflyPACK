#pragma once
// Device construction of the compression-only H2 (butterfly_compression.hpp,
// precon=2 with H2_use_gpu).  Per level:
//   ids     the ID of every local box: its two-hop ring (all points of every
//           ring box, as gather_id_workspace forms it when no box is
//           eliminated) and its static training rows (H2_ID_radius > 2,
//           H2_ID_proxy 1), sketched on the device with the host's sparse sign
//           draws (compute_id_sparse_sketch: std::mt19937_64 seeded with
//           morton + 1, 4 per row, sketch size n) and factored by the device
//           pivoted QR.  Skeleton, redundant and T go to the box as
//           h2_skeletonize_box leaves them.  A box the device cannot sketch
//           (sketch too large for the packed lists) is left to the caller.
//   blocks  the coupling block K(skeleton t, skeleton s) of every local
//           target t and source s of its interaction list, and at the leaf
//           level the near block K(points t, points s), evaluated on the
//           device and copied into the boxes' h2 block lists in the host's
//           order (as build_h2_blocks_for_level forms them).
// Kind 3 (EMSURF EFIE) evaluates the ring rows and the blocks by triangle
// pairs (emsurf_plan.hpp).  Empty boxes (unstructured trees) have no ID and
// no blocks.  The host keeps the metadata exchanges and the parent boxes.

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "device_kernels.hpp"
#include "emsurf_plan.hpp"
#include "gpu_runtime.hpp"
#include "h2_matvec_store.hpp"
#include "host_copier.hpp"
#include "kernel_tables.hpp"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fmm {
namespace gpu {

// Totals of the device construction (seconds and bytes, this rank).
struct CompressionStats {
    double ids = 0.0, id_plan = 0.0, id_device = 0.0, id_store = 0.0;
    double blocks = 0.0, blocks_plan = 0.0, blocks_device = 0.0, blocks_wait = 0.0;
    double bytes_up = 0.0, bytes_down = 0.0;
    int64_t id_boxes = 0, host_id_boxes = 0, interaction_blocks = 0, near_blocks = 0;
    size_t heap_peak = 0;
};

inline CompressionStats& compression_stats() {
    static CompressionStats stats;
    return stats;
}

// The device covers the compression of this tree: a registered device form
// of the kernel for the data type, and 3D.  Its IDs also need
// device_sketch_supported (no adaptive training rows).
template<typename CoordType, typename DataType, typename KernelType>
bool compression_supported(const ParallelTree<CoordType, DataType>* tree, const KernelType* kernel,
                           std::string* reason) {
    auto fail = [&](const char* why) {
        if (reason) *reason = why;
        return false;
    };
    if constexpr (!gpu_data_type<DataType>) {
        return fail("single precision is not run on the GPU");
    } else {
        if (tree->dimension != 3) return fail("the device kernels are 3D");
        const char* why = nullptr;
        if (!device_kernel_registered<DataType>(kernel->gpu_spec, &why)) return fail(why);
        return true;
    }
}

// A new compression: device data of earlier factorizations or compressions
// is gone (the caller has released their stores).  The blocks are kept for
// the device matvec unless H2_GPU_MATVEC=0.
inline void begin_device_compression(int num_levels) {
    Context::instance().activate();
    DeviceMatvecStore& store = device_matvec_store();
    store.release();
    DeviceHeap& heap = DeviceHeap::instance();
    heap.ensure_initialized();
    heap.reset();
    heap.reset_peak();
    compression_stats() = CompressionStats{};
    store.building = device_matvec_enabled();
    store.levels.assign(static_cast<size_t>(std::max(num_levels, 0)), H2MatvecLevel{});
}

namespace compression_detail {

// Points of a level on the device: every local box, then every remote box
// with point data (the assisting records the metadata exchanges filled),
// then the static training points of the level's IDs (global indices, with
// the tree's retained global coordinates when it has them: kind 3 reads only
// the ids).  The ids stay on the host too, for the kind-3 blocks.
template<typename CoordType, typename DataType>
class LevelPoints {
public:
    struct Entry {
        int slot = 0;
        int n = 0;
        const std::vector<int64_t>* skeleton = nullptr;  // positions in the box
    };

    LevelPoints() = default;
    LevelPoints(const LevelPoints&) = delete;
    LevelPoints& operator=(const LevelPoints&) = delete;
    ~LevelPoints() { release(); }

    void build(TreeLevel<CoordType, DataType>& level, cudaStream_t stream, const std::vector<int64_t>& training = {},
               const std::vector<CoordType>* global_coords = nullptr) {  // global_coords: 3 per point, or none
        release();
        entries_.clear();
        int64_t slots = 0;
        for (const auto& box : level.local_boxes) {
            entries_[box.morton_index] = Entry{static_cast<int>(slots), static_cast<int>(box.num_points),
                                               &box.skeleton_indices};
            slots += box.num_points;
        }
        std::vector<std::pair<int64_t, int64_t>> remote(level.assisting_box_points_for_kernel_evaluation.begin(),
                                                         level.assisting_box_points_for_kernel_evaluation.end());
        std::sort(remote.begin(), remote.end());
        for (const auto& [morton, index] : remote) {
            if (entries_.count(morton)) continue;
            const auto& assist = level.assisting_boxes[static_cast<size_t>(index)];
            if (assist.indices.empty()) continue;
            if (assist.coords.size() != 3 * assist.indices.size()) {
                throw std::runtime_error("device compression: remote box " + std::to_string(morton) +
                                         " has no coordinates");
            }
            entries_[morton] = Entry{static_cast<int>(slots), static_cast<int>(assist.indices.size()),
                                     &assist.skel_indices};
            slots += static_cast<int64_t>(assist.indices.size());
        }
        training_slot0_ = static_cast<int>(slots);
        slots += static_cast<int64_t>(training.size());
        if (slots >= std::numeric_limits<int>::max()) throw std::runtime_error("device compression: too many points");
        std::vector<double> xyz(static_cast<size_t>(3 * slots));
        std::vector<int64_t>& ids = ids_;
        ids.assign(static_cast<size_t>(slots), 0);
        for (const auto& box : level.local_boxes) {
            const size_t at = static_cast<size_t>(entries_.at(box.morton_index).slot);
            for (int64_t i = 0; i < box.num_points; ++i) {
                for (int d = 0; d < 3; ++d) xyz[3 * (at + i) + d] = static_cast<double>(box.point_coords[3 * i + d]);
                ids[at + i] = box.point_indices[static_cast<size_t>(i)];
            }
        }
        for (const auto& [morton, index] : remote) {
            auto it = entries_.find(morton);
            if (it == entries_.end() || level.find_local_box(morton) != nullptr) continue;
            const auto& assist = level.assisting_boxes[static_cast<size_t>(index)];
            const size_t at = static_cast<size_t>(it->second.slot);
            for (size_t i = 0; i < assist.indices.size(); ++i) {
                for (int d = 0; d < 3; ++d) xyz[3 * (at + i) + d] = static_cast<double>(assist.coords[3 * i + d]);
                ids[at + i] = assist.indices[i];
            }
        }
        for (size_t t = 0; t < training.size(); ++t) {
            const size_t at = static_cast<size_t>(training_slot0_) + t;
            ids[at] = training[t];
            for (int d = 0; d < 3; ++d) {
                xyz[3 * at + d] =
                    global_coords != nullptr ? static_cast<double>((*global_coords)[3 * static_cast<size_t>(training[t]) + d])
                                             : 0.0;
            }
        }
        DeviceHeap& heap = DeviceHeap::instance();
        xyz_ = heap.alloc_resident<double>(std::max<size_t>(xyz.size(), 1) * sizeof(double));
        d_ids_ = heap.alloc_resident<int64_t>(std::max<size_t>(ids.size(), 1) * sizeof(int64_t));
        if (!xyz.empty()) {
            check_cuda(cudaMemcpyAsync(xyz_, xyz.data(), xyz.size() * sizeof(double), cudaMemcpyHostToDevice, stream),
                       "compression points");
            check_cuda(cudaMemcpyAsync(d_ids_, ids.data(), ids.size() * sizeof(int64_t), cudaMemcpyHostToDevice, stream),
                       "compression ids");
            // pageable copies: staged before the calls return, so the host
            // vectors may go
        }
        bytes_ = static_cast<double>(xyz.size() * sizeof(double) + ids.size() * sizeof(int64_t));
    }

    const Entry* find(int64_t morton) const {
        auto it = entries_.find(morton);
        return it == entries_.end() ? nullptr : &it->second;
    }
    PointTable table() const { return PointTable{xyz_, d_ids_}; }
    int training_slot0() const { return training_slot0_; }
    int64_t host_id(int slot) const { return ids_[static_cast<size_t>(slot)]; }
    double bytes() const { return bytes_; }

    void release() {
        DeviceHeap& heap = DeviceHeap::instance();
        heap.free(xyz_);
        heap.free(d_ids_);
        xyz_ = nullptr;
        d_ids_ = nullptr;
    }

private:
    std::unordered_map<int64_t, Entry> entries_;
    double* xyz_ = nullptr;
    int64_t* d_ids_ = nullptr;
    std::vector<int64_t> ids_;
    int training_slot0_ = 0;
    double bytes_ = 0.0;
};

// Device memory a chunk may take: the free heap is shared with the level's
// other allocations and fragmented.
inline size_t chunk_budget() {
    DeviceHeap& heap = DeviceHeap::instance();
    return std::max<size_t>(size_t{64} << 20,
                            std::min({size_t{2} << 30, (heap.capacity() - heap.used()) / 4, heap.largest_free() / 2}));
}

}  // namespace compression_detail

// ---------------------------------------------------------------------------
// IDs of the local boxes of a level.  Returns the boxes left to the host
// (local indices).
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
std::vector<int64_t> compress_level_ids(ParallelTree<CoordType, DataType>* tree, int level_number,
                                        KernelType* kernel, double tolerance) {
    using S = typename DeviceScalar<DataType>::type;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    auto& stats = compression_stats();
    auto& level = tree->levels[static_cast<size_t>(level_number)];
    std::vector<int64_t> host_boxes;
    if (!level.is_process_active) return host_boxes;
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    DeviceHeap& heap = DeviceHeap::instance();
    const KernelSpec spec = device_kernel_spec(kernel->gpu_spec);
    const bool trace = h2_id_trace_enabled();

    // ---- static training rows (H2_ID_radius > 2, H2_ID_proxy 1), as the
    // host selects them: their points after the level's, each once, and a
    // list of their slots per box
    const size_t nb = level.local_boxes.size();
    std::vector<std::vector<int64_t>> training(nb);
    std::vector<int64_t> training_points;  // global indices, in slot order
    std::vector<int> training_lists;       // per box, slots relative to the first training slot
    std::vector<size_t> training_at(nb, 0);
    if (tree->id_neighborhood_radius > 2 || tree->id_proxy_mode == 1) {
        std::exception_ptr failure;
        std::mutex failure_mutex;
        #pragma omp parallel for schedule(dynamic)
        for (int64_t b = 0; b < static_cast<int64_t>(nb); ++b) {
            const auto& box = level.local_boxes[static_cast<size_t>(b)];
            if (box.num_points <= 0) continue;
            try {
                training[static_cast<size_t>(b)] = select_static_id_training_indices(tree, &box);
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
        std::unordered_map<int64_t, int> slot_of;
        for (size_t b = 0; b < nb; ++b) {
            training_at[b] = training_lists.size();
            for (int64_t index : training[b]) {
                auto it = slot_of.emplace(index, static_cast<int>(training_points.size())).first;
                if (it->second == static_cast<int>(training_points.size())) training_points.push_back(index);
                training_lists.push_back(it->second);
            }
        }
    }
    const bool have_coords = tree->id_source_point_coords.size() == 3 * static_cast<size_t>(tree->num_points);
    compression_detail::LevelPoints<CoordType, DataType> points;
    points.build(level, stream, training_points, have_coords ? &tree->id_source_point_coords : nullptr);
    stats.bytes_up += points.bytes();
    int* d_training = nullptr;
    if (!training_lists.empty()) {
        d_training = heap.alloc_resident<int>(training_lists.size() * sizeof(int));
        check_cuda(cudaMemcpyAsync(d_training, training_lists.data(), training_lists.size() * sizeof(int),
                                   cudaMemcpyHostToDevice, stream), "compression training rows");
        stats.bytes_up += static_cast<double>(training_lists.size() * sizeof(int));
    }

    // ---- per box: its ring blocks (all points of each two-hop box, in
    // two_hop order), then its training rows, and the sketch shape
    // (compute_id_sparse_sketch's); kind 3: the rows' edges for the evaluator
    const bool efie = spec.kind == 3;
    struct Plan {
        int64_t box = -1;
        int n = 0, d = 0, sk = 0, slot = 0;
        int64_t rows = 0;
        double scale = 0.0;
        std::vector<RowBlockDesc> blocks;
        std::vector<int> test_edges, row_edges;
    };
    std::vector<Plan> plans;
    for (size_t b = 0; b < nb; ++b) {
        auto& box = level.local_boxes[b];
        box.h2_interaction_blocks.clear();
        box.h2_near_blocks.clear();
        if (box.num_points <= 0) {  // no columns: no skeleton (h2_skeletonize_box)
            box.skeleton_indices.clear();
            box.redundant_indices.clear();
            box.interpolation_matrix.allocate(0, 0, MatrixStorage<DataType>::NONE);
            continue;
        }
        Plan p;
        p.box = static_cast<int64_t>(b);
        p.n = static_cast<int>(box.num_points);
        p.slot = points.find(box.morton_index)->slot;
        int row_base = 0;
        for (int64_t rm : box.two_hop) {
            const auto* e = points.find(rm);
            if (e == nullptr || e->n == 0) {
                throw std::runtime_error("device compression: 2-hop neighbor box " + std::to_string(rm) +
                                         " not found or has no points");
            }
            p.blocks.push_back(RowBlockDesc{row_base, e->n, e->slot, nullptr});
            row_base += e->n;
            if (efie) {
                for (int i = 0; i < e->n; ++i) p.row_edges.push_back(static_cast<int>(points.host_id(e->slot + i)));
            }
        }
        // static training rows after the ring (their draws continue the box's stream, as on the host)
        if (!training[b].empty()) {
            const int nt = static_cast<int>(training[b].size());
            p.blocks.push_back(RowBlockDesc{row_base, nt, points.training_slot0(), d_training + training_at[b]});
            row_base += nt;
            if (efie) {
                for (int64_t index : training[b]) p.row_edges.push_back(static_cast<int>(index));
            }
        }
        p.rows = row_base;
        if (p.rows == 0) {  // no ring rows: every point is a skeleton point (h2_skeletonize_box)
            box.skeleton_indices.resize(static_cast<size_t>(box.num_points));
            std::iota(box.skeleton_indices.begin(), box.skeleton_indices.end(), int64_t{0});
            box.redundant_indices.clear();
            box.interpolation_matrix.allocate(box.num_points, 0);
            if (trace) {
                h2_id_trace_write("L" + std::to_string(box.level) + " m=" + std::to_string(box.morton_index) +
                                  " local=1 ob=" + std::to_string(box.on_boundary) + " n=" +
                                  std::to_string(box.num_points) + " k=" + std::to_string(box.num_points) +
                                  " wave=-1 | rows=0 d=0 sketch_norm=0");
            }
            continue;
        }
        p.d = static_cast<int>(std::max<int64_t>(std::min<int64_t>(p.n, p.rows), p.n));
        p.sk = std::min(4, p.d);
        p.scale = 1.0 / std::sqrt(static_cast<double>(p.sk));
        if (p.rows > kEntryRowMask || p.d > kEntryDestMask) {
            host_boxes.push_back(p.box);  // too large for the packed list entries
            continue;
        }
        if (efie) {
            for (int i = 0; i < p.n; ++i) p.test_edges.push_back(static_cast<int>(points.host_id(p.slot + i)));
        }
        plans.push_back(std::move(p));
    }
    stats.id_plan += std::chrono::duration<double>(clock::now() - t0).count();

    // ---- chunks of boxes: sketch lists, sketch, pivoted QR, then the
    // skeletons and T to the host
    const size_t D = sizeof(S), I = sizeof(int);
    const size_t hist_ints = 32 * static_cast<size_t>(kSketchOwners);
    auto lists_bytes = [&](const Plan& p) {
        const size_t draws = static_cast<size_t>(p.rows) * p.sk;
        return align_up(draws * I) + align_up(static_cast<size_t>(p.rows) * I) + align_up((kSketchOwners + 1) * I) +
               align_up(draws * I) + align_up(hist_ints * I);
    };
    auto y_bytes = [&](const Plan& p) { return align_up(static_cast<size_t>(p.d) * p.n * D); };
    MetaBuilder& meta = pinned_pool().meta;
    DeviceBuffer meta_device;
    const size_t budget = compression_detail::chunk_budget();
    for (size_t c0 = 0; c0 < plans.size();) {
        const auto tc = clock::now();
        size_t c1 = c0, bytes = 0;
        while (c1 < plans.size() && (c1 == c0 || bytes + y_bytes(plans[c1]) + lists_bytes(plans[c1]) <= budget)) {
            bytes += y_bytes(plans[c1]) + lists_bytes(plans[c1]);
            ++c1;
        }
        const int count = static_cast<int>(c1 - c0);
        // one block: Y of every box, then its lists; id results in a small one
        std::vector<size_t> off_y(static_cast<size_t>(count)), off_lists(static_cast<size_t>(count));
        size_t total = 0;
        for (int i = 0; i < count; ++i) {
            off_y[static_cast<size_t>(i)] = total;
            total += y_bytes(plans[c0 + static_cast<size_t>(i)]);
        }
        for (int i = 0; i < count; ++i) {
            off_lists[static_cast<size_t>(i)] = total;
            total += lists_bytes(plans[c0 + static_cast<size_t>(i)]);
        }
        char* block = heap.alloc(std::max<size_t>(total, 1));
        size_t id_bytes = 0;
        const size_t id_rank = id_bytes; id_bytes = align_up(id_bytes + static_cast<size_t>(count) * I);
        const size_t id_flag = id_bytes; id_bytes = align_up(id_bytes + static_cast<size_t>(count) * I);
        const size_t id_norm = id_bytes; id_bytes = align_up(id_bytes + static_cast<size_t>(count) * sizeof(double));
        std::vector<size_t> id_jpvt(static_cast<size_t>(count));
        int max_n = 0, max_d = 0, max_blocks = 0;
        for (int i = 0; i < count; ++i) {
            const Plan& p = plans[c0 + static_cast<size_t>(i)];
            id_jpvt[static_cast<size_t>(i)] = id_bytes;
            id_bytes = align_up(id_bytes + static_cast<size_t>(p.n) * I);
            max_n = std::max(max_n, p.n);
            max_d = std::max(max_d, p.d);
            max_blocks = std::max(max_blocks, static_cast<int>(p.blocks.size()));
        }
        char* d_id = heap.alloc(std::max<size_t>(id_bytes, 1));

        meta.clear();
        std::vector<SketchListsItem> list_items(static_cast<size_t>(count));
        std::vector<OrderedSketchItemT<S>> row_items(static_cast<size_t>(count));
        std::vector<QrcpItemT<S>> id_items(static_cast<size_t>(count));
        for (int i = 0; i < count; ++i) {
            const Plan& p = plans[c0 + static_cast<size_t>(i)];
            const auto& box = level.local_boxes[static_cast<size_t>(p.box)];
            char* at = block + off_lists[static_cast<size_t>(i)];
            const size_t draws = static_cast<size_t>(p.rows) * p.sk;
            int* d_draws = reinterpret_cast<int*>(at);    at += align_up(draws * I);
            int* d_rows = reinterpret_cast<int*>(at);     at += align_up(static_cast<size_t>(p.rows) * I);
            int* d_ptr = reinterpret_cast<int*>(at);      at += align_up((kSketchOwners + 1) * I);
            int* d_list = reinterpret_cast<int*>(at);     at += align_up(draws * I);
            int* d_hist = reinterpret_cast<int*>(at);
            SketchListsItem& li = list_items[static_cast<size_t>(i)];
            li.seed = static_cast<uint64_t>(box.morton_index + 1);
            li.d = p.d;
            li.sk = p.sk;
            li.rows = static_cast<int>(p.rows);
            li.num_blocks = static_cast<int>(p.blocks.size());
            li.blocks_offset = static_cast<int64_t>(meta.append(p.blocks));
            li.draws = d_draws;
            li.rows_out = d_rows;
            li.ptr = d_ptr;
            li.entries = d_list;
            li.hist = d_hist;
            S* Y = reinterpret_cast<S*>(block + off_y[static_cast<size_t>(i)]);
            OrderedSketchItemT<S> item{};
            item.out = Y;
            item.ldo = p.d;
            item.d = p.d;
            item.ncols = p.n;
            item.rows = static_cast<int>(p.rows);
            item.ptr = d_ptr;
            item.entries = d_list;
            item.scale = p.scale;
            item.row_slots = d_rows;
            item.col_base = p.slot;
            row_items[static_cast<size_t>(i)] = item;
            id_items[static_cast<size_t>(i)] =
                QrcpItemT<S>{Y, p.d, p.n, p.d, reinterpret_cast<int*>(d_id + id_jpvt[static_cast<size_t>(i)]),
                             reinterpret_cast<int*>(d_id + id_rank) + i, reinterpret_cast<double*>(d_id + id_norm) + i,
                             reinterpret_cast<int*>(d_id + id_flag) + i};
        }
        EfieBlockSet& efie_rows = efie_block_sets().compression;
        if (efie) {  // row c (ring or training edge), column r (box edge) at out[c * n + r]
            efie_rows.clear();
            for (int i = 0; i < count; ++i) {
                Plan& p = plans[c0 + static_cast<size_t>(i)];
                efie_rows.add(std::move(p.test_edges), std::move(p.row_edges), nullptr, 1, p.n);
            }
            efie_rows.plan(kernel->gpu_spec.table_int, compression_detail::chunk_budget());
            for (int i = 0; i < count; ++i) {
                OrderedSketchItemT<S>& item = row_items[static_cast<size_t>(i)];
                item.src = reinterpret_cast<const S*>(efie_rows.out(static_cast<size_t>(i)));
                item.row_stride = item.ncols;
                item.row_index = nullptr;
            }
        }
        const size_t off_list_items = meta.append(list_items);
        const size_t off_rows = meta.append(row_items);
        const size_t off_id = meta.append(id_items);
        stats.id_plan += std::chrono::duration<double>(clock::now() - tc).count();

        const auto td = clock::now();
        char* md = meta.upload(meta_device, stream);
        stats.bytes_up += static_cast<double>(meta.size());
        launch_sketch_lists(reinterpret_cast<const SketchListsItem*>(md + off_list_items), count, max_blocks, nullptr, 0,
                            md, stream);
        if (efie) {
            if constexpr (std::is_same_v<S, dcomplex>) {
                for (size_t g = 0; g < efie_rows.groups(); ++g) {
                    efie_rows.launch(g, spec, stream);
                    launch_ordered_sketch(reinterpret_cast<const OrderedSketchItemT<S>*>(md + off_rows) +
                                              efie_rows.group_begin(g),
                                          static_cast<int>(efie_rows.group_end(g) - efie_rows.group_begin(g)), max_d,
                                          max_n, false, spec, points.table(), stream);
                }
            }
            efie_rows.release();  // later launches are ordered after the sketch
        } else {
            launch_ordered_sketch(reinterpret_cast<const OrderedSketchItemT<S>*>(md + off_rows), count, max_d, max_n,
                                  true, spec, points.table(), stream);
        }
        char* d_qrcp_work = nullptr;
        if (const size_t w = qrcp_work_bytes<S>(count, max_n)) d_qrcp_work = heap.alloc(w);
        launch_qrcp(reinterpret_cast<const QrcpItemT<S>*>(md + off_id), count, max_n, tolerance, d_qrcp_work, stream);
        heap.free(d_qrcp_work);
        std::vector<char> h_id(id_bytes);
        check_cuda(cudaMemcpyAsync(h_id.data(), d_id, id_bytes, cudaMemcpyDeviceToHost, stream), "compression ids");
        check_cuda(cudaStreamSynchronize(stream), "compression ids");
        const int* ranks = reinterpret_cast<const int*>(h_id.data() + id_rank);
        const int* flags = reinterpret_cast<const int*>(h_id.data() + id_flag);
        const double* norms = reinterpret_cast<const double*>(h_id.data() + id_norm);

        // T (K x (n - K), rows 0..K-1 of columns K.. of Y) of every box, packed
        std::vector<size_t> off_t(static_cast<size_t>(count));
        size_t t_bytes = 0;
        std::vector<GatherItemT<S>> gathers;
        int max_tm = 0, max_tn = 0;
        for (int i = 0; i < count; ++i) {
            const Plan& p = plans[c0 + static_cast<size_t>(i)];
            if (flags[i] == 1) {
                throw std::runtime_error("ID input contains NaN or Inf (box " +
                                         std::to_string(level.local_boxes[static_cast<size_t>(p.box)].morton_index) + ")");
            }
            if (flags[i] == 2) {
                throw std::runtime_error("Triangular solve produced a non-finite interpolation matrix (box " +
                                         std::to_string(level.local_boxes[static_cast<size_t>(p.box)].morton_index) + ")");
            }
            const int K = ranks[i] == 0 ? 1 : ranks[i];  // rank 0 keeps column 0 with a zero row of T
            const int R = p.n - K;
            off_t[static_cast<size_t>(i)] = t_bytes;
            t_bytes = align_up(t_bytes + static_cast<size_t>(K) * R * D);
            if (R > 0) {
                S* Y = reinterpret_cast<S*>(block + off_y[static_cast<size_t>(i)]);
                gathers.push_back(GatherItemT<S>{nullptr, K, K, R, Y + static_cast<size_t>(K) * p.d, 1, p.d,
                                                 IndexList{}, IndexList{}});
                max_tm = std::max(max_tm, K);
                max_tn = std::max(max_tn, R);
            }
        }
        char* d_t = heap.alloc(std::max<size_t>(t_bytes, 1));
        {
            size_t g = 0;
            for (int i = 0; i < count; ++i) {
                const Plan& p = plans[c0 + static_cast<size_t>(i)];
                const int K = ranks[i] == 0 ? 1 : ranks[i];
                if (p.n - K > 0) gathers[g++].out = reinterpret_cast<S*>(d_t + off_t[static_cast<size_t>(i)]);
            }
        }
        meta.clear();
        const size_t off_g = meta.append(gathers);
        md = meta.upload(meta_device, stream);
        launch_gather(reinterpret_cast<const GatherItemT<S>*>(md + off_g), static_cast<int>(gathers.size()), max_tm, max_tn,
                      md, stream);
        std::vector<char> h_t(std::max<size_t>(t_bytes, 1));
        check_cuda(cudaMemcpyAsync(h_t.data(), d_t, t_bytes, cudaMemcpyDeviceToHost, stream), "compression T");
        check_cuda(cudaStreamSynchronize(stream), "compression T");
        stats.bytes_down += static_cast<double>(id_bytes + t_bytes);
        heap.free(d_t);
        heap.free(d_id);
        heap.free(block);
        stats.heap_peak = std::max(stats.heap_peak, heap.peak());
        stats.id_device += std::chrono::duration<double>(clock::now() - td).count();

        // ---- the boxes (compute_id_complex's conventions: full rank keeps the
        // natural order; rank 0 keeps column 0 with a zero T)
        const auto ts = clock::now();
        #pragma omp parallel for schedule(dynamic, 8)
        for (int i = 0; i < count; ++i) {
            const Plan& p = plans[c0 + static_cast<size_t>(i)];
            auto& box = level.local_boxes[static_cast<size_t>(p.box)];
            const int n = p.n;
            const int rank = ranks[i];
            const int* jpvt = reinterpret_cast<const int*>(h_id.data() + id_jpvt[static_cast<size_t>(i)]);
            box.skeleton_indices.clear();
            box.redundant_indices.clear();
            const int K = rank == 0 ? 1 : rank;
            if (rank == n) {
                for (int j = 0; j < n; ++j) box.skeleton_indices.push_back(j);
            } else if (rank == 0) {
                box.skeleton_indices.push_back(0);
                for (int j = 1; j < n; ++j) box.redundant_indices.push_back(j);
            } else {
                box.skeleton_indices.assign(jpvt, jpvt + K);
                box.redundant_indices.assign(jpvt + K, jpvt + n);
            }
            const int R = n - K;
            std::vector<DataType> t(static_cast<size_t>(K) * R);
            if (!t.empty()) std::memcpy(t.data(), h_t.data() + off_t[static_cast<size_t>(i)], t.size() * sizeof(DataType));
            box.interpolation_matrix.set_owned(K, R, std::move(t), MatrixStorage<DataType>::FULL);
            if (trace) {
                std::ostringstream tail;
                tail << std::setprecision(10) << " | rows=" << p.rows << " d=" << p.d << " sketch_norm=" << norms[i];
                h2_id_trace_write("L" + std::to_string(box.level) + " m=" + std::to_string(box.morton_index) +
                                  " local=1 ob=" + std::to_string(box.on_boundary) + " n=" + std::to_string(n) +
                                  " k=" + std::to_string(K) + " wave=-1" + tail.str());
            }
        }
        stats.id_store += std::chrono::duration<double>(clock::now() - ts).count();
        stats.id_boxes += count;
        c0 = c1;
    }
    heap.free(d_training);  // later launches are ordered after the sketches
    stats.host_id_boxes += static_cast<int64_t>(host_boxes.size());
    stats.ids += std::chrono::duration<double>(clock::now() - t0).count();
    return host_boxes;
}

// ---------------------------------------------------------------------------
// The tables of a level's device matvec (kept blocks in their host order per
// box): vector and compact skeleton offsets of the local boxes, their
// skeleton and redundant lists and T, and the sources of every block (local
// boxes, or remote ones in the ghost areas).
// ---------------------------------------------------------------------------
template<typename S, typename CoordType, typename DataType, typename BoxOf, typename NearOf, typename SourceOf,
         typename ColsOf>
void register_matvec_level(TreeLevel<CoordType, DataType>& level, H2MatvecLevel& lv, BoxOf box_of, NearOf near_of,
                           SourceOf source_of, ColsOf cols_of, const std::vector<const void*>& block_ptr,
                           bool symmetric, cudaStream_t stream) {
    const size_t nb = level.local_boxes.size();
    lv.active = true;
    lv.boxes.assign(nb, H2Box{});
    // offsets, and the lists and T packed into one block
    size_t bytes = 0;
    std::vector<size_t> off_skel(nb), off_red(nb), off_t(nb);
    for (size_t b = 0; b < nb; ++b) {
        const auto& box = level.local_boxes[b];
        H2Box& hb = lv.boxes[b];
        hb.vec = lv.points;
        hb.q = lv.q_points;
        hb.n = static_cast<int>(box.num_points);
        hb.k = static_cast<int>(box.skeleton_indices.size());
        hb.r = static_cast<int>(box.redundant_indices.size());
        lv.points += hb.n;
        lv.q_points += hb.k;
        lv.max_n = std::max(lv.max_n, hb.n);
        lv.max_k = std::max(lv.max_k, hb.k);
        lv.max_r = std::max(lv.max_r, hb.r);
        off_skel[b] = bytes; bytes = align_up(bytes + static_cast<size_t>(hb.k) * sizeof(int));
        off_red[b] = bytes;  bytes = align_up(bytes + static_cast<size_t>(hb.r) * sizeof(int));
        off_t[b] = bytes;
        if (hb.r > 0) {
            const auto& T = box.interpolation_matrix;
            if (T.rows != hb.k || T.cols != hb.r || T.lda != hb.k) {
                throw std::runtime_error("device matvec: unexpected T shape for box " + std::to_string(box.morton_index));
            }
            bytes = align_up(bytes + static_cast<size_t>(hb.k) * hb.r * sizeof(S));
        }
    }
    // blocks: coupling blocks of every box in order, then near blocks (they
    // were planned box by box, so each box's blocks are contiguous).  With a
    // symmetric kernel, the near block of local boxes t > s is kept once and
    // serves both through the partial vectors of their pair.
    std::unordered_map<int64_t, size_t> local_index;
    for (size_t b = 0; b < nb; ++b) local_index.emplace(level.local_boxes[b].morton_index, b);
    struct PairKey {
        int64_t t, s;
        bool operator==(const PairKey& o) const { return t == o.t && s == o.s; }
    };
    struct PairHash {
        size_t operator()(const PairKey& k) const { return std::hash<int64_t>()(k.t * 1000003 + k.s); }
    };
    std::unordered_map<PairKey, size_t, PairHash> pair_of;
    auto local_pair = [&](size_t i) {
        return symmetric && near_of(i) && local_index.count(source_of(i)) &&
               source_of(i) != level.local_boxes[box_of(i)].morton_index;
    };
    if (symmetric) {  // the pairs, from their kept blocks (t > s)
        for (size_t i = 0; i < block_ptr.size(); ++i) {
            if (!local_pair(i)) continue;
            const auto& t = level.local_boxes[box_of(i)];
            const int64_t sm = source_of(i);
            if (sm > t.morton_index) continue;
            const H2Box& tb = lv.boxes[box_of(i)];
            const H2Box& sb = lv.boxes[local_index.at(sm)];
            H2Pair pr{block_ptr[i], tb.n, sb.n, tb.vec, sb.vec, 0, 0};
            if (pr.K == nullptr) throw std::runtime_error("device matvec: kept near block missing");
            pr.pa = lv.partial_points;
            lv.partial_points += pr.rows;
            pr.pb = lv.partial_points;
            lv.partial_points += pr.cols;
            lv.max_pair_cols = std::max(lv.max_pair_cols, pr.cols);
            pair_of.emplace(PairKey{t.morton_index, sm}, lv.pairs.size());
            lv.pairs.push_back(pr);
        }
    }
    std::vector<int> ncoupling(nb, 0), nnear(nb, 0);
    for (int pass = 0; pass < 2; ++pass) {
        for (size_t i = 0; i < block_ptr.size(); ++i) {
            if (near_of(i) != (pass == 1)) continue;
            int64_t partial = -1;
            const void* K = block_ptr[i];
            if (local_pair(i)) {
                const int64_t tm = level.local_boxes[box_of(i)].morton_index, sm = source_of(i);
                const H2Pair& pr = lv.pairs.at(tm > sm ? pair_of.at(PairKey{tm, sm}) : pair_of.at(PairKey{sm, tm}));
                partial = tm > sm ? pr.pa : pr.pb;
                K = nullptr;
            } else if (K == nullptr) {
                throw std::runtime_error("device matvec: kept block missing");
            }
            lv.block_K.push_back(K);
            lv.block_cols.push_back(cols_of(i));
            lv.block_source.push_back(source_of(i));
            lv.block_partial.push_back(partial);
            ++(pass == 0 ? ncoupling : nnear)[box_of(i)];
            if (pass == 0) {
                lv.max_q_cols = std::max(lv.max_q_cols, cols_of(i));
            } else if (partial < 0) {
                lv.max_near_cols = std::max(lv.max_near_cols, cols_of(i));
            }
        }
        if (pass == 0) lv.ncoupling = static_cast<int>(lv.block_K.size());
    }
    int at_c = 0, at_n = lv.ncoupling;
    for (size_t b = 0; b < nb; ++b) {
        lv.boxes[b].block0 = at_c;
        lv.boxes[b].nblocks = ncoupling[b];
        lv.boxes[b].near0 = at_n;
        lv.boxes[b].nnear = nnear[b];
        at_c += ncoupling[b];
        at_n += nnear[b];
    }
    const size_t off_boxes = bytes;  bytes = align_up(bytes + nb * sizeof(H2Box));
    char* d = DeviceHeap::instance().alloc_resident(std::max<size_t>(bytes, 1));
    lv.allocations.push_back(d);
    std::vector<char> image(bytes, 0);
    for (size_t b = 0; b < nb; ++b) {
        const auto& box = level.local_boxes[b];
        H2Box& hb = lv.boxes[b];
        int* skel = reinterpret_cast<int*>(image.data() + off_skel[b]);
        for (int i = 0; i < hb.k; ++i) skel[i] = static_cast<int>(box.skeleton_indices[static_cast<size_t>(i)]);
        int* red = reinterpret_cast<int*>(image.data() + off_red[b]);
        for (int i = 0; i < hb.r; ++i) red[i] = static_cast<int>(box.redundant_indices[static_cast<size_t>(i)]);
        hb.skel = reinterpret_cast<const int*>(d + off_skel[b]);
        hb.red = reinterpret_cast<const int*>(d + off_red[b]);
        hb.T = nullptr;
        if (hb.r > 0) {
            std::memcpy(image.data() + off_t[b], box.interpolation_matrix.data.data(),
                        static_cast<size_t>(hb.k) * hb.r * sizeof(S));
            hb.T = d + off_t[b];
        }
    }
    std::memcpy(image.data() + off_boxes, lv.boxes.data(), nb * sizeof(H2Box));
    check_cuda(cudaMemcpyAsync(d, image.data(), bytes, cudaMemcpyHostToDevice, stream), "matvec tables");
    check_cuda(cudaStreamSynchronize(stream), "matvec tables");
    lv.d_boxes = reinterpret_cast<const H2Box*>(d + off_boxes);
    lv.bytes += static_cast<double>(bytes);
}

// ---------------------------------------------------------------------------
// Coupling blocks (sources[b]: the interaction list of local box b, as
// h2_interaction_list forms it) and, when `near`, the leaf near blocks of
// the local boxes of a level.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void build_level_blocks(ParallelTree<CoordType, DataType>* tree, int level_number, KernelType* kernel,
                        const std::vector<std::vector<int64_t>>& sources, bool near) {
    using S = typename DeviceScalar<DataType>::type;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    auto& stats = compression_stats();
    auto& level = tree->levels[static_cast<size_t>(level_number)];
    if (!level.is_process_active) return;
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    DeviceHeap& heap = DeviceHeap::instance();
    const KernelSpec spec = device_kernel_spec(kernel->gpu_spec);
    const size_t nb = level.local_boxes.size();
    if (sources.size() != nb) throw std::runtime_error("build_level_blocks: one interaction list per local box expected");

    compression_detail::LevelPoints<CoordType, DataType> points;
    points.build(level, stream);
    stats.bytes_up += points.bytes();

    // Every block: its box, rows (skeleton or all points of the target) and
    // columns (of the source), in the host's order per box.
    MetaBuilder& meta = pinned_pool().meta;
    meta.clear();
    std::unordered_map<int64_t, int64_t> skeleton_offset;  // morton -> int32 list in the image
    auto skeleton_list = [&](int64_t morton, const std::vector<int64_t>& skeleton) {
        auto it = skeleton_offset.find(morton);
        if (it != skeleton_offset.end()) return it->second;
        std::vector<int> list(skeleton.begin(), skeleton.end());
        const int64_t offset = static_cast<int64_t>(meta.append(list));
        skeleton_offset.emplace(morton, offset);
        return offset;
    };
    struct Block {
        size_t box;
        bool near;
        size_t index;  // in the box's block list
        int m, n;
        IndexList rows, cols;
        int64_t source;
        // the same rows and columns for the host (kind 3): slots, and the
        // skeleton positions in the box (none: all its points)
        int row_slot, col_slot;
        const std::vector<int64_t>* row_skeleton;
        const std::vector<int64_t>* col_skeleton;
    };
    std::vector<Block> blocks;
    for (size_t b = 0; b < nb; ++b) {
        auto& target = level.local_boxes[b];
        target.h2_interaction_blocks.clear();
        target.h2_near_blocks.clear();
        const auto* te = points.find(target.morton_index);
        const int kt = static_cast<int>(target.skeleton_indices.size());
        if (kt > 0) {
            std::vector<int64_t> kept;
            for (int64_t sm : sources[b]) {
                const auto* se = points.find(sm);
                if (se == nullptr) {
                    throw std::runtime_error("device compression: missing remote box " + std::to_string(sm));
                }
                if (se->skeleton->empty()) continue;  // (h2_box_global_indices is empty)
                kept.push_back(sm);
                blocks.push_back(Block{b, false, kept.size() - 1, kt, static_cast<int>(se->skeleton->size()),
                                       IndexList{skeleton_list(target.morton_index, target.skeleton_indices), te->slot},
                                       IndexList{skeleton_list(sm, *se->skeleton), se->slot}, sm, te->slot, se->slot,
                                       &target.skeleton_indices, se->skeleton});
            }
            target.h2_interaction_blocks.resize(kept.size());
            for (size_t j = 0; j < kept.size(); ++j) target.h2_interaction_blocks[j].source_morton = kept[j];
        }
        if (near && target.num_points > 0) {
            std::vector<int64_t> list = target.one_hop;
            list.push_back(target.morton_index);
            std::sort(list.begin(), list.end());
            list.erase(std::unique(list.begin(), list.end()), list.end());
            std::vector<int64_t> kept;
            for (int64_t sm : list) {
                const auto* se = points.find(sm);
                if (se == nullptr) throw std::runtime_error("device compression: missing near box " + std::to_string(sm));
                if (se->n == 0) continue;
                kept.push_back(sm);
                blocks.push_back(Block{b, true, kept.size() - 1, static_cast<int>(target.num_points), se->n,
                                       IndexList{-1, te->slot, nullptr}, IndexList{-1, se->slot, nullptr}, sm, te->slot,
                                       se->slot, nullptr, nullptr});
            }
            target.h2_near_blocks.resize(kept.size());
            for (size_t j = 0; j < kept.size(); ++j) target.h2_near_blocks[j].source_morton = kept[j];
        }
    }
    const size_t lists_end = meta.size();  // the skeleton lists stay in every chunk's image
    stats.blocks_plan += std::chrono::duration<double>(clock::now() - t0).count();

    // ---- chunks: evaluate on the device, copy to the host in the background;
    // for the device matvec, the chunks stay while the heap has room
    const auto td = clock::now();
    const size_t D = sizeof(S);
    auto copier = std::make_unique<HostCopier>(ctx.device());
    std::vector<char*> in_flight;
    DeviceMatvecStore& store = device_matvec_store();
    bool keep = store.building && !store.failed && level_number < static_cast<int>(store.levels.size());
    std::unordered_set<char*> keep_set;   // chunks retained for the matvec
    std::vector<char*> kept_done;         // ... whose host copies are complete
    std::vector<const void*> block_ptr(blocks.size(), nullptr);
    auto reclaim = [&](std::vector<char*>&& done) {
        for (char* p : done) {
            in_flight.erase(std::find(in_flight.begin(), in_flight.end(), p));
            if (keep_set.count(p)) {
                kept_done.push_back(p);
            } else {
                heap.free(p);
            }
        }
    };
    // this rank cannot keep everything: nothing is kept (the matvec runs on
    // the host on every rank, decided at the end of the compression)
    auto give_up = [&] {
        keep = false;
        store.failed = true;
        for (char* p : kept_done) heap.free(p);
        kept_done.clear();
        keep_set.clear();  // chunks still being copied are freed when released
        for (auto& lv : store.levels) {
            for (char* p : lv.allocations) heap.free(p);
            lv = H2MatvecLevel{};
        }
    };
    const double keep_limit = device_matvec_keep_fraction() * static_cast<double>(heap.capacity());
    // With a symmetric kernel, the near block of local boxes t < s is not
    // kept (the block of s, t serves both): such blocks go in chunks of
    // their own after the others, and those chunks only go to the host.
    // Kind 3 keeps every block: the EMSURF near field is symmetric only to
    // about 1e-6 (a self-triangle pair integrates the test side by Gauss
    // points and the source side analytically), and the host matvec applies
    // each block as evaluated.
    const bool symmetric = keep && near && device_matvec_symmetric() && spec.kind != 3;
    std::vector<size_t> order(blocks.size());
    std::iota(order.begin(), order.end(), size_t{0});
    size_t split = blocks.size();  // order[0, split): blocks the matvec keeps
    if (symmetric) {
        auto kept = [&](const Block& bl) {
            if (!bl.near) return true;
            const int64_t tm = level.local_boxes[bl.box].morton_index;
            return bl.source == tm || level.find_local_box(bl.source) == nullptr || bl.source < tm;
        };
        auto mid = std::stable_partition(order.begin(), order.end(), [&](size_t i) { return kept(blocks[i]); });
        split = static_cast<size_t>(mid - order.begin());
        for (size_t i = split; i < order.size(); ++i) {
            if (blocks[order[i]].n > kH2PairMaxCols) throw std::runtime_error("device matvec: near block too wide");
        }
    }
    const size_t budget = compression_detail::chunk_budget();
    DeviceBuffer meta_device;
    std::vector<char> lists_image(meta.host(0), meta.host(0) + lists_end);
    for (size_t c0 = 0; c0 < blocks.size();) {
        size_t c1 = c0, bytes = 0;
        const size_t limit = c0 < split ? split : blocks.size();  // a chunk is kept or not as a whole
        while (c1 < limit && (c1 == c0 || bytes + align_up(static_cast<size_t>(blocks[order[c1]].m) *
                                                           blocks[order[c1]].n * D) <= budget)) {
            bytes += align_up(static_cast<size_t>(blocks[order[c1]].m) * blocks[order[c1]].n * D);
            ++c1;
        }
        const bool keep_chunk = keep && c0 < split;
        reclaim(copier->take_finished());
        if (keep_chunk && static_cast<double>(heap.used() + bytes) > keep_limit) give_up();
        char* out = heap.try_alloc(std::max<size_t>(bytes, 1));
        while (out == nullptr && !in_flight.empty()) {  // wait for earlier chunks' host copies
            const auto tw = clock::now();
            reclaim(copier->wait_released());
            stats.blocks_wait += std::chrono::duration<double>(clock::now() - tw).count();
            out = heap.try_alloc(std::max<size_t>(bytes, 1));
        }
        if (out == nullptr && keep) {  // the kept blocks take the room: keep none
            give_up();
            out = heap.try_alloc(std::max<size_t>(bytes, 1));
        }
        if (out == nullptr) out = heap.alloc(std::max<size_t>(bytes, 1));
        in_flight.push_back(out);
        if (keep && keep_chunk) keep_set.insert(out);

        meta.clear();
        (void)meta.append(lists_image.data(), lists_image.size());
        std::vector<EvalItemT<S>> evals;
        HostCopier::Job job;
        job.device = out;
        job.bytes = bytes;
        struct Shape { MatrixStorage<DataType>* m; int rows, cols; };
        auto shapes = std::make_shared<std::vector<Shape>>();
        int max_m = 0, max_n = 0;
        size_t at = 0;
        for (size_t o = c0; o < c1; ++o) {
            const size_t i = order[o];
            const Block& bl = blocks[i];
            auto& target = level.local_boxes[bl.box];
            auto& storage = bl.near ? target.h2_near_blocks[bl.index].matrix : target.h2_interaction_blocks[bl.index].matrix;
            evals.push_back(EvalItemT<S>{reinterpret_cast<S*>(out + at), bl.m, bl.m, bl.n, bl.rows, bl.cols});
            if (keep_chunk) block_ptr[i] = out + at;
            job.segments.push_back(HostCopier::Segment{&storage.data, nullptr, at, static_cast<size_t>(bl.m) * bl.n * D});
            shapes->push_back(Shape{&storage, bl.m, bl.n});
            max_m = std::max(max_m, bl.m);
            max_n = std::max(max_n, bl.n);
            at += align_up(static_cast<size_t>(bl.m) * bl.n * D);
            (bl.near ? stats.near_blocks : stats.interaction_blocks) += 1;
        }
        if (spec.kind == 3) {  // by triangle pairs, into the chunk (out(r, c) = out[r + c * m])
            if constexpr (std::is_same_v<S, dcomplex>) {
                auto edges = [&](int slot, const std::vector<int64_t>* skeleton, int count) {
                    std::vector<int> e(static_cast<size_t>(count));
                    for (int i = 0; i < count; ++i) {
                        const int at = skeleton != nullptr ? static_cast<int>((*skeleton)[static_cast<size_t>(i)]) : i;
                        e[static_cast<size_t>(i)] = static_cast<int>(points.host_id(slot + at));
                    }
                    return e;
                };
                EfieBlockSet& set = efie_block_sets().compression;
                set.clear();
                for (size_t o = c0; o < c1; ++o) {
                    const Block& bl = blocks[order[o]];
                    set.add(edges(bl.row_slot, bl.row_skeleton, bl.m), edges(bl.col_slot, bl.col_skeleton, bl.n),
                            reinterpret_cast<dcomplex*>(evals[o - c0].out), 1, bl.m);
                }
                set.plan(kernel->gpu_spec.table_int, compression_detail::chunk_budget());
                set.launch_all(spec, stream);
                set.release();  // later launches are ordered after these
            }
        } else {
            const size_t off_evals = meta.append(evals);
            char* md = meta.upload(meta_device, stream);
            stats.bytes_up += static_cast<double>(meta.size());
            launch_eval(reinterpret_cast<const EvalItemT<S>*>(md + off_evals), static_cast<int>(evals.size()), max_m,
                        max_n, md, spec, points.table(), stream);
        }
        check_cuda(cudaEventCreateWithFlags(&job.ready, cudaEventDisableTiming), "cudaEventCreate");
        check_cuda(cudaEventRecord(job.ready, stream), "cudaEventRecord");
        job.finalize = [shapes] {
            for (const Shape& s : *shapes) {
                s.m->rows = s.rows;
                s.m->cols = s.cols;
                s.m->lda = s.rows;
                s.m->format = MatrixStorage<DataType>::FULL;
            }
        };
        copier->submit(std::move(job));
        stats.heap_peak = std::max(stats.heap_peak, heap.peak());
        c0 = c1;
    }
    {
        const auto tw = clock::now();
        reclaim(copier->wait_all());
        stats.blocks_wait += std::chrono::duration<double>(clock::now() - tw).count();
    }
    stats.bytes_down += copier->bytes();
    copier.reset();
    for (char* p : in_flight) heap.free(p);
    if (keep) {
        H2MatvecLevel& lv = store.levels[static_cast<size_t>(level_number)];
        lv.allocations = std::move(kept_done);
        for (size_t o = 0; o < split; ++o) lv.bytes += static_cast<double>(blocks[order[o]].m) * blocks[order[o]].n * D;
        register_matvec_level<S>(level, lv, [&](size_t i) { return blocks[i].box; },
                                 [&](size_t i) { return blocks[i].near; }, [&](size_t i) { return blocks[i].source; },
                                 [&](size_t i) { return blocks[i].n; }, block_ptr, symmetric, stream);
        stats.heap_peak = std::max(stats.heap_peak, heap.peak());
    }
    stats.blocks_device += std::chrono::duration<double>(clock::now() - td).count();
    stats.blocks += std::chrono::duration<double>(clock::now() - t0).count();
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
