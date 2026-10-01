#pragma once
// The application's evaluator of the matrix's entries on the device
// (GPU_INTERFACE/bpack_gpu.h; user guide doc/gpu_kernels.md), one per matrix
// (H2: H2Kernel::gpu_evaluator; HODLR: its GpuState), in one of three forms:
//
//  - entry: the library's entry and sketch kernels for the application's
//    entry evaluator, instantiated in the application (bpack_gpu_entry.cuh)
//    or compiled here with NVRTC from its source text;
//  - block: the application's function filling a batch of blocks;
//  - list: the application's function filling the values of a list of
//    (row, column) pairs.
//
// The library describes the blocks it needs as eval items (index lists into
// the point table of the level, possibly built on the device); for the block
// and list forms the ids (and coordinates) of each block's points are
// gathered on the device from those lists, so the host need not know them.
// The sketch of kernel rows runs fused with an entry evaluator, else its
// rows are evaluated into scratch first and sketched from there.

#ifdef H2_HAVE_GPU

#include "../../GPU_INTERFACE/bpack_gpu.h"
#include "device_heap.hpp"
#include "device_kernels.hpp"
#include "evaluator_device.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fmm {
namespace gpu {

// ---- the evaluator

class Evaluator {
public:
    enum class Form { entry, block, list };

    Evaluator(const Evaluator&) = delete;
    Evaluator& operator=(const Evaluator&) = delete;
    ~Evaluator() {
        if (params_ != nullptr) cudaFree(params_);
    }

    // Entry kernels instantiated in the application (bpack_gpu_entry.cuh).
    static std::shared_ptr<Evaluator> from_launchers(const bpack_gpu_entry_launchers& l, int scalar_bytes) {
        if (l.version != BPACK_GPU_INTERFACE_VERSION) {
            throw std::invalid_argument("GPU entry evaluator compiled against version " + std::to_string(l.version) +
                                        " of the ButterflyPACK GPU headers, the library is version " +
                                        std::to_string(BPACK_GPU_INTERFACE_VERSION) + ": rebuild it");
        }
        if (l.scalar_bytes != scalar_bytes) {
            throw std::invalid_argument("GPU entry evaluator of the wrong value type for this matrix (" +
                                        std::to_string(l.scalar_bytes) + " bytes, the matrix's " +
                                        std::to_string(scalar_bytes) + ")");
        }
        if (l.eval == nullptr || l.sketch == nullptr || l.entry == nullptr) {
            throw std::invalid_argument("GPU entry evaluator without its launches");
        }
        std::shared_ptr<Evaluator> e(new Evaluator(Form::entry, l.flags, scalar_bytes));
        e->launchers_ = l;
        e->entry_.assign(static_cast<const char*>(l.entry), static_cast<const char*>(l.entry) + l.entry_bytes);
        e->launchers_.entry = nullptr;
        e->origin_ = "entry evaluator";
        return e;
    }
    // Entry kernels compiled from source text (NVRTC, at the first use).
    static std::shared_ptr<Evaluator> from_source(const char* source, const double* params, int nparams, int flags,
                                                  int scalar_bytes) {
        if (source == nullptr) throw std::invalid_argument("GPU entry evaluator: no source");
        std::shared_ptr<Evaluator> e(new Evaluator(Form::entry, flags, scalar_bytes));
        e->source_ = source;
        e->host_params_.assign(params, params + std::max(nparams, 0));
        e->origin_ = "entry evaluator (source, NVRTC)";
        return e;
    }
    static std::shared_ptr<Evaluator> from_block(bpack_gpu_block_evaluator f, void* user, int flags,
                                                 int scalar_bytes) {
        if (f == nullptr) throw std::invalid_argument("GPU block evaluator: null function");
        std::shared_ptr<Evaluator> e(new Evaluator(Form::block, flags, scalar_bytes));
        e->block_ = f;
        e->user_ = user;
        e->origin_ = "block evaluator";
        return e;
    }
    static std::shared_ptr<Evaluator> from_list(bpack_gpu_list_evaluator f, void* user, int flags,
                                                int scalar_bytes) {
        if (f == nullptr) throw std::invalid_argument("GPU entry-list evaluator: null function");
        std::shared_ptr<Evaluator> e(new Evaluator(Form::list, flags, scalar_bytes));
        e->list_ = f;
        e->user_ = user;
        e->origin_ = "entry-list evaluator";
        return e;
    }

    Form form() const { return form_; }
    int scalar_bytes() const { return scalar_bytes_; }
    bool symmetric() const { return (flags_ & BPACK_GPU_SYMMETRIC) != 0; }
    bool needs_coordinates() const { return (flags_ & BPACK_GPU_COORDINATES) != 0; }
    // The sketch's kernel rows can be evaluated inside the sketch.
    bool inline_entries() const { return form_ == Form::entry; }
    const std::string& origin() const { return origin_; }
    // Seconds of the NVRTC compile of an entry evaluator given as source (0
    // before its first use, and for the other forms)
    double compile_seconds() const { return compile_seconds_; }

    // out = K for the eval items (on the device at `items`; host_items: the
    // same items, for their sizes and outputs).
    template<typename T>
    void eval(const std::vector<EvalItemT<T>>& host_items, const EvalItemT<T>* items, int max_m, int max_n,
              const char* meta, const PointTable& points, cudaStream_t stream) const {
        if (host_items.empty() || max_m <= 0 || max_n <= 0) return;
        check_type<T>();
        if (form_ == Form::entry) {
            if (!source_.empty()) {
                const NvrtcEntryKernels& kernels = nvrtc();  // (sets params_ at the first call)
                nvrtc_launch_eval(kernels, items, static_cast<int>(host_items.size()), max_m, max_n, meta, points,
                                  params_, stream);
            } else {
                launchers_.eval(items, static_cast<int>(host_items.size()), max_m, max_n, meta, &points,
                                entry_.data(), stream);
            }
            check_launch("the entry evaluator");
            return;
        }
        eval_blocks(host_items, items, meta, points, stream);
    }

    // The ordered sketch of kernel rows (OrderedSketchItemT, value(row, c) =
    // K(point col_base + c, point row_slots[row])).
    template<typename T>
    void sketch_rows(const std::vector<OrderedSketchItemT<T>>& host_items, const OrderedSketchItemT<T>* items,
                     int max_d, int max_cols, const PointTable& points, cudaStream_t stream) const {
        if (host_items.empty() || max_d <= 0 || max_cols <= 0) return;
        check_type<T>();
        if (form_ == Form::entry) {
            if (!source_.empty()) {
                const NvrtcEntryKernels& kernels = nvrtc();  // (sets params_ at the first call)
                nvrtc_launch_sketch(kernels, items, static_cast<int>(host_items.size()), max_d, max_cols, points,
                                    params_, stream);
            } else {
                launchers_.sketch(items, static_cast<int>(host_items.size()), max_d, max_cols, &points, entry_.data(),
                                  stream);
            }
            check_launch("the entry evaluator's sketch");
            return;
        }
        // the rows into scratch (ncols x rows per item, row by row), then
        // their sketch, in groups within a budget
        DeviceHeap& heap = DeviceHeap::instance();
        const size_t budget = scratch_budget(heap);
        size_t i = 0;
        while (i < host_items.size()) {
            size_t j = i, bytes = 0;
            while (j < host_items.size()) {
                const OrderedSketchItemT<T>& it = host_items[j];
                const size_t b = align_up(static_cast<size_t>(std::max(it.ncols, 0)) * std::max(it.rows, 0) * sizeof(T));
                if (j > i && bytes + b > budget) break;
                bytes += b;
                ++j;
            }
            char* scratch = bytes > 0 ? heap.alloc(bytes) : nullptr;
            std::vector<EvalItemT<T>> evals;
            std::vector<OrderedSketchItemT<T>> stored;
            int max_m = 0, max_n = 0, md = 0, mc = 0;
            size_t at = 0;
            for (size_t k = i; k < j; ++k) {
                OrderedSketchItemT<T> s = host_items[k];
                T* block = reinterpret_cast<T*>(scratch + at);
                at += align_up(static_cast<size_t>(std::max(s.ncols, 0)) * std::max(s.rows, 0) * sizeof(T));
                if (s.ncols > 0 && s.rows > 0) {
                    // block(c, row) = K(point col_base + c, point row_slots[row])
                    EvalItemT<T> e;
                    e.out = block;
                    e.ld = s.ncols;
                    e.m = s.ncols;
                    e.n = s.rows;
                    e.rows = IndexList{-1, s.col_base, nullptr};
                    e.cols = IndexList{-1, 0, s.row_slots};
                    evals.push_back(e);
                    max_m = std::max(max_m, e.m);
                    max_n = std::max(max_n, e.n);
                }
                s.src = block;
                s.row_stride = s.ncols;
                s.row_index = nullptr;
                stored.push_back(s);
                md = std::max(md, s.d);
                mc = std::max(mc, s.ncols);
            }
            if (!evals.empty()) {
                EvalItemT<T>* d_evals = upload(heap, evals, stream);
                eval_blocks(evals, d_evals, nullptr, points, stream);
                heap.free(d_evals);
            }
            OrderedSketchItemT<T>* d_stored = upload(heap, stored, stream);
            launch_ordered_sketch_stored(d_stored, static_cast<int>(stored.size()), md, mc, stream);
            heap.free(d_stored);
            if (scratch != nullptr) heap.free(scratch);  // later launches are ordered after the sketch
            i = j;
        }
    }

    // Before the evaluator's first use: K on the block of points 0 .. count -
    // 1 of `points` with themselves, and the sketch of their kernel rows, so
    // that its first-use costs (an NVRTC compile, the loading of the kernels
    // instantiated for it, the application's own, e.g. CuPy compiling its
    // kernels) are not paid inside a level.  Once per evaluator; returns its
    // seconds (0 when it already ran or count is 0).
    template<typename T>
    double warm_up(const PointTable& points, int count, cudaStream_t stream) const {
        if (warmed_ || count <= 0) return 0.0;
        warmed_ = true;
        check_type<T>();
        const auto t0 = std::chrono::steady_clock::now();
        DeviceHeap& heap = DeviceHeap::instance();
        const size_t n = static_cast<size_t>(count);
        T* block = heap.alloc<T>(n * n * sizeof(T));
        T* sketch = heap.alloc<T>(n * sizeof(T));
        std::vector<EvalItemT<T>> evals{EvalItemT<T>{block, count, count, count, IndexList{-1, 0, nullptr},
                                                     IndexList{-1, 0, nullptr}}};
        EvalItemT<T>* d_evals = upload(heap, evals, stream);
        eval<T>(evals, d_evals, count, count, nullptr, points, stream);
        // the kernel rows of points 0 .. count - 1 at columns 0 .. count - 1,
        // one draw (row 0, destination 0, owner warp 0)
        std::vector<int> lists(kSketchOwners + 2 + n, 1);
        lists[0] = 0;               // ptr: warp 0 owns entry 0, the others none
        lists[kSketchOwners + 1] = 0;  // the entry: row 0, destination 0, +
        for (size_t r = 0; r < n; ++r) lists[kSketchOwners + 2 + r] = static_cast<int>(r);
        int* d_lists = upload(heap, lists, stream);
        OrderedSketchItemT<T> s{};
        s.out = sketch;
        s.ldo = 1;
        s.d = 1;
        s.ncols = count;
        s.rows = count;
        s.ptr = d_lists;
        s.entries = d_lists + kSketchOwners + 1;
        s.scale = 1.0;
        s.row_slots = d_lists + kSketchOwners + 2;
        s.col_base = 0;
        std::vector<OrderedSketchItemT<T>> sketches{s};
        OrderedSketchItemT<T>* d_sketches = upload(heap, sketches, stream);
        sketch_rows<T>(sketches, d_sketches, 1, count, points, stream);
        check_cuda(cudaStreamSynchronize(stream), "GPU evaluator warm-up");
        for (void* p : {static_cast<void*>(d_sketches), static_cast<void*>(d_lists), static_cast<void*>(d_evals),
                        static_cast<void*>(sketch), static_cast<void*>(block)}) {
            heap.free(p);
        }
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }

private:
    Evaluator(Form form, int flags, int scalar_bytes) : form_(form), flags_(flags), scalar_bytes_(scalar_bytes) {}

    template<typename T>
    void check_type() const {
        if (sizeof(T) != static_cast<size_t>(scalar_bytes_)) {
            throw std::logic_error("GPU evaluator: value type does not match the matrix");
        }
    }
    static void check_launch(const char* what) {
        const cudaError_t status = cudaGetLastError();
        if (status != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA error after ") + what + ": " + cudaGetErrorString(status));
        }
    }
    static size_t scratch_budget(const DeviceHeap& heap) {
        return std::max<size_t>(size_t{64} << 20, std::min<size_t>(size_t{1} << 30, heap.largest_free() / 4));
    }
    template<typename V>
    static V* upload(DeviceHeap& heap, const std::vector<V>& v, cudaStream_t stream) {
        V* d = heap.alloc<V>(std::max<size_t>(v.size(), 1) * sizeof(V));
        if (!v.empty()) {
            // (pageable: the host data is staged before the call returns)
            check_cuda(cudaMemcpyAsync(d, v.data(), v.size() * sizeof(V), cudaMemcpyHostToDevice, stream),
                       "GPU evaluator upload");
        }
        return d;
    }
    // the scratch of the application's block evaluator: from the heap
    static void* allocate_scratch(void*, size_t bytes) {
        try {
            return DeviceHeap::instance().alloc(std::max<size_t>(bytes, 1));
        } catch (const std::exception&) {
            return nullptr;
        }
    }
    static void release_scratch(void*, void* ptr) {
        if (ptr != nullptr) DeviceHeap::instance().free(ptr);
    }

    // The block and list forms: the items' points gathered, then the
    // application's function, in groups within a budget.
    template<typename T>
    void eval_blocks(const std::vector<EvalItemT<T>>& host_items, const EvalItemT<T>*, const char* meta,
                     const PointTable& points, cudaStream_t stream) const {
        DeviceHeap& heap = DeviceHeap::instance();
        const size_t budget = scratch_budget(heap);
        const bool coords = needs_coordinates() && points.xyz != nullptr;
        const bool host_ids = (flags_ & BPACK_GPU_HOST_IDS) != 0;
        const size_t per_point = sizeof(int64_t) + (coords ? sizeof(double) * static_cast<size_t>(points.dim) : 0);
        const size_t per_entry = form_ == Form::list ? 2 * sizeof(int64_t) + sizeof(T) : 0;
        size_t i = 0;
        while (i < host_items.size()) {
            size_t j = i, bytes = 0;
            while (j < host_items.size()) {
                const EvalItemT<T>& it = host_items[j];
                const size_t b = static_cast<size_t>(it.m + it.n) * per_point +
                                 static_cast<size_t>(it.m) * static_cast<size_t>(it.n) * per_entry;
                if (j > i && bytes + b > budget) break;
                bytes += b;
                ++j;
            }
            run_group<T>(host_items, i, j, meta, points, coords, host_ids, stream);
            i = j;
        }
    }

    template<typename T>
    void run_group(const std::vector<EvalItemT<T>>& host_items, size_t i0, size_t i1, const char* meta,
                   const PointTable& points, bool coords, bool host_ids, cudaStream_t stream) const {
        DeviceHeap& heap = DeviceHeap::instance();
        std::vector<BlockGatherItem> gather;
        int64_t points_total = 0, entries = 0;
        int max_mn = 0, max_m = 0, max_n = 0;
        for (size_t k = i0; k < i1; ++k) {
            const EvalItemT<T>& it = host_items[k];
            if (it.m <= 0 || it.n <= 0) continue;
            gather.push_back(BlockGatherItem{it.rows, it.cols, it.m, it.n, points_total, points_total + it.m});
            points_total += it.m + it.n;
            entries += static_cast<int64_t>(it.m) * it.n;
            max_mn = std::max(max_mn, std::max(it.m, it.n));
            max_m = std::max(max_m, it.m);
            max_n = std::max(max_n, it.n);
        }
        if (gather.empty()) return;
        int64_t* ids = heap.alloc<int64_t>(static_cast<size_t>(points_total) * sizeof(int64_t));
        double* xyz = coords ? heap.alloc<double>(static_cast<size_t>(points_total) * points.dim * sizeof(double))
                             : nullptr;
        BlockGatherItem* d_gather = upload(heap, gather, stream);
        launch_gather_block_points(d_gather, static_cast<int>(gather.size()), max_mn, meta, points, ids, xyz, stream);
        check_launch("the gather of a batch's points");
        std::vector<int64_t> h_ids;
        if (host_ids) {
            h_ids.resize(static_cast<size_t>(points_total));
            check_cuda(cudaMemcpyAsync(h_ids.data(), ids, h_ids.size() * sizeof(int64_t), cudaMemcpyDeviceToHost,
                                       stream),
                       "GPU evaluator ids");
            check_cuda(cudaStreamSynchronize(stream), "GPU evaluator ids");
        }
        std::vector<bpack_gpu_block> blocks;
        std::vector<FlatBlockItem> flat;
        size_t g = 0;
        int64_t at = 0;
        for (size_t k = i0; k < i1; ++k) {
            const EvalItemT<T>& it = host_items[k];
            if (it.m <= 0 || it.n <= 0) continue;
            const BlockGatherItem& gi = gather[g++];
            bpack_gpu_block b;
            b.m = it.m;
            b.n = it.n;
            b.ld = it.ld;
            b.out = it.out;
            b.row_ids = ids + gi.row_at;
            b.col_ids = ids + gi.col_at;
            b.row_ids_host = host_ids ? h_ids.data() + gi.row_at : nullptr;
            b.col_ids_host = host_ids ? h_ids.data() + gi.col_at : nullptr;
            b.row_coords = xyz != nullptr ? xyz + gi.row_at * points.dim : nullptr;
            b.col_coords = xyz != nullptr ? xyz + gi.col_at * points.dim : nullptr;
            blocks.push_back(b);
            flat.push_back(FlatBlockItem{it.m, it.n, it.ld, it.out, b.row_ids, b.col_ids, at});
            at += static_cast<int64_t>(it.m) * it.n;
        }
        if (form_ == Form::block) {
            bpack_gpu_block* d_blocks = upload(heap, blocks, stream);
            bpack_gpu_block_batch batch;
            batch.count = static_cast<int>(blocks.size());
            batch.blocks = blocks.data();
            batch.blocks_device = d_blocks;
            batch.max_m = max_m;
            batch.max_n = max_n;
            batch.dim = points.dim;
            batch.stream = stream;
            batch.allocate = &Evaluator::allocate_scratch;
            batch.release = &Evaluator::release_scratch;
            batch.allocator = nullptr;
            block_(&batch, user_);
            check_launch("the block evaluator");
            heap.free(d_blocks);
        } else {
            FlatBlockItem* d_flat = upload(heap, flat, stream);
            int64_t* rows = heap.alloc<int64_t>(static_cast<size_t>(entries) * sizeof(int64_t));
            int64_t* cols = heap.alloc<int64_t>(static_cast<size_t>(entries) * sizeof(int64_t));
            T* values = heap.alloc<T>(static_cast<size_t>(entries) * sizeof(T));
            launch_flatten_entries(d_flat, static_cast<int>(flat.size()), max_m, max_n, rows, cols, stream);
            check_launch("the entry list of a batch");
            list_(entries, rows, cols, values, stream, user_);
            check_launch("the entry-list evaluator");
            launch_scatter_entries<T>(d_flat, static_cast<int>(flat.size()), max_m, max_n, values, stream);
            check_launch("the entry list's values into their blocks");
            heap.free(values);
            heap.free(cols);
            heap.free(rows);
            heap.free(d_flat);
        }
        heap.free(d_gather);
        if (xyz != nullptr) heap.free(xyz);
        heap.free(ids);  // later launches are ordered after the evaluator's
    }

    const NvrtcEntryKernels& nvrtc() const {
        if (!nvrtc_) {
            Context::instance().activate();
            const auto t0 = std::chrono::steady_clock::now();
            nvrtc_ = compile_entry_source(source_, scalar_bytes_ == static_cast<int>(sizeof(dcomplex)));
            compile_seconds_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            const size_t bytes = std::max<size_t>(host_params_.size(), 1) * sizeof(double);
            check_cuda(cudaMalloc(reinterpret_cast<void**>(&params_), bytes), "GPU entry evaluator parameters");
            if (!host_params_.empty()) {
                check_cuda(cudaMemcpy(params_, host_params_.data(), host_params_.size() * sizeof(double),
                                      cudaMemcpyHostToDevice),
                           "GPU entry evaluator parameters");
            }
        }
        return *nvrtc_;
    }

    Form form_;
    int flags_ = 0;
    int scalar_bytes_ = 0;
    std::string origin_;
    // entry, C++
    bpack_gpu_entry_launchers launchers_{};
    std::vector<char> entry_;
    // entry, source
    std::string source_;
    std::vector<double> host_params_;
    mutable std::shared_ptr<NvrtcEntryKernels> nvrtc_;
    mutable double* params_ = nullptr;
    mutable bool warmed_ = false;  // (warm_up ran)
    mutable double compile_seconds_ = 0.0;
    // block, list
    bpack_gpu_block_evaluator block_ = nullptr;
    bpack_gpu_list_evaluator list_ = nullptr;
    void* user_ = nullptr;
};

// The evaluator a matrix registered (H2Kernel::gpu_evaluator, opaque outside
// the GPU backend), or null.
inline const Evaluator* evaluator_of(const std::shared_ptr<void>& registered) {
    return static_cast<const Evaluator*>(registered.get());
}

// The registered evaluator if the H2 backend can use it for DataType (double
// or complex<double>), else null and why: none, of the other value type, or
// not symmetric (H2 factors symmetric matrices).
template<typename DataType>
const Evaluator* registered_evaluator(const std::shared_ptr<void>& registered, const char** why) {
    const Evaluator* e = evaluator_of(registered);
    if (e == nullptr) {
        *why = "no GPU evaluator registered (doc/gpu_kernels.md)";
        return nullptr;
    }
    if (e->scalar_bytes() != static_cast<int>(sizeof(DataType))) {
        *why = "the GPU evaluator's value type does not match the matrix";
        return nullptr;
    }
    if (!e->symmetric()) {
        *why = "the GPU evaluator is not marked symmetric (BPACK_GPU_SYMMETRIC), which the H2 format requires";
        return nullptr;
    }
    return e;
}

// The device sketch covers a level's ID targets (factorization and
// compression).  Static training rows (H2_ID_radius > 2, H2_ID_proxy 1:
// points anywhere in the tree) join the level's point table: an evaluator
// that reads coordinates needs their global coordinates (kept for
// H2_ID_proxy 1), others only their ids.  Adaptive rows (H2_ID_proxy 2) are
// not streamed.
template<typename Tree>
bool device_sketch_supported(const Tree* tree, const Evaluator& evaluator) {
    if (tree->id_proxy_mode == 2) return false;
    if (tree->id_neighborhood_radius <= 2 && tree->id_proxy_mode != 1) return true;
    const size_t coordinates = static_cast<size_t>(tree->num_points) * static_cast<size_t>(tree->dimension);
    return !evaluator.needs_coordinates() || tree->id_source_point_coords.size() == coordinates;
}

// The warm-up of the evaluator an H2 matrix registered (Evaluator::warm_up)
// before a compression or factorization uses it: on up to 32 points of the
// rank's first non-empty leaf box.  Its seconds (0: no usable evaluator, it
// is warm already, or the rank has no points).
template<typename DataType, typename Tree>
double warm_up_registered_evaluator(const Tree* tree, const std::shared_ptr<void>& registered) {
    const char* why = nullptr;
    const Evaluator* evaluator = registered_evaluator<DataType>(registered, &why);
    if (evaluator == nullptr || tree->num_levels < 1) return 0.0;
    const size_t dim = static_cast<size_t>(tree->dimension);
    std::vector<int64_t> ids;
    std::vector<double> xyz;
    for (const auto& box : tree->levels[static_cast<size_t>(tree->num_levels - 1)].local_boxes) {
        if (box.num_points <= 0) continue;
        const int64_t n = std::min<int64_t>(box.num_points, 32);
        for (int64_t i = 0; i < n; ++i) {
            ids.push_back(box.point_indices[static_cast<size_t>(i)]);
            for (size_t d = 0; d < dim; ++d) xyz.push_back(static_cast<double>(box.point_coords[dim * i + d]));
        }
        break;
    }
    if (ids.empty()) return 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    Context& ctx = Context::instance();
    ctx.activate();
    DeviceHeap& heap = DeviceHeap::instance();
    int64_t* d_ids = heap.alloc<int64_t>(ids.size() * sizeof(int64_t));
    double* d_xyz = heap.alloc<double>(std::max<size_t>(xyz.size(), 1) * sizeof(double));
    check_cuda(cudaMemcpy(d_ids, ids.data(), ids.size() * sizeof(int64_t), cudaMemcpyHostToDevice), "warm-up points");
    check_cuda(cudaMemcpy(d_xyz, xyz.data(), xyz.size() * sizeof(double), cudaMemcpyHostToDevice), "warm-up points");
    PointTable points;
    points.xyz = d_xyz;
    points.ids = d_ids;
    points.dim = tree->dimension;
    const double evaluated = evaluator->template warm_up<typename DeviceScalar<DataType>::type>(
        points, static_cast<int>(ids.size()), ctx.stream());
    heap.free(d_xyz);
    heap.free(d_ids);
    return evaluated > 0.0 ? std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() : 0.0;
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
