#pragma once
// Device-resident elimination of the Color waves of one level (H2 GPU backend).
//
// For the whole level the GPU holds the matrix blocks that later waves modify:
// the near-field block of every one-hop pair (one physical matrix per edge,
// rows = the higher Morton box, columns = the lower one) and the self (Schur)
// block of every box.  Blocks always have their current shape: a box's rows
// or columns are its full point set until it is eliminated and its skeleton
// afterwards, so a block never needs slicing when it is read.
//
// Per wave, the host computes the sketch and the ID of every box (the CPU
// code), plans all device work, uploads one metadata image, and the device
// runs the elimination batched over the wave's boxes: X_RR, X_SR, X_NR and
// the LU of X_RR, temp1 / temp2, the Schur block, the near-block update of
// step 5, and the deferred X_NN owner pass.  The host then receives the
// factors each box keeps for the solve and for the sketches of later waves
// (LU and pivots, temp1, temp2, X_RR_full, X_RS_entry).  At level end the
// blocks come back to the host so the CPU level transition runs unchanged.
//
// Scope (checked by level_eliminator_supported): real or complex symmetric
// kernels with a registered device form, streamed sketches with lazy far
// fill, owner-deferred X_NN updates, LU of X_RR.  Complex data is symmetric,
// not Hermitian: every transpose is plain.  The arithmetic matches the CPU
// path up to rounding (GEMM summation order differs).

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "device_kernels.hpp"
#include "device_solve.hpp"
#include "solve_kernels.hpp"
#include "gpu_runtime.hpp"
#include "host_copier.hpp"
#include "kernel_tables.hpp"
#include "emsurf_blocks.hpp"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <deque>
#include <map>
#include <exception>
#include <cmath>
#include <iomanip>
#include <limits>
#include <random>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fmm {
namespace gpu {

// Per-level totals of the device box path (seconds and bytes, this rank).
struct EliminatorStats {
    double begin = 0.0;       // point table and block upload
    double sketch = 0.0;      // device sketch (plan, launches, Y download)
    double sketch_plan = 0.0; //   host index lists
    double sketch_gpu = 0.0;  //   metadata build and upload, kernels, Y download
    double sk_wait = 0.0;     // host blocked on the sketch (and the previous wave's device work)
    double sk_meta = 0.0, sk_upload = 0.0, sk_rows = 0.0, sk_stored = 0.0, sk_p = 0.0, sk_fill = 0.0,
           sk_id = 0.0, sk_download = 0.0;  //   ... its parts (device times from events)
    double el[8] = {0, 0, 0, 0, 0, 0, 0, 0};  // device: fills, X_RR/X_SR, LU, X_NR, solves, Schur+step 5, owner targets, owner GEMMs
    double finish_sources = 0.0;  // level end: wait for the background copies
    double finish_store = 0.0;    // background copier: busy time (overlapped)
    double copier_wait = 0.0;     //   ... of which waiting for device chunks
    double id = 0.0;          // host ID (and host sketch without the device sketch)
    double plan = 0.0;        // host planning and metadata
    double plan_boxes = 0.0, plan_owner = 0.0;  //   ... box region, owner pass
    double reclaim_wait = 0.0;    // full heap: waits for background copies to free blocks
    double exchange_wait = 0.0;   // full exchange arena: waits for receive buffers' host copies
    double launch = 0.0;                        // host time issuing a wave's device work
    double device = 0.0;      // device elimination and owner pass (device time)
    double download = 0.0;    // factor download
    double store = 0.0;       // host BoxData update
    double finish = 0.0;      // level end (wait for the background copies)
    double transition = 0.0;  // device transition to the parent level
    double tr_plan = 0.0, tr_blocks = 0.0, tr_p = 0.0, tr_fill = 0.0;  // its parts
    double tr_restore = 0.0, tr_restore_bytes = 0.0;  // restored fill sources (host gather, in plan)
    int64_t tr_chunks = 0;
    int64_t transition_fill_gemms = 0;
    double owner_flops = 0.0, transition_flops = 0.0, solve_flops = 0.0;
    double bytes_up = 0.0;
    double bytes_down = 0.0;
    int64_t boxes = 0;
    int64_t owner_gemms = 0;
    int64_t owner_batches = 0;
    int64_t new_targets = 0;
    int64_t heap_reclaims = 0;  // allocations that waited for the background copies
    double remote = 0.0;          // host time of the remote-generator updates (multi-rank)
    double emit = 0.0, exchange = 0.0;  // generator emission, transport between ranks
    bool device_exchange = false;       // generators sent from device memory
    int64_t exchange_fallbacks = 0;     // MPI buffers the exchange arena could not hold
    int64_t remote_generators = 0;
    size_t heap_peak = 0;
};

inline EliminatorStats& eliminator_stats() {
    static EliminatorStats stats;
    return stats;
}

template<typename T>
struct DeviceMatrixT {
    T* ptr = nullptr;
    int rows = 0;
    int cols = 0;
};

// Device blocks of one level, in the eliminator's layout: the Schur block of
// every box (by Morton index) and the near-field block of every one-hop pair
// (by (lo, hi) key: rows = hi's points, columns = lo's).  Built by the device
// transition and adopted by the next level's eliminator, or copied to the
// host (download_level_blocks) when that level runs on the host.
template<typename T>
struct DeviceLevelBlocksT {
    std::unordered_map<int64_t, DeviceMatrixT<T>> schur;
    std::unordered_map<uint64_t, DeviceMatrixT<T>> edges;
};
// The device blocks of a factorization's DataType (double or dcomplex elements).
template<typename DataType>
using DeviceLevelBlocks = DeviceLevelBlocksT<typename DeviceScalar<DataType>::type>;

// Handle of a level's eliminator, independent of the kernel type.
template<typename CoordType, typename DataType>
class LevelEliminatorBase {
public:
    virtual ~LevelEliminatorBase() = default;
    // Eliminates one wave; returns the number of boundary boxes in it.
    virtual int eliminate_wave(const std::vector<int64_t>& wave, int wave_index) = 0;
    // Waits for the factor copies; the level's blocks stay on the device.
    virtual void finish() = 0;
    // Near-field and Schur blocks into the host BoxData (host transition).
    virtual void download_blocks() = 0;
    // The parent level's blocks, assembled on the device from this level's
    // (the device transition); this level's device data is released.
    virtual bool can_build_parent() const = 0;
    // A level that is not eliminated (level 1): its boxes keep all points as
    // skeleton, so the adopted blocks go straight into the transition.
    virtual void adopt_without_elimination() = 0;
    // Multi-rank levels, around each transport of the Color loop: the
    // generators of the last wave's boxes with remote neighbors (their host
    // copies are complete on return), then the device updates from the
    // installed remote generators.
    virtual void emit_generators(PendingFactorUpdates<DataType>& pending) = 0;
    virtual void receive_remote(const std::vector<int64_t>& installed) = 0;
    // Multi-rank levels with CUDA-aware MPI: the whole transport instead,
    // generators sent from and received into device memory.  Returns the
    // MPI time.
    virtual bool device_exchange() const = 0;
    virtual std::chrono::high_resolution_clock::duration exchange() = 0;
    // The solve's factors copied to a device block as each wave produces
    // them (device_solve.hpp), while the heap has room: whether this rank
    // still holds all of the level's, then the verdict all ranks agreed on.
    virtual bool solve_factors_kept() const = 0;
    virtual void commit_solve_factors(bool keep) = 0;
    virtual std::unique_ptr<DeviceLevelBlocks<DataType>> build_parent(
        std::vector<BoxData<CoordType, DataType>>& parents) = 0;
};

// The device box path covers this level (the caller has already selected a
// streamed, lazy, owner-deferred Color level).
template<typename CoordType, typename DataType, typename KernelType>
bool level_eliminator_supported(const TreeLevel<CoordType, DataType>& level, const KernelType* kernel,
                                int dimension, FactorizationMethod method, std::string* reason) {
    auto fail = [&](const char* why) {
        if (reason) *reason = why;
        return false;
    };
    static_assert(std::is_same_v<DataType, double> || std::is_same_v<DataType, std::complex<double>>,
                  "the GPU box path supports double and complex<double>");
    if (dimension != 3) return fail("the device kernels are 3D");
    const char* why = nullptr;
    if (!device_kernel_registered<DataType>(kernel->gpu_spec, &why)) return fail(why);
    if (method != FactorizationMethod::LU) return fail("X_RR is factored by LU on the GPU (use H2_XRR_factor=1)");
    if (level.num_active_processes != 1) {
        // other ranks exchange lazy generators only (generated near blocks)
        if (!level.ghost_boxes.empty()) return fail("ghost boxes (CA levels) are not ported");
        if (!(generator_near_enabled() && lazy_far_field_mode() == LazyFarFieldMode::LAZY)) {
            return fail("multi-rank levels need H2_lazy_schur=2");
        }
    }
    return true;
}

// The device sketch covers the level's ID targets, unless H2_GPU_SKETCH=0
// keeps the host sketch (for comparisons).  Static training rows
// (H2_ID_radius > 2, H2_ID_proxy 1: points anywhere in the tree) join the
// level's point table: a kernel of coordinates needs their global
// coordinates (kept for H2_ID_proxy 1), kind 3 only their ids.  Adaptive rows
// (H2_ID_proxy 2) are not streamed.
template<typename CoordType, typename DataType>
bool device_sketch_supported(const ParallelTree<CoordType, DataType>* tree, int kernel_kind) {
    static const bool enabled = [] {
        const char* v = std::getenv("H2_GPU_SKETCH");
        return v == nullptr || std::atoi(v) != 0;
    }();
    if (!enabled || tree->id_proxy_mode == 2) return false;
    if (tree->id_neighborhood_radius <= 2 && tree->id_proxy_mode != 1) return true;
    const size_t coordinates = static_cast<size_t>(tree->num_points) * static_cast<size_t>(tree->dimension);
    return kernel_kind == 3 || tree->id_source_point_coords.size() == coordinates;
}

template<typename CoordType, typename DataType, typename KernelType>
class LevelEliminator : public LevelEliminatorBase<CoordType, DataType> {
    // Device element type (double, or dcomplex for std::complex<double>)
    // and the launch items and blocks over it.
    using S = typename DeviceScalar<DataType>::type;
    using DeviceMatrix = DeviceMatrixT<S>;
    using Blocks = DeviceLevelBlocksT<S>;
    using EvalItem = EvalItemT<S>;
    using GatherItem = GatherItemT<S>;
    using AddStoreItem = AddStoreItemT<S>;
    using IdentityItem = IdentityItemT<S>;
    using SymAddItem = SymAddItemT<S>;
    using ColumnSwapItem = ColumnSwapItemT<S>;
    using OrderedSketchItem = OrderedSketchItemT<S>;
    using TransposeItem = TransposeItemT<S>;
    using QrcpItem = QrcpItemT<S>;
    using SumAddItem = SumAddItemT<S>;
    // flops per multiply-add of S, in units of real ones (statistics)
    static constexpr double kFlopScale = is_complex_scalar<S> ? 4.0 : 1.0;

public:
    using Box = BoxData<CoordType, DataType>;
    using Level = TreeLevel<CoordType, DataType>;
    using Tree = ParallelTree<CoordType, DataType>;

    // occupancy: a tree of occupied boxes (color_unstructured), whose local
    // slabs also hold empty boxes; they are never in a wave nor one hop from
    // a box, and stay eliminated with an empty skeleton.  Its transports
    // follow the unstructured host transport's protocol (exchange()).
    LevelEliminator(Tree* tree, int level_index, KernelType* kernel, double tolerance,
                    std::unique_ptr<Blocks> adopt = nullptr, bool occupancy = false, bool device_sketch = true)
        : tree_(tree), level_(tree->levels[static_cast<size_t>(level_index)]), kernel_(kernel),
          tolerance_(tolerance), heap_(DeviceHeap::instance()), adopt_(std::move(adopt)), occupancy_(occupancy) {
        spec_ = device_kernel_spec(kernel->gpu_spec);
        gpu_sketch_ = device_sketch && device_sketch_supported(tree, spec_.kind);
        early_free_sources_ = level_.num_active_processes > 1;
        level_index_ = level_index;
        // (the host sketch would read remote generators' X_NR during the level)
        device_exchange_ = level_.num_active_processes > 1 && device_exchange_enabled() && gpu_sketch_;
    }

    bool device_sketch() const { return gpu_sketch_; }

    LevelEliminator(const LevelEliminator&) = delete;
    LevelEliminator& operator=(const LevelEliminator&) = delete;

    ~LevelEliminator() override {
        if (heap_initialized_) heap_.set_reclaimer(nullptr);
        copier_.reset();  // joins the worker
        exchange_copier_.reset();
        drop_kept();
        release_level_data();
    }

    void begin();
    int eliminate_wave(const std::vector<int64_t>& wave, int wave_index) override;
    void finish() override;
    void download_blocks() override { download_blocks_to_host(); }
    // needs the device fill sources (freed ones are restored from their
    // host copies), unless nothing was eliminated
    bool can_build_parent() const override { return gpu_sketch_ || !eliminated_any_; }
    void adopt_without_elimination() override {
        for (BoxState& st : states_) {
            st.eliminated = true;
            st.skeleton.resize(static_cast<size_t>(st.n));
            for (int i = 0; i < st.n; ++i) st.skeleton[static_cast<size_t>(i)] = i;
            st.skeleton_identity = true;
        }
    }
    std::unique_ptr<Blocks> build_parent(std::vector<Box>& parents) override;
    void emit_generators(PendingFactorUpdates<DataType>& pending) override;
    void receive_remote(const std::vector<int64_t>& installed) override;
    bool device_exchange() const override { return device_exchange_; }
    std::chrono::high_resolution_clock::duration exchange() override;
    bool solve_factors_kept() const override { return keep_solve_ && !keep_failed_; }
    void commit_solve_factors(bool keep) override {
        if (keep && keep_solve_ && !keep_failed_) {
            device_solve_store().kept[level_index_] = std::move(kept_);
            kept_ = KeptSolveLevel{};
        } else {
            drop_kept();
        }
    }

private:
    struct BoxState {
        int slot = 0;                 // first point slot of the box
        int n = 0;
        bool eliminated = false;
        bool skeleton_identity = true;
        std::vector<int> skeleton;    // local positions, when eliminated
        DeviceMatrix schur;           // current x current, or none
        int64_t skeleton_stamp = -1;  // wave whose metadata holds the skeleton list
        int64_t skeleton_offset = -1;
        const int* d_skeleton = nullptr;  // device skeleton positions (device ID), unless identity
        // lazy fill source (eliminated, device sketch): temp2 stored row-major
        // (r x ntot column-major, i.e. temp2^T) and X_RR_full (r x r)
        const S* temp2t = nullptr;
        const S* xrr_full = nullptr;
        int ntot = 0;
        int r = 0;
        // boxes of other ranks (assisting boxes): points appended to the
        // table; eliminated once their skeleton has arrived
        bool remote = false;
        bool on_boundary = false;
        // multi-rank levels: the fill-source block (temp2^T, X_RR_full) is
        // freed once every local neighbor, its only reader, is eliminated
        char* persist = nullptr;
        int readers = 0;
        // static ID training rows of the box: slots after training_slot0_
        const int* d_training = nullptr;
        int ntraining = 0;
        std::vector<int64_t> training_ids;  // their global ids (kind 3 sketch)
    };

    struct WaveBox {
        Box* box = nullptr;
        BoxState* state = nullptr;
        int n = 0, k = 0, r = 0, ntot = 0;
        std::vector<int64_t> counts;   // rows per one-hop neighbor (current sizes)
        std::vector<int> row0;         // first X_NR row of each neighbor
        int64_t skeleton_offset = -1;  // metadata offsets of the local index lists
        int64_t redundant_offset = -1;
        // device buffers
        S* T = nullptr;         // interpolation matrix (k x r), leading dimension ldt
        int ldt = 0;
        S* xrr = nullptr;       // A_RR -> X_RR -> LU(X_RR)
        S* xrr_full = nullptr;  // X_RR before factorization
        S* xsr = nullptr;       // A_SR -> X_SR
        S* s = nullptr;         // A_SS -> new Schur block
        S* tmp1 = nullptr;      // -A_SR^T T
        S* tmp2 = nullptr;      // A_SS T
        S* ans = nullptr;       // A_NS of all neighbors (ntot x k)
        S* xnr = nullptr;       // A_NR -> X_NR (ntot x r)
        S* xns = nullptr;       // temp2 X_SR^T (ntot x k)
        S* temp1 = nullptr;
        S* temp2 = nullptr;
        int* piv = nullptr;
        S* persist_xrr = nullptr;     // X_RR_full, kept for later sketches
        S* persist_temp2t = nullptr;  // temp2^T (r x ntot), kept for later sketches
        // offsets inside the result, source and persist blocks
        size_t off_xrr = 0, off_piv = 0, off_temp1 = 0, off_xrs = 0, off_temp2 = 0, off_xrr_full = 0;
        size_t off_t = 0;            // compact T for the host (device ID)
        int source_group = 0;
        bool has_remote = false;     // a neighbor on another rank: its generator is sent
        size_t off_xnr_orig = 0;     // original X_NR in the source group (generator boxes)
    };

    // (lo, hi) Morton pair of a near-field edge.
    static uint64_t edge_key(int64_t a, int64_t b) {
        const int64_t lo = std::min(a, b), hi = std::max(a, b);
        return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
    }

    bool is_local(int64_t morton) const {
        const int64_t idx = morton - level_.local_morton_start;
        return idx >= 0 && idx < static_cast<int64_t>(states_.size());
    }
    BoxState& state_of(int64_t morton) {
        if (is_local(morton)) return states_[static_cast<size_t>(morton - level_.local_morton_start)];
        auto it = remote_index_.find(morton);
        if (it == remote_index_.end()) {
            throw std::runtime_error("LevelEliminator: box " + std::to_string(morton) +
                                     " is neither local nor an assisting box");
        }
        return remote_states_[it->second];
    }
    // Rows of a ring box of the sketch: its full points or its skeleton.
    struct RingInfo {
        int64_t full = 0;
        const std::vector<int64_t>* skeleton = nullptr;
        bool on_boundary = false;
        bool remote = false;
    };
    RingInfo ring_info(int64_t morton) {
        RingInfo info;
        if (Box* nb = level_.find_local_box(morton)) {
            info.full = nb->num_points;
            info.skeleton = &nb->skeleton_indices;
            info.on_boundary = nb->on_boundary;
            return info;
        }
        auto it = level_.assisting_box_points_for_kernel_evaluation.find(morton);
        if (it == level_.assisting_box_points_for_kernel_evaluation.end()) {
            throw std::runtime_error("LevelEliminator: ring box " + std::to_string(morton) + " has no point data");
        }
        const auto& assist = level_.assisting_boxes[static_cast<size_t>(it->second)];
        info.full = static_cast<int64_t>(assist.indices.size());
        info.skeleton = &assist.skel_indices;
        info.on_boundary = assist.on_boundary;
        info.remote = true;
        return info;
    }
    // A fill-source block of box `st` whose local neighbors among `hop` are
    // its readers (multi-rank levels), or a level allocation.
    void keep_fill_source(BoxState& st, char* block, const std::vector<int64_t>& hop) {
        if (!early_free_sources_) {
            level_allocs_.push_back(block);
            return;
        }
        st.persist = block;
        st.readers = 0;
        for (int64_t m : hop) {
            if (is_local(m) && !state_of(m).eliminated) ++st.readers;
        }
    }
    void set_remote_skeletons(const std::vector<std::pair<BoxState*, const std::vector<int64_t>*>>& skeletons);
    void sync_remote_boxes();
    // a generator received into device memory (device exchange)
    struct Incoming { int r = 0, k = 0, ntot = 0, nx = 0; const char* bulk = nullptr; };
    void apply_remote_generators(const std::vector<int64_t>& installed,
                                 const std::unordered_map<int64_t, Incoming>* incoming = nullptr);
    // Bulk of a generator sent between devices: temp2 (ntot x r), the rows
    // of the original X_NR that belong to the destination's boxes (nx x r,
    // its one-hop slots in order), X_RS (r x k) and X_RR_full (r x r).
    struct GenLayout {
        size_t xnr, xrs, xrr, bytes;
        GenLayout(size_t ntot, size_t r, size_t k, size_t nx) {
            xnr = align_up(ntot * r * sizeof(S));
            xrs = xnr + align_up(nx * r * sizeof(S));
            xrr = xrs + align_up(r * k * sizeof(S));
            bytes = xrr + align_up(r * r * sizeof(S));
        }
    };
    void pack_generators(const std::vector<WaveBox>& boxes, const std::vector<char*>& d_sources, const char* d_xnr,
                         char* d_result, cudaStream_t stream);
    // MPI buffers: from the exchange arena, or the main one when it is full
    char* exchange_alloc(size_t bytes) {
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        if (char* p = arena.try_alloc(bytes)) return p;
        // receive buffers the copier still holds, as their host copies finish
        const auto tw = std::chrono::steady_clock::now();
        auto waited = [&] {
            eliminator_stats().exchange_wait +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - tw).count();
        };
        while (exchange_copier_ && arena.initialized() && arena.capacity() >= bytes) {
            std::vector<char*> done = exchange_copier_->wait_released();
            if (done.empty()) break;
            free_copied_blocks(std::move(done));
            if (char* p = arena.try_alloc(bytes)) {
                waited();
                return p;
            }
        }
        waited();
        ++eliminator_stats().exchange_fallbacks;
        return heap_.alloc(bytes);
    }
    void free_block(const void* ptr) {
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        if (arena.owns(ptr)) {
            arena.free(ptr);
        } else {
            heap_.free(ptr);
        }
    }
    void keep_solve_factors(const std::vector<WaveBox>& boxes, char* d_result, const std::vector<char*>& d_sources,
                            cudaStream_t stream);
    void drop_kept() {
        for (char* b : kept_.blocks) heap_.free(b);
        kept_ = KeptSolveLevel{};
    }
    int owner_of(int64_t morton) {
        auto it = owner_cache_.find(morton);
        if (it != owner_cache_.end()) return it->second;
        const std::vector<uint64_t> single = {static_cast<uint64_t>(morton)};
        const std::vector<uint32_t> region = morton::assign_to_processes_nd(
            tree_->dimension, single, level_.num_active_processes, 1u << level_index_);
        const int rank = level_.morton_to_rank.at(static_cast<int>(region[0]));
        owner_cache_.emplace(morton, rank);
        return rank;
    }
    static int current_size(const BoxState& s) {
        return s.eliminated ? static_cast<int>(s.skeleton.size()) : s.n;
    }
    // Point slots of a box's current rows or columns.
    IndexList current_slots(BoxState& s) { return current_slots_in(s, meta_, wave_stamp_); }
    // Skeleton lists without a device copy go into the image `meta`, once
    // per `stamp` (which names the image being built).
    IndexList current_slots_in(BoxState& s, MetaBuilder& meta, uint64_t stamp) {
        IndexList list;
        list.base = s.slot;
        if (s.eliminated && !s.skeleton_identity && s.d_skeleton != nullptr) {
            list.ptr = s.d_skeleton;
        } else if (s.eliminated && !s.skeleton_identity) {
            if (s.skeleton_stamp != stamp) {
                s.skeleton_offset = static_cast<int64_t>(meta.append(s.skeleton));
                s.skeleton_stamp = stamp;
            }
            list.offset = s.skeleton_offset;
        }
        return list;
    }
    // Device times of the last wave's elimination, once its work is done.
    void collect_elimination_marks();

    void run_ids(const std::vector<int64_t>& wave, std::vector<WaveBox>& boxes, int& boundary_count);
    void sketch_wave(const std::vector<int64_t>& wave, std::vector<WaveBox>& boxes, int& boundary_count);
    void free_copied_blocks(std::vector<char*>&& blocks);
    void download_blocks_to_host();
    void release_level_data();
    void check_device_kernel();
    // Kernel entries of `items` (in the metadata image `meta`, on the device
    // at md + offset): launch_eval, or for kind 3 the triangle-pair
    // evaluator (emsurf_blocks.hpp), whose scratch comes from the heap and is
    // written at once: no block freed ahead of a pending read may be live.
    void eval_blocks(const std::vector<EvalItem>& items, MetaBuilder& meta, const char* md, size_t offset, int max_m,
                     int max_n, cudaStream_t stream);

    Tree* tree_;
    Level& level_;
    KernelType* kernel_;
    double tolerance_;
    DeviceHeap& heap_;
    KernelSpec spec_;
    bool heap_initialized_ = false;
    std::unique_ptr<Blocks> adopt_;  // blocks from the device transition, if any
    std::vector<char*> level_allocs_;            // level-lifetime allocations besides the blocks
    bool gpu_sketch_ = false;
    bool occupancy_ = false;           // empty boxes in the local slab (see the constructor)
    bool early_free_sources_ = false;  // multi-rank level: fill sources released once read
    // Fill sources whose local readers have all sketched, oldest first.  They
    // stay on the device for the transition (which otherwise restores them
    // from their host copies) until the heap needs their memory.
    std::deque<BoxState*> released_sources_;
    bool reclaim_released_source() {
        while (!released_sources_.empty()) {
            BoxState* es = released_sources_.front();
            released_sources_.pop_front();
            if (es->persist == nullptr || es->readers > 0) continue;  // gone, or kept again with readers
            heap_.free(es->persist);  // its readers' launches are all ordered before later ones
            es->persist = nullptr;
            es->temp2t = nullptr;
            es->xrr_full = nullptr;
            return true;
        }
        return false;
    }
    bool eliminated_any_ = false;  // a wave ran (the transition may need fill sources)
    char* wave_sketch_ = nullptr;  // sketches of the current wave; T in place after the device ID
    std::unique_ptr<HostCopier> copier_;  // background copies of the factors into BoxData
    // device exchange: host copies of the received generators, apart from the
    // wave downloads so the receive buffers return to the exchange arena early
    std::unique_ptr<HostCopier> exchange_copier_;

    std::vector<BoxState> states_;
    std::deque<BoxState> remote_states_;               // stable addresses (plans keep pointers)
    std::unordered_map<int64_t, size_t> remote_index_;  // Morton -> remote_states_
    std::unordered_map<uint64_t, DeviceMatrix> edges_;
    PointTable points_;
    int64_t num_slots_ = 0;                            // point slots: local, training, then remote
    int training_slot0_ = 0;                           // first slot of the static ID training points
    double* d_xyz_ = nullptr;
    int64_t* d_ids_ = nullptr;
    std::vector<int64_t> last_wave_;                   // boxes of the last wave (their generators)
    int last_wave_index_ = -1;

    // device exchange (multi-rank levels, CUDA-aware MPI)
    bool device_exchange_ = false;
    int level_index_ = 0;
    struct PeerOut { std::vector<int64_t> header; size_t offset = 0, bytes = 0; };
    char* outbox_ = nullptr;                           // generator bulk of the last wave, a segment per rank
    std::map<int, PeerOut> peer_out_;                  // rank -> its generators (header, segment)
    std::vector<int64_t> since_transport_;             // local boxes eliminated since the last transport
    std::unordered_map<int64_t, std::vector<int>> requesters_;  // local box -> ranks holding it (assisting)
    bool registered_ = false;                          // requests exchanged (first transport of the level)
    struct RemoteHost { std::vector<DataType> temp2, xrr; int ntot = 0, r = 0; };
    std::unordered_map<int64_t, std::unique_ptr<RemoteHost>> remote_host_;  // background host copies
    std::unordered_map<int64_t, int> owner_cache_;

    // the solve's factors kept on the device (device_solve.hpp)
    bool keep_solve_ = false;
    bool keep_failed_ = false;
    KeptSolveLevel kept_;
    DeviceBuffer solve_meta_device_;

    int64_t wave_stamp_ = 0;
    MetaBuilder& meta_ = pinned_pool().meta;
    DeviceBuffer meta_device_;
    MetaBuilder& owner_meta_ = pinned_pool().owner_meta;  // a wave's owner pass
    // kind 3 blocks of eval_blocks: the host ids of the point slots, the host
    // copies of the device skeleton lists, and the lists' image
    std::vector<int64_t> host_ids_;
    std::unordered_map<const int*, const std::vector<int>*> host_lists_;
    MetaBuilder efie_meta_;
    DeviceBuffer efie_meta_device_;
    DeviceBuffer owner_meta_device_;
    std::unique_ptr<StreamMarks> elim_marks_;  // of the last wave, read at the next synchronization
    DeviceBuffer getrf_work_;
};

// ---------------------------------------------------------------------------
// Level start: point table, and the Schur and near-field blocks the level
// transition assembled on the host.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::begin() {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context& ctx = Context::instance();
    ctx.activate();
    if (device_exchange_) {
        // before the main arena takes its share of the free memory
        size_t mb = 2560;
        if (const char* env = std::getenv("H2_GPU_EXCHANGE_MB")) mb = static_cast<size_t>(std::atoll(env));
        DeviceHeap::exchange_arena().initialize_fixed(mb << 20);
    }
    heap_.ensure_initialized();
    // a fresh factorization (unless earlier levels' solve factors are kept);
    // otherwise the parent blocks live on
    if (!adopt_ && device_solve_store().empty()) heap_.reset();
    keep_solve_ = device_solve_enabled() && device_solve_keep() && gpu_sketch_;
    if (keep_solve_) solve_meta_device_.reserve(size_t{1} << 20);  // span lists of the waves
    heap_.reset_peak();
    heap_initialized_ = true;
    copier_ = std::make_unique<HostCopier>(ctx.device());
    if (device_exchange_) exchange_copier_ = std::make_unique<HostCopier>(ctx.device(), 1, 8);
    // a full heap waits for the next background copy, whose block then frees
    heap_.set_reclaimer([this] {
        const auto tw = std::chrono::steady_clock::now();
        bool freed = false;
        for (HostCopier* c : {copier_.get(), exchange_copier_.get()}) {
            if (c == nullptr) continue;
            std::vector<char*> done = c->take_finished();
            freed = freed || !done.empty();
            free_copied_blocks(std::move(done));
        }
        if (!freed) freed = reclaim_released_source();
        if (!freed && copier_) {
            std::vector<char*> done = copier_->wait_released();
            freed = !done.empty();
            free_copied_blocks(std::move(done));
        }
        eliminator_stats().reclaim_wait += std::chrono::duration<double>(std::chrono::steady_clock::now() - tw).count();
        return freed;
    });
    cudaStream_t stream = ctx.stream();
    auto& stats = eliminator_stats();

    // ---- point table: slots are the local boxes' points in box order
    const size_t num_boxes = level_.local_boxes.size();
    states_.assign(num_boxes, BoxState{});
    int64_t num_slots = 0;
    for (size_t b = 0; b < num_boxes; ++b) {
        const Box& box = level_.local_boxes[b];
        if (box.morton_index >= (int64_t{1} << 32)) throw std::runtime_error("LevelEliminator: Morton index exceeds 32 bits");
        if (box.num_points <= 0) {
            if (!occupancy_) throw std::runtime_error("LevelEliminator: empty box " + std::to_string(box.morton_index));
            states_[b].slot = static_cast<int>(num_slots);
            states_[b].eliminated = true;  // empty skeleton
            continue;
        }
        if (!box.far_field_modified_interactions.empty()) {
            throw std::runtime_error("LevelEliminator: stored far-field blocks are not supported (lazy mode expected)");
        }
        states_[b].slot = static_cast<int>(num_slots);
        states_[b].n = static_cast<int>(box.num_points);
        num_slots += box.num_points;
    }
    // static ID training rows of the device sketch (H2_ID_radius > 2,
    // H2_ID_proxy 1), as the host selects them: their points after the local
    // ones, each once, and a list of their slots per box
    std::vector<std::vector<int64_t>> training(num_boxes);
    std::vector<int64_t> training_points;  // global indices, in slot order
    std::vector<int> training_lists;       // per box, slots relative to training_slot0_
    std::vector<size_t> training_at(num_boxes, 0);
    if (gpu_sketch_ && (tree_->id_neighborhood_radius > 2 || tree_->id_proxy_mode == 1)) {
        std::exception_ptr failure;
        std::mutex failure_mutex;
        #pragma omp parallel for schedule(dynamic)
        for (int64_t b = 0; b < static_cast<int64_t>(num_boxes); ++b) {
            if (states_[static_cast<size_t>(b)].n == 0) continue;
            try {
                training[static_cast<size_t>(b)] =
                    select_static_id_training_indices(tree_, &level_.local_boxes[static_cast<size_t>(b)]);
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
        std::unordered_map<int64_t, int> slot_of;
        for (size_t b = 0; b < num_boxes; ++b) {
            training_at[b] = training_lists.size();
            for (int64_t index : training[b]) {
                auto it = slot_of.emplace(index, static_cast<int>(training_points.size())).first;
                if (it->second == static_cast<int>(training_points.size())) training_points.push_back(index);
                training_lists.push_back(it->second);
            }
        }
    }
    training_slot0_ = static_cast<int>(num_slots);
    num_slots += static_cast<int64_t>(training_points.size());
    if (num_slots >= std::numeric_limits<int>::max()) throw std::runtime_error("LevelEliminator: too many points per rank");
    std::vector<double> xyz(static_cast<size_t>(3 * num_slots));
    std::vector<int64_t> ids(static_cast<size_t>(num_slots));
    for (size_t b = 0; b < num_boxes; ++b) {
        const Box& box = level_.local_boxes[b];
        const size_t slot = static_cast<size_t>(states_[b].slot);
        for (int64_t i = 0; i < box.num_points; ++i) {
            for (int d = 0; d < 3; ++d) xyz[3 * (slot + i) + d] = static_cast<double>(box.point_coords[3 * i + d]);
            ids[slot + i] = box.point_indices[static_cast<size_t>(i)];
        }
    }
    {
        // (kind 3 reads only the ids; the global coordinates exist for H2_ID_proxy 1)
        const auto& coords = tree_->id_source_point_coords;
        const bool have_coords = coords.size() == static_cast<size_t>(tree_->num_points) * tree_->dimension;
        for (size_t t = 0; t < training_points.size(); ++t) {
            const size_t slot = static_cast<size_t>(training_slot0_) + t;
            const int64_t index = training_points[t];
            ids[slot] = index;
            for (int d = 0; d < 3; ++d) {
                xyz[3 * slot + d] = have_coords && d < tree_->dimension
                    ? static_cast<double>(coords[static_cast<size_t>(index * tree_->dimension + d)]) : 0.0;
            }
        }
    }
    double* d_xyz = heap_.alloc_resident<double>(xyz.size() * sizeof(double));
    int64_t* d_ids = heap_.alloc_resident<int64_t>(ids.size() * sizeof(int64_t));
    check_cuda(cudaMemcpyAsync(d_xyz, xyz.data(), xyz.size() * sizeof(double), cudaMemcpyHostToDevice, stream), "upload points");
    check_cuda(cudaMemcpyAsync(d_ids, ids.data(), ids.size() * sizeof(int64_t), cudaMemcpyHostToDevice, stream), "upload ids");
    points_.xyz = d_xyz;
    points_.ids = d_ids;
    d_xyz_ = d_xyz;  // the table grows when remote boxes arrive (sync_remote_boxes)
    d_ids_ = d_ids;
    num_slots_ = num_slots;
    if (spec_.kind == 3) host_ids_ = ids;
    stats.bytes_up += static_cast<double>(xyz.size() * sizeof(double) + ids.size() * sizeof(int64_t));
    if (!training_lists.empty()) {
        int* d_lists = heap_.alloc_resident<int>(training_lists.size() * sizeof(int));
        level_allocs_.push_back(reinterpret_cast<char*>(d_lists));
        check_cuda(cudaMemcpyAsync(d_lists, training_lists.data(), training_lists.size() * sizeof(int),
                                   cudaMemcpyHostToDevice, stream), "upload training lists");
        for (size_t b = 0; b < num_boxes; ++b) {
            states_[b].d_training = d_lists + training_at[b];
            states_[b].ntraining = static_cast<int>(training[b].size());
            if (spec_.kind == 3) states_[b].training_ids = training[b];
        }
        stats.bytes_up += static_cast<double>(training_lists.size() * sizeof(int));
    }

    // H2_GPU_KERNEL_CHECK=1: the device kernel against the host's entries
    static const bool kernel_check = [] {
        const char* v = std::getenv("H2_GPU_KERNEL_CHECK");
        return v != nullptr && std::atoi(v) != 0;
    }();
    if (kernel_check) check_device_kernel();

    // ---- Schur and near-field blocks: adopted from the device transition,
    //      or packed on the host (in parallel), uploaded once, and scattered
    //      into their own heap allocations with one batched copy
    if (adopt_) {
        for (size_t b = 0; b < num_boxes; ++b) {
            const Box& box = level_.local_boxes[b];
            if (states_[b].n == 0) continue;  // an empty box (occupancy)
            auto it = adopt_->schur.find(box.morton_index);
            if (it == adopt_->schur.end() || it->second.rows != states_[b].n || it->second.cols != states_[b].n) {
                throw std::runtime_error("LevelEliminator::begin: adopted Schur block of box " +
                                         std::to_string(box.morton_index) + " is missing or mis-shaped");
            }
            states_[b].schur = it->second;
        }
        for (const auto& kv : adopt_->edges) {
            const int64_t lo = static_cast<int64_t>(kv.first >> 32), hi = static_cast<int64_t>(kv.first & 0xffffffffu);
            // a box of another rank is known only from its block here
            if ((is_local(hi) && kv.second.rows != state_of(hi).n) || (is_local(lo) && kv.second.cols != state_of(lo).n) ||
                (!is_local(lo) && !is_local(hi))) {
                throw std::runtime_error("LevelEliminator::begin: adopted near block is mis-shaped");
            }
        }
        edges_ = std::move(adopt_->edges);
        adopt_.reset();
        check_cuda(cudaStreamSynchronize(stream), "level upload");
        stats.begin += std::chrono::duration<double>(clock::now() - t0).count();
        return;
    }
    struct Placement {
        const MatrixStorage<DataType>* schur = nullptr;      // or
        const ModifiedBlock<DataType>* block = nullptr;      // its view, transposed when `transpose`
        bool transpose = false;
        DeviceMatrix* dst = nullptr;
        size_t offset = 0;                                    // in elements
    };
    std::vector<Placement> placements;
    size_t image_elems = 0;
    auto place = [&](Placement p, int rows, int cols) {
        p.dst->rows = rows;
        p.dst->cols = cols;
        p.dst->ptr = heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * sizeof(S));
        p.offset = image_elems;
        image_elems += align_up(static_cast<size_t>(rows) * cols, 32);
        placements.push_back(p);
    };
    for (size_t b = 0; b < num_boxes; ++b) {
        Box& box = level_.local_boxes[b];
        BoxState& st = states_[b];
        if (st.n == 0) continue;  // an empty box (occupancy): no blocks
        if (box.schur_complement.is_allocated()) {
            if (box.schur_complement.rows != st.n || box.schur_complement.cols != st.n ||
                box.schur_complement.lda != st.n) {
                throw std::runtime_error("LevelEliminator::begin: Schur block of box " +
                                         std::to_string(box.morton_index) + " is not n x n");
            }
            Placement p;
            p.schur = &box.schur_complement;
            p.dst = &st.schur;
            place(p, st.n, st.n);
        }
        for (const auto& entry : box.near_field_interaction_map) {
            const int64_t other = entry.first;
            const auto& block = box.near_field_modified_interactions[static_cast<size_t>(entry.second)];
            if (!block.a_ns_is_allocated()) continue;
            const uint64_t key = edge_key(box.morton_index, other);
            if (edges_.count(key)) continue;
            // view A_NS of this box for `other`: rows = other's points, cols =
            // this box's (a remote box's size is known only from its block)
            const int other_n = is_local(other) ? state_of(other).n : static_cast<int>(block.a_ns_rows());
            if (block.a_ns_rows() != other_n || block.a_ns_cols() != st.n) {
                std::ostringstream oss;
                oss << "LevelEliminator::begin: near block (" << box.morton_index << ", " << other
                    << ") is " << block.a_ns_rows() << " x " << block.a_ns_cols() << ", expected "
                    << other_n << " x " << st.n;
                throw std::runtime_error(oss.str());
            }
            Placement p;
            p.block = &block;
            p.transpose = box.morton_index > other;  // physical rows = the higher Morton box
            p.dst = &edges_[key];
            if (p.transpose) place(p, st.n, other_n); else place(p, other_n, st.n);
        }
    }
    if (!placements.empty()) {
        // Pageable image: pinning a multi-GB buffer costs more than the
        // slower pageable copy of a once-per-level upload.
        std::unique_ptr<DataType[]> image_owner(new DataType[image_elems]);
        DataType* h_image = image_owner.get();
        #pragma omp parallel for schedule(dynamic)
        for (int64_t i = 0; i < static_cast<int64_t>(placements.size()); ++i) {
            const Placement& p = placements[static_cast<size_t>(i)];
            DataType* out = h_image + p.offset;
            const int rows = p.dst->rows, cols = p.dst->cols;
            if (p.schur != nullptr) {
                std::memcpy(out, p.schur->data.data(), static_cast<size_t>(rows) * cols * sizeof(DataType));
            } else if (!p.transpose) {
                for (int j = 0; j < cols; ++j)
                    for (int r = 0; r < rows; ++r) out[r + static_cast<size_t>(j) * rows] = p.block->a_ns(r, j);
            } else {
                for (int j = 0; j < cols; ++j)
                    for (int r = 0; r < rows; ++r) out[r + static_cast<size_t>(j) * rows] = p.block->a_ns(j, r);
            }
        }
        S* d_image = heap_.alloc<S>(image_elems * sizeof(S));
        check_cuda(cudaMemcpyAsync(d_image, h_image, image_elems * sizeof(S), cudaMemcpyHostToDevice, stream),
                   "upload blocks");
        std::vector<GatherItem> copies;
        int max_m = 0, max_n = 0;
        for (const Placement& p : placements) {
            copies.push_back(GatherItem{p.dst->ptr, p.dst->rows, p.dst->rows, p.dst->cols, d_image + p.offset, 1,
                                        p.dst->rows, IndexList{}, IndexList{}});
            max_m = std::max(max_m, p.dst->rows);
            max_n = std::max(max_n, p.dst->cols);
        }
        meta_.clear();
        const size_t off = meta_.append(copies);
        char* md = meta_.upload(meta_device_, stream);
        launch_gather(reinterpret_cast<const GatherItem*>(md + off), static_cast<int>(copies.size()), max_m, max_n,
                      md, stream);
        check_cuda(cudaStreamSynchronize(stream), "block upload");
        heap_.free(d_image);
        stats.bytes_up += static_cast<double>(image_elems * sizeof(S));
    }
    check_cuda(cudaStreamSynchronize(stream), "level upload");

    // The device copies are now authoritative; drop the host ones so a stale
    // read fails loudly instead of using outdated values.
    for (auto& box : level_.local_boxes) {
        box.schur_complement = MatrixStorage<DataType>{};
        std::vector<ModifiedBlock<DataType>>().swap(box.near_field_modified_interactions);
        box.near_field_interaction_map.clear();
    }
    stats.begin += std::chrono::duration<double>(clock::now() - t0).count();
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::eval_blocks(const std::vector<EvalItem>& items, MetaBuilder& meta,
                                                                    const char* md, size_t offset, int max_m, int max_n,
                                                                    cudaStream_t stream) {
    if (items.empty()) return;
    if (spec_.kind != 3) {
        launch_eval(reinterpret_cast<const EvalItem*>(md + offset), static_cast<int>(items.size()), max_m, max_n, md,
                    spec_, points_, stream);
        return;
    }
    if constexpr (std::is_same_v<S, dcomplex>) {
        const size_t D = sizeof(S);
        const std::vector<int>& mesh_ints = kernel_->gpu_spec.table_int;
        // the point ids of an index list (as index_at on the device)
        auto ids_of = [&](const IndexList& list, int count, std::vector<int>& out) {
            out.resize(static_cast<size_t>(count));
            const int* host_list = nullptr;
            if (list.ptr != nullptr) {
                auto it = host_lists_.find(list.ptr);
                if (it == host_lists_.end()) throw std::runtime_error("LevelEliminator: device list without a host copy");
                host_list = it->second->data();
            } else if (list.offset >= 0) {
                host_list = reinterpret_cast<const int*>(meta.host(static_cast<size_t>(list.offset)));
            }
            for (int i = 0; i < count; ++i) {
                const int slot = list.base + (host_list != nullptr ? host_list[i] : i);
                out[static_cast<size_t>(i)] = static_cast<int>(host_ids_[static_cast<size_t>(slot)]);
            }
        };
        auto side = [&](const std::vector<int>& edges, std::vector<int>& tri_of, std::vector<int>& tris) {
            tris.clear();
            for (int e : edges) {
                for (int a = 0; a < 2; ++a) {
                    const int t = mesh_ints[6 * static_cast<size_t>(e) + 2 + static_cast<size_t>(a)];
                    if (t >= 0) tris.push_back(t);
                }
            }
            std::sort(tris.begin(), tris.end());
            tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
            tri_of.resize(2 * edges.size());
            for (size_t i = 0; i < edges.size(); ++i) {
                for (int a = 0; a < 2; ++a) {
                    const int t = mesh_ints[6 * static_cast<size_t>(edges[i]) + 2 + static_cast<size_t>(a)];
                    tri_of[2 * i + static_cast<size_t>(a)] =
                        t < 0 ? -1 : static_cast<int>(std::lower_bound(tris.begin(), tris.end(), t) - tris.begin());
                }
            }
        };
        // the blocks' lists (in parallel), then groups of pair sums within a budget
        struct Lists { std::vector<int> rows, cols, row_tri, col_tri, tr, tc; };
        std::vector<Lists> lists(items.size());
        #pragma omp parallel for schedule(dynamic)
        for (int64_t b = 0; b < static_cast<int64_t>(items.size()); ++b) {
            const EvalItem& it = items[static_cast<size_t>(b)];
            Lists& l = lists[static_cast<size_t>(b)];
            if (it.m <= 0 || it.n <= 0) continue;
            ids_of(it.rows, it.m, l.rows);
            ids_of(it.cols, it.n, l.cols);
            side(l.rows, l.row_tri, l.tr);
            side(l.cols, l.col_tri, l.tc);
        }
        efie_meta_.clear();
        std::vector<EfieBlockItem> blocks;
        std::vector<size_t> sums_at;
        std::vector<size_t> group{0};
        std::vector<int64_t> group_pairs, group_entries;
        const size_t budget = std::max<size_t>(size_t{256} << 20,
                                               std::min<size_t>(size_t{2} << 30, heap_.largest_free() / 4));
        size_t group_bytes = 0, scratch_bytes = 0;
        int64_t max_pairs = 0, max_entries = 0;
        for (size_t b = 0; b < items.size(); ++b) {
            const EvalItem& it = items[b];
            const Lists& l = lists[b];
            if (it.m <= 0 || it.n <= 0) continue;
            const size_t sums = align_up(l.tr.size() * l.tc.size() * kEfieSums * D);
            if (group_bytes > 0 && group_bytes + sums > budget) {
                group.push_back(blocks.size());
                group_pairs.push_back(max_pairs);
                group_entries.push_back(max_entries);
                group_bytes = 0;
                max_pairs = max_entries = 0;
            }
            EfieBlockItem e;
            e.nrow = it.m;
            e.ncol = it.n;
            e.ntr = static_cast<int>(l.tr.size());
            e.ntc = static_cast<int>(l.tc.size());
            e.row_edges = static_cast<int64_t>(efie_meta_.append(l.rows));
            e.col_edges = static_cast<int64_t>(efie_meta_.append(l.cols));
            e.row_tri = static_cast<int64_t>(efie_meta_.append(l.row_tri));
            e.col_tri = static_cast<int64_t>(efie_meta_.append(l.col_tri));
            e.tr = static_cast<int64_t>(efie_meta_.append(l.tr));
            e.tc = static_cast<int64_t>(efie_meta_.append(l.tc));
            e.out = reinterpret_cast<dcomplex*>(it.out);  // out(r, c) = out[r + c * ld]
            e.rs = 1;
            e.cs = it.ld;
            sums_at.push_back(group_bytes);
            group_bytes += sums;
            scratch_bytes = std::max(scratch_bytes, group_bytes);
            max_pairs = std::max<int64_t>(max_pairs, static_cast<int64_t>(l.tr.size() * l.tc.size()));
            max_entries = std::max<int64_t>(max_entries, static_cast<int64_t>(it.m) * it.n);
            blocks.push_back(e);
        }
        group.push_back(blocks.size());
        group_pairs.push_back(max_pairs);
        group_entries.push_back(max_entries);
        if (blocks.empty()) return;
        char* scratch = heap_.alloc(std::max<size_t>(scratch_bytes, 1));
        for (size_t b = 0; b < blocks.size(); ++b) blocks[b].M = reinterpret_cast<dcomplex*>(scratch + sums_at[b]);
        const size_t off_blocks = efie_meta_.append(blocks);
        char* emd = efie_meta_.upload(efie_meta_device_, stream);
        for (size_t g = 0; g + 1 < group.size(); ++g) {
            launch_efie_blocks(reinterpret_cast<const EfieBlockItem*>(emd + off_blocks) + group[g],
                               static_cast<int>(group[g + 1] - group[g]), group_pairs[g], group_entries[g], emd, spec_,
                               stream);
        }
        heap_.free(scratch);  // later launches are ordered after these
    } else {
        (void)meta;
    }
}

// Entries of the device kernel against the host kernel's (the application's
// callback), up to 32 x 32 of each of three blocks: a box with itself (the
// singular self terms), with a local neighbor, and with the last local box.
// Prints the largest difference relative to the block's largest entry.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::check_device_kernel() {
    cudaStream_t stream = Context::instance().stream();
    int b0 = -1, b2 = -1;
    for (size_t b = 0; b < states_.size(); ++b) {
        if (states_[b].n == 0) continue;
        if (b0 < 0) b0 = static_cast<int>(b);
        b2 = static_cast<int>(b);
    }
    if (b0 < 0) return;
    int b1 = b0;
    for (int64_t m : level_.local_boxes[static_cast<size_t>(b0)].one_hop) {
        if (m != level_.local_boxes[static_cast<size_t>(b0)].morton_index && is_local(m) && state_of(m).n > 0) {
            b1 = static_cast<int>(m - level_.local_morton_start);
            break;
        }
    }
    const int pairs[3][2] = {{b0, b0}, {b0, b1}, {b0, b2}};
    const char* names[3] = {"self", "near", "far"};
    std::vector<EvalItem> evals;
    std::vector<size_t> offset(4, 0);
    int max_m = 0, max_n = 0;
    for (int p = 0; p < 3; ++p) {
        const BoxState& a = states_[static_cast<size_t>(pairs[p][0])];
        const BoxState& c = states_[static_cast<size_t>(pairs[p][1])];
        const int m = std::min(a.n, 32), n = std::min(c.n, 32);
        offset[static_cast<size_t>(p) + 1] = offset[static_cast<size_t>(p)] + static_cast<size_t>(m) * n;
        max_m = std::max(max_m, m);
        max_n = std::max(max_n, n);
    }
    S* d_out = heap_.alloc<S>(offset[3] * sizeof(S));
    for (int p = 0; p < 3; ++p) {
        const BoxState& a = states_[static_cast<size_t>(pairs[p][0])];
        const BoxState& c = states_[static_cast<size_t>(pairs[p][1])];
        evals.push_back(EvalItem{d_out + offset[static_cast<size_t>(p)], std::min(a.n, 32), std::min(a.n, 32),
                                 std::min(c.n, 32), IndexList{-1, a.slot, nullptr}, IndexList{-1, c.slot, nullptr}});
    }
    meta_.clear();
    const size_t off = meta_.append(evals);
    char* md = meta_.upload(meta_device_, stream);
    launch_eval(reinterpret_cast<const EvalItem*>(md + off), 3, max_m, max_n, md, spec_, points_, stream);
    std::vector<DataType> dev(offset[3]);
    check_cuda(cudaMemcpyAsync(dev.data(), d_out, offset[3] * sizeof(S), cudaMemcpyDeviceToHost, stream), "kernel check");
    check_cuda(cudaStreamSynchronize(stream), "kernel check");
    // kind 3: the same blocks by the triangle-pair evaluator of the sketch
    std::vector<DataType> dev_pairs;
    if constexpr (std::is_same_v<S, dcomplex>) {
        if (spec_.kind == 3) {
            const std::vector<int>& mesh_ints = kernel_->gpu_spec.table_int;
            auto side = [&](const std::vector<int>& edges, std::vector<int>& tri_of, std::vector<int>& tris) {
                for (int e : edges) {
                    for (int q = 0; q < 2; ++q) {
                        const int t = mesh_ints[6 * static_cast<size_t>(e) + 2 + static_cast<size_t>(q)];
                        if (t >= 0) tris.push_back(t);
                    }
                }
                std::sort(tris.begin(), tris.end());
                tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
                for (int e : edges) {
                    for (int q = 0; q < 2; ++q) {
                        const int t = mesh_ints[6 * static_cast<size_t>(e) + 2 + static_cast<size_t>(q)];
                        tri_of.push_back(t < 0 ? -1
                                               : static_cast<int>(std::lower_bound(tris.begin(), tris.end(), t) -
                                                                  tris.begin()));
                    }
                }
            };
            meta_.clear();
            std::vector<EfieBlockItem> items;
            std::vector<size_t> sums_at;
            size_t sums_bytes = 0;
            int64_t max_pairs = 0, max_entries = 0;
            for (int p = 0; p < 3; ++p) {
                const Box& a = level_.local_boxes[static_cast<size_t>(pairs[p][0])];
                const Box& c = level_.local_boxes[static_cast<size_t>(pairs[p][1])];
                const int m = std::min(static_cast<int>(a.num_points), 32), n = std::min(static_cast<int>(c.num_points), 32);
                std::vector<int> test(a.point_indices.begin(), a.point_indices.begin() + m);
                std::vector<int> src(c.point_indices.begin(), c.point_indices.begin() + n);
                std::vector<int> test_tri, src_tri, tt, st;
                side(test, test_tri, tt);
                side(src, src_tri, st);
                EfieBlockItem e;
                e.nrow = m;
                e.ncol = n;
                e.ntr = static_cast<int>(tt.size());
                e.ntc = static_cast<int>(st.size());
                e.row_edges = static_cast<int64_t>(meta_.append(test));
                e.col_edges = static_cast<int64_t>(meta_.append(src));
                e.row_tri = static_cast<int64_t>(meta_.append(test_tri));
                e.col_tri = static_cast<int64_t>(meta_.append(src_tri));
                e.tr = static_cast<int64_t>(meta_.append(tt));
                e.tc = static_cast<int64_t>(meta_.append(st));
                e.rs = 1;
                e.cs = m;
                sums_at.push_back(sums_bytes);
                sums_bytes = align_up(sums_bytes + tt.size() * st.size() * kEfieSums * sizeof(S));
                max_pairs = std::max<int64_t>(max_pairs, static_cast<int64_t>(tt.size() * st.size()));
                max_entries = std::max<int64_t>(max_entries, static_cast<int64_t>(m) * n);
                items.push_back(e);
            }
            char* scratch = heap_.alloc(std::max<size_t>(sums_bytes, 1) + offset[3] * sizeof(S));
            for (int p = 0; p < 3; ++p) {
                items[static_cast<size_t>(p)].M = reinterpret_cast<dcomplex*>(scratch + sums_at[static_cast<size_t>(p)]);
                items[static_cast<size_t>(p)].out =
                    reinterpret_cast<dcomplex*>(scratch + std::max<size_t>(sums_bytes, 1)) + offset[static_cast<size_t>(p)];
            }
            const size_t off_items = meta_.append(items);
            char* md2 = meta_.upload(meta_device_, stream);
            launch_efie_blocks(reinterpret_cast<const EfieBlockItem*>(md2 + off_items), 3, max_pairs, max_entries, md2,
                               spec_, stream);
            dev_pairs.resize(offset[3]);
            check_cuda(cudaMemcpyAsync(dev_pairs.data(), scratch + std::max<size_t>(sums_bytes, 1),
                                       offset[3] * sizeof(S), cudaMemcpyDeviceToHost, stream), "kernel check");
            check_cuda(cudaStreamSynchronize(stream), "kernel check");
            heap_.free(scratch);
        }
    }
    heap_.free(d_out);
    std::ostringstream oss;
    oss << "  [gpu] level " << level_index_ << " kernel check (rank " << tree_->mpi_rank << ", kind " << spec_.kind
        << "): max |K_gpu - K_host| / max |K_host|";
    for (int p = 0; p < 3; ++p) {
        const Box& a = level_.local_boxes[static_cast<size_t>(pairs[p][0])];
        const Box& c = level_.local_boxes[static_cast<size_t>(pairs[p][1])];
        const int m = std::min(static_cast<int>(a.num_points), 32), n = std::min(static_cast<int>(c.num_points), 32);
        std::vector<DataType> host(static_cast<size_t>(m) * n);
        kernel_->evaluate_block_by_index(a.point_indices.data(), m, c.point_indices.data(), n, host.data(), m);
        double diff = 0.0, diff_pairs = 0.0, scale = 0.0;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                const DataType h = host[static_cast<size_t>(i + j * m)];
                const size_t at = offset[static_cast<size_t>(p)] + static_cast<size_t>(i + j * m);
                diff = std::max(diff, static_cast<double>(std::abs(dev[at] - h)));
                if (!dev_pairs.empty()) diff_pairs = std::max(diff_pairs, static_cast<double>(std::abs(dev_pairs[at] - h)));
                scale = std::max(scale, static_cast<double>(std::abs(h)));
            }
        }
        oss << (p == 0 ? " " : ", ") << names[p] << " " << std::scientific << std::setprecision(2)
            << (scale > 0.0 ? diff / scale : diff);
        if (!dev_pairs.empty()) oss << " (triangle pairs " << (scale > 0.0 ? diff_pairs / scale : diff_pairs) << ")";
        oss << " (" << m << "x" << n << ")";
    }
    std::printf("%s\n", oss.str().c_str());
    std::fflush(stdout);
}

// ---------------------------------------------------------------------------
// Host sketch and ID of every box of the wave (the CPU code).
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::run_ids(
    const std::vector<int64_t>& wave, std::vector<WaveBox>& boxes, int& boundary_count) {
    std::exception_ptr failure;
    std::mutex failure_mutex;
    int boundary = 0;
    const int team = std::max(1, omp_get_max_threads());
    const int split = split_threads_for(static_cast<int64_t>(wave.size()), team);
    #pragma omp parallel default(shared) reduction(+ : boundary)
    {
        FactorizationThreadScratch<CoordType, DataType> scratch;
        scratch.split_threads = split;
        #pragma omp for schedule(dynamic)
        for (int64_t bi = 0; bi < static_cast<int64_t>(wave.size()); ++bi) {
            try {
                Box* box = level_.find_local_box(wave[static_cast<size_t>(bi)]);
                if (box == nullptr) throw std::runtime_error("LevelEliminator: wave box is not local");
                gather_id_target_streamed(tree_, box, level_, kernel_, scratch, box->on_boundary);
                if (!scratch.streamed_sketch_valid) {
                    throw std::runtime_error("LevelEliminator: streamed sketch unavailable for box " +
                                             std::to_string(box->morton_index));
                }
                boundary += box->on_boundary ? 1 : 0;
                IDResult<DataType> id = fmm::compute_id_complex(
                    scratch.sketch_storage.data(), scratch.streamed_sketch_rows, scratch.workspace_cols,
                    scratch.streamed_sketch_rows, tolerance_, 0);
                ensure_nonempty_box_id_rank(id, scratch.workspace_cols);
                box->skeleton_indices = id.skeleton_indices;
                box->redundant_indices = id.redundant_indices;
                box->interpolation_matrix = std::move(id.interpolation);
                if (h2_id_trace_enabled()) {
                    auto wave_it = level_.elimination_wave.find(box->morton_index);
                    h2_id_trace_write(
                        "L" + std::to_string(level_.level) + " m=" + std::to_string(box->morton_index) +
                        " local=1 ob=" + std::to_string(box->on_boundary) +
                        " n=" + std::to_string(box->num_points) +
                        " k=" + std::to_string(box->skeleton_indices.size()) +
                        " wave=" + std::to_string(wave_it == level_.elimination_wave.end() ? -1 : wave_it->second) +
                        " |" + scratch.id_trace);
                }
                WaveBox& wb = boxes[static_cast<size_t>(bi)];
                wb.box = box;
                wb.n = static_cast<int>(box->num_points);
                wb.k = static_cast<int>(box->skeleton_indices.size());
                wb.r = static_cast<int>(box->redundant_indices.size());
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        }
    }
    if (failure) std::rethrow_exception(failure);
    boundary_count = boundary;
}

// ---------------------------------------------------------------------------
// Device sketch of the wave's ID targets, then the host ID.
//
// Mirrors gather_id_target_streamed: the rows are the distance-2 ring (full
// or skeleton rows by the same policy), sketched with the same sparse sign
// draws (4 per row, stream seeded with morton + 1).  Every destination row
// sums its contributions in the CPU order (row, then draw), so the kernel
// part of Y matches the host sketch up to kernel rounding.  The lazy fill is
// applied in sketch space as on the host: Y -= W_E P_E per fill source E, in
// the host's source order, with P_E = X_RR_full_E temp2_E[B rows]^T and
// W_E the sketch of E's temp2 rows of the ring blocks that list E.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::sketch_wave(
    const std::vector<int64_t>& wave, std::vector<WaveBox>& boxes, int& boundary_count) {
    using clock = std::chrono::steady_clock;
    auto& stats = eliminator_stats();
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    auto t0 = clock::now();

    struct Source {
        int64_t morton = -1;
        const BoxState* state = nullptr;
        int slot_b = 0;             // first temp2 row of this box
        std::vector<RunDesc> runs;  // its stored rows, by ring block
        int64_t entries = 0;        // sketched stored rows * sk
        int offset = 0;             // first column of W / row of P in the box's fill buffers
        int* ptr = nullptr;         // device lists (built on the device)
        int* list = nullptr;
        int* row_index = nullptr;   // source row of each sketched row
        int rows = 0;
    };
    struct Plan {
        int n = 0, d = 0, sk = 0;
        double scale = 0.0;
        int64_t total_rows = 0;
        std::vector<RowBlockDesc> blocks;  // ring blocks: rows and their point slots
        std::vector<Source> sources;
        // device lists (built on the device): draws, row slots, destination lists
        int* draws = nullptr;
        int* rows = nullptr;
        int* ptr = nullptr;
        int* list = nullptr;
        int* hist = nullptr;
        char* lists_block = nullptr;
        std::string trace;
        // kind 3: the kernel rows as a block of the triangle-pair evaluator
        // (emsurf_blocks.hpp): test edges = the box's points, source edges =
        // the rows, and each side's triangles
        std::vector<int> efie_test, efie_src, efie_test_tri, efie_src_tri, efie_tt, efie_st;
        S* Y = nullptr;        // d x n
        size_t off_y = 0;
        int fill_rank = 0;          // sum of the sources' r
        S* W = nullptr;        // d x fill_rank: [W_E1 W_E2 ...]
        S* P = nullptr;        // fill_rank x n: [P_E1; P_E2; ...]
    };
    std::vector<Plan> plans(wave.size());
    const bool trace = h2_id_trace_enabled();
    // kind 3: the rows' kernel values by triangle pairs (efie_* of the plans)
    const bool efie = spec_.kind == 3;
    const std::vector<int>& mesh_ints = kernel_->gpu_spec.table_int;
    auto efie_side = [&](const std::vector<int>& edges, std::vector<int>& tri_of, std::vector<int>& tris) {
        tris.clear();
        for (int e : edges) {
            for (int a = 0; a < 2; ++a) {
                const int t = mesh_ints[6 * static_cast<size_t>(e) + 2 + static_cast<size_t>(a)];
                if (t >= 0) tris.push_back(t);
            }
        }
        std::sort(tris.begin(), tris.end());
        tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
        tri_of.resize(2 * edges.size());
        for (size_t i = 0; i < edges.size(); ++i) {
            for (int a = 0; a < 2; ++a) {
                const int t = mesh_ints[6 * static_cast<size_t>(edges[i]) + 2 + static_cast<size_t>(a)];
                tri_of[2 * i + static_cast<size_t>(a)] =
                    t < 0 ? -1 : static_cast<int>(std::lower_bound(tris.begin(), tris.end(), t) - tris.begin());
            }
        }
    };

    // ---- host plans (index lists only)
    {
        std::exception_ptr failure;
        std::mutex failure_mutex;
        int boundary = 0;
        #pragma omp parallel for schedule(dynamic) reduction(+ : boundary)
        for (int64_t bi = 0; bi < static_cast<int64_t>(wave.size()); ++bi) {
            try {
                Box* box = level_.find_local_box(wave[static_cast<size_t>(bi)]);
                if (box == nullptr) throw std::runtime_error("LevelEliminator: wave box is not local");
                boundary += box->on_boundary ? 1 : 0;
                Plan& p = plans[static_cast<size_t>(bi)];
                p.n = static_cast<int>(box->num_points);
                for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
                    if (sketch_box_eliminated(level_, box->one_hop[idx])) box->use_full_set[idx] = 0;
                }

                // ring blocks: row policy and pair fill sources (a box of
                // another rank is an assisting box: skeleton rows once
                // eliminated and its skeleton has arrived, else all points)
                struct Block { int64_t morton; RingInfo info; int row_base; int a; bool skeleton_rows; std::vector<int64_t> sources; };
                std::vector<Block> blocks(box->two_hop.size());
                std::vector<LazyFarSource> pair_sources;
                bool any_sources = false;
                int row_base = 0;
                for (size_t j = 0; j < box->two_hop.size(); ++j) {
                    const int64_t rm = box->two_hop[j];
                    const RingInfo ri = ring_info(rm);
                    collect_lazy_far_sources(level_, box, rm, level_.dimension, pair_sources);
                    Block& blk = blocks[j];
                    blk.morton = rm;
                    blk.info = ri;
                    for (const auto& src : pair_sources) blk.sources.push_back(src.morton);
                    any_sources = any_sources || !pair_sources.empty();
                    const bool both_on_boundary = box->on_boundary && ri.on_boundary;
                    if (ri.remote) {
                        blk.skeleton_rows = sketch_box_eliminated(level_, rm) && !ri.skeleton->empty();
                    } else {
                        blk.skeleton_rows = !(both_on_boundary && pair_sources.empty()) &&
                                            sketch_box_eliminated(level_, rm);
                    }
                    blk.a = blk.skeleton_rows ? static_cast<int>(ri.skeleton->size()) : static_cast<int>(ri.full);
                    if (blk.a == 0) throw std::runtime_error("LevelEliminator: empty ring block");
                    blk.row_base = row_base;
                    const BoxState& rs = state_of(rm);
                    if (blk.skeleton_rows && !rs.skeleton_identity && rs.d_skeleton == nullptr) {
                        throw std::runtime_error("LevelEliminator: ring skeleton not on the device");
                    }
                    p.blocks.push_back(RowBlockDesc{row_base, blk.a, rs.slot,
                                                    blk.skeleton_rows && !rs.skeleton_identity ? rs.d_skeleton : nullptr});
                    row_base += blk.a;
                    if (efie) {  // the block's edges, in the device's row order
                        const Box* nb = level_.find_local_box(rm);
                        const std::vector<int64_t>& pts = nb != nullptr
                            ? nb->point_indices
                            : level_.assisting_boxes[static_cast<size_t>(
                                  level_.assisting_box_points_for_kernel_evaluation.at(rm))].indices;
                        for (int i = 0; i < blk.a; ++i) {
                            p.efie_src.push_back(static_cast<int>(
                                pts[static_cast<size_t>(blk.skeleton_rows ? (*ri.skeleton)[static_cast<size_t>(i)] : i)]));
                        }
                    }
                    if (trace) {
                        char policy = ri.remote ? (blk.skeleton_rows ? 's' : 'a')
                                                : (blk.skeleton_rows ? 'S' : ((both_on_boundary && pair_sources.empty()) ? 'B' : 'F'));
                        p.trace += ' ' + std::to_string(rm) + ':' + std::to_string(blk.a) + policy;
                        for (int64_t m : blk.sources) p.trace += '+' + std::to_string(m);
                    }
                }
                // static training rows after the ring: kernel rows without fill
                // (their draws continue the box's stream, as on the host)
                if (const BoxState& bs = state_of(box->morton_index); bs.ntraining > 0) {
                    p.blocks.push_back(RowBlockDesc{row_base, bs.ntraining, training_slot0_, bs.d_training});
                    row_base += bs.ntraining;
                    if (efie) {
                        for (int64_t id : bs.training_ids) p.efie_src.push_back(static_cast<int>(id));
                    }
                }
                p.total_rows = row_base;
                if (efie) {
                    if (static_cast<int64_t>(p.efie_src.size()) != p.total_rows) {
                        throw std::runtime_error("LevelEliminator: kind 3 rows do not match the sketch rows");
                    }
                    p.efie_test.assign(box->point_indices.begin(), box->point_indices.end());
                    efie_side(p.efie_test, p.efie_test_tri, p.efie_tt);
                    efie_side(p.efie_src, p.efie_src_tri, p.efie_st);
                }

                // sketch parameters (compute_id_sparse_sketch's); the draws
                // themselves are generated on the device
                int64_t d = static_cast<int64_t>(std::ceil(1.0 * p.n));
                d = std::min<int64_t>(d, p.total_rows);
                d = std::max<int64_t>(d, p.n);
                p.d = static_cast<int>(d);
                p.sk = std::min(4, p.d);
                p.scale = 1.0 / std::sqrt(static_cast<double>(p.sk));

                // fill sources of the box, in the host order (one_hop order),
                // with the runs of ring rows each one contributes to
                if (any_sources) {
                    LazyFarEndpoint<CoordType> col_end;
                    col_end.morton = box->morton_index;
                    col_end.full_size = p.n;
                    col_end.skeleton = &box->skeleton_indices;
                    col_end.wanted = nullptr;
                    col_end.wanted_count = p.n;
                    std::vector<int64_t> positions;
                    for (int64_t cm : box->one_hop) {
                        if (!sketch_box_eliminated(level_, cm)) continue;
                        if (level_.elimination_wave.find(cm) == level_.elimination_wave.end()) continue;
                        const BoxState& ss = state_of(cm);
                        Box* src_box = ss.remote ? level_.find_generator_box(cm) : level_.find_local_box(cm);
                        if (src_box == nullptr) {
                            throw std::runtime_error("LevelEliminator: fill source " + std::to_string(cm) + " not found");
                        }
                        if (ss.r == 0 || ss.temp2t == nullptr) continue;  // no Schur contribution
                        int64_t slot_offset = 0, slot_count = 0;
                        if (!lazy_far_locate_endpoint_rows(src_box, col_end, slot_offset, slot_count, positions)) continue;
                        if (slot_count != p.n) throw std::runtime_error("LevelEliminator: fill source holds a sliced box");
                        Source src;
                        src.morton = cm;
                        src.state = &ss;
                        src.slot_b = static_cast<int>(slot_offset);
                        for (const Block& blk : blocks) {
                            if (std::find(blk.sources.begin(), blk.sources.end(), cm) == blk.sources.end()) continue;
                            LazyFarEndpoint<CoordType> row_end;
                            row_end.morton = blk.morton;
                            row_end.full_size = blk.info.full;
                            row_end.skeleton = blk.info.skeleton;
                            row_end.wanted = blk.skeleton_rows ? blk.info.skeleton : nullptr;
                            row_end.wanted_count = blk.a;
                            int64_t row_offset = 0, row_count = 0;
                            if (!lazy_far_locate_endpoint_rows(src_box, row_end, row_offset, row_count, positions)) continue;
                            // rows of a full-size slot are the wanted positions
                            // (the skeleton for skeleton rows); a skeleton-size
                            // slot is already the skeleton, in order
                            const BoxState& rs = state_of(blk.morton);
                            const bool by_skeleton = row_count == blk.info.full && blk.skeleton_rows &&
                                                     !rs.skeleton_identity;
                            src.runs.push_back(RunDesc{blk.row_base, blk.a, static_cast<int>(row_offset),
                                                       by_skeleton ? rs.d_skeleton : nullptr});
                            src.entries += static_cast<int64_t>(blk.a) * p.sk;
                        }
                        if (src.runs.empty()) continue;
                        if (src.runs.size() > static_cast<size_t>(kMaxSourceRuns)) {
                            throw std::runtime_error("LevelEliminator: too many ring blocks per fill source");
                        }
                        p.sources.push_back(std::move(src));
                    }
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
        boundary_count = boundary;
    }

    stats.sketch_plan += std::chrono::duration<double>(clock::now() - t0).count();
    const auto t_gpu = clock::now();

    // ---- device buffers and launch items
    ++wave_stamp_;
    meta_.clear();
    const size_t D = sizeof(S);
    size_t y_bytes = 0;
    for (Plan& p : plans) {
        p.off_y = y_bytes;
        y_bytes = align_up(y_bytes + static_cast<size_t>(p.d) * p.n * D);
        for (Source& src : p.sources) {
            src.offset = p.fill_rank;
            p.fill_rank += src.state->r;
        }
    }
    char* d_y = heap_.alloc(std::max<size_t>(y_bytes, 1));
    std::vector<char*> fill_blocks;  // per box [P; W] and lists, freed after the sketch
    std::vector<SketchListsItem> list_items;
    std::vector<SourceListsItem> source_list_items;
    int max_row_blocks = 0;
    for (size_t bi = 0; bi < plans.size(); ++bi) {
        Plan& p = plans[bi];
        const size_t I = sizeof(int);
        const size_t hist_ints = 32 * static_cast<size_t>(kSketchOwners);
        if (p.total_rows > kEntryRowMask || p.d > kEntryDestMask) {
            throw std::runtime_error("LevelEliminator: sketch too large for the packed list entries");
        }
        const size_t draws_n = static_cast<size_t>(p.total_rows) * p.sk;
        size_t bytes = 0;
        const size_t o_draws = bytes; bytes = align_up(bytes + draws_n * I);
        const size_t o_rows = bytes;  bytes = align_up(bytes + static_cast<size_t>(p.total_rows) * I);
        const size_t o_ptr = bytes;   bytes = align_up(bytes + (kSketchOwners + 1) * I);
        const size_t o_list = bytes;  bytes = align_up(bytes + draws_n * I);
        const size_t o_hist = bytes;  bytes = align_up(bytes + hist_ints * I);
        std::vector<size_t> o_src(p.sources.size());
        for (size_t j = 0; j < p.sources.size(); ++j) {
            o_src[j] = bytes;
            bytes = align_up(bytes + (kSketchOwners + 1) * I);
            bytes = align_up(bytes + static_cast<size_t>(p.sources[j].entries) * I);
            bytes = align_up(bytes + hist_ints * I);
            bytes = align_up(bytes + static_cast<size_t>(p.sources[j].entries / p.sk) * I);  // row_index
        }
        char* block = heap_.alloc(bytes);
        fill_blocks.push_back(block);
        p.draws = reinterpret_cast<int*>(block + o_draws);
        p.rows = reinterpret_cast<int*>(block + o_rows);
        p.ptr = reinterpret_cast<int*>(block + o_ptr);
        p.list = reinterpret_cast<int*>(block + o_list);
        p.hist = reinterpret_cast<int*>(block + o_hist);
        SketchListsItem li;
        li.seed = static_cast<uint64_t>(wave[bi] + 1);
        li.d = p.d;
        li.sk = p.sk;
        li.rows = static_cast<int>(p.total_rows);
        li.num_blocks = static_cast<int>(p.blocks.size());
        li.blocks_offset = static_cast<int64_t>(meta_.append(p.blocks));
        li.draws = p.draws;
        li.rows_out = p.rows;
        li.ptr = p.ptr;
        li.entries = p.list;
        li.hist = p.hist;
        list_items.push_back(li);
        max_row_blocks = std::max(max_row_blocks, li.num_blocks);
        for (size_t j = 0; j < p.sources.size(); ++j) {
            Source& src = p.sources[j];
            char* at = block + o_src[j];
            src.ptr = reinterpret_cast<int*>(at);
            at += align_up((kSketchOwners + 1) * I);
            src.list = reinterpret_cast<int*>(at);
            at += align_up(static_cast<size_t>(src.entries) * I);
            int* hist = reinterpret_cast<int*>(at);
            at += align_up(hist_ints * I);
            src.row_index = reinterpret_cast<int*>(at);
            src.rows = static_cast<int>(src.entries / p.sk);
            SourceListsItem si;
            si.box = static_cast<int>(bi);
            si.num_runs = static_cast<int>(src.runs.size());
            si.runs_offset = static_cast<int64_t>(meta_.append(src.runs));
            si.ptr = src.ptr;
            si.entries = src.list;
            si.hist = hist;
            si.row_index = src.row_index;
            source_list_items.push_back(si);
        }
    }
    std::vector<OrderedSketchItem> row_items;     // kernel rows -> Y
    std::vector<OrderedSketchItem> stored_items;  // stored rows of the fill sources -> W
    VBatch<S> p_batch;
    VBatch<S> fill_batch;
    int max_d = 0, max_n = 0, max_r = 0;
    {
        for (Plan& p : plans) {
            p.Y = reinterpret_cast<S*>(d_y + p.off_y);
            OrderedSketchItem item{};
            item.out = p.Y;
            item.ldo = p.d;
            item.d = p.d;
            item.ncols = p.n;
            item.rows = static_cast<int>(p.total_rows);
            item.ptr = p.ptr;
            item.entries = p.list;
            item.scale = p.scale;
            item.row_slots = p.rows;
            item.col_base = state_of(wave[static_cast<size_t>(&p - plans.data())]).slot;
            if (efie) {  // the rows are stored by the evaluator (row-major, set below)
                item.row_stride = p.n;
                item.row_index = nullptr;
            }
            row_items.push_back(item);
            max_d = std::max(max_d, p.d);
            max_n = std::max(max_n, p.n);
            if (p.fill_rank == 0) continue;
            const size_t p_bytes = align_up(static_cast<size_t>(p.fill_rank) * p.n * D);
            char* fill = heap_.alloc(p_bytes + static_cast<size_t>(p.d) * p.fill_rank * D);
            fill_blocks.push_back(fill);
            p.P = reinterpret_cast<S*>(fill);
            p.W = reinterpret_cast<S*>(fill + p_bytes);
            for (Source& src : p.sources) {
                const BoxState& ss = *src.state;
                const int r = ss.r;
                // P_E = X_RR_full temp2[B rows]^T; temp2 is resident row-major (NN)
                p_batch.entries.push_back({ss.xrr_full, ss.temp2t + static_cast<size_t>(src.slot_b) * r,
                                           p.P + src.offset, r, p.n, r, r, r, p.fill_rank});
                OrderedSketchItem w{};
                w.out = p.W + static_cast<size_t>(src.offset) * p.d;
                w.ldo = p.d;
                w.d = p.d;
                w.ncols = r;
                w.rows = src.rows;
                w.ptr = src.ptr;
                w.entries = src.list;
                w.scale = p.scale;
                w.src = ss.temp2t;
                w.row_stride = r;
                w.row_index = src.row_index;
                stored_items.push_back(w);
                max_r = std::max(max_r, r);
            }
            // Y -= [W_E1 W_E2 ...] [P_E1; P_E2; ...]: all sources of the box
            // in one product (the host applies them one by one; the sum is
            // the same up to rounding)
            fill_batch.entries.push_back({p.W, p.P, p.Y, p.d, p.n, p.fill_rank, p.d, p.fill_rank, p.d});
        }
    }
    // device ID: rank, flag, norm and pivots of every box, in one small block
    size_t id_bytes = 0;
    const size_t id_rank = id_bytes; id_bytes = align_up(id_bytes + plans.size() * sizeof(int));
    const size_t id_flag = id_bytes; id_bytes = align_up(id_bytes + plans.size() * sizeof(int));
    const size_t id_norm = id_bytes; id_bytes = align_up(id_bytes + plans.size() * sizeof(double));
    std::vector<size_t> id_jpvt(plans.size());
    int max_id_n = 0;
    for (size_t bi = 0; bi < plans.size(); ++bi) {
        id_jpvt[bi] = id_bytes;
        id_bytes = align_up(id_bytes + static_cast<size_t>(plans[bi].n) * sizeof(int));
        max_id_n = std::max(max_id_n, plans[bi].n);
    }
    char* d_id = heap_.alloc_resident(std::max<size_t>(id_bytes, 1));  // pivots: the level's skeleton lists
    level_allocs_.push_back(d_id);  // its pivots are the device skeleton lists
    std::vector<QrcpItem> id_items(plans.size());
    for (size_t bi = 0; bi < plans.size(); ++bi) {
        const Plan& p = plans[bi];
        id_items[bi] = QrcpItem{p.Y, p.d, p.n, p.d, reinterpret_cast<int*>(d_id + id_jpvt[bi]),
                                reinterpret_cast<int*>(d_id + id_rank) + bi, reinterpret_cast<double*>(d_id + id_norm) + bi,
                                reinterpret_cast<int*>(d_id + id_flag) + bi};
    }
    // kind 3: the rows of groups of boxes, one group at a time in one
    // scratch buffer (triangle-pair sums, then the rows) within a memory
    // budget; each group's rows are sketched before the next group runs
    std::vector<EfieBlockItem> efie_items;
    std::vector<size_t> efie_group{0};  // group g: boxes [efie_group[g], efie_group[g + 1])
    std::vector<int64_t> efie_max_pairs, efie_max_entries;
    char* efie_scratch = nullptr;
    if (efie) {
        const size_t budget = std::max<size_t>(size_t{256} << 20,
                                               std::min<size_t>(size_t{4} << 30, heap_.largest_free() / 4));
        std::vector<size_t> sums_off(plans.size()), rows_off(plans.size());
        size_t group_bytes = 0, scratch_bytes = 0;
        int64_t max_pairs = 0, max_entries = 0;
        for (size_t bi = 0; bi < plans.size(); ++bi) {
            const Plan& p = plans[bi];
            const size_t sums = align_up(p.efie_tt.size() * p.efie_st.size() * kEfieSums * D);
            const size_t rows = align_up(static_cast<size_t>(p.total_rows) * p.n * D);
            if (group_bytes > 0 && group_bytes + sums + rows > budget) {
                efie_group.push_back(bi);
                efie_max_pairs.push_back(max_pairs);
                efie_max_entries.push_back(max_entries);
                group_bytes = 0;
                max_pairs = max_entries = 0;
            }
            sums_off[bi] = group_bytes;
            rows_off[bi] = group_bytes + sums;
            group_bytes += sums + rows;
            scratch_bytes = std::max(scratch_bytes, group_bytes);
            max_pairs = std::max<int64_t>(max_pairs, static_cast<int64_t>(p.efie_tt.size() * p.efie_st.size()));
            max_entries = std::max<int64_t>(max_entries, p.total_rows * p.n);
        }
        efie_group.push_back(plans.size());
        efie_max_pairs.push_back(max_pairs);
        efie_max_entries.push_back(max_entries);
        efie_scratch = heap_.alloc(std::max<size_t>(scratch_bytes, 1));
        for (size_t bi = 0; bi < plans.size(); ++bi) {
            Plan& p = plans[bi];
            EfieBlockItem e;
            e.nrow = p.n;  // test edges: the box's points
            e.ncol = static_cast<int>(p.total_rows);
            e.ntr = static_cast<int>(p.efie_tt.size());
            e.ntc = static_cast<int>(p.efie_st.size());
            e.row_edges = static_cast<int64_t>(meta_.append(p.efie_test));
            e.col_edges = static_cast<int64_t>(meta_.append(p.efie_src));
            e.row_tri = static_cast<int64_t>(meta_.append(p.efie_test_tri));
            e.col_tri = static_cast<int64_t>(meta_.append(p.efie_src_tri));
            e.tr = static_cast<int64_t>(meta_.append(p.efie_tt));
            e.tc = static_cast<int64_t>(meta_.append(p.efie_st));
            e.M = reinterpret_cast<dcomplex*>(efie_scratch + sums_off[bi]);
            e.out = reinterpret_cast<dcomplex*>(efie_scratch + rows_off[bi]);
            e.rs = 1;  // stored row c (source edge), column r (test edge): out[c * n + r]
            e.cs = p.n;
            efie_items.push_back(e);
            row_items[bi].src = reinterpret_cast<const S*>(e.out);
        }
    }
    const size_t off_efie = meta_.append(efie_items);
    const size_t off_list_items = meta_.append(list_items);
    const size_t off_source_list_items = meta_.append(source_list_items);
    const size_t off_rows = meta_.append(row_items);
    const size_t off_stored = meta_.append(stored_items);
    const size_t off_id = meta_.append(id_items);
    p_batch.stage(meta_);
    fill_batch.stage(meta_);

    // ---- launches, then Y to the host
    stats.sk_meta += std::chrono::duration<double>(clock::now() - t_gpu).count();
    StreamMarks marks;
    marks.mark(stream);
    char* md = meta_.upload(meta_device_, stream);
    stats.bytes_up += static_cast<double>(meta_.size());
    launch_sketch_lists(reinterpret_cast<const SketchListsItem*>(md + off_list_items),
                        static_cast<int>(list_items.size()), max_row_blocks,
                        reinterpret_cast<const SourceListsItem*>(md + off_source_list_items),
                        static_cast<int>(source_list_items.size()), md, stream);
    marks.mark(stream);
    if (efie) {
        if constexpr (std::is_same_v<S, dcomplex>) {
            for (size_t g = 0; g + 1 < efie_group.size(); ++g) {
                const size_t b0 = efie_group[g], count = efie_group[g + 1] - b0;
                launch_efie_blocks(reinterpret_cast<const EfieBlockItem*>(md + off_efie) + b0, static_cast<int>(count),
                                   efie_max_pairs[g], efie_max_entries[g], md, spec_, stream);
                launch_ordered_sketch(reinterpret_cast<const OrderedSketchItem*>(md + off_rows) + b0,
                                      static_cast<int>(count), max_d, max_n, false, spec_, points_, stream);
            }
        }
        heap_.free(efie_scratch);  // later launches are ordered after the sketch
    } else {
        launch_ordered_sketch(reinterpret_cast<const OrderedSketchItem*>(md + off_rows),
                              static_cast<int>(row_items.size()), max_d, max_n, true, spec_, points_, stream);
    }
    marks.mark(stream);
    launch_ordered_sketch(reinterpret_cast<const OrderedSketchItem*>(md + off_stored),
                          static_cast<int>(stored_items.size()), max_d, max_r, false, spec_, points_, stream);
    marks.mark(stream);
    p_batch.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
    marks.mark(stream);
    fill_batch.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 1.0, queue);
    marks.mark(stream);
    char* d_qrcp_work = nullptr;
    if (const size_t w = qrcp_work_bytes<S>(static_cast<int>(id_items.size()), max_id_n)) d_qrcp_work = heap_.alloc(w);
    launch_qrcp(reinterpret_cast<const QrcpItem*>(md + off_id), static_cast<int>(id_items.size()), max_id_n,
                tolerance_, d_qrcp_work, stream);
    heap_.free(d_qrcp_work);  // later launches are ordered after the ID
    marks.mark(stream);
    // written by the SMs into mapped pinned memory: a copy-engine transfer
    // would queue behind the background downloads of the previous wave
    char* h_id = static_cast<char*>(pinned_pool().result.reserve(std::max<size_t>(id_bytes, 256)));
    launch_copy_bytes(h_id, d_id, id_bytes, stream);
    marks.mark(stream);
    {
        // also waits for the previous wave's elimination, queued before
        const auto tw = clock::now();
        check_cuda(cudaStreamSynchronize(stream), "device sketch");
        stats.sk_wait += std::chrono::duration<double>(clock::now() - tw).count();
    }
    stats.sk_upload += marks.seconds(0);
    stats.sk_rows += marks.seconds(1);
    stats.sk_stored += marks.seconds(2);
    stats.sk_p += marks.seconds(3);
    stats.sk_fill += marks.seconds(4);
    stats.sk_id += marks.seconds(5);
    stats.sk_download += marks.seconds(6);
    stats.bytes_down += static_cast<double>(id_bytes);
    for (char* fill : fill_blocks) heap_.free(fill);
    // d_id stays: its pivot arrays are the device skeleton lists of the boxes
    wave_sketch_ = d_y;  // T lives here until the elimination has run
    stats.sketch_gpu += std::chrono::duration<double>(clock::now() - t_gpu).count();
    stats.sketch += std::chrono::duration<double>(clock::now() - t0).count();

    // ---- skeletons from the device ID (compute_id_complex's conventions:
    // full rank keeps the natural order; rank 0 keeps column 0 with a zero T)
    t0 = clock::now();
    std::exception_ptr failure;
    std::mutex failure_mutex;
    const int* ranks = reinterpret_cast<const int*>(h_id + id_rank);
    const int* flags = reinterpret_cast<const int*>(h_id + id_flag);
    const double* norms = reinterpret_cast<const double*>(h_id + id_norm);
    #pragma omp parallel for schedule(dynamic, 8)
    for (int64_t bi = 0; bi < static_cast<int64_t>(wave.size()); ++bi) {
        try {
            const Plan& p = plans[static_cast<size_t>(bi)];
            Box* box = level_.find_local_box(wave[static_cast<size_t>(bi)]);
            if (flags[bi] == 1) throw std::runtime_error("ID input contains NaN or Inf (box " + std::to_string(box->morton_index) + ")");
            if (flags[bi] == 2) throw std::runtime_error("Triangular solve produced a non-finite interpolation matrix (box " + std::to_string(box->morton_index) + ")");
            const int n = p.n;
            const int K = ranks[bi];
            const int* jpvt = reinterpret_cast<const int*>(h_id + id_jpvt[static_cast<size_t>(bi)]);
            WaveBox& wb = boxes[static_cast<size_t>(bi)];
            box->skeleton_indices.clear();
            box->redundant_indices.clear();
            box->interpolation_matrix = MatrixStorage<DataType>{};  // arrives with the factors
            if (K == n) {
                for (int i = 0; i < n; ++i) box->skeleton_indices.push_back(i);
                box->interpolation_matrix.allocate(n, 0);
            } else if (K == 0) {
                box->skeleton_indices.push_back(0);
                for (int i = 1; i < n; ++i) box->redundant_indices.push_back(i);
            } else {
                box->skeleton_indices.assign(jpvt, jpvt + K);
                box->redundant_indices.assign(jpvt + K, jpvt + n);
            }
            wb.box = box;
            wb.n = n;
            wb.k = static_cast<int>(box->skeleton_indices.size());
            wb.r = static_cast<int>(box->redundant_indices.size());
            wb.T = p.Y + static_cast<size_t>(wb.k) * p.d;  // R12 columns, solved in place
            wb.state = &state_of(box->morton_index);
            wb.state->d_skeleton = K == n ? nullptr : reinterpret_cast<const int*>(d_id + id_jpvt[static_cast<size_t>(bi)]);
            wb.ldt = p.d;
            if (trace) {
                std::ostringstream tail;
                tail << std::setprecision(10) << " | rows=" << p.total_rows << " d=" << p.d
                     << " sketch_norm=" << norms[bi] << " box_sources=";
                for (const Source& src : p.sources) tail << src.morton << ',';
                h2_id_trace_write(
                    "L" + std::to_string(level_.level) + " m=" + std::to_string(box->morton_index) +
                    " local=1 ob=" + std::to_string(box->on_boundary) + " n=" + std::to_string(box->num_points) +
                    " k=" + std::to_string(wb.k) + " wave=-1 |" + p.trace + tail.str());
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(failure_mutex);
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) std::rethrow_exception(failure);
    // after the parallel loop: the map is not thread safe
    for (const WaveBox& wb : boxes) {
        if (wb.state->d_skeleton != nullptr) host_lists_[wb.state->d_skeleton] = &wb.state->skeleton;
    }
    stats.id += std::chrono::duration<double>(clock::now() - t0).count();
}

// ---------------------------------------------------------------------------
// One wave.  Returns the number of boundary boxes in it.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
int LevelEliminator<CoordType, DataType, KernelType>::eliminate_wave(
    const std::vector<int64_t>& wave, int wave_index) {
    using clock = std::chrono::steady_clock;
    auto& stats = eliminator_stats();
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    eliminated_any_ = true;
    // a transport must have sent the generators of the previous wave's
    // boxes with remote neighbors (the transport schedule of the loop)
    for (int64_t m : last_wave_) {
        const Box* box = level_.find_local_box(m);
        if (box == nullptr) continue;
        for (int64_t nb : box->one_hop) {
            if (!is_local(nb)) {
                throw std::runtime_error("LevelEliminator: wave " + std::to_string(wave_index) +
                                         " starts before the generators of wave " +
                                         std::to_string(last_wave_index_) + " were sent");
            }
        }
    }
    last_wave_ = wave;  // their generators go out at the next transport
    last_wave_index_ = wave_index;
    since_transport_.insert(since_transport_.end(), wave.begin(), wave.end());

    // ---- 1. sketch and ID (device, or host with H2_GPU_SKETCH=0)
    auto t0 = clock::now();
    std::vector<WaveBox> boxes(wave.size());
    int boundary_count = 0;
    if (gpu_sketch_) {
        sketch_wave(wave, boxes, boundary_count);
    } else {
        run_ids(wave, boxes, boundary_count);
        stats.id += std::chrono::duration<double>(clock::now() - t0).count();
    }
    collect_elimination_marks();
    free_copied_blocks(copier_->take_finished());  // copies that finished during the sketch
    if (exchange_copier_) free_copied_blocks(exchange_copier_->take_finished());
    stats.boxes += static_cast<int64_t>(wave.size());

    // ---- 2. plan: shapes, buffers, launch items
    t0 = clock::now();
    ++wave_stamp_;
    meta_.clear();
    std::unordered_map<int64_t, size_t> wave_index_of;  // morton -> index in boxes
    for (size_t i = 0; i < boxes.size(); ++i) {
        WaveBox& wb = boxes[i];
        wb.state = &state_of(wb.box->morton_index);
        wave_index_of.emplace(wb.box->morton_index, i);
        if (wb.state->eliminated) throw std::runtime_error("LevelEliminator: box eliminated twice");
    }

    // Current neighbor sizes, X_NR row layout and per-box buffer sizes.
    // Buffers of the wave: result (LU, pivots, temp1, X_SR, T; downloaded
    // now), source (temp2 and X_RR_full for the host, in groups of at most
    // kSourceGroup bytes; downloaded while the next wave runs), persist
    // (X_RR_full and temp2^T per box, fill sources of later device sketches;
    // kept to level end), work (per box, freed at wave end).  Only the small
    // result block is one piece: multi-GB blocks would fragment the heap.
    constexpr size_t kSourceGroup = size_t{256} << 20;
    size_t result_bytes = 0, t_bytes = 0;
    std::vector<size_t> source_group_bytes;
    size_t xnr_bytes = 0;  // original X_NR of generator boxes, device exchange (packed, not downloaded)
    const size_t D = sizeof(S);
    std::vector<int> active;  // boxes with a redundant part
    for (size_t i = 0; i < boxes.size(); ++i) {
        WaveBox& wb = boxes[i];
        Box* box = wb.box;
        if (wb.r == 0) {
            // Full rank: no elimination; owner payloads still carry one
            // (empty) row count per one-hop slot.
            box->deferred_xnn_neighbor_point_counts.assign(box->one_hop.size(), 0);
            continue;
        }
        active.push_back(static_cast<int>(i));
        wb.counts.resize(box->one_hop.size());
        wb.row0.resize(box->one_hop.size());
        int rows = 0;
        for (size_t a = 0; a < box->one_hop.size(); ++a) {
            const BoxState& ns = state_of(box->one_hop[a]);
            wb.row0[a] = rows;
            wb.counts[a] = current_size(ns);
            rows += current_size(ns);
            wb.has_remote = wb.has_remote || ns.remote;
        }
        wb.ntot = rows;
        const size_t k = static_cast<size_t>(wb.k), r = static_cast<size_t>(wb.r), nt = static_cast<size_t>(wb.ntot);
        wb.off_xrr = result_bytes;      result_bytes = align_up(result_bytes + r * r * D);
        wb.off_piv = result_bytes;      result_bytes = align_up(result_bytes + r * sizeof(int));
        wb.off_temp1 = result_bytes;    result_bytes = align_up(result_bytes + k * r * D);
        wb.off_xrs = result_bytes;      result_bytes = align_up(result_bytes + r * k * D);
        {
            const size_t bytes = align_up(nt * r * D) + align_up(r * r * D) +
                                 (wb.has_remote && !device_exchange_ ? align_up(nt * r * D) : 0);
            if (source_group_bytes.empty() || (source_group_bytes.back() > 0 &&
                                               source_group_bytes.back() + bytes > kSourceGroup)) {
                source_group_bytes.push_back(0);
            }
            wb.source_group = static_cast<int>(source_group_bytes.size()) - 1;
            size_t& group = source_group_bytes.back();
            wb.off_temp2 = group;    group = align_up(group + nt * r * D);
            wb.off_xrr_full = group; group = align_up(group + r * r * D);
            if (wb.has_remote && !device_exchange_) {
                wb.off_xnr_orig = group;
                group = align_up(group + nt * r * D);
            } else if (wb.has_remote) {
                wb.off_xnr_orig = xnr_bytes;
                xnr_bytes = align_up(xnr_bytes + nt * r * D);
            }
        }
        if (gpu_sketch_) {
            wb.off_t = result_bytes;    result_bytes = align_up(result_bytes + k * r * D);
        } else {
            t_bytes = align_up(t_bytes + k * r * D);
        }
    }
    const size_t info_offset = result_bytes;
    result_bytes = align_up(result_bytes + std::max<size_t>(active.size(), 1) * sizeof(int));

    char* d_result = heap_.alloc(std::max<size_t>(result_bytes, 1));
    std::vector<char*> d_sources(source_group_bytes.size());
    for (size_t g = 0; g < d_sources.size(); ++g) d_sources[g] = heap_.alloc(std::max<size_t>(source_group_bytes[g], 1));
    char* d_xnr = xnr_bytes > 0 ? heap_.alloc(xnr_bytes) : nullptr;
    std::vector<char*> work_blocks;
    char* d_t = heap_.alloc(std::max<size_t>(t_bytes, 1));
    int* d_info = reinterpret_cast<int*>(d_result + info_offset);

    // Host image of all interpolation matrices, laid out as d_t.
    std::vector<char> t_image(t_bytes);
    {
        size_t t_off = 0;
        for (int i : active) {
            WaveBox& wb = boxes[static_cast<size_t>(i)];
            const size_t k = static_cast<size_t>(wb.k), r = static_cast<size_t>(wb.r), nt = static_cast<size_t>(wb.ntot);
            if (!gpu_sketch_) {
                const auto& T = wb.box->interpolation_matrix;
                if (T.rows != wb.k || T.cols != wb.r || T.lda != wb.k) throw std::runtime_error("LevelEliminator: unexpected T shape");
                std::memcpy(t_image.data() + t_off, T.data.data(), k * r * D);
                wb.T = reinterpret_cast<S*>(d_t + t_off);
                wb.ldt = wb.k;
                t_off = align_up(t_off + k * r * D);
            }
            wb.xrr = reinterpret_cast<S*>(d_result + wb.off_xrr);
            wb.piv = reinterpret_cast<int*>(d_result + wb.off_piv);
            wb.temp1 = reinterpret_cast<S*>(d_result + wb.off_temp1);

            char* source = d_sources[static_cast<size_t>(wb.source_group)];
            wb.temp2 = reinterpret_cast<S*>(source + wb.off_temp2);
            wb.xrr_full = reinterpret_cast<S*>(source + wb.off_xrr_full);
            if (gpu_sketch_) {
                char* persist = heap_.alloc_resident(align_up(r * r * D) + std::max<size_t>(nt * r * D, 1));
                keep_fill_source(*wb.state, persist, wb.box->one_hop);
                wb.persist_xrr = reinterpret_cast<S*>(persist);
                wb.persist_temp2t = reinterpret_cast<S*>(persist + align_up(r * r * D));
            }
            size_t w_off = 0;
            const size_t o_tmp1 = w_off; w_off = align_up(w_off + r * r * D);
            const size_t o_tmp2 = w_off; w_off = align_up(w_off + k * r * D);
            const size_t o_xsr = w_off;  w_off = align_up(w_off + k * r * D);
            const size_t o_ans = w_off;  w_off = align_up(w_off + nt * k * D);
            const size_t o_xnr = w_off;  w_off = align_up(w_off + nt * r * D);
            const size_t o_xns = w_off;  w_off = align_up(w_off + nt * k * D);
            char* work = heap_.alloc(w_off);
            work_blocks.push_back(work);
            wb.tmp1 = reinterpret_cast<S*>(work + o_tmp1);
            wb.tmp2 = reinterpret_cast<S*>(work + o_tmp2);
            wb.xsr = reinterpret_cast<S*>(work + o_xsr);
            wb.ans = reinterpret_cast<S*>(work + o_ans);
            wb.xnr = reinterpret_cast<S*>(work + o_xnr);
            wb.xns = reinterpret_cast<S*>(work + o_xns);
            wb.s = heap_.alloc_resident<S>(std::max<size_t>(k * k * D, 1));  // becomes the Schur block
        }
    }

    // Launch items of the box region.
    std::vector<EvalItem> evals;
    std::vector<GatherItem> fills;
    std::vector<GatherItem> copy_full;   // X_RR_full <- X_RR
    std::vector<IdentityItem> identities;  // X_RR^{-1} = I U^{-1} L^{-1}, before the row swaps
    std::vector<IdentityItem> symmetrizes;  // X_RR := its symmetric part, before the LU (as the host)
    std::vector<GatherItem> copy_xnr;      // original X_NR of generator boxes, for the host
    int max_xnr_m = 0, max_xnr_n = 0;
    std::vector<SymAddItem> sym_adds;
    std::vector<ColumnSwapItem> swaps;
    std::vector<AddStoreItem> stores;
    std::vector<TransposeItem> transposes;  // temp2^T into the persist block
    int max_tr_m = 0, max_tr_n = 0;
    VBatch<S> g1, g2, g3, g4, g5, g6, g7, trsm, getrf, gsolve;
    int max_eval_m = 0, max_eval_n = 0, max_fill_m = 0, max_fill_n = 0, max_copy_m = 0, max_copy_n = 0;
    int max_r = 0, max_store_m = 0, max_store_n = 0;
    std::vector<DeviceMatrix> retired;  // blocks replaced by this wave

    auto add_eval = [&](S* out, int ld, int m, int n, IndexList rows, IndexList cols) {
        if (m <= 0 || n <= 0) return;
        evals.push_back(EvalItem{out, ld, m, n, rows, cols});
        max_eval_m = std::max(max_eval_m, m);
        max_eval_n = std::max(max_eval_n, n);
    };
    auto add_gather = [&](std::vector<GatherItem>& list, int& mm, int& mn, S* out, int ld, int m, int n,
                          const S* src, int64_t rs, int64_t cs, IndexList rows, IndexList cols) {
        if (m <= 0 || n <= 0) return;
        list.push_back(GatherItem{out, ld, m, n, src, rs, cs, rows, cols});
        mm = std::max(mm, m);
        mn = std::max(mn, n);
    };
    auto push = [](VBatch<S>& batch, const S* a, int lda, const S* b, int ldb, S* c, int ldc,
                   int m, int n, int k) {
        if (m <= 0 || n <= 0 || k <= 0) return;
        batch.entries.push_back({a, b, c, m, n, k, std::max(lda, 1), std::max(ldb, 1), std::max(ldc, 1)});
    };
    const IndexList identity{};

    for (int i : active) {
        WaveBox& wb = boxes[static_cast<size_t>(i)];
        Box* box = wb.box;
        BoxState& st = *wb.state;
        const int n = wb.n, k = wb.k, r = wb.r, nt = wb.ntot;
        std::vector<int> skel(box->skeleton_indices.begin(), box->skeleton_indices.end());
        std::vector<int> red(box->redundant_indices.begin(), box->redundant_indices.end());
        wb.skeleton_offset = static_cast<int64_t>(meta_.append(skel));
        wb.redundant_offset = static_cast<int64_t>(meta_.append(red));
        const IndexList skel_local{wb.skeleton_offset, 0}, red_local{wb.redundant_offset, 0};
        const IndexList skel_slots{wb.skeleton_offset, st.slot}, red_slots{wb.redundant_offset, st.slot};

        // A_RR, A_SR, A_SS from the Schur block or the kernel
        if (st.schur.ptr != nullptr) {
            if (st.schur.rows != n || st.schur.cols != n) throw std::runtime_error("LevelEliminator: Schur block shape");
            add_gather(fills, max_fill_m, max_fill_n, wb.xrr, r, r, r, st.schur.ptr, 1, n, red_local, red_local);
            add_gather(fills, max_fill_m, max_fill_n, wb.xsr, k, k, r, st.schur.ptr, 1, n, skel_local, red_local);
            add_gather(fills, max_fill_m, max_fill_n, wb.s, k, k, k, st.schur.ptr, 1, n, skel_local, skel_local);
            retired.push_back(st.schur);
        } else {
            add_eval(wb.xrr, r, r, r, red_slots, red_slots);
            add_eval(wb.xsr, k, k, r, skel_slots, red_slots);
            add_eval(wb.s, k, k, k, skel_slots, skel_slots);
        }

        // A_NS / A_NR of every neighbor: stored block or kernel
        for (size_t a = 0; a < box->one_hop.size(); ++a) {
            const int64_t nm = box->one_hop[a];
            BoxState& ns = state_of(nm);
            const int c = static_cast<int>(wb.counts[a]);
            if (c == 0) continue;
            S* ans_rows = wb.ans + wb.row0[a];
            S* anr_rows = wb.xnr + wb.row0[a];
            auto it = edges_.find(edge_key(box->morton_index, nm));
            if (it != edges_.end()) {
                const DeviceMatrix& e = it->second;
                int64_t rs = 1, cs = e.rows;  // view rows = neighbor, cols = this box
                if (box->morton_index < nm) {
                    if (e.rows != c || e.cols != n) throw std::runtime_error("LevelEliminator: near block shape (lo)");
                } else {
                    if (e.rows != n || e.cols != c) throw std::runtime_error("LevelEliminator: near block shape (hi)");
                    rs = e.rows;
                    cs = 1;
                }
                add_gather(fills, max_fill_m, max_fill_n, ans_rows, nt, c, k, e.ptr, rs, cs, identity, skel_local);
                add_gather(fills, max_fill_m, max_fill_n, anr_rows, nt, c, r, e.ptr, rs, cs, identity, red_local);
                retired.push_back(e);
            } else {
                const IndexList rows = current_slots(ns);
                add_eval(ans_rows, nt, c, k, rows, skel_slots);
                add_eval(anr_rows, nt, c, r, rows, red_slots);
            }
        }

        // X_RR = A_RR - A_SR^T T - T^T A_SR + T^T A_SS T,  X_SR = A_SR - A_SS T
        push(g1, wb.xsr, k, wb.T, wb.ldt, wb.tmp1, r, r, r, k);      // tmp1 = -A_SR^T T   (TN)
        sym_adds.push_back(SymAddItem{wb.xrr, r, wb.tmp1, r, r});
        push(g2, wb.s, k, wb.T, wb.ldt, wb.tmp2, k, k, r, k);        // tmp2 = A_SS T       (NN)
        push(g3, wb.T, wb.ldt, wb.tmp2, k, wb.xrr, r, r, r, k);      // X_RR += T^T tmp2    (TN)
        push(g4, wb.s, k, wb.T, wb.ldt, wb.xsr, k, k, r, k);         // X_SR -= A_SS T      (NN)
        // the LU reads both triangles; the symmetric path assumes X_RR^T = X_RR
        // (compute_and_modify factors the symmetric part too)
        symmetrizes.push_back(IdentityItem{wb.xrr, r, r});
        add_gather(copy_full, max_copy_m, max_copy_n, wb.xrr_full, r, r, r, wb.xrr, 1, r, identity, identity);
        transposes.push_back(TransposeItem{wb.xsr, k, reinterpret_cast<S*>(d_result + wb.off_xrs), r, k, r});
        max_tr_m = std::max(max_tr_m, k);
        max_tr_n = std::max(max_tr_n, r);
        if (gpu_sketch_) {
            add_gather(copy_full, max_copy_m, max_copy_n, reinterpret_cast<S*>(d_result + wb.off_t), k, k, r,
                       wb.T, 1, wb.ldt, identity, identity);
            add_gather(copy_full, max_copy_m, max_copy_n, wb.persist_xrr, r, r, r, wb.xrr, 1, r, identity, identity);
            if (nt > 0) {
                transposes.push_back(TransposeItem{wb.temp2, nt, wb.persist_temp2t, r, nt, r});
                max_tr_m = std::max(max_tr_m, nt);
                max_tr_n = std::max(max_tr_n, r);
            }
        }
        push(getrf, wb.xrr, r, nullptr, 1, nullptr, 1, r, r, 1);
        // X_NR = A_NR - A_NS T
        push(g5, wb.ans, nt, wb.T, wb.ldt, wb.xnr, nt, nt, r, k);
        if (wb.has_remote && nt > 0) {
            char* source = d_xnr != nullptr ? d_xnr : d_sources[static_cast<size_t>(wb.source_group)];
            add_gather(copy_xnr, max_xnr_m, max_xnr_n, reinterpret_cast<S*>(source + wb.off_xnr_orig), nt, nt, r,
                       wb.xnr, 1, nt, identity, identity);
        }
        // temp1 = -X_SR X_RR^{-1},  temp2 = -X_NR X_RR^{-1}: with X_RR = P L U,
        // W = U^{-1} L^{-1} from two small solves (tmp1 is free after the
        // X_RR update), then one product each, then the row swaps of P
        identities.push_back(IdentityItem{wb.tmp1, r, r});
        push(trsm, wb.xrr, r, nullptr, 1, wb.tmp1, r, r, r, 1);
        push(gsolve, wb.xsr, k, wb.tmp1, r, wb.temp1, k, k, r, r);
        push(gsolve, wb.xnr, nt, wb.tmp1, r, wb.temp2, nt, nt, r, r);
        swaps.push_back(ColumnSwapItem{wb.temp1, k, k, wb.piv, r});
        if (nt > 0) swaps.push_back(ColumnSwapItem{wb.temp2, nt, nt, wb.piv, r});
        max_r = std::max(max_r, r);
        // S = A_SS + temp1 X_SR^T,  X_NS update = temp2 X_SR^T
        push(g6, wb.temp1, k, wb.xsr, k, wb.s, k, k, k, r);
        push(g7, wb.temp2, nt, wb.xsr, k, wb.xns, nt, nt, k, r);

        // Step 5: the near block of every neighbor becomes A_NS + update
        // (current neighbor rows x this box's skeleton), stored in the
        // canonical orientation.
        for (size_t a = 0; a < box->one_hop.size(); ++a) {
            const int64_t nm = box->one_hop[a];
            const int c = static_cast<int>(wb.counts[a]);
            if (c == 0) continue;
            DeviceMatrix fresh;
            fresh.ptr = heap_.alloc_resident<S>(static_cast<size_t>(c) * k * D);
            int64_t ors = 1, ocs = c;
            if (box->morton_index < nm) {
                fresh.rows = c;  // rows = neighbor (hi)
                fresh.cols = k;
            } else {
                fresh.rows = k;  // rows = this box (hi)
                fresh.cols = c;
                ors = k;
                ocs = 1;
            }
            stores.push_back(AddStoreItem{fresh.ptr, ors, ocs, c, k, wb.ans + wb.row0[a], nt, wb.xns + wb.row0[a], nt});
            max_store_m = std::max(max_store_m, c);
            max_store_n = std::max(max_store_n, k);
            edges_[edge_key(box->morton_index, nm)] = fresh;
        }
        st.schur = DeviceMatrix{wb.s, k, k};
    }

    stats.plan_boxes += std::chrono::duration<double>(clock::now() - t0).count();
    // ---- metadata image: item arrays and MAGMA arrays
    const size_t off_evals = meta_.append(evals);
    const size_t off_fills = meta_.append(fills);
    const size_t off_copy_full = meta_.append(copy_full);
    const size_t off_identities = meta_.append(identities);
    const size_t off_symmetrizes = meta_.append(symmetrizes);
    const size_t off_copy_xnr = meta_.append(copy_xnr);
    const size_t off_sym = meta_.append(sym_adds);
    const size_t off_swaps = meta_.append(swaps);
    const size_t off_stores = meta_.append(stores);
    const size_t off_transposes = meta_.append(transposes);
    for (VBatch<S>* b : {&g1, &g2, &g3, &g4, &g5, &g6, &g7, &trsm, &getrf, &gsolve}) b->stage(meta_);
    // getrf pivot pointers: one per matrix, in getrf order
    std::vector<int*> piv_ptrs;
    for (int i : active) piv_ptrs.push_back(boxes[static_cast<size_t>(i)].piv);
    const size_t off_piv_ptrs = meta_.append(piv_ptrs);
    stats.plan += std::chrono::duration<double>(clock::now() - t0).count();

    // ---- 3. device work, in stream order: the box region, then the owner
    // pass, planned on the host while the box region runs.  Nothing waits
    // for the device here: the next wave's sketch synchronizes.
    t0 = clock::now();
    char* md = meta_.upload(meta_device_, stream);
    stats.bytes_up += static_cast<double>(meta_.size() + t_bytes);
    if (t_bytes > 0) {
        check_cuda(cudaMemcpyAsync(d_t, t_image.data(), t_bytes, cudaMemcpyHostToDevice, stream), "upload T");
    }
    auto items = [&](size_t offset) { return md + offset; };
    elim_marks_ = std::make_unique<StreamMarks>();
    StreamMarks& marks = *elim_marks_;
    marks.mark(stream);
    eval_blocks(evals, meta_, md, off_evals, max_eval_m, max_eval_n, stream);
    launch_gather(reinterpret_cast<const GatherItem*>(items(off_fills)), static_cast<int>(fills.size()),
                  max_fill_m, max_fill_n, md, stream);
    // The fills were the last reads of the blocks this wave replaces: later
    // allocations may reuse them, every later launch being ordered after the
    // reads on the stream (not before: eval_blocks allocates and writes).
    for (const DeviceMatrix& m : retired) heap_.free(m.ptr);
    marks.mark(stream);  // 0: fills
    g1.gemm(md, MagmaTrans, MagmaNoTrans, -1.0, 0.0, queue);
    launch_sym_add(reinterpret_cast<const SymAddItem*>(items(off_sym)), static_cast<int>(sym_adds.size()), max_r, stream);
    g2.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
    g3.gemm(md, MagmaTrans, MagmaNoTrans, 1.0, 1.0, queue);
    g4.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 1.0, queue);
    launch_symmetrize(reinterpret_cast<const IdentityItem*>(items(off_symmetrizes)), static_cast<int>(symmetrizes.size()),
                      max_r, stream);
    marks.mark(stream);  // 1: X_RR, X_SR
    launch_gather(reinterpret_cast<const GatherItem*>(items(off_copy_full)), static_cast<int>(copy_full.size()),
                  max_copy_m, max_copy_n, md, stream);
    getrf_vbatched<S>(getrf.max_n, getrf.size_array(md, 0), getrf.size_array(md, 1),
                           getrf.template pointer_array<S*>(md, 0), getrf.size_array(md, 3),
                           reinterpret_cast<magma_int_t**>(md + off_piv_ptrs), d_info,
                           static_cast<magma_int_t>(getrf.count()), getrf_work_, queue);
    marks.mark(stream);  // 2: copies, LU
    g5.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 1.0, queue);
    launch_gather(reinterpret_cast<const GatherItem*>(items(off_copy_xnr)), static_cast<int>(copy_xnr.size()),
                  max_xnr_m, max_xnr_n, md, stream);
    marks.mark(stream);  // 3: X_NR
    launch_identity(reinterpret_cast<const IdentityItem*>(items(off_identities)), static_cast<int>(identities.size()),
                    max_r, stream);
    trsm_vbatched<S>(MagmaRight, MagmaUpper, MagmaNoTrans, MagmaNonUnit, trsm.max_m, trsm.max_n,
                          trsm.size_array(md, 0), trsm.size_array(md, 1), 1.0,
                          trsm.template pointer_array<S*>(md, 0), trsm.size_array(md, 3),
                          trsm.template pointer_array<S*>(md, 2), trsm.size_array(md, 5),
                          static_cast<magma_int_t>(trsm.count()), queue);
    trsm_vbatched<S>(MagmaRight, MagmaLower, MagmaNoTrans, MagmaUnit, trsm.max_m, trsm.max_n,
                          trsm.size_array(md, 0), trsm.size_array(md, 1), 1.0,
                          trsm.template pointer_array<S*>(md, 0), trsm.size_array(md, 3),
                          trsm.template pointer_array<S*>(md, 2), trsm.size_array(md, 5),
                          static_cast<magma_int_t>(trsm.count()), queue);
    gsolve.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 0.0, queue);
    launch_column_swaps(reinterpret_cast<const ColumnSwapItem*>(items(off_swaps)), static_cast<int>(swaps.size()), stream);
    marks.mark(stream);  // 4: temp1, temp2
    for (const auto& e : trsm.entries) stats.solve_flops += kFlopScale * 2.0 * e.m * e.n * e.n;  // two triangular solves
    for (const auto& e : gsolve.entries) stats.solve_flops += kFlopScale * 2.0 * e.m * e.n * e.k;
    launch_transpose(reinterpret_cast<const TransposeItem*>(items(off_transposes)), static_cast<int>(transposes.size()),
                     max_tr_m, max_tr_n, stream);
    g6.gemm(md, MagmaNoTrans, MagmaTrans, 1.0, 1.0, queue);
    g7.gemm(md, MagmaNoTrans, MagmaTrans, 1.0, 0.0, queue);
    launch_add_store(reinterpret_cast<const AddStoreItem*>(items(off_stores)), static_cast<int>(stores.size()),
                     max_store_m, max_store_n, stream);
    marks.mark(stream);  // 5: transposes, Schur, step 5
    stats.launch += std::chrono::duration<double>(clock::now() - t0).count();

    const auto t_owner = clock::now();
    owner_meta_.clear();
    const uint64_t owner_stamp = ~wave_stamp_;  // skeleton lists of this image
    // ---- deferred X_NN owner pass: for every source E of the wave and every
    // pair of its neighbors, target(lo, hi) += temp2_E[hi rows] X_NR_E[lo rows]^T,
    // and Schur(C) += temp2_E[C rows] X_NR_E[C rows]^T.  Contributions to a
    // target arrive in the CPU order: sources in the owner's one_hop order.
    struct Target { DeviceMatrix* block; std::vector<typename VBatch<S>::Entry> tasks; };
    std::vector<Target> targets;
    std::unordered_map<uint64_t, size_t> target_index;  // edge key, or ~morton for a Schur block
    std::vector<int64_t> candidates;
    {
        std::unordered_set<int64_t> seen;
        for (int i : active) {
            const WaveBox& wb = boxes[static_cast<size_t>(i)];
            if (wb.ntot == 0) continue;
            for (int64_t c : wb.box->one_hop) {
                if (is_local(c) && wave_index_of.count(c) == 0 && seen.insert(c).second) candidates.push_back(c);
            }
        }
        std::sort(candidates.begin(), candidates.end());
    }
    std::vector<EvalItem> owner_evals;
    int max_owner_m = 0, max_owner_n = 0;
    auto target_for = [&](uint64_t key, DeviceMatrix& block, int rows, int cols, BoxState& row_box,
                          BoxState& col_box) -> Target& {
        auto found = target_index.find(key);
        if (found != target_index.end()) return targets[found->second];
        if (block.ptr == nullptr) {
            block.ptr = heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * D);
            block.rows = rows;
            block.cols = cols;
            owner_evals.push_back(EvalItem{block.ptr, rows, rows, cols, current_slots_in(row_box, owner_meta_, owner_stamp),
                                             current_slots_in(col_box, owner_meta_, owner_stamp)});
            max_owner_m = std::max(max_owner_m, rows);
            max_owner_n = std::max(max_owner_n, cols);
            ++stats.new_targets;
        } else if (block.rows != rows || block.cols != cols) {
            std::ostringstream oss;
            oss << "LevelEliminator: owner target is " << block.rows << " x " << block.cols
                << ", expected " << rows << " x " << cols;
            throw std::runtime_error(oss.str());
        }
        target_index.emplace(key, targets.size());
        targets.push_back(Target{&block, {}});
        return targets.back();
    };
    for (int64_t cm : candidates) {
        Box* cbox = level_.find_local_box(cm);
        BoxState& cs = state_of(cm);
        for (int64_t em : cbox->one_hop) {
            auto wit = wave_index_of.find(em);
            if (wit == wave_index_of.end()) continue;
            const WaveBox& src = boxes[wit->second];
            if (src.r == 0 || src.ntot == 0) continue;
            const auto& src_hop = src.box->one_hop;
            const size_t c_pos = static_cast<size_t>(std::find(src_hop.begin(), src_hop.end(), cm) - src_hop.begin());
            if (c_pos == src_hop.size()) throw std::runtime_error("LevelEliminator: asymmetric one_hop lists");
            const int n_c = static_cast<int>(src.counts[c_pos]);
            if (n_c == 0) continue;
            const S* x_nr_c = src.xnr + src.row0[c_pos];
            for (size_t a = 0; a < src_hop.size(); ++a) {
                const int64_t am = src_hop[a];
                const int n_a = static_cast<int>(src.counts[a]);
                if (n_a == 0) continue;
                Target* target = nullptr;
                if (am == cm) {
                    target = &target_for(~static_cast<uint64_t>(cm), cs.schur, n_c, n_c, cs, cs);
                } else if (cm < am && deferred_xnn_boxes_are_one_hop(level_.dimension, cm, am)) {
                    BoxState& as = state_of(am);
                    target = &target_for(edge_key(cm, am), edges_[edge_key(cm, am)], n_a, n_c, as, cs);
                } else if (!is_local(am) && deferred_xnn_boxes_are_one_hop(level_.dimension, cm, am)) {
                    // A pair with a lower box of another rank: this rank's
                    // copy (the other rank's comes from the generator), rows =
                    // the candidate, so the product is the transpose,
                    // X_NR_E[c rows] temp2_E[a rows]^T, with the same terms.
                    BoxState& as = state_of(am);
                    target = &target_for(edge_key(cm, am), edges_[edge_key(cm, am)], n_c, n_a, cs, as);
                    target->tasks.push_back({x_nr_c, src.temp2 + src.row0[a], target->block->ptr,
                                             n_c, n_a, src.r, src.ntot, src.ntot, n_c});
                    continue;
                } else {
                    continue;  // owned by the other endpoint, or a lazy far pair
                }
                target->tasks.push_back({src.temp2 + src.row0[a], x_nr_c, target->block->ptr,
                                         n_a, n_c, src.r, src.ntot, src.ntot, n_a});
            }
        }
    }
    // Sub-batch j holds the j-th contribution of every target, so no launch
    // writes a block twice and each block sees its contributions in order.
    std::vector<VBatch<S>> owner_batches;
    {
        size_t depth = 0;
        for (const Target& t : targets) depth = std::max(depth, t.tasks.size());
        owner_batches.resize(depth);
        for (const Target& t : targets)
            for (size_t j = 0; j < t.tasks.size(); ++j) owner_batches[j].entries.push_back(t.tasks[j]);
    }

    stats.plan_owner += std::chrono::duration<double>(clock::now() - t_owner).count();
    const size_t off_owner_evals = owner_meta_.append(owner_evals);
    for (auto& b : owner_batches) b.stage(owner_meta_);
    stats.plan += std::chrono::duration<double>(clock::now() - t_owner).count();
    t0 = clock::now();
    char* omd = owner_meta_.upload(owner_meta_device_, stream);
    stats.bytes_up += static_cast<double>(owner_meta_.size());
    eval_blocks(owner_evals, owner_meta_, omd, off_owner_evals, max_owner_m, max_owner_n, stream);
    marks.mark(stream);  // 6: owner targets
    for (const auto& b : owner_batches) {
        b.gemm(omd, MagmaNoTrans, MagmaTrans, 1.0, 1.0, queue);
        stats.owner_gemms += static_cast<int64_t>(b.count());
        for (const auto& e : b.entries) stats.owner_flops += kFlopScale * 2.0 * e.m * e.n * e.k;
        ++stats.owner_batches;
    }
    marks.mark(stream);  // 7: owner GEMMs
    if (device_exchange_) {
        pack_generators(boxes, d_sources, d_xnr, d_result, stream);
        heap_.free(d_xnr);  // read by the pack only; later launches are ordered after it
    }
    if (keep_solve_ && !keep_failed_) keep_solve_factors(boxes, d_result, d_sources, stream);
    // ---- 4. factors to the host, in the background: the result block (LU,
    // pivots, temp1, X_RS, T; LU status checked on arrival) and the source
    // groups (temp2 as X_NR, X_RR_full).
    {
        auto ready = [&] {
            cudaEvent_t ev = nullptr;
            check_cuda(cudaEventCreateWithFlags(&ev, cudaEventDisableTiming), "cudaEventCreate");
            check_cuda(cudaEventRecord(ev, stream), "cudaEventRecord");
            return ev;
        };
        using Seg = HostCopier::Segment;
        const bool with_t = gpu_sketch_;
        HostCopier::Job result;
        result.device = d_result;
        result.bytes = result_bytes;
        result.ready = ready();
        auto info_host = std::make_shared<std::vector<int>>();
        result.segments.push_back(Seg{nullptr, info_host.get(), info_offset, active.size() * sizeof(int)});
        struct Dims { Box* box; int k, r; };
        auto dims = std::make_shared<std::vector<Dims>>();
        for (int i : active) {
            const WaveBox& wb = boxes[static_cast<size_t>(i)];
            Box* box = wb.box;
            const size_t k = static_cast<size_t>(wb.k), r = static_cast<size_t>(wb.r);
            box->X_RR = MatrixStorage<DataType>{};
            box->X_SR = MatrixStorage<DataType>{};
            box->X_RS_entry = MatrixStorage<DataType>{};
            if (with_t) box->interpolation_matrix = MatrixStorage<DataType>{};
            result.segments.push_back(Seg{&box->X_RR.data, nullptr, wb.off_xrr, r * r * D});
            result.segments.push_back(Seg{nullptr, &box->X_RR_pivots, wb.off_piv, r * sizeof(int)});
            result.segments.push_back(Seg{&box->X_SR.data, nullptr, wb.off_temp1, k * r * D});
            result.segments.push_back(Seg{&box->X_RS_entry.data, nullptr, wb.off_xrs, r * k * D});
            if (with_t) result.segments.push_back(Seg{&box->interpolation_matrix.data, nullptr, wb.off_t, k * r * D});
            dims->push_back({box, wb.k, wb.r});
        }
        result.finalize = [info_host, dims, with_t] {
            for (size_t j = 0; j < dims->size(); ++j) {
                const Dims& d = (*dims)[j];
                if ((*info_host)[j] != 0) {
                    throw std::runtime_error("LevelEliminator: LU of X_RR failed with INFO = " +
                                             std::to_string((*info_host)[j]) + " for box " +
                                             std::to_string(d.box->morton_index));
                }
                auto shape = [](MatrixStorage<DataType>& m, int rows, int cols, typename MatrixStorage<DataType>::Format f) {
                    m.rows = rows;
                    m.cols = cols;
                    m.lda = rows;
                    m.format = f;
                };
                shape(d.box->X_RR, d.r, d.r, MatrixStorage<DataType>::LU_FACTORED);
                shape(d.box->X_SR, d.k, d.r, MatrixStorage<DataType>::FULL);
                shape(d.box->X_RS_entry, d.r, d.k, MatrixStorage<DataType>::FULL);
                if (with_t) shape(d.box->interpolation_matrix, d.k, d.r, MatrixStorage<DataType>::FULL);
            }
        };
        copier_->submit(std::move(result));
        for (size_t g = 0; g < d_sources.size(); ++g) {
            HostCopier::Job job;
            job.device = d_sources[g];
            job.bytes = source_group_bytes[g];
            job.ready = ready();
            struct Shape { Box* box; int ntot, r; };
            auto shapes = std::make_shared<std::vector<Shape>>();
            for (int i : active) {
                const WaveBox& wb = boxes[static_cast<size_t>(i)];
                if (wb.source_group != static_cast<int>(g)) continue;
                const size_t nt = static_cast<size_t>(wb.ntot), r = static_cast<size_t>(wb.r);
                wb.box->X_NR = MatrixStorage<DataType>{};
                wb.box->X_RR_full = MatrixStorage<DataType>{};
                if (nt > 0) job.segments.push_back(Seg{&wb.box->X_NR.data, nullptr, wb.off_temp2, nt * r * D});
                job.segments.push_back(Seg{&wb.box->X_RR_full.data, nullptr, wb.off_xrr_full, r * r * D});
                if (wb.has_remote && nt > 0 && !device_exchange_) {
                    job.segments.push_back(Seg{&wb.box->lazy_original_x_nr, nullptr, wb.off_xnr_orig, nt * r * D});
                }
                shapes->push_back({wb.box, wb.ntot, wb.r});
            }
            job.finalize = [shapes] {
                for (const Shape& sh : *shapes) {
                    if (sh.ntot > 0) {
                        sh.box->X_NR.rows = sh.ntot;
                        sh.box->X_NR.cols = sh.r;
                        sh.box->X_NR.lda = sh.ntot;
                        sh.box->X_NR.format = MatrixStorage<DataType>::FULL;
                    }
                    sh.box->X_RR_full.rows = sh.r;
                    sh.box->X_RR_full.cols = sh.r;
                    sh.box->X_RR_full.lda = sh.r;
                    sh.box->X_RR_full.format = MatrixStorage<DataType>::FULL;
                }
            };
            copier_->submit(std::move(job));
        }
    }
    for (int i : active) {
        const WaveBox& wb = boxes[static_cast<size_t>(i)];
        Box* box = wb.box;
        box->deferred_xnn_neighbor_point_counts = wb.counts;
        std::vector<DataType>().swap(box->deferred_xnn_temp2);
        level_.solve_neighbor_size[static_cast<size_t>(box->morton_index - level_.local_morton_start)] = wb.counts;
        if (gpu_sketch_) {
            BoxState& st = *wb.state;
            st.temp2t = wb.persist_temp2t;
            st.xrr_full = wb.persist_xrr;
            st.ntot = wb.ntot;
            st.r = wb.r;
        }
    }
    stats.launch += std::chrono::duration<double>(clock::now() - t0).count();
    t0 = clock::now();
    if (!gpu_sketch_) {
        free_copied_blocks(copier_->wait_all());  // the host sketch of the next wave reads them
    } else {
        free_copied_blocks(copier_->take_finished());
    }
    for (WaveBox& wb : boxes) {
        BoxState& st = *wb.state;
        st.eliminated = true;
        st.skeleton.assign(wb.box->skeleton_indices.begin(), wb.box->skeleton_indices.end());
        st.skeleton_identity = true;
        for (size_t i = 0; i < st.skeleton.size(); ++i) {
            if (st.skeleton[i] != static_cast<int>(i)) {
                st.skeleton_identity = false;
                break;
            }
        }
        if (st.skeleton_identity && static_cast<int>(st.skeleton.size()) != st.n) st.skeleton_identity = false;
    }
    for (char* work : work_blocks) heap_.free(work);
    heap_.free(d_t);
    if (early_free_sources_) {
        // the wave's sketches were the last reads of sources whose local
        // neighbors are now all eliminated: the heap may take them back
        for (const WaveBox& wb : boxes) {
            for (int64_t m : wb.box->one_hop) {
                BoxState& es = state_of(m);
                if (es.persist != nullptr && --es.readers == 0) released_sources_.push_back(&es);
            }
        }
    }
    if (wave_sketch_ != nullptr) {
        heap_.free(wave_sketch_);
        wave_sketch_ = nullptr;
    }
    stats.store += std::chrono::duration<double>(clock::now() - t0).count();
    stats.heap_peak = std::max(stats.heap_peak, heap_.peak());
    return boundary_count;
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::free_copied_blocks(std::vector<char*>&& blocks) {
    for (char* block : blocks) free_block(block);
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::collect_elimination_marks() {
    if (!elim_marks_) return;
    elim_marks_->synchronize();
    auto& stats = eliminator_stats();
    for (size_t i = 0; i < 8; ++i) {
        const double t = elim_marks_->seconds(i);
        stats.el[i] += t;
        stats.device += t;
    }
    elim_marks_.reset();
}

// ---------------------------------------------------------------------------
// Level end: Schur and near-field blocks back into the host BoxData.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::download_blocks_to_host() {
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    auto& stats = eliminator_stats();

    struct Piece { const DeviceMatrix* m; size_t offset; };
    std::vector<Piece> pieces;
    size_t total = 0;
    for (const BoxState& st : states_) {
        if (st.schur.ptr) { pieces.push_back({&st.schur, total}); total += static_cast<size_t>(st.schur.rows) * st.schur.cols; }
    }
    std::vector<std::pair<uint64_t, const DeviceMatrix*>> edge_list(edges_.size());
    {
        size_t e = 0;
        for (const auto& kv : edges_) edge_list[e++] = {kv.first, &kv.second};
        std::sort(edge_list.begin(), edge_list.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
    }
    for (const auto& kv : edge_list) {
        if (kv.second->ptr == nullptr) throw std::runtime_error("LevelEliminator: unallocated edge at level end");
        pieces.push_back({kv.second, total});
        total += static_cast<size_t>(kv.second->rows) * kv.second->cols;
    }
    if (total == 0) return;

    // pack on the device with one batched copy, then one download
    S* d_pack = heap_.alloc<S>(total * sizeof(S));
    {
        std::vector<GatherItem> copies;
        int max_m = 0, max_n = 0;
        for (const Piece& p : pieces) {
            copies.push_back(GatherItem{d_pack + p.offset, p.m->rows, p.m->rows, p.m->cols, p.m->ptr, 1, p.m->rows,
                                        IndexList{}, IndexList{}});
            max_m = std::max(max_m, p.m->rows);
            max_n = std::max(max_n, p.m->cols);
        }
        meta_.clear();
        const size_t off = meta_.append(copies);
        char* md = meta_.upload(meta_device_, stream);
        launch_gather(reinterpret_cast<const GatherItem*>(md + off), static_cast<int>(copies.size()), max_m, max_n,
                      md, stream);
    }
    DataType* h_pack = static_cast<DataType*>(pinned_pool().result.reserve(total * sizeof(S)));
    check_cuda(cudaMemcpyAsync(h_pack, d_pack, total * sizeof(S), cudaMemcpyDeviceToHost, stream), "download blocks");
    check_cuda(cudaStreamSynchronize(stream), "download blocks");
    stats.bytes_down += static_cast<double>(total * sizeof(S));

    // host copies in parallel, then the (serial) map updates
    std::vector<std::vector<DataType>> data(pieces.size());
    #pragma omp parallel for schedule(dynamic, 64)
    for (int64_t i = 0; i < static_cast<int64_t>(pieces.size()); ++i) {
        const Piece& p = pieces[static_cast<size_t>(i)];
        data[static_cast<size_t>(i)].assign(h_pack + p.offset,
                                            h_pack + p.offset + static_cast<size_t>(p.m->rows) * p.m->cols);
    }
    size_t piece = 0;
    for (size_t b = 0; b < states_.size(); ++b) {
        const BoxState& st = states_[b];
        if (!st.schur.ptr) continue;
        level_.local_boxes[b].schur_complement.set_owned(st.schur.rows, st.schur.cols, std::move(data[piece++]),
                                                         MatrixStorage<DataType>::FULL);
    }
    for (const auto& kv : edge_list) {
        std::vector<DataType>& block_data = data[piece++];
        const int64_t lo = static_cast<int64_t>(kv.first >> 32);
        const int64_t hi = static_cast<int64_t>(kv.first & 0xffffffffu);
        const DeviceMatrix& m = *kv.second;
        Box* lo_box = level_.find_local_box(lo);
        Box* hi_box = level_.find_local_box(hi);
        if (lo_box == nullptr || hi_box == nullptr) {
            // a pair with a box of another rank: this rank's copy, as the
            // local box's view (rows = the other box's points)
            if (lo_box == nullptr && hi_box == nullptr) {
                throw std::runtime_error("LevelEliminator: edge without a local box at level end");
            }
            ModifiedBlock<DataType> view;
            if (lo_box != nullptr) {
                view.neighbor_morton = hi;
                view.set_a_ns_owned(m.rows, m.cols, std::move(block_data), MatrixStorage<DataType>::FULL);
                lo_box->near_field_interaction_map[hi] = static_cast<int64_t>(lo_box->near_field_modified_interactions.size());
                lo_box->near_field_modified_interactions.push_back(std::move(view));
            } else {
                std::vector<DataType> transposed(block_data.size());
                for (int j = 0; j < m.cols; ++j)
                    for (int i = 0; i < m.rows; ++i)
                        transposed[j + static_cast<size_t>(i) * m.cols] = block_data[i + static_cast<size_t>(j) * m.rows];
                view.neighbor_morton = lo;
                view.set_a_ns_owned(m.cols, m.rows, std::move(transposed), MatrixStorage<DataType>::FULL);
                hi_box->near_field_interaction_map[lo] = static_cast<int64_t>(hi_box->near_field_modified_interactions.size());
                hi_box->near_field_modified_interactions.push_back(std::move(view));
            }
            continue;
        }
        // lo's view (rows = hi, cols = lo) is the physical block; hi's view
        // shares it transposed (share_symmetric_level_edges below).
        ModifiedBlock<DataType> lo_block;
        lo_block.neighbor_morton = hi;
        lo_block.set_a_ns_owned(m.rows, m.cols, std::move(block_data), MatrixStorage<DataType>::FULL);
        lo_box->near_field_interaction_map[hi] = static_cast<int64_t>(lo_box->near_field_modified_interactions.size());
        lo_box->near_field_modified_interactions.push_back(std::move(lo_block));
        ModifiedBlock<DataType> hi_block;
        hi_block.neighbor_morton = lo;
        hi_box->near_field_interaction_map[lo] = static_cast<int64_t>(hi_box->near_field_modified_interactions.size());
        hi_box->near_field_modified_interactions.push_back(std::move(hi_block));
    }
    share_symmetric_level_edges(level_);
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::finish() {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context::instance().activate();
    for (const BoxState& st : states_) {
        if (!st.eliminated) throw std::runtime_error("LevelEliminator::finish: a box was not eliminated");
    }
    auto& stats = eliminator_stats();
    {
        const auto t1 = clock::now();
        free_copied_blocks(copier_->wait_all());
        if (exchange_copier_) free_copied_blocks(exchange_copier_->wait_all());
        collect_elimination_marks();
        // host copies of the remote generators (restored fill sources, host transition)
        for (auto& [m, rh] : remote_host_) {
            Box* gb = level_.find_generator_box(m);
            if (gb == nullptr || rh->temp2.size() != static_cast<size_t>(rh->ntot) * rh->r ||
                rh->xrr.size() != static_cast<size_t>(rh->r) * rh->r) {
                throw std::runtime_error("LevelEliminator::finish: host copy of generator " + std::to_string(m) +
                                         " incomplete");
            }
            gb->X_NR.set_owned(rh->ntot, rh->r, std::move(rh->temp2), MatrixStorage<DataType>::FULL);
            gb->X_RR_full.set_owned(rh->r, rh->r, std::move(rh->xrr), MatrixStorage<DataType>::FULL);
        }
        remote_host_.clear();
        stats.finish_sources += std::chrono::duration<double>(clock::now() - t1).count();
        stats.finish_store = copier_->busy_seconds();
        stats.copier_wait = copier_->wait_seconds();
        stats.bytes_down += copier_->bytes();
        if (exchange_copier_) stats.bytes_down += exchange_copier_->bytes();
        heap_.set_reclaimer(nullptr);
        stats.heap_reclaims = heap_.reclaims();
        copier_.reset();
        exchange_copier_.reset();
    }
    stats.heap_peak = std::max(stats.heap_peak, heap_.peak());
    stats.finish += std::chrono::duration<double>(clock::now() - t0).count();
}

// ---------------------------------------------------------------------------
// Multi-rank levels.  Boxes of other ranks within two hops (the level's
// assisting boxes) get states and point slots; their skeleton, and so their
// current rows, arrive with the assisting data of a transport.  Remote
// sources arrive as generators (lazy Schur mode 2): their updates of this
// rank's blocks are formed here on the device (the host's
// form_near_updates_from_generators + apply_updates_with_kernel_symmetric),
// and they stay resident as fill sources of later sketches.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::set_remote_skeletons(
    const std::vector<std::pair<BoxState*, const std::vector<int64_t>*>>& skeletons) {
    std::vector<int> all;
    std::vector<std::pair<BoxState*, size_t>> listed;
    for (const auto& [st, skeleton] : skeletons) {
        if (st->eliminated) continue;
        st->eliminated = true;
        st->skeleton.assign(skeleton->begin(), skeleton->end());
        st->skeleton_identity = static_cast<int>(st->skeleton.size()) == st->n;
        for (size_t i = 0; st->skeleton_identity && i < st->skeleton.size(); ++i) {
            st->skeleton_identity = st->skeleton[i] == static_cast<int>(i);
        }
        if (st->skeleton_identity) continue;
        listed.emplace_back(st, all.size());
        all.insert(all.end(), st->skeleton.begin(), st->skeleton.end());
    }
    if (all.empty()) return;
    int* d_lists = heap_.alloc_resident<int>(all.size() * sizeof(int));
    level_allocs_.push_back(reinterpret_cast<char*>(d_lists));
    check_cuda(cudaMemcpyAsync(d_lists, all.data(), all.size() * sizeof(int), cudaMemcpyHostToDevice,
                               Context::instance().stream()),
               "upload remote skeletons");
    for (const auto& [st, start] : listed) {
        st->d_skeleton = d_lists + start;
        host_lists_[st->d_skeleton] = &st->skeleton;
    }
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::sync_remote_boxes() {
    cudaStream_t stream = Context::instance().stream();
    auto& stats = eliminator_stats();
    std::vector<std::pair<int64_t, int64_t>> fresh;  // Morton, assisting index
    for (const auto& kv : level_.assisting_box_points_for_kernel_evaluation) {
        if (remote_index_.count(kv.first) == 0) fresh.emplace_back(kv.first, kv.second);
    }
    std::sort(fresh.begin(), fresh.end());
    if (!fresh.empty()) {
        int64_t extra = 0;
        for (const auto& [morton, index] : fresh) {
            const auto& assist = level_.assisting_boxes[static_cast<size_t>(index)];
            if (assist.indices.empty() || assist.coords.size() != 3 * assist.indices.size()) {
                throw std::runtime_error("LevelEliminator: assisting box " + std::to_string(morton) +
                                         " has no point data");
            }
            extra += static_cast<int64_t>(assist.indices.size());
        }
        if (num_slots_ + extra >= std::numeric_limits<int>::max()) {
            throw std::runtime_error("LevelEliminator: too many points per rank");
        }
        std::vector<double> xyz(static_cast<size_t>(3 * extra));
        std::vector<int64_t> ids(static_cast<size_t>(extra));
        int64_t at = 0;
        for (const auto& [morton, index] : fresh) {
            const auto& assist = level_.assisting_boxes[static_cast<size_t>(index)];
            BoxState st;
            st.remote = true;
            st.slot = static_cast<int>(num_slots_ + at);
            st.n = static_cast<int>(assist.indices.size());
            st.on_boundary = assist.on_boundary;
            for (size_t i = 0; i < assist.indices.size(); ++i) {
                for (int d = 0; d < 3; ++d) {
                    xyz[static_cast<size_t>(3 * (at + static_cast<int64_t>(i)) + d)] =
                        static_cast<double>(assist.coords[3 * i + static_cast<size_t>(d)]);
                }
                ids[static_cast<size_t>(at) + i] = assist.indices[i];
            }
            remote_index_[morton] = remote_states_.size();
            remote_states_.push_back(std::move(st));
            at += static_cast<int64_t>(assist.indices.size());
        }
        // a larger table: the local slots and earlier remote ones keep their places
        const int64_t total = num_slots_ + extra;
        double* xyz_new = heap_.alloc_resident<double>(static_cast<size_t>(3 * total) * sizeof(double));
        int64_t* ids_new = heap_.alloc_resident<int64_t>(static_cast<size_t>(total) * sizeof(int64_t));
        check_cuda(cudaMemcpyAsync(xyz_new, d_xyz_, static_cast<size_t>(3 * num_slots_) * sizeof(double),
                                   cudaMemcpyDeviceToDevice, stream), "point table");
        check_cuda(cudaMemcpyAsync(ids_new, d_ids_, static_cast<size_t>(num_slots_) * sizeof(int64_t),
                                   cudaMemcpyDeviceToDevice, stream), "point table");
        check_cuda(cudaMemcpyAsync(xyz_new + 3 * num_slots_, xyz.data(), xyz.size() * sizeof(double),
                                   cudaMemcpyHostToDevice, stream), "remote points");
        check_cuda(cudaMemcpyAsync(ids_new + num_slots_, ids.data(), ids.size() * sizeof(int64_t),
                                   cudaMemcpyHostToDevice, stream), "remote ids");
        heap_.free(d_xyz_);  // later launches are ordered after the copies
        heap_.free(d_ids_);
        d_xyz_ = xyz_new;
        d_ids_ = ids_new;
        points_.xyz = xyz_new;
        points_.ids = ids_new;
        num_slots_ = total;
        if (spec_.kind == 3) host_ids_.insert(host_ids_.end(), ids.begin(), ids.end());
        stats.bytes_up += static_cast<double>(xyz.size() * sizeof(double) + ids.size() * sizeof(int64_t));
    }
    // boxes eliminated since the last transport (their skeleton came with it)
    std::vector<std::pair<BoxState*, const std::vector<int64_t>*>> skeletons;
    for (const auto& [morton, index] : remote_index_) {
        BoxState& st = remote_states_[index];
        if (st.eliminated) continue;
        const auto& assist = level_.assisting_boxes[static_cast<size_t>(
            level_.assisting_box_points_for_kernel_evaluation.at(morton))];
        if (!assist.skel_indices.empty()) skeletons.emplace_back(&st, &assist.skel_indices);
    }
    set_remote_skeletons(skeletons);
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::receive_remote(const std::vector<int64_t>& installed) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context::instance().activate();
    sync_remote_boxes();
    apply_remote_generators(installed);
    eliminator_stats().remote += std::chrono::duration<double>(clock::now() - t0).count();
}

// For each installed source E (wave, then Morton order) and each local
// neighbor L of E (rows of E's temp2 at the counts of E's elimination):
//   Schur(L)  += temp2_L X_NR_L^T
//   edge(L,E) := A(E skeleton, L) + X_RS^T temp2_L^T        (E's rows shrink)
//   edge(L,N) += temp2_N X_NR_L^T  for N one hop from L (a local N once, from
//                                  its lower endpoint)
// with X_NR the original coupling.  Deltas of a block from several sources
// are summed in source order before they are added, as on the host.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::apply_remote_generators(
    const std::vector<int64_t>& installed, const std::unordered_map<int64_t, Incoming>* incoming) {
    if (installed.empty()) return;
    auto& stats = eliminator_stats();
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    const size_t D = sizeof(S);

    std::vector<int64_t> order(installed);
    std::sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
        const int32_t wa = level_.elimination_wave.at(a), wb = level_.elimination_wave.at(b);
        return wa != wb ? wa < wb : a < b;
    });

    struct Gen {
        int64_t morton = -1;
        Box* box = nullptr;
        BoxState* st = nullptr;
        int r = 0, k = 0, ntot = 0;
        std::vector<int> offset;  // first temp2 row of each one-hop neighbor
        std::vector<int> xnr_offset;  // first X_NR row of each local one-hop neighbor
        int xnr_ld = 0;
        const S* temp2 = nullptr;
        const S* xnr = nullptr;
        const S* xrs = nullptr;
        const char* bulk = nullptr;  // received into device memory (GenLayout)
    };
    std::vector<Gen> gens;
    std::vector<std::pair<BoxState*, const std::vector<int64_t>*>> skeletons;
    for (int64_t m : order) {
        Box* gb = level_.find_generator_box(m);
        if (gb == nullptr) throw std::runtime_error("LevelEliminator: generator " + std::to_string(m) + " missing");
        BoxState& st = state_of(m);
        if (!st.remote) throw std::runtime_error("LevelEliminator: generator of local box " + std::to_string(m));
        if (!gb->skeleton_indices.empty()) skeletons.emplace_back(&st, &gb->skeleton_indices);
        Gen g;
        g.morton = m;
        g.box = gb;
        g.st = &st;
        ++stats.remote_generators;
        const auto& counts = gb->deferred_xnn_neighbor_point_counts;
        if (incoming != nullptr) {  // dimensions from the header, data in the receive buffer
            auto it = incoming->find(m);
            if (it == incoming->end()) continue;  // full rank: nothing to add, no fill
            g.r = it->second.r;
            g.k = it->second.k;
            g.ntot = it->second.ntot;
            g.xnr_ld = it->second.nx;
            g.bulk = it->second.bulk;
            if (counts.size() != gb->one_hop.size() || g.k != static_cast<int>(gb->skeleton_indices.size())) {
                throw std::runtime_error("LevelEliminator: incomplete generator " + std::to_string(m));
            }
        } else {
            g.r = gb->X_NR.is_allocated() ? static_cast<int>(gb->X_NR.cols) : 0;
            if (g.r == 0) continue;  // full rank: nothing to add, no fill
            g.ntot = static_cast<int>(gb->X_NR.rows);
            g.k = gb->X_RS.is_allocated() ? static_cast<int>(gb->X_RS.cols) : 0;
            if (counts.size() != gb->one_hop.size() || gb->X_NR.lda != gb->X_NR.rows ||
                gb->lazy_original_x_nr.size() != static_cast<size_t>(g.ntot) * g.r ||
                !gb->X_RR_full.is_allocated() || gb->X_RR_full.rows != g.r || gb->X_RR_full.cols != g.r ||
                (g.k > 0 && (gb->X_RS.rows != g.r || gb->X_RS.lda != g.r)) ||
                g.k != static_cast<int>(gb->skeleton_indices.size())) {
                throw std::runtime_error("LevelEliminator: incomplete generator " + std::to_string(m));
            }
        }
        g.offset.resize(counts.size());
        int rows = 0;
        for (size_t i = 0; i < counts.size(); ++i) {
            g.offset[i] = rows;
            rows += static_cast<int>(counts[i]);
        }
        if (rows != g.ntot) throw std::runtime_error("LevelEliminator: generator row counts " + std::to_string(m));
        if (g.bulk != nullptr) {  // X_NR: rows of the local neighbors only
            g.xnr_offset.assign(counts.size(), -1);
            int local_rows = 0;
            for (size_t a = 0; a < counts.size(); ++a) {
                if (counts[a] == 0 || !is_local(gb->one_hop[a])) continue;
                g.xnr_offset[a] = local_rows;
                local_rows += static_cast<int>(counts[a]);
            }
            if (local_rows != g.xnr_ld) {
                throw std::runtime_error("LevelEliminator: generator " + std::to_string(m) + " X_NR rows mismatch");
            }
        } else {
            g.xnr_offset = g.offset;
            g.xnr_ld = g.ntot;
        }
        gens.push_back(std::move(g));
    }
    set_remote_skeletons(skeletons);
    if (gens.empty()) return;

    // uploads (unless received into device memory): [temp2 | original X_NR
    // | X_RS] for this pass, and [temp2^T | X_RR_full], kept as a fill
    // source of later sketches
    size_t tmp_bytes = 0;
    std::vector<size_t> tmp_off(gens.size());
    for (size_t i = 0; i < gens.size(); ++i) {
        if (gens[i].bulk != nullptr) continue;
        const size_t nr = static_cast<size_t>(gens[i].ntot) * gens[i].r * D;
        tmp_off[i] = tmp_bytes;
        tmp_bytes = align_up(align_up(tmp_bytes + nr) + nr) + align_up(static_cast<size_t>(gens[i].r) * gens[i].k * D);
    }
    char* d_tmp = tmp_bytes > 0 ? heap_.alloc(tmp_bytes) : nullptr;
    std::vector<TransposeItem> transposes;
    int max_tr_m = 0, max_tr_n = 0;
    for (size_t i = 0; i < gens.size(); ++i) {
        Gen& g = gens[i];
        Box* gb = g.box;
        const size_t nr = static_cast<size_t>(g.ntot) * g.r * D;
        const char* xrr_device = nullptr;
        if (g.bulk != nullptr) {
            const GenLayout lay(static_cast<size_t>(g.ntot), static_cast<size_t>(g.r), static_cast<size_t>(g.k),
                                static_cast<size_t>(g.xnr_ld));
            g.temp2 = reinterpret_cast<const S*>(g.bulk);
            g.xnr = reinterpret_cast<const S*>(g.bulk + lay.xnr);
            g.xrs = reinterpret_cast<const S*>(g.bulk + lay.xrs);
            xrr_device = g.bulk + lay.xrr;
        } else {
            char* base = d_tmp + tmp_off[i];
            g.temp2 = reinterpret_cast<const S*>(base);
            g.xnr = reinterpret_cast<const S*>(base + align_up(nr));
            g.xrs = reinterpret_cast<const S*>(base + align_up(align_up(nr) + nr));
            // pageable copies: the host data is staged before each call returns
            check_cuda(cudaMemcpyAsync(base, gb->X_NR.data.data(), nr, cudaMemcpyHostToDevice, stream),
                       "generator temp2");
            check_cuda(cudaMemcpyAsync(base + align_up(nr), gb->lazy_original_x_nr.data(), nr, cudaMemcpyHostToDevice,
                                       stream), "generator X_NR");
            if (g.k > 0) {
                check_cuda(cudaMemcpyAsync(base + align_up(align_up(nr) + nr), gb->X_RS.data.data(),
                                           static_cast<size_t>(g.r) * g.k * D, cudaMemcpyHostToDevice, stream),
                           "generator X_RS");
            }
            std::vector<DataType>().swap(gb->lazy_original_x_nr);  // needed only here (as on the host)
            stats.bytes_up += static_cast<double>(2 * nr + static_cast<size_t>(g.r) * (g.r + g.k) * D);
        }
        const size_t t_bytes = align_up(nr);
        char* persist = heap_.alloc_resident(t_bytes + static_cast<size_t>(g.r) * g.r * D);
        keep_fill_source(*g.st, persist, gb->one_hop);
        g.st->temp2t = reinterpret_cast<const S*>(persist);
        g.st->xrr_full = reinterpret_cast<const S*>(persist + t_bytes);
        g.st->ntot = g.ntot;
        g.st->r = g.r;
        if (xrr_device != nullptr) {
            check_cuda(cudaMemcpyAsync(persist + t_bytes, xrr_device, static_cast<size_t>(g.r) * g.r * D,
                                       cudaMemcpyDeviceToDevice, stream), "generator X_RR_full");
        } else {
            check_cuda(cudaMemcpyAsync(persist + t_bytes, gb->X_RR_full.data.data(),
                                       static_cast<size_t>(g.r) * g.r * D, cudaMemcpyHostToDevice, stream),
                       "generator X_RR_full");
        }
        transposes.push_back(TransposeItem{g.temp2, g.ntot, reinterpret_cast<S*>(persist), g.r, g.ntot, g.r});
        max_tr_m = std::max(max_tr_m, g.ntot);
        max_tr_n = std::max(max_tr_n, g.r);
    }

    ++wave_stamp_;
    meta_.clear();
    std::vector<EvalItem> evals;
    std::vector<GatherItem> gathers;
    int max_ev_m = 0, max_ev_n = 0, max_ga_m = 0, max_ga_n = 0;
    auto add_eval = [&](S* out, int ld, int m, int n, IndexList rows, IndexList cols) {
        evals.push_back(EvalItem{out, ld, m, n, rows, cols});
        max_ev_m = std::max(max_ev_m, m);
        max_ev_n = std::max(max_ev_n, n);
    };
    VBatch<S> rep_tt, rep_nn;  // replaced blocks: + X_RS^T temp2_L^T (or its transpose)
    struct Product { size_t off; const S* a; int lda; const S* b; int ldb; int m, n, k; };
    std::vector<Product> products;  // a b^T into the work block
    struct Accum { S* target; int ld, rows, cols; std::vector<size_t> parts; };
    std::vector<Accum> accums;
    std::unordered_map<uint64_t, size_t> accum_of;  // edge key, or ~Morton for a Schur block
    std::unordered_set<uint64_t> replaced;
    std::vector<char*> retired;
    size_t work_bytes = 0;

    auto ensure_schur = [&](BoxState& st) -> DeviceMatrix& {
        if (st.schur.ptr == nullptr) {
            const int c = current_size(st);
            st.schur = DeviceMatrix{heap_.alloc_resident<S>(static_cast<size_t>(c) * c * D), c, c};
            add_eval(st.schur.ptr, c, c, c, current_slots(st), current_slots(st));
        }
        return st.schur;
    };
    auto ensure_edge = [&](int64_t lo, int64_t hi) -> DeviceMatrix& {
        const uint64_t key = edge_key(lo, hi);
        auto it = edges_.find(key);
        if (it != edges_.end()) return it->second;
        BoxState& ls = state_of(lo);
        BoxState& hs = state_of(hi);
        const int rows = current_size(hs), cols = current_size(ls);
        DeviceMatrix& e = edges_[key];
        e = DeviceMatrix{heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * D), rows, cols};
        add_eval(e.ptr, rows, rows, cols, current_slots(hs), current_slots(ls));
        return e;
    };
    auto accumulate = [&](uint64_t key, const DeviceMatrix& target, int rows, int cols, const S* a, int lda,
                          const S* b, int ldb, int k) {
        if (target.rows != rows || target.cols != cols) {
            throw std::runtime_error("LevelEliminator: remote update does not match its block");
        }
        const size_t off = work_bytes;
        work_bytes = align_up(work_bytes + static_cast<size_t>(rows) * cols * D);
        products.push_back(Product{off, a, lda, b, ldb, rows, cols, k});
        auto it = accum_of.find(key);
        if (it == accum_of.end()) {
            it = accum_of.emplace(key, accums.size()).first;
            accums.push_back(Accum{target.ptr, target.rows, rows, cols, {}});
        }
        accums[it->second].parts.push_back(off);
    };

    for (const Gen& g : gens) {
        const auto& hop = g.box->one_hop;
        const auto& counts = g.box->deferred_xnn_neighbor_point_counts;
        for (size_t li = 0; li < hop.size(); ++li) {
            const int64_t lm = hop[li];
            const int c_l = static_cast<int>(counts[li]);
            if (c_l == 0 || !is_local(lm)) continue;
            BoxState& ls = state_of(lm);
            if (current_size(ls) != c_l) {
                throw std::runtime_error("LevelEliminator: generator " + std::to_string(g.morton) +
                                         " sees box " + std::to_string(lm) + " at another size");
            }
            const S* t2_l = g.temp2 + g.offset[li];
            const S* x_l = g.xnr + g.xnr_offset[li];
            accumulate(~static_cast<uint64_t>(lm), ensure_schur(ls), c_l, c_l, t2_l, g.ntot, x_l, g.xnr_ld, g.r);

            if (g.k > 0) {  // edge (L, E): E's skeleton rows
                const uint64_t key = edge_key(lm, g.morton);
                const bool e_hi = g.morton > lm;
                const int rows = e_hi ? g.k : c_l, cols = e_hi ? c_l : g.k;
                S* fresh = heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * D);
                IndexList skeleton = current_slots(*g.st);
                skeleton.base = 0;  // positions within E
                auto it = edges_.find(key);
                if (it != edges_.end()) {
                    const DeviceMatrix& old = it->second;
                    if (old.rows != (e_hi ? g.st->n : c_l) || old.cols != (e_hi ? c_l : g.st->n)) {
                        throw std::runtime_error("LevelEliminator: stored block of a remote source has another shape");
                    }
                    gathers.push_back(e_hi ? GatherItem{fresh, rows, rows, cols, old.ptr, 1, old.rows, skeleton, IndexList{}}
                                           : GatherItem{fresh, rows, rows, cols, old.ptr, 1, old.rows, IndexList{}, skeleton});
                    max_ga_m = std::max(max_ga_m, rows);
                    max_ga_n = std::max(max_ga_n, cols);
                    retired.push_back(reinterpret_cast<char*>(old.ptr));
                } else if (e_hi) {
                    add_eval(fresh, rows, rows, cols, current_slots(*g.st), current_slots(ls));
                } else {
                    add_eval(fresh, rows, rows, cols, current_slots(ls), current_slots(*g.st));
                }
                if (e_hi) {
                    rep_tt.entries.push_back({g.xrs, t2_l, fresh, g.k, c_l, g.r, g.r, g.ntot, rows});
                } else {
                    rep_nn.entries.push_back({t2_l, g.xrs, fresh, c_l, g.k, g.r, g.ntot, g.r, rows});
                }
                edges_[key] = DeviceMatrix{fresh, rows, cols};
                replaced.insert(key);
            }

            for (size_t ni = 0; ni < hop.size(); ++ni) {
                if (ni == li) continue;
                const int64_t nm = hop[ni];
                const int c_n = static_cast<int>(counts[ni]);
                if (c_n == 0 || !deferred_xnn_boxes_are_one_hop(level_.dimension, lm, nm)) continue;
                if (is_local(nm) && nm < lm) continue;  // from the lower local endpoint
                if (current_size(state_of(nm)) != c_n) {
                    throw std::runtime_error("LevelEliminator: generator " + std::to_string(g.morton) +
                                             " sees box " + std::to_string(nm) + " at another size");
                }
                const uint64_t key = edge_key(lm, nm);
                if (replaced.count(key)) throw std::runtime_error("LevelEliminator: a replaced block is also updated");
                const DeviceMatrix& target = ensure_edge(std::min(lm, nm), std::max(lm, nm));
                const S* t2_n = g.temp2 + g.offset[ni];
                if (lm < nm) {
                    accumulate(key, target, c_n, c_l, t2_n, g.ntot, x_l, g.xnr_ld, g.r);  // rows: N
                } else {
                    accumulate(key, target, c_l, c_n, x_l, g.xnr_ld, t2_n, g.ntot, g.r);  // rows: L
                }
            }
        }
    }

    char* d_work = heap_.alloc(std::max<size_t>(work_bytes, 1));
    VBatch<S> prod;
    for (const Product& p : products) {
        prod.entries.push_back({p.a, p.b, reinterpret_cast<S*>(d_work + p.off), p.m, p.n, p.k, p.lda, p.ldb, p.m});
    }
    std::vector<const S*> part_ptrs;
    std::vector<size_t> part_start(accums.size());
    for (size_t a = 0; a < accums.size(); ++a) {
        part_start[a] = part_ptrs.size();
        for (size_t off : accums[a].parts) part_ptrs.push_back(reinterpret_cast<const S*>(d_work + off));
    }
    const size_t off_parts = meta_.append(part_ptrs);
    std::vector<SumAddItem> sums;
    int max_sum_m = 0, max_sum_n = 0;
    for (size_t a = 0; a < accums.size(); ++a) {
        const Accum& acc = accums[a];
        sums.push_back(SumAddItem{acc.target, acc.ld, acc.rows, acc.cols,
                                  static_cast<int64_t>(off_parts + part_start[a] * sizeof(const S*)),
                                  static_cast<int>(acc.parts.size())});
        max_sum_m = std::max(max_sum_m, acc.rows);
        max_sum_n = std::max(max_sum_n, acc.cols);
    }
    const size_t off_transposes = meta_.append(transposes);
    const size_t off_evals = meta_.append(evals);
    const size_t off_gathers = meta_.append(gathers);
    const size_t off_sums = meta_.append(sums);
    prod.stage(meta_);
    rep_tt.stage(meta_);
    rep_nn.stage(meta_);
    char* md = meta_.upload(meta_device_, stream);
    stats.bytes_up += static_cast<double>(meta_.size());
    launch_transpose(reinterpret_cast<const TransposeItem*>(md + off_transposes), static_cast<int>(transposes.size()),
                     max_tr_m, max_tr_n, stream);
    eval_blocks(evals, meta_, md, off_evals, max_ev_m, max_ev_n, stream);
    launch_gather(reinterpret_cast<const GatherItem*>(md + off_gathers), static_cast<int>(gathers.size()), max_ga_m,
                  max_ga_n, md, stream);
    prod.gemm(md, MagmaNoTrans, MagmaTrans, 1.0, 0.0, queue);
    rep_tt.gemm(md, MagmaTrans, MagmaTrans, 1.0, 1.0, queue);
    rep_nn.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 1.0, queue);
    launch_sum_add(reinterpret_cast<const SumAddItem*>(md + off_sums), static_cast<int>(sums.size()), max_sum_m,
                   max_sum_n, md, stream);
    // later launches are ordered after these reads
    heap_.free(d_work);
    heap_.free(d_tmp);
    for (char* p : retired) heap_.free(p);
    stats.heap_peak = std::max(stats.heap_peak, heap_.peak());
}

// Generators of the last wave's boxes with a neighbor on another rank, as
// the host's emit_lazy_generators packages them (X_NR here is the original
// coupling, downloaded for these boxes; the host X_NR holds temp2).
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::emit_generators(PendingFactorUpdates<DataType>& pending) {
    if (last_wave_.empty()) return;
    if (copier_) free_copied_blocks(copier_->wait_all());  // the host copies of the wave
    for (int64_t m : last_wave_) {
        Box* box = level_.find_local_box(m);
        if (box == nullptr) continue;
        bool has_remote = false;
        for (int64_t nb : box->one_hop) has_remote = has_remote || !is_local(nb);
        if (!has_remote) continue;
        auto generator = std::make_shared<GeneratorPayload<DataType>>();
        generator->wave = static_cast<int32_t>(last_wave_index_);
        generator->one_hop = box->one_hop;
        const bool has_fill = box->X_NR.is_allocated() && box->X_NR.cols > 0;
        if (has_fill) {
            const int64_t r = box->X_NR.cols, rows = box->X_NR.rows;
            const int64_t k = static_cast<int64_t>(box->skeleton_indices.size());
            if (box->deferred_xnn_neighbor_point_counts.size() != box->one_hop.size() ||
                static_cast<int64_t>(box->lazy_original_x_nr.size()) != rows * r ||
                box->X_RR_full.rows != r || box->X_RS_entry.rows != r || box->X_RS_entry.cols != k) {
                throw std::runtime_error("LevelEliminator::emit_generators: incomplete state for source " +
                                         std::to_string(m));
            }
            generator->r = r;
            generator->total_rows = rows;
            generator->neighbor_point_counts = box->deferred_xnn_neighbor_point_counts;
            generator->temp2.assign(box->X_NR.data.begin(), box->X_NR.data.begin() + rows * r);
            generator->original_x_nr = std::move(box->lazy_original_x_nr);
            box->lazy_original_x_nr.clear();
            generator->x_rr_full.assign(box->X_RR_full.data.begin(), box->X_RR_full.data.begin() + r * r);
            generator->k = k;
            generator->skeleton_indices = box->skeleton_indices;
            generator->x_rs.assign(box->X_RS_entry.data.begin(), box->X_RS_entry.data.begin() + r * k);
        }
        pending.generators[m] = std::move(generator);
    }
    last_wave_.clear();
}

// ---------------------------------------------------------------------------
// The solve's factors of a wave (T, X_SR, X_NR, the LU of X_RR and its
// pivots, as the host keeps them) copied into a device block of their own,
// while the heap stays under device_solve_keep_fraction() of its capacity;
// past that, this rank drops the level's copies.  Launched before the
// wave's result and source blocks can be released.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::keep_solve_factors(
    const std::vector<WaveBox>& boxes, char* d_result, const std::vector<char*>& d_sources, cudaStream_t stream) {
    const size_t D = sizeof(S), I = sizeof(int);
    struct Span { size_t off; const char* src; size_t bytes; };
    std::vector<Span> spans;
    struct Made { int64_t morton; size_t T, xsr, xnr, lu, piv; int k, r, ntot; };
    std::vector<Made> made;
    size_t total = 0;
    for (const WaveBox& wb : boxes) {
        if (wb.k == 0 || wb.r == 0) continue;
        const size_t k = static_cast<size_t>(wb.k), r = static_cast<size_t>(wb.r), nt = static_cast<size_t>(wb.ntot);
        const char* source = d_sources[static_cast<size_t>(wb.source_group)];
        Made m{wb.box->morton_index, 0, 0, 0, 0, 0, wb.k, wb.r, wb.ntot};
        m.T = total;   spans.push_back({total, d_result + wb.off_t, k * r * D});        total = align_up(total + k * r * D);
        m.xsr = total; spans.push_back({total, d_result + wb.off_temp1, k * r * D});    total = align_up(total + k * r * D);
        m.xnr = total; if (nt > 0) spans.push_back({total, source + wb.off_temp2, nt * r * D});
        total = align_up(total + nt * r * D);
        m.lu = total;  spans.push_back({total, d_result + wb.off_xrr, r * r * D});      total = align_up(total + r * r * D);
        m.piv = total; spans.push_back({total, d_result + wb.off_piv, r * I});          total = align_up(total + r * I);
        made.push_back(m);
    }
    if (total == 0) return;
    const double limit = device_solve_keep_fraction() * static_cast<double>(heap_.capacity());
    if (static_cast<double>(heap_.used() + total) > limit || heap_.largest_free() < total) {
        drop_kept();
        keep_failed_ = true;
        return;
    }
    char* block = heap_.alloc_resident(total);
    kept_.blocks.push_back(block);
    kept_.bytes += static_cast<double>(total);
    // the span list is a few KB: a pageable copy (staged before the call
    // returns) into a buffer reserved at level start, so no allocation lands
    // in the wave's pipeline
    std::vector<CopySpan> items;
    int64_t max_words = 0;
    for (const Span& sp : spans) {
        items.push_back(CopySpan{sp.src, block + sp.off, static_cast<int64_t>(sp.bytes / 4)});
        max_words = std::max(max_words, static_cast<int64_t>(sp.bytes / 4));
    }
    const size_t list_bytes = items.size() * sizeof(CopySpan);
    if (list_bytes > solve_meta_device_.capacity()) {
        check_cuda(cudaStreamSynchronize(stream), "kept factors");  // the old list may still be read
        solve_meta_device_.reserve(list_bytes);
    }
    auto* d_items = static_cast<CopySpan*>(solve_meta_device_.reserve(list_bytes));
    check_cuda(cudaMemcpyAsync(d_items, items.data(), list_bytes, cudaMemcpyHostToDevice, stream), "kept factors");
    launch_copy_spans(d_items, static_cast<int>(items.size()), max_words, stream);
    for (const Made& m : made) {
        KeptSolveBox kb;
        kb.T = block + m.T;
        kb.xsr = block + m.xsr;
        kb.xnr = block + m.xnr;
        kb.lu = block + m.lu;
        kb.ipiv = reinterpret_cast<const int*>(block + m.piv);
        kb.k = m.k;
        kb.r = m.r;
        kb.ntot = m.ntot;
        kept_.boxes[m.morton] = kb;
    }
}

// ---------------------------------------------------------------------------
// Device exchange (multi-rank levels, CUDA-aware MPI).
//
// Sender: the generators of a wave's boxes with a neighbor on another rank
// go into the outbox, one segment per destination rank (GenLayout each),
// copied from the wave's source groups and result block before the copier
// releases them; host headers describe them.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::pack_generators(
    const std::vector<WaveBox>& boxes, const std::vector<char*>& d_sources, const char* d_xnr, char* d_result,
    cudaStream_t stream) {
    if (outbox_ != nullptr || !peer_out_.empty()) {
        throw std::runtime_error("LevelEliminator: generators of two waves wait for a transport");
    }
    const size_t D = sizeof(S);
    struct Run { int row, len; };  // rows of the original X_NR for the destination
    struct Piece { int peer; size_t off; const WaveBox* wb; int nx; std::vector<Run> runs; };
    std::vector<Piece> pieces;
    std::vector<int> peers;
    for (const WaveBox& wb : boxes) {
        const Box* box = wb.box;
        peers.clear();
        for (int64_t nb : box->one_hop) {
            if (is_local(nb)) continue;
            const int p = owner_of(nb);
            if (std::find(peers.begin(), peers.end(), p) == peers.end()) peers.push_back(p);
        }
        // occupancy: also the ranks holding the box as an assisting box (the
        // unstructured host transport's lazy_generator_requesters)
        auto rq = occupancy_ && !peers.empty() ? requesters_.find(box->morton_index) : requesters_.end();
        if (rq != requesters_.end()) {
            for (int p : rq->second) {
                if (std::find(peers.begin(), peers.end(), p) == peers.end()) peers.push_back(p);
            }
        }
        if (peers.empty()) continue;
        const bool fill = wb.r > 0;
        if (fill && !wb.has_remote) throw std::runtime_error("LevelEliminator: generator without its X_NR copy");
        const auto& skeleton = box->skeleton_indices;
        for (int p : peers) {
            Piece piece{p, 0, &wb, 0, {}};
            if (fill) {
                for (size_t a = 0; a < box->one_hop.size(); ++a) {
                    const int len = wb.counts[a];
                    if (len == 0 || is_local(box->one_hop[a]) || owner_of(box->one_hop[a]) != p) continue;
                    if (!piece.runs.empty() && piece.runs.back().row + piece.runs.back().len == wb.row0[a]) {
                        piece.runs.back().len += len;
                    } else {
                        piece.runs.push_back(Run{wb.row0[a], len});
                    }
                    piece.nx += len;
                }
            }
            const GenLayout lay(static_cast<size_t>(wb.ntot), static_cast<size_t>(wb.r), static_cast<size_t>(wb.k),
                                static_cast<size_t>(piece.nx));
            PeerOut& po = peer_out_[p];
            std::vector<int64_t>& h = po.header;
            h.push_back(box->morton_index);
            h.push_back(last_wave_index_);
            h.push_back(fill ? wb.r : 0);
            h.push_back(fill ? wb.k : 0);
            h.push_back(fill ? wb.ntot : 0);
            h.push_back(piece.nx);
            h.push_back(fill ? static_cast<int64_t>(po.bytes) : -1);
            h.push_back(static_cast<int64_t>(box->one_hop.size()));
            h.insert(h.end(), box->one_hop.begin(), box->one_hop.end());
            if (fill) {
                h.insert(h.end(), wb.counts.begin(), wb.counts.end());
                h.push_back(static_cast<int64_t>(skeleton.size()));
                h.insert(h.end(), skeleton.begin(), skeleton.end());
                piece.off = po.bytes;
                pieces.push_back(std::move(piece));
                po.bytes += lay.bytes;
            }
        }
    }
    size_t total = 0;
    for (auto& kv : peer_out_) {
        kv.second.offset = total;
        total = align_up(total + kv.second.bytes);
    }
    if (total == 0) return;
    outbox_ = exchange_alloc(total);
    auto copy = [&](char* to, const char* from, size_t bytes, const char* what) {
        if (bytes > 0) check_cuda(cudaMemcpyAsync(to, from, bytes, cudaMemcpyDeviceToDevice, stream), what);
    };
    for (const Piece& pc : pieces) {
        const WaveBox& wb = *pc.wb;
        const size_t r = static_cast<size_t>(wb.r), k = static_cast<size_t>(wb.k), nt = static_cast<size_t>(wb.ntot);
        const size_t nx = static_cast<size_t>(pc.nx);
        const GenLayout lay(nt, r, k, nx);
        char* dst = outbox_ + peer_out_.at(pc.peer).offset + pc.off;
        const char* source = d_sources[static_cast<size_t>(wb.source_group)];
        const char* xnr = d_xnr != nullptr ? d_xnr : source;
        copy(dst, source + wb.off_temp2, nt * r * D, "pack temp2");
        size_t row = 0;
        for (const Run& run : pc.runs) {  // the destination's rows of every column
            check_cuda(cudaMemcpy2DAsync(dst + lay.xnr + row * D, nx * D, xnr + wb.off_xnr_orig + run.row * D, nt * D,
                                         static_cast<size_t>(run.len) * D, r, cudaMemcpyDeviceToDevice, stream),
                       "pack X_NR");
            row += static_cast<size_t>(run.len);
        }
        copy(dst + lay.xrs, d_result + wb.off_xrs, r * k * D, "pack X_RS");
        copy(dst + lay.xrr, source + wb.off_xrr_full, r * r * D, "pack X_RR_full");
    }
}

// The transport of the Color loop.  Each rank sends each neighbor rank, in
// two rounds (sizes, then all messages at once), a host header -- on the
// first transport the assisting boxes it holds of that rank (whose owner
// then sends their skeletons as they are eliminated), the generators of the
// last wave, and the skeletons of the rank's assisting boxes eliminated since
// the last transport -- and the generators' bulk, device to device.  The
// receiver installs the generator boxes (metadata), refreshes the assisting
// skeletons, forms the updates of its blocks on the device straight from the
// receive buffer, and copies the generators' temp2 and X_RR_full to the host
// in the background (attached at finish()).  Point data of the assisting
// boxes does not change within a level: it comes once, with the first
// transport.
template<typename CoordType, typename DataType, typename KernelType>
std::chrono::high_resolution_clock::duration LevelEliminator<CoordType, DataType, KernelType>::exchange() {
    using hclock = std::chrono::high_resolution_clock;
    const auto t_start = hclock::now();
    auto& stats = eliminator_stats();
    auto& timers = transport_timers();
    auto seconds_since = [](hclock::time_point t) { return std::chrono::duration<double>(hclock::now() - t).count(); };
    const size_t D = sizeof(S);
    hclock::duration comm_time{};
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    // the outbox is packed, and no queued launch still uses heap memory the
    // receive buffer may take
    check_cuda(cudaStreamSynchronize(stream), "exchange: last wave");

    MPI_Comm comm = tree_->comm;
    const std::vector<int> peers = compute_one_hop_neighbor_ranks(tree_, level_, level_index_);
    const size_t np = peers.size();
    std::unordered_map<int, size_t> peer_index;
    for (size_t i = 0; i < np; ++i) peer_index.emplace(peers[i], i);
    auto index_of = [&](int rank, const char* what) {
        auto it = peer_index.find(rank);
        if (it == peer_index.end()) {
            throw std::runtime_error(std::string("LevelEliminator::exchange: ") + what + " on rank " +
                                     std::to_string(rank) + ", not a neighbor rank");
        }
        return it->second;
    };

    // headers: [requests] [generators] [skeletons], each after its length
    std::vector<std::vector<int64_t>> requests(np), skeletons(np);
    const bool first = !registered_;
    if (first) {
        if (!since_transport_.empty()) {
            throw std::runtime_error("LevelEliminator::exchange: boxes eliminated before the first transport");
        }
        for (const auto& kv : level_.assisting_box_points_for_kernel_evaluation) {
            requests[index_of(owner_of(kv.first), "assisting box")].push_back(kv.first);
        }
        for (auto& v : requests) std::sort(v.begin(), v.end());
    }
    for (int64_t m : since_transport_) {
        auto it = requesters_.find(m);
        if (it == requesters_.end()) continue;
        const Box* box = level_.find_local_box(m);
        for (int p : it->second) {
            auto& sk = skeletons[index_of(p, "requester")];
            sk.push_back(m);
            sk.push_back(static_cast<int64_t>(box->skeleton_indices.size()));
            sk.insert(sk.end(), box->skeleton_indices.begin(), box->skeleton_indices.end());
        }
    }
    for (const auto& kv : peer_out_) index_of(kv.first, "generator destination");
    std::vector<std::vector<int64_t>> head_out(np), head_in(np);
    std::vector<uint64_t> sizes_out(2 * np, 0), sizes_in(2 * np, 0);
    std::vector<size_t> bulk_off(np, 0);
    const std::vector<int64_t> none;
    for (size_t i = 0; i < np; ++i) {
        auto& h = head_out[i];
        auto append = [&](const std::vector<int64_t>& part) {
            h.push_back(static_cast<int64_t>(part.size()));
            h.insert(h.end(), part.begin(), part.end());
        };
        auto po = peer_out_.find(peers[i]);
        append(requests[i]);
        append(po != peer_out_.end() ? po->second.header : none);
        append(skeletons[i]);
        sizes_out[2 * i] = h.size();
        if (po != peer_out_.end()) {
            sizes_out[2 * i + 1] = po->second.bytes;
            bulk_off[i] = po->second.offset;
        }
    }

    auto check_mpi = [](int err, const char* what) {
        if (err != MPI_SUCCESS) throw std::runtime_error(std::string("LevelEliminator::exchange: ") + what + " failed");
    };
    constexpr int kTagSizes = 720, kTagHead = 721, kTagBulk = 722;
    std::vector<MPI_Request> reqs;
    // round 1: sizes
    static const bool trace = [] {
        const char* v = std::getenv("H2_GPU_EXCHANGE_TRACE");
        return v != nullptr && std::atoi(v) != 0;
    }();
    const int64_t fallbacks0 = stats.exchange_fallbacks;
    const bool outbox_in_arena = outbox_ == nullptr || DeviceHeap::exchange_arena().owns(outbox_);
    auto t_comm = hclock::now();
    for (size_t i = 0; i < np; ++i) {
        reqs.emplace_back();
        check_mpi(MPI_Irecv(&sizes_in[2 * i], 2, MPI_UINT64_T, peers[i], kTagSizes, comm, &reqs.back()), "size Irecv");
    }
    for (size_t i = 0; i < np; ++i) {
        reqs.emplace_back();
        check_mpi(MPI_Isend(&sizes_out[2 * i], 2, MPI_UINT64_T, peers[i], kTagSizes, comm, &reqs.back()), "size Isend");
    }
    check_mpi(MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE), "size Waitall");
    reqs.clear();
    comm_time += hclock::now() - t_comm;
    const double t_sizes = seconds_since(t_comm);
    timers.sizes += t_sizes;

    // round 2: headers (host) and generator bulk (device to device)
    size_t recv_total = 0;
    std::vector<size_t> recv_off(np, 0);
    for (size_t i = 0; i < np; ++i) {
        recv_off[i] = recv_total;
        recv_total = align_up(recv_total + sizes_in[2 * i + 1]);
        head_in[i].resize(sizes_in[2 * i]);
        timers.bytes_sent += static_cast<double>(sizes_out[2 * i] * sizeof(int64_t) + sizes_out[2 * i + 1]);
        timers.bytes_received += static_cast<double>(sizes_in[2 * i] * sizeof(int64_t) + sizes_in[2 * i + 1]);
    }
    char* d_recv = recv_total > 0 ? exchange_alloc(recv_total) : nullptr;
    t_comm = hclock::now();
    for (size_t i = 0; i < np; ++i) {
        if (!head_in[i].empty()) {
            reqs.emplace_back();
            check_mpi(MPI_Irecv(head_in[i].data(), static_cast<int>(head_in[i].size()), MPI_INT64_T, peers[i], kTagHead,
                                comm, &reqs.back()), "header Irecv");
        }
        if (sizes_in[2 * i + 1] > 0) {
            check_mpi(MPI_Irecv_large(d_recv + recv_off[i], sizes_in[2 * i + 1], MPI_BYTE, peers[i], kTagBulk, comm, reqs),
                      "bulk Irecv");
        }
    }
    for (size_t i = 0; i < np; ++i) {
        if (!head_out[i].empty()) {
            reqs.emplace_back();
            check_mpi(MPI_Isend(head_out[i].data(), static_cast<int>(head_out[i].size()), MPI_INT64_T, peers[i], kTagHead,
                                comm, &reqs.back()), "header Isend");
        }
        if (sizes_out[2 * i + 1] > 0) {
            check_mpi(MPI_Isend_large(outbox_ + bulk_off[i], sizes_out[2 * i + 1], MPI_BYTE, peers[i], kTagBulk, comm,
                                      reqs), "bulk Isend");
        }
    }
    check_mpi(MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE), "payload Waitall");
    reqs.clear();
    comm_time += hclock::now() - t_comm;
    const double t_payload = seconds_since(t_comm);
    timers.payload += t_payload;
    if (trace) {
        uint64_t out = 0, in = 0;
        for (size_t i = 0; i < np; ++i) {
            out += sizes_out[2 * i + 1];
            in += sizes_in[2 * i + 1];
        }
        const DeviceHeap& arena = DeviceHeap::exchange_arena();
        std::printf("[exchange] rank %d level %d wave %d: sizes %.3f s, payload %.3f s, out %.3f GB%s, in %.3f GB%s, "
                    "arena used %.2f of %.2f GB\n",
                    tree_->mpi_rank, level_index_, last_wave_index_, t_sizes, t_payload, out / 1e9,
                    outbox_in_arena ? "" : " (main heap)", in / 1e9,
                    stats.exchange_fallbacks > fallbacks0 ? " (main heap)" : "", arena.used() / 1e9,
                    arena.capacity() / 1e9);
        std::fflush(stdout);
    }
    free_block(outbox_);
    outbox_ = nullptr;
    peer_out_.clear();

    // headers: requesters, generator boxes, assisting skeletons
    const auto t_parse = hclock::now();
    std::vector<int64_t> installed;
    std::unordered_map<int64_t, Incoming> incoming;
    std::vector<std::pair<int64_t, std::vector<int64_t>>> fresh;  // assisting skeletons
    for (size_t i = 0; i < np; ++i) {
        const std::vector<int64_t>& h = head_in[i];
        size_t at = 0;
        auto next = [&]() -> int64_t {
            if (at >= h.size()) throw std::runtime_error("LevelEliminator::exchange: truncated header");
            return h[at++];
        };
        const int64_t num_requests = next();
        for (int64_t q = 0; q < num_requests; ++q) {
            const int64_t m = next();
            if (!is_local(m)) throw std::runtime_error("LevelEliminator::exchange: request for a box of another rank");
            requesters_[m].push_back(peers[i]);
        }
        const size_t generators_end = static_cast<size_t>(next()) + at;
        while (at < generators_end) {
            const int64_t m = next();
            const int32_t wave = static_cast<int32_t>(next());
            const int r = static_cast<int>(next()), k = static_cast<int>(next()), ntot = static_cast<int>(next());
            const int nx = static_cast<int>(next());
            const int64_t off = next();
            Box box;
            box.morton_index = m;
            box.one_hop.resize(static_cast<size_t>(next()));
            for (auto& x : box.one_hop) x = next();
            std::vector<int64_t> skeleton;
            if (r > 0) {
                box.deferred_xnn_neighbor_point_counts.resize(box.one_hop.size());
                for (auto& c : box.deferred_xnn_neighbor_point_counts) c = next();
                skeleton.resize(static_cast<size_t>(next()));
                for (auto& x : skeleton) x = next();
                box.skeleton_indices.assign(skeleton.begin(), skeleton.end());
            }
            if (level_.generator_id_to_index.count(m) != 0) continue;
            level_.generator_id_to_index[m] = static_cast<int64_t>(level_.generator_boxes.size());
            level_.generator_boxes.push_back(std::move(box));
            level_.elimination_wave[m] = wave;
            installed.push_back(m);
            if (r > 0) {
                const GenLayout lay(static_cast<size_t>(ntot), static_cast<size_t>(r), static_cast<size_t>(k),
                                    static_cast<size_t>(nx));
                if (off < 0 || static_cast<uint64_t>(off) + lay.bytes > sizes_in[2 * i + 1]) {
                    throw std::runtime_error("LevelEliminator::exchange: generator outside its bulk message");
                }
                incoming.emplace(m, Incoming{r, k, ntot, nx, d_recv + recv_off[i] + static_cast<size_t>(off)});
                fresh.emplace_back(m, std::move(skeleton));
            }
        }
        if (at != generators_end) throw std::runtime_error("LevelEliminator::exchange: malformed generator header");
        const size_t skeletons_end = static_cast<size_t>(next()) + at;
        while (at < skeletons_end) {
            const int64_t m = next();
            std::vector<int64_t> skeleton(static_cast<size_t>(next()));
            for (auto& x : skeleton) x = next();
            fresh.emplace_back(m, std::move(skeleton));
        }
        if (at != h.size()) throw std::runtime_error("LevelEliminator::exchange: malformed header");
    }
    for (auto& [m, skeleton] : fresh) {
        auto it = level_.assisting_box_points_for_kernel_evaluation.find(m);
        if (it == level_.assisting_box_points_for_kernel_evaluation.end()) continue;
        auto& assist = level_.assisting_boxes[static_cast<size_t>(it->second)];
        assist.skel_indices.assign(skeleton.begin(), skeleton.end());
    }
    timers.deserialize += seconds_since(t_parse);
    if (first) {
        // point data of the assisting boxes, once per level
        const auto t_assisting = hclock::now();
        std::vector<int64_t> need_assist;
        for (const auto& kv : level_.assisting_box_points_for_kernel_evaluation) need_assist.push_back(kv.first);
        if (occupancy_) {
            // every rank enters (an empty process region has no requests but
            // serves its neighbors'), as the unstructured host transport does
            comm_time += exchange_assisting_for_mortons_onehop(tree_, level_, level_index_, peers, need_assist, {},
                                                               true);
        } else if (!need_assist.empty()) {
            comm_time += exchange_assisting_for_mortons_onehop(tree_, level_, level_index_, peers, need_assist);
        }
        kernel_->register_level_coordinates(level_);
        timers.assisting += seconds_since(t_assisting);
    }
    registered_ = true;

    // device updates, straight from the receive buffer
    const auto t_apply = hclock::now();
    sync_remote_boxes();
    apply_remote_generators(installed, &incoming);
    if (d_recv != nullptr) {
        // host copies of temp2 and X_RR_full; the copier releases the buffer
        HostCopier::Job job;
        job.device = d_recv;
        job.bytes = recv_total;
        for (int64_t m : installed) {
            auto it = incoming.find(m);
            if (it == incoming.end()) continue;
            const Incoming& in = it->second;
            const GenLayout lay(static_cast<size_t>(in.ntot), static_cast<size_t>(in.r), static_cast<size_t>(in.k),
                                static_cast<size_t>(in.nx));
            auto rh = std::make_unique<RemoteHost>();
            rh->ntot = in.ntot;
            rh->r = in.r;
            const size_t base = static_cast<size_t>(in.bulk - d_recv);
            const size_t nr = static_cast<size_t>(in.ntot) * in.r * D;
            if (nr > 0) job.segments.push_back(HostCopier::Segment{&rh->temp2, nullptr, base, nr});
            job.segments.push_back(HostCopier::Segment{&rh->xrr, nullptr, base + lay.xrr, static_cast<size_t>(in.r) * in.r * D});
            remote_host_[m] = std::move(rh);
        }
        if (job.segments.empty()) {
            free_block(d_recv);  // later launches are ordered after its reads
        } else {
            check_cuda(cudaEventCreateWithFlags(&job.ready, cudaEventDisableTiming), "cudaEventCreate");
            check_cuda(cudaEventRecord(job.ready, stream), "cudaEventRecord");
            exchange_copier_->submit(std::move(job));
        }
    }
    since_transport_.clear();
    last_wave_.clear();
    stats.remote += seconds_since(t_apply);
    stats.exchange += std::chrono::duration<double>(t_apply - t_start).count();
    stats.device_exchange = true;
    return comm_time;
}

// ---------------------------------------------------------------------------
// Device transition: the parent level's Schur and near-field blocks, from
// this level's blocks, fill sources and kernel (build_parent_level_
// interactions of the host).  A parent block is the 8 x 8 array of child
// blocks C(i, j) (skeleton x skeleton):
//   i == j            the child's Schur block (kernel if it has none)
//   one-hop           the stored near-field block
//   two-hop, sources  K - sum_E temp2_E[i] X_RR_full_E temp2_E[j]^T, in the
//                     host's source order (lazy far fill)
//   otherwise         K
// The parent Schur block holds C(i, j) of the parent's own children; the
// near-field block of parents P < Q holds C(i, j)^T for i in P, j in Q (rows
// = Q's points).  The lazy fill of C^T is P_{E,j}^T G_{E,i} with
// G_{E,x} = temp2_E[x rows]^T (columns of the resident temp2^T) and
// P_{E,j} = X_RR_full_E G_{E,j}.  Parents run in chunks that bound the
// memory of G and P.  A tree of occupied boxes (occupancy) builds only the
// occupied parents (num_children > 0) and their occupied neighbors
// (one_hop), and skips the pairs with an empty child (no skeleton), as the
// host's build_parent_level_interactions_unstructured does.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
std::unique_ptr<DeviceLevelBlocks<DataType>> LevelEliminator<CoordType, DataType, KernelType>::build_parent(
    std::vector<Box>& parents) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    auto& stats = eliminator_stats();
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    if (!can_build_parent()) throw std::runtime_error("LevelEliminator::build_parent: needs the resident fill sources");
    // Released fill sources still on the device are read in place, unless
    // keeping them would leave less than half the heap for the parents'
    // blocks and the chunk buffers (the parents' blocks are at most the size
    // of this level's).  None is reclaimed during the transition: its plan
    // holds their pointers.
    while (heap_.used() > heap_.capacity() / 2 && reclaim_released_source()) {
    }
    released_sources_.clear();
    const int dim = level_.dimension;
    const int nc = morton::children_per_box(dim);
    const size_t np = parents.size();
    if (np * static_cast<size_t>(nc) != level_.local_boxes.size()) {
        throw std::runtime_error("LevelEliminator::build_parent: parent/child count mismatch");
    }
    const size_t D = sizeof(S);
    const int64_t parent_start = parents.front().morton_index;
    const uint32_t grid = 1u << (level_.level - 1);
    auto is_local_parent = [&](int64_t pm) { return pm >= parent_start && pm < parent_start + static_cast<int64_t>(np); };
    auto built = [&](size_t p) { return !occupancy_ || parents[p].num_children > 0; };
    // with occupancy, an empty child: a local box without points, or a box
    // of another rank without point data (its record is the zero-sized slot
    // of add_empty_parent_transition_assisting_slots)
    auto empty_child = [&](int64_t cm) {
        if (!occupancy_) return false;
        if (is_local(cm)) return state_of(cm).n == 0;
        if (remote_index_.count(cm) != 0) return false;
        auto it = level_.assisting_box_points_for_kernel_evaluation.find(cm);
        if (it != level_.assisting_box_points_for_kernel_evaluation.end() &&
            !level_.assisting_boxes[static_cast<size_t>(it->second)].indices.empty()) {
            throw std::runtime_error("LevelEliminator::build_parent: child " + std::to_string(cm) +
                                     " has point data but no device state");
        }
        return true;
    };

    // Parents: the local ones and their neighbors on other ranks (whose
    // children are assisting boxes, all eliminated by now), with the parent
    // point offsets of their children (skeletons, in child order).
    struct Parent {
        int64_t morton = -1;
        std::vector<int64_t> child;
        std::vector<int> offset;  // nc + 1
        int n = 0;
    };
    std::unordered_map<int64_t, Parent> parent_of;
    auto parent = [&](int64_t pm) -> Parent& {
        auto it = parent_of.find(pm);
        if (it != parent_of.end()) return it->second;
        Parent& pa = parent_of[pm];
        pa.morton = pm;
        pa.offset.assign(static_cast<size_t>(nc) + 1, 0);
        for (int c = 0; c < nc; ++c) {
            const int64_t cm = pm * nc + c;
            pa.child.push_back(cm);
            if (empty_child(cm)) {
                pa.offset[static_cast<size_t>(c) + 1] = pa.offset[static_cast<size_t>(c)];
                continue;
            }
            const BoxState& cs = state_of(cm);
            if (!cs.eliminated) throw std::runtime_error("LevelEliminator::build_parent: child " + std::to_string(cm) + " not eliminated");
            pa.offset[static_cast<size_t>(c) + 1] = pa.offset[static_cast<size_t>(c)] + static_cast<int>(cs.skeleton.size());
        }
        pa.n = pa.offset[static_cast<size_t>(nc)];
        return pa;
    };
    for (size_t p = 0; p < np; ++p) {
        if (!built(p)) continue;
        if (parent(parents[p].morton_index).n != parents[p].num_points) {
            throw std::runtime_error("LevelEliminator::build_parent: parent point count mismatch");
        }
    }
    // host skeleton (positions) of a child, for the fill-row lookups
    auto host_skeleton = [&](int64_t cm) -> const std::vector<int64_t>* { return ring_info(cm).skeleton; };
    auto host_full = [&](int64_t cm) -> int64_t { return ring_info(cm).full; };

    // Parent blocks: the Schur block of each local parent, and a near-field
    // block for each neighbor parent (a local pair once, from its lower
    // parent; a pair with a parent of another rank always: this rank's copy).
    auto out = std::make_unique<Blocks>();
    struct Pair { int64_t pm, qm; DeviceMatrix* block; };  // qm == pm: Schur block
    std::vector<std::vector<Pair>> pairs_of(np);
    for (size_t p = 0; p < np; ++p) {
        if (!built(p) || (occupancy_ && parents[p].num_points == 0)) continue;
        const int64_t pm = parents[p].morton_index;
        const int n = parents[p].num_points;
        DeviceMatrix& sm = out->schur[pm];
        sm = DeviceMatrix{heap_.alloc_resident<S>(static_cast<size_t>(n) * n * D), n, n};
        pairs_of[p].push_back(Pair{pm, pm, &sm});
    }
    for (size_t p = 0; p < np; ++p) {
        if (!built(p) || (occupancy_ && parents[p].num_points == 0)) continue;
        const int64_t pm = parents[p].morton_index;
        std::vector<int64_t> neighbors;
        if (occupancy_) {
            neighbors = parents[p].one_hop;  // the occupied ones
        } else {
            for (uint64_t qu : morton::neighbors_nd(dim, static_cast<uint64_t>(pm), grid)) {
                neighbors.push_back(static_cast<int64_t>(qu));
            }
        }
        for (const int64_t qm : neighbors) {
            if (is_local_parent(qm) && qm <= pm) continue;
            const Parent& q = parent(qm);
            if (occupancy_ && q.n == 0) continue;
            const int64_t lo = std::min(pm, qm), hi = std::max(pm, qm);
            const int rows = hi == qm ? q.n : parents[p].num_points, cols = hi == qm ? parents[p].num_points : q.n;
            DeviceMatrix& em = out->edges[edge_key(lo, hi)];
            em = DeviceMatrix{heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * D), rows, cols};
            pairs_of[p].push_back(Pair{pm, qm, &em});
        }
    }

    // Fill sources freed during the waves (multi-rank levels) are restored
    // once, in the chunk of their first reader, as [G_all | X_RR_full]
    // (ld r): G_all holds the columns G_{E,x} = temp2_E(rows of x's
    // skeleton, :)^T of every neighbor x.  The host gathers the rows (and
    // X_RR_full^T) column by column into pinned slices, uploaded as they
    // fill, and the device transposes them into the source's block.  Each
    // block is freed after the chunk of its last reader (the parent of the
    // source's last local neighbor).
    struct Restored {
        const S* gall = nullptr;
        const S* xrr = nullptr;
        std::unordered_map<int64_t, int> col_of;  // neighbor -> first column of its G
        char* block = nullptr;
        size_t last_parent = 0;
    };
    std::unordered_map<int64_t, Restored> restored;
    cudaEvent_t staged[2] = {nullptr, nullptr};
    int stage_turn = 0;
    auto last_reader = [&](const Box* sb) {
        size_t last = 0;
        for (int64_t m : sb->one_hop) {
            if (is_local(m)) last = std::max(last, static_cast<size_t>(m - level_.local_morton_start) / static_cast<size_t>(nc));
        }
        return last;
    };

    size_t p0 = 0;
    std::vector<LazyFarSource> sources;
    std::vector<int64_t> positions;
    int64_t fill_gemms = 0;
    // Chunks alternate between two metadata images and are not waited for:
    // the host plans chunk c + 1 while the device runs chunk c.
    std::vector<std::unique_ptr<StreamMarks>> chunk_marks;
    while (p0 < np) {
        const auto t_plan = clock::now();
        ++wave_stamp_;
        const bool second = chunk_marks.size() % 2 == 1;
        MetaBuilder& meta = second ? owner_meta_ : meta_;
        DeviceBuffer& meta_device = second ? owner_meta_device_ : meta_device_;
        meta.clear();
        std::vector<EvalItem> evals;
        std::vector<GatherItem> gathers, g_gathers;
        int max_eval_m = 0, max_eval_n = 0, max_gather_m = 0, max_gather_n = 0, max_g_m = 0, max_g_n = 0;
        // chunk budget for G, P and restored fill sources (the free space is
        // fragmented by the level's blocks)
        const size_t budget = std::min({size_t{2} << 30, (heap_.capacity() - heap_.used()) / 4, heap_.largest_free() / 2});
        size_t buffer_bytes = 0;
        // A fill source: its resident temp2^T and X_RR_full, or its restored
        // G columns (restored in an earlier chunk, or in this one: offsets in
        // this chunk's restore block).
        struct Source {
            const S* t2t = nullptr;   // resident temp2^T
            const S* gall = nullptr;  // restored G columns
            const S* xrr = nullptr;
            size_t gall_off = SIZE_MAX, xrr_off = SIZE_MAX;  // in this chunk's restore block
            const std::unordered_map<int64_t, int>* col_of = nullptr;
            int r = 0;
            Box* box = nullptr;
        };
        std::unordered_map<int64_t, Source> source_of;
        // restored in this chunk, at increasing offsets of the staging buffer
        struct Pending { int64_t morton; Box* box; int r; size_t gall_off, xrr_off, end; size_t last_parent; std::vector<int64_t> rows; };
        std::vector<Pending> pending;
        std::unordered_map<int64_t, std::unordered_map<int64_t, int>> pending_cols;
        size_t restore_bytes = 0;
        std::vector<int64_t> rows_of;
        auto source = [&](int64_t em) -> Source& {
            auto it = source_of.find(em);
            if (it != source_of.end()) return it->second;
            const BoxState& es = state_of(em);
            Source& s = source_of[em];
            s.r = es.r;
            s.box = es.remote ? level_.find_generator_box(em) : level_.find_local_box(em);
            if (s.box == nullptr) throw std::runtime_error("LevelEliminator::build_parent: fill source " + std::to_string(em) + " not found");
            if (es.temp2t != nullptr) {
                s.t2t = es.temp2t;
                s.xrr = es.xrr_full;
                return s;
            }
            auto rit = restored.find(em);
            if (rit != restored.end()) {
                s.gall = rit->second.gall;
                s.xrr = rit->second.xrr;
                s.col_of = &rit->second.col_of;
                return s;
            }
            // restore in this chunk: the skeleton rows of each neighbor with
            // point data here (the only ones a pair of this rank can read)
            Box* sb = s.box;
            if (sb->X_NR.rows != es.ntot || sb->X_NR.cols != s.r || sb->X_RR_full.rows != s.r) {
                throw std::runtime_error("LevelEliminator::build_parent: host copy of fill source " +
                                         std::to_string(em) + " missing");
            }
            Pending pd;
            pd.morton = em;
            pd.box = sb;
            pd.r = s.r;
            auto& cols = pending_cols[em];
            for (int64_t xm : sb->one_hop) {
                if (!is_local(xm) && remote_index_.count(xm) == 0) continue;
                LazyFarEndpoint<CoordType> endpoint;
                endpoint.morton = xm;
                endpoint.full_size = host_full(xm);
                endpoint.skeleton = host_skeleton(xm);
                endpoint.wanted = endpoint.skeleton;
                endpoint.wanted_count = static_cast<int64_t>(endpoint.skeleton->size());
                int64_t slot_offset = 0, slot_count = 0;
                if (endpoint.wanted_count == 0 ||
                    !lazy_far_locate_endpoint_rows(sb, endpoint, slot_offset, slot_count, rows_of)) {
                    continue;
                }
                cols[xm] = static_cast<int>(pd.rows.size());
                for (int64_t q : rows_of) pd.rows.push_back(slot_offset + q);
            }
            pd.gall_off = restore_bytes;
            pd.xrr_off = restore_bytes + static_cast<size_t>(s.r) * pd.rows.size() * D;
            restore_bytes = align_up(pd.xrr_off + static_cast<size_t>(s.r) * s.r * D);
            pd.end = restore_bytes;
            pd.last_parent = last_reader(sb);
            s.gall_off = pd.gall_off;
            s.xrr_off = pd.xrr_off;
            s.col_of = &cols;
            pending.push_back(std::move(pd));
            return s;
        };
        // G_{E,x} and P_{E,j}, deduplicated per chunk; a pointer is absolute,
        // or an offset in the chunk's buffers (G gathers, P) or in its
        // restore block
        struct GInfo { const S* ptr; int ld; size_t buffer = SIZE_MAX; size_t restore = SIZE_MAX; };
        std::unordered_map<uint64_t, GInfo> g_of;       // (source Morton << 32) | child Morton
        std::unordered_map<uint64_t, size_t> p_of;      // same key -> P buffer offset
        struct GPending { size_t offset; const S* src; int r; int k; IndexList cols; };
        std::vector<GPending> g_pending;
        struct PPending { size_t offset; uint64_t g_key; const S* xrr; size_t xrr_restore; int r; int k; };
        std::vector<PPending> p_pending;
        // out -= P_{E,j}^T G_{E,i} (rows j) or G_{E,i}^T P_{E,j} (rows i)
        struct FillTask { uint64_t g_key_i; uint64_t p_key_j; S* out; int ld; int m; int n; int r; bool rows_i; };
        std::vector<std::vector<FillTask>> fills;  // per target child pair, in source order

        // G_{E,x}: temp2_E rows of child x (its skeleton), as columns of temp2_E^T
        auto g_for = [&](int64_t em, int64_t xm) -> uint64_t {
            const uint64_t key = (static_cast<uint64_t>(em) << 32) | static_cast<uint64_t>(xm);
            if (g_of.count(key)) return key;
            Source& s = source(em);
            BoxState& xs = state_of(xm);
            if (s.t2t == nullptr) {  // restored: contiguous columns
                auto cit = s.col_of->find(xm);
                if (cit == s.col_of->end()) {
                    throw std::runtime_error("LevelEliminator::build_parent: restored fill source without rows for its endpoint");
                }
                GInfo info;
                info.ld = s.r;
                if (s.gall != nullptr) {
                    info.ptr = s.gall + static_cast<size_t>(cit->second) * s.r;
                } else {
                    info.ptr = nullptr;
                    info.restore = s.gall_off + static_cast<size_t>(cit->second) * s.r * D;
                }
                g_of.emplace(key, info);
                return key;
            }
            LazyFarEndpoint<CoordType> endpoint;
            endpoint.morton = xm;
            endpoint.full_size = host_full(xm);
            endpoint.skeleton = host_skeleton(xm);
            endpoint.wanted = endpoint.skeleton;
            endpoint.wanted_count = static_cast<int64_t>(endpoint.skeleton->size());
            int64_t slot_offset = 0, slot_count = 0;
            if (!lazy_far_locate_endpoint_rows(s.box, endpoint, slot_offset, slot_count, positions)) {
                throw std::runtime_error("LevelEliminator::build_parent: fill source without rows for its endpoint");
            }
            GInfo info;
            info.ld = s.r;
            const int k = static_cast<int>(endpoint.skeleton->size());
            const bool by_skeleton = slot_count == endpoint.full_size && !xs.skeleton_identity;
            if (!by_skeleton) {  // contiguous columns
                info.ptr = s.t2t + static_cast<size_t>(slot_offset) * s.r;
            } else {
                IndexList cols = current_slots_in(xs, meta, wave_stamp_);
                cols.base = static_cast<int>(slot_offset);
                info.ptr = nullptr;
                info.buffer = buffer_bytes;
                g_pending.push_back(GPending{buffer_bytes, s.t2t, s.r, k, cols});
                buffer_bytes = align_up(buffer_bytes + static_cast<size_t>(s.r) * k * D);
            }
            g_of.emplace(key, info);
            return key;
        };
        auto p_for = [&](int64_t em, int64_t jm) -> uint64_t {
            const uint64_t key = (static_cast<uint64_t>(em) << 32) | static_cast<uint64_t>(jm);
            if (p_of.count(key)) return key;
            const uint64_t gk = g_for(em, jm);
            const Source& s = source(em);
            const int k = static_cast<int>(state_of(jm).skeleton.size());
            p_of.emplace(key, buffer_bytes);
            p_pending.push_back(PPending{buffer_bytes, gk, s.xrr, s.xrr == nullptr ? s.xrr_off : SIZE_MAX, s.r, k});
            buffer_bytes = align_up(buffer_bytes + static_cast<size_t>(s.r) * k * D);
            return key;
        };

        // plan the chunk's parents
        size_t p1 = p0;
        while (p1 < np && (p1 == p0 || buffer_bytes + 2 * restore_bytes < budget)) {
            const size_t p = p1++;
            for (const Pair& pr : pairs_of[p]) {
                DeviceMatrix& M = *pr.block;
                const bool diag = pr.qm == pr.pm;
                const Parent& pa = parent(pr.pm);
                const Parent& qa = parent(pr.qm);
                // rows of M: the higher parent (Q for a local pair)
                const bool rows_q = !diag && pr.qm > pr.pm;
                for (int ci = 0; ci < nc; ++ci) {
                    const int64_t im = pa.child[static_cast<size_t>(ci)];
                    for (int cj = 0; cj < nc; ++cj) {
                        const int64_t jm = qa.child[static_cast<size_t>(cj)];
                        // target view: rows R, columns S, at (row0, col0) of M
                        const int64_t Rm = rows_q ? jm : im, Sm = rows_q ? im : jm;
                        if (empty_child(Rm) || empty_child(Sm)) continue;
                        const int row0 = rows_q ? qa.offset[static_cast<size_t>(cj)] : pa.offset[static_cast<size_t>(ci)];
                        const int col0 = rows_q ? pa.offset[static_cast<size_t>(ci)] : qa.offset[static_cast<size_t>(cj)];
                        BoxState& rs = state_of(Rm);
                        BoxState& ss = state_of(Sm);
                        const int kr = static_cast<int>(rs.skeleton.size()), ks = static_cast<int>(ss.skeleton.size());
                        if (kr == 0 || ks == 0) continue;  // (an empty block)
                        S* target = M.ptr + row0 + static_cast<size_t>(col0) * M.rows;
                        auto eval = [&] {
                            evals.push_back(EvalItem{target, M.rows, kr, ks, current_slots_in(rs, meta, wave_stamp_),
                                                     current_slots_in(ss, meta, wave_stamp_)});
                            max_eval_m = std::max(max_eval_m, kr);
                            max_eval_n = std::max(max_eval_n, ks);
                        };
                        auto gather = [&](const S* src, int64_t rsd, int64_t csd) {
                            gathers.push_back(GatherItem{target, M.rows, kr, ks, src, rsd, csd, IndexList{}, IndexList{}});
                            max_gather_m = std::max(max_gather_m, kr);
                            max_gather_n = std::max(max_gather_n, ks);
                        };
                        if (im == jm) {
                            if (rs.schur.ptr != nullptr) {
                                if (rs.schur.rows != kr || rs.schur.cols != kr) {
                                    throw std::runtime_error("LevelEliminator::build_parent: child Schur block shape");
                                }
                                gather(rs.schur.ptr, 1, kr);
                            } else {
                                eval();
                            }
                            continue;
                        }
                        int64_t chebyshev = 0;
                        bool reaches_two = false;
                        {
                            int64_t ci3[3], cj3[3];
                            lazy_far_decode_coords(dim, im, ci3);
                            lazy_far_decode_coords(dim, jm, cj3);
                            for (int d = 0; d < dim; ++d) {
                                const int64_t dist = std::abs(ci3[d] - cj3[d]);
                                chebyshev = std::max(chebyshev, dist);
                                reaches_two = reaches_two || dist == 2;
                            }
                        }
                        if (chebyshev <= 1) {
                            auto it = edges_.find(edge_key(im, jm));
                            if (it == edges_.end()) {
                                // a pair no update reached (e.g. of full-rank
                                // boxes, or with a box of another rank) keeps the
                                // kernel, as on the host (extract_child_interaction)
                                eval();
                                continue;
                            }
                            const DeviceMatrix& E = it->second;  // rows = hi, columns = lo
                            if (Rm > Sm) {
                                if (E.rows != kr || E.cols != ks) throw std::runtime_error("LevelEliminator::build_parent: child near block shape");
                                gather(E.ptr, 1, E.rows);
                            } else {
                                if (E.rows != ks || E.cols != kr) throw std::runtime_error("LevelEliminator::build_parent: child near block shape");
                                gather(E.ptr, E.rows, 1);
                            }
                            continue;
                        }
                        eval();
                        if (chebyshev != 2 || !reaches_two) continue;  // beyond two hops: kernel only
                        // two-hop: lazy fill of C(i, j) from i's eliminated neighbors next to j
                        collect_lazy_far_sources(level_, level_.find_local_box(im), jm, dim, sources);
                        size_t t = 0;
                        for (const LazyFarSource& src : sources) {
                            const BoxState& es = state_of(src.morton);
                            if (es.r == 0) continue;  // no Schur contribution
                            if (fills.size() <= t) fills.resize(t + 1);
                            fills[t].push_back(FillTask{g_for(src.morton, im), p_for(src.morton, jm), target, M.rows, kr, ks,
                                                        es.r, !rows_q});
                            ++t;
                        }
                    }
                }
            }
        }

        // the chunk's restored sources: [temp2 rows; X_RR_full^T] of each
        // (L = columns + r rows, r columns) gathered into pinned slices,
        // uploaded into a raw buffer and transposed into the source's block
        std::vector<char*> d_restore(pending.size());
        auto restore_ptr = [&](size_t off) -> const S* {  // staging offset -> device
            size_t lo = 0, hi = pending.size();
            while (hi - lo > 1) {
                const size_t mid = (lo + hi) / 2;
                (pending[mid].gall_off <= off ? lo : hi) = mid;
            }
            return reinterpret_cast<const S*>(d_restore[lo] + (off - pending[lo].gall_off));
        };
        char* d_raw = nullptr;
        std::vector<TransposeItem> restores;
        int max_rs_m = 0, max_rs_n = 0;
        if (!pending.empty()) {
            const auto t_restore = clock::now();
            d_raw = heap_.alloc(restore_bytes);
            for (size_t q = 0; q < pending.size(); ++q) {
                const Pending& pd = pending[q];
                d_restore[q] = heap_.alloc_resident(pd.end - pd.gall_off);
                const int L = static_cast<int>(pd.rows.size()) + pd.r;
                restores.push_back(TransposeItem{reinterpret_cast<const S*>(d_raw + pd.gall_off), L,
                                                 reinterpret_cast<S*>(d_restore[q]), pd.r, L, pd.r});
                max_rs_m = std::max(max_rs_m, L);
                max_rs_n = std::max(max_rs_n, pd.r);
                Restored& rs = restored[pd.morton];
                rs.gall = restore_ptr(pd.gall_off);
                rs.xrr = restore_ptr(pd.xrr_off);
                rs.col_of = std::move(pending_cols[pd.morton]);
                rs.block = d_restore[q];
                rs.last_parent = pd.last_parent;
            }
            constexpr size_t kSlice = size_t{128} << 20;
            for (size_t q0 = 0; q0 < pending.size();) {
                const size_t base = pending[q0].gall_off;
                size_t q1 = q0 + 1;
                while (q1 < pending.size() && pending[q1].end - base <= kSlice) ++q1;
                const size_t bytes = pending[q1 - 1].end - base;
                if (staged[stage_turn] == nullptr) {
                    check_cuda(cudaEventCreateWithFlags(&staged[stage_turn], cudaEventDisableTiming), "cudaEventCreate");
                } else {
                    check_cuda(cudaEventSynchronize(staged[stage_turn]), "restore staging");  // its last upload is done
                }
                char* h_stage = static_cast<char*>(pinned_pool().staging[stage_turn].reserve(bytes));
                #pragma omp parallel for schedule(dynamic)
                for (int64_t q = static_cast<int64_t>(q0); q < static_cast<int64_t>(q1); ++q) {
                    const Pending& pd = pending[static_cast<size_t>(q)];
                    const DataType* temp2 = pd.box->X_NR.data.data();
                    const int64_t ld = pd.box->X_NR.lda;
                    const DataType* xrr = pd.box->X_RR_full.data.data();
                    const int64_t ldx = pd.box->X_RR_full.lda;
                    const size_t nrows = pd.rows.size();
                    const size_t L = nrows + static_cast<size_t>(pd.r);
                    DataType* h = reinterpret_cast<DataType*>(h_stage + (pd.gall_off - base));
                    for (int t = 0; t < pd.r; ++t) {
                        DataType* col = h + static_cast<size_t>(t) * L;
                        const DataType* src = temp2 + static_cast<size_t>(t) * ld;
                        for (size_t c = 0; c < nrows; ++c) col[c] = src[pd.rows[c]];
                        for (int u = 0; u < pd.r; ++u) col[nrows + u] = xrr[t + static_cast<size_t>(u) * ldx];
                    }
                }
                check_cuda(cudaMemcpyAsync(d_raw + base, h_stage, bytes, cudaMemcpyHostToDevice, stream),
                           "restore fill sources");
                check_cuda(cudaEventRecord(staged[stage_turn], stream), "cudaEventRecord");
                stage_turn ^= 1;
                q0 = q1;
            }
            stats.bytes_up += static_cast<double>(restore_bytes);
            stats.tr_restore_bytes += static_cast<double>(restore_bytes);
            stats.tr_restore += std::chrono::duration<double>(clock::now() - t_restore).count();
        }

        // buffers of G and P for the chunk
        char* d_buffers = heap_.alloc(std::max<size_t>(buffer_bytes, 1));
        auto at = [&](size_t off) { return reinterpret_cast<S*>(d_buffers + off); };
        auto g_ptr = [&](uint64_t key) -> const S* {
            const GInfo& g = g_of.at(key);
            if (g.restore != SIZE_MAX) return restore_ptr(g.restore);
            return g.buffer == SIZE_MAX ? g.ptr : at(g.buffer);
        };
        for (const GPending& g : g_pending) {
            g_gathers.push_back(GatherItem{at(g.offset), g.r, g.r, g.k, g.src, 1, g.r, IndexList{}, g.cols});
            max_g_m = std::max(max_g_m, g.r);
            max_g_n = std::max(max_g_n, g.k);
        }
        VBatch<S> p_batch;
        for (const PPending& pp : p_pending) {
            const S* xrr = pp.xrr_restore == SIZE_MAX ? pp.xrr : restore_ptr(pp.xrr_restore);
            p_batch.entries.push_back({xrr, g_ptr(pp.g_key), at(pp.offset), pp.r, pp.k, pp.r, pp.r,
                                       g_of.at(pp.g_key).ld, pp.r});
        }
        // two batches per source rank: fills with rows j and fills with rows i
        std::vector<VBatch<S>> fill_j(fills.size()), fill_i(fills.size());
        for (size_t t = 0; t < fills.size(); ++t) {
            for (const FillTask& f : fills[t]) {
                const S* P = reinterpret_cast<const S*>(d_buffers + p_of.at(f.p_key_j));
                if (!f.rows_i) {
                    // out (m x n) -= P_{E,j}^T (m x r) G_{E,i} (r x n)
                    fill_j[t].entries.push_back({P, g_ptr(f.g_key_i), f.out, f.m, f.n, f.r, f.r, g_of.at(f.g_key_i).ld, f.ld});
                } else {
                    // out (m x n) -= G_{E,i}^T (m x r) P_{E,j} (r x n): the same terms, transposed
                    fill_i[t].entries.push_back({g_ptr(f.g_key_i), P, f.out, f.m, f.n, f.r, g_of.at(f.g_key_i).ld, f.r, f.ld});
                }
                ++fill_gemms;
                stats.transition_flops += kFlopScale * 2.0 * f.m * f.n * f.r;
            }
        }
        const size_t off_restores = meta.append(restores);
        const size_t off_evals = meta.append(evals);
        const size_t off_gathers = meta.append(gathers);
        const size_t off_g = meta.append(g_gathers);
        p_batch.stage(meta);
        for (auto& b : fill_j) b.stage(meta);
        for (auto& b : fill_i) b.stage(meta);
        stats.tr_plan += std::chrono::duration<double>(clock::now() - t_plan).count();
        ++stats.tr_chunks;
        char* md = meta.upload(meta_device, stream);
        stats.bytes_up += static_cast<double>(meta.size());
        chunk_marks.push_back(std::make_unique<StreamMarks>());
        StreamMarks& tmarks = *chunk_marks.back();
        tmarks.mark(stream);
        launch_transpose(reinterpret_cast<const TransposeItem*>(md + off_restores), static_cast<int>(restores.size()),
                         max_rs_m, max_rs_n, stream);
        eval_blocks(evals, meta, md, off_evals, max_eval_m, max_eval_n, stream);
        launch_gather(reinterpret_cast<const GatherItem*>(md + off_gathers), static_cast<int>(gathers.size()),
                      max_gather_m, max_gather_n, md, stream);
        launch_gather(reinterpret_cast<const GatherItem*>(md + off_g), static_cast<int>(g_gathers.size()), max_g_m,
                      max_g_n, md, stream);
        tmarks.mark(stream);
        p_batch.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
        tmarks.mark(stream);
        for (size_t t = 0; t < fills.size(); ++t) {
            fill_j[t].gemm(md, MagmaTrans, MagmaNoTrans, -1.0, 1.0, queue);
            fill_i[t].gemm(md, MagmaTrans, MagmaNoTrans, -1.0, 1.0, queue);
        }
        tmarks.mark(stream);
        heap_.free(d_buffers);  // later launches are ordered after this chunk
        heap_.free(d_raw);
        // restored sources whose last reader was in this chunk or before
        for (auto it = restored.begin(); it != restored.end();) {
            if (it->second.last_parent < p1) {
                heap_.free(it->second.block);
                it = restored.erase(it);
            } else {
                ++it;
            }
        }
        stats.heap_peak = std::max(stats.heap_peak, heap_.peak());
        p0 = p1;
    }
    for (auto& kv : restored) heap_.free(kv.second.block);
    for (cudaEvent_t e : staged) {
        if (e != nullptr) cudaEventDestroy(e);
    }
    stats.transition_fill_gemms += fill_gemms;
    check_cuda(cudaStreamSynchronize(stream), "device transition");
    for (const auto& m : chunk_marks) {
        stats.tr_blocks += m->seconds(0);
        stats.tr_p += m->seconds(1);
        stats.tr_fill += m->seconds(2);
    }

    // the child level is done: release its device data (not the new blocks)
    release_level_data();
    stats.transition += std::chrono::duration<double>(clock::now() - t0).count();
    return out;
}

// Free this level's device data (blocks, fill sources, point table).  Blocks
// handed to the parent level are not the level's any more.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::release_level_data() {
    released_sources_.clear();  // their blocks are freed with the states below
    if (!heap_initialized_) return;
    if (outbox_ != nullptr) free_block(outbox_);
    outbox_ = nullptr;
    for (auto& kv : edges_) heap_.free(kv.second.ptr);
    edges_.clear();
    for (BoxState& st : states_) {
        heap_.free(st.schur.ptr);
        st.schur = DeviceMatrix{};
    }
    for (BoxState& st : remote_states_) {
        heap_.free(st.schur.ptr);
        st.schur = DeviceMatrix{};
    }
    for (BoxState& st : states_) {
        heap_.free(st.persist);
        st.persist = nullptr;
    }
    for (BoxState& st : remote_states_) {
        heap_.free(st.persist);
        st.persist = nullptr;
    }
    for (char* p : level_allocs_) heap_.free(p);
    level_allocs_.clear();
    heap_.free(d_xyz_);
    heap_.free(d_ids_);
    d_xyz_ = nullptr;
    d_ids_ = nullptr;
    if (wave_sketch_ != nullptr) {
        heap_.free(wave_sketch_);
        wave_sketch_ = nullptr;
    }
    if (adopt_) {
        for (auto& kv : adopt_->schur) heap_.free(kv.second.ptr);
        for (auto& kv : adopt_->edges) heap_.free(kv.second.ptr);
        adopt_.reset();
    }
    heap_initialized_ = false;
}

// Device blocks of a level into its host boxes (for a level that runs on
// the host, e.g. level 1 before the root): Schur blocks, and each near-field
// block as the lower box's view (rows = the higher box) with an empty
// reciprocal for share_symmetric_level_edges.  The device copies are freed.
template<typename CoordType, typename DataType>
void download_level_blocks(DeviceLevelBlocks<DataType>& blocks, std::vector<BoxData<CoordType, DataType>>& boxes) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)blocks; (void)boxes;
    throw std::runtime_error("download_level_blocks: data type not run on the GPU");
  } else {
    using S = typename DeviceScalar<DataType>::type;
    static_assert(sizeof(S) == sizeof(DataType), "device and host elements share their layout");
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    DeviceHeap& heap = DeviceHeap::instance();
    if (boxes.empty()) return;
    const int64_t start = boxes.front().morton_index;
    auto box_of = [&](int64_t morton) -> BoxData<CoordType, DataType>& {
        const int64_t idx = morton - start;
        if (idx < 0 || idx >= static_cast<int64_t>(boxes.size())) {
            throw std::runtime_error("download_level_blocks: box is not local");
        }
        return boxes[static_cast<size_t>(idx)];
    };
    auto fetch = [&](const DeviceMatrixT<S>& m) {
        std::vector<DataType> data(static_cast<size_t>(m.rows) * m.cols);
        check_cuda(cudaMemcpy(data.data(), m.ptr, data.size() * sizeof(S), cudaMemcpyDeviceToHost),
                   "download level blocks");
        return data;
    };
    check_cuda(cudaStreamSynchronize(stream), "download level blocks");
    for (auto& kv : blocks.schur) {
        box_of(kv.first).schur_complement.set_owned(kv.second.rows, kv.second.cols, fetch(kv.second),
                                                    MatrixStorage<DataType>::FULL);
        heap.free(kv.second.ptr);
    }
    auto local_box = [&](int64_t morton) -> BoxData<CoordType, DataType>* {
        const int64_t idx = morton - start;
        return idx >= 0 && idx < static_cast<int64_t>(boxes.size()) ? &boxes[static_cast<size_t>(idx)] : nullptr;
    };
    for (auto& kv : blocks.edges) {
        const int64_t lo = static_cast<int64_t>(kv.first >> 32), hi = static_cast<int64_t>(kv.first & 0xffffffffu);
        if (local_box(lo) == nullptr || local_box(hi) == nullptr) {
            // a pair with a box of another rank: the local box's view
            // (rows = the other box's points)
            std::vector<DataType> data = fetch(kv.second);
            ModifiedBlock<DataType> view;
            if (auto* lb = local_box(lo)) {
                view.neighbor_morton = hi;
                view.set_a_ns_owned(kv.second.rows, kv.second.cols, std::move(data), MatrixStorage<DataType>::FULL);
                lb->near_field_interaction_map[hi] = static_cast<int64_t>(lb->near_field_modified_interactions.size());
                lb->near_field_modified_interactions.push_back(std::move(view));
            } else if (auto* hb = local_box(hi)) {
                std::vector<DataType> transposed(data.size());
                for (int j = 0; j < kv.second.cols; ++j)
                    for (int i = 0; i < kv.second.rows; ++i)
                        transposed[j + static_cast<size_t>(i) * kv.second.cols] = data[i + static_cast<size_t>(j) * kv.second.rows];
                view.neighbor_morton = lo;
                view.set_a_ns_owned(kv.second.cols, kv.second.rows, std::move(transposed), MatrixStorage<DataType>::FULL);
                hb->near_field_interaction_map[lo] = static_cast<int64_t>(hb->near_field_modified_interactions.size());
                hb->near_field_modified_interactions.push_back(std::move(view));
            } else {
                throw std::runtime_error("download_level_blocks: edge without a local box");
            }
            heap.free(kv.second.ptr);
            continue;
        }
        auto& lo_box = box_of(lo);
        auto& hi_box = box_of(hi);
        ModifiedBlock<DataType> lo_block;
        lo_block.neighbor_morton = hi;
        lo_block.set_a_ns_owned(kv.second.rows, kv.second.cols, fetch(kv.second), MatrixStorage<DataType>::FULL);
        lo_box.near_field_interaction_map[hi] = static_cast<int64_t>(lo_box.near_field_modified_interactions.size());
        lo_box.near_field_modified_interactions.push_back(std::move(lo_block));
        ModifiedBlock<DataType> hi_block;
        hi_block.neighbor_morton = lo;
        hi_box.near_field_interaction_map[lo] = static_cast<int64_t>(hi_box.near_field_modified_interactions.size());
        hi_box.near_field_modified_interactions.push_back(std::move(hi_block));
        heap.free(kv.second.ptr);
    }
    blocks.schur.clear();
    blocks.edges.clear();
  }
}

// LU of the root block, left on the device by the device transition of
// level 1, with the factors copied to the root box (as the host does:
// X_RR = LU, 1-based pivots; of the block's symmetric part when `symmetric`).
// The device block is freed.
template<typename CoordType, typename DataType>
void factor_root_on_device(DeviceLevelBlocks<DataType>& blocks, BoxData<CoordType, DataType>& root, bool symmetric) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)blocks; (void)root; (void)symmetric;
    throw std::runtime_error("factor_root_on_device: data type not run on the GPU");
  } else {
    using S = typename DeviceScalar<DataType>::type;
    auto it = blocks.schur.find(root.morton_index);
    if (it == blocks.schur.end() || blocks.schur.size() != 1 || !blocks.edges.empty()) {
        throw std::runtime_error("factor_root_on_device: expected the root block alone");
    }
    const DeviceMatrixT<S> m = it->second;
    if (m.rows != root.num_points || m.cols != root.num_points) {
        throw std::runtime_error("factor_root_on_device: root block shape");
    }
    Context& ctx = Context::instance();
    ctx.activate();
    const magma_int_t n = m.rows;
    if (symmetric) launch_symmetrize(m.ptr, static_cast<int>(n), static_cast<int>(n), ctx.stream());
    check_cuda(cudaStreamSynchronize(ctx.stream()), "root block");
    std::vector<magma_int_t> piv(static_cast<size_t>(n));
    magma_int_t info = 0;
    if constexpr (std::is_same_v<DataType, double>) {
        magma_dgetrf_gpu(n, n, m.ptr, n, piv.data(), &info);
    } else {
        magma_zgetrf_gpu(n, n, reinterpret_cast<magmaDoubleComplex*>(m.ptr), n, piv.data(), &info);
    }
    if (info != 0) {
        throw std::runtime_error("factor_root_on_device: LU factorization of root failed with INFO = " +
                                 std::to_string(info));
    }
    std::vector<DataType> lu(static_cast<size_t>(n) * static_cast<size_t>(n));
    check_cuda(cudaMemcpy(lu.data(), m.ptr, lu.size() * sizeof(S), cudaMemcpyDeviceToHost), "root factors");
    DeviceHeap::instance().free(m.ptr);
    blocks.schur.clear();
    root.X_RR.set_owned(n, n, std::move(lu), MatrixStorage<DataType>::LU_FACTORED);
    root.X_RR_pivots.assign(piv.begin(), piv.end());
  }
}

// LU of a root block assembled on the host (multi-rank runs, whose level-1
// blocks pass through the host for the process reduction) on the device,
// when the heap has room for it: the root box gets the factors as the host
// LU leaves them (X_RR = LU, 1-based pivots; of the block's symmetric part
// when `symmetric`).  Returns false, with nothing changed, otherwise.
template<typename CoordType, typename DataType>
bool factor_host_root_on_device(BoxData<CoordType, DataType>& root, bool symmetric) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)root; (void)symmetric;
    return false;
  } else {
    using S = typename DeviceScalar<DataType>::type;
    const MatrixStorage<DataType>& A = root.schur_complement;
    if (!A.is_allocated() || A.rows <= 0 || A.cols != A.rows || A.lda != A.rows) return false;
    DeviceHeap& heap = DeviceHeap::instance();
    if (!heap.initialized()) return false;
    const magma_int_t n = static_cast<magma_int_t>(A.rows);
    const size_t bytes = static_cast<size_t>(n) * static_cast<size_t>(n) * sizeof(S);
    Context& ctx = Context::instance();
    ctx.activate();
    check_cuda(cudaStreamSynchronize(ctx.stream()), "root block");
    char* d = heap.try_alloc(bytes);
    if (d == nullptr) return false;
    check_cuda(cudaMemcpy(d, A.data.data(), bytes, cudaMemcpyHostToDevice), "root block");
    if (symmetric) {
        launch_symmetrize(reinterpret_cast<S*>(d), static_cast<int>(n), static_cast<int>(n), ctx.stream());
        check_cuda(cudaStreamSynchronize(ctx.stream()), "root block");
    }
    std::vector<magma_int_t> piv(static_cast<size_t>(n));
    magma_int_t info = 0;
    if constexpr (std::is_same_v<DataType, double>) {
        magma_dgetrf_gpu(n, n, reinterpret_cast<double*>(d), n, piv.data(), &info);
    } else {
        magma_zgetrf_gpu(n, n, reinterpret_cast<magmaDoubleComplex*>(d), n, piv.data(), &info);
    }
    if (info != 0) {
        heap.free(d);
        throw std::runtime_error("factor_host_root_on_device: LU factorization of root failed with INFO = " +
                                 std::to_string(info));
    }
    std::vector<DataType> lu(static_cast<size_t>(n) * static_cast<size_t>(n));
    check_cuda(cudaMemcpy(lu.data(), d, bytes, cudaMemcpyDeviceToHost), "root factors");
    heap.free(d);
    root.X_RR.set_owned(n, n, std::move(lu), MatrixStorage<DataType>::LU_FACTORED);
    root.X_RR_pivots.assign(piv.begin(), piv.end());
    return true;
  }
}

// Whether the device box path would take level `level_index` (same test as
// make_level_eliminator), for deciding where the parent blocks go.
template<typename CoordType, typename DataType, typename KernelType>
bool level_eliminator_would_run(const ParallelTree<CoordType, DataType>* tree, int level_index, const KernelType* kernel,
                                FactorizationMethod method) {
    if constexpr (gpu_data_type<DataType>) {
        return level_eliminator_supported(tree->levels[static_cast<size_t>(level_index)], kernel, tree->dimension,
                                          method, nullptr);
    } else {
        (void)tree; (void)level_index; (void)kernel; (void)method;
        return false;
    }
}

// The device box path for this level, uploaded and ready, or nullptr (with
// the reason) when the level stays on the host.
template<typename CoordType, typename DataType, typename KernelType>
std::unique_ptr<LevelEliminatorBase<CoordType, DataType>> make_level_eliminator(
    ParallelTree<CoordType, DataType>* tree, int level_index, KernelType* kernel, double tolerance,
    FactorizationMethod method, std::string* reason, std::unique_ptr<DeviceLevelBlocks<DataType>> adopt = nullptr,
    bool occupancy = false, bool device_sketch = true) {
    if constexpr (gpu_data_type<DataType>) {
        const auto& level = tree->levels[static_cast<size_t>(level_index)];
        if (!level_eliminator_supported(level, kernel, tree->dimension, method, reason)) {
            if (adopt) throw std::runtime_error("make_level_eliminator: device blocks for a level the device cannot run");
            return nullptr;
        }
        auto eliminator = std::make_unique<LevelEliminator<CoordType, DataType, KernelType>>(
            tree, level_index, kernel, tolerance, std::move(adopt), occupancy, device_sketch);
        eliminator->begin();
        return eliminator;
    } else {
        (void)tree; (void)level_index; (void)kernel; (void)tolerance; (void)method; (void)adopt; (void)occupancy;
        (void)device_sketch;
        if (reason) *reason = "single precision is not run on the GPU";
        return nullptr;
    }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
