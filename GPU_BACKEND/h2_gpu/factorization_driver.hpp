#pragma once
// The GPU side of the Color factorization loop (H2_use_gpu), shared by the
// full-grid loop (butterfly_factorization.hpp) and the unstructured one
// (unstructured/factorization_impl.hpp).  Per level: the device box
// path (level_eliminator.hpp) and its transports, the level end, the level
// transition, and the root.  The loops keep their host code and their own
// wave schedules; each call here replaces the host step it names.
//
// Replicated CA levels (start_ca_level) run in the same eliminator with their
// ghost boxes; the loop's CA schedule calls eliminate_wave, and nothing is
// transported during the level (h2_gpu/GPU_CA_PLAN.md, M1).
//
// A tree of occupied boxes (`occupancy`, the unstructured backend) keeps
// empty boxes in the local slabs; the eliminator skips them, and its
// transports follow the unstructured host transport's protocol (every rank
// in the assisting exchange, generators also to the assisting boxes'
// requesters).  The loop's host_transport callback sets the host
// transport's flags when the device exchange is off.

#ifdef H2_HAVE_GPU

#include "ca_halo.hpp"
#include "collective_evaluator.hpp"
#include "device_reduction.hpp"
#include "level_eliminator.hpp"
#include "owner_plan.hpp"
#include "wave_trace.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fmm {
namespace gpu {

template<typename CoordType, typename DataType, typename KernelType>
class ColorGpuDriver {
public:
    using Tree = ParallelTree<CoordType, DataType>;
    using Level = TreeLevel<CoordType, DataType>;
    using Box = BoxData<CoordType, DataType>;
    using Duration = std::chrono::high_resolution_clock::duration;

    struct Options {
        double tolerance = 0.0;
        FactorizationMethod method = FactorizationMethod::LU;
        int use_sketch = 0;
        int lazy_schur = 0;
        bool is_symmetric = true;
        bool is_hermitian = false;
        bool occupancy = false;
        int ca_owner_component = 0;  // H2_CA_owner_component
        int ca_pack_shared = 0;      // H2_CA_pack_shared: 0 auto, 1 always, 2 never
        int ca_solve = 0;            // H2_CA_solve: 0 replicated, 1 home
        int compress_factors = 0;    // H2_compress_factors (tolerance: `tolerance`)
    };

    ColorGpuDriver(Tree* tree, KernelType* kernel, const Options& options)
        : tree_(tree), kernel_(kernel), opt_(options) {}

    // Before the first level of the process's first factorization (of this
    // data type): the batched LU, the two solves of X_RR^{-1} and the three
    // GEMM forms once per size range, so no level pays the first-use loading
    // of their kernels (lazy module loading; up to ~0.3 s inside a level).
    // Returns its time (0 when it already ran), which the factor time leaves
    // out (butterfly_factorization.hpp).
    double warm_up() {
        static bool done = false;
        if (done || !color_gpu_enabled()) return 0.0;
        done = true;
        if constexpr (!gpu_data_type<DataType>) {
            return 0.0;
        } else {
            using S = typename DeviceScalar<DataType>::type;
            const auto t0 = std::chrono::steady_clock::now();
            Context& ctx = Context::instance();
            ctx.activate();
            cudaStream_t stream = ctx.stream();
            magma_queue_t queue = ctx.queue();
            MetaBuilder& meta = pinned_pool().meta;
            DeviceBuffer meta_device, work;
            for (int n : {8, 16, 32, 64, 128, 256, 512, 1024}) {
                const size_t bytes = static_cast<size_t>(n) * n * sizeof(S);
                std::vector<S> identity(static_cast<size_t>(n) * n, S(0.0));
                for (int i = 0; i < n; ++i) identity[static_cast<size_t>(i) * n + i] = S(1.0);
                S *d_a = nullptr, *d_b = nullptr, *d_c = nullptr;
                magma_int_t *d_piv = nullptr, *d_info = nullptr;
                magma_int_t** d_piv_ptr = nullptr;
                check_cuda(cudaMalloc(&d_a, bytes), "warm-up");
                check_cuda(cudaMalloc(&d_b, bytes), "warm-up");
                check_cuda(cudaMalloc(&d_c, bytes), "warm-up");
                check_cuda(cudaMalloc(&d_piv, n * sizeof(magma_int_t)), "warm-up");
                check_cuda(cudaMalloc(&d_info, sizeof(magma_int_t)), "warm-up");
                check_cuda(cudaMalloc(&d_piv_ptr, sizeof(magma_int_t*)), "warm-up");
                check_cuda(cudaMemcpy(d_a, identity.data(), bytes, cudaMemcpyHostToDevice), "warm-up");
                check_cuda(cudaMemcpy(d_b, identity.data(), bytes, cudaMemcpyHostToDevice), "warm-up");
                check_cuda(cudaMemcpy(d_piv_ptr, &d_piv, sizeof(magma_int_t*), cudaMemcpyHostToDevice), "warm-up");
                // (a, lda) the LU, (c, ldc) the solves' right side: as the eliminator's batches
                VBatch<S> lu;
                lu.entries.push_back({d_a, nullptr, d_b, n, n, 1, n, 1, n});
                VBatch<S> g;
                g.entries.push_back({d_a, d_b, d_c, n, n, n, n, n, n});
                meta.clear();
                lu.stage(meta);
                g.stage(meta);
                char* md = meta.upload(meta_device, stream);
                getrf_vbatched<S>(n, lu.size_array(md, 0), lu.size_array(md, 1), lu.template pointer_array<S*>(md, 0),
                                  lu.size_array(md, 3), d_piv_ptr, d_info, 1, work, queue);
                for (magma_uplo_t uplo : {MagmaUpper, MagmaLower}) {
                    trsm_vbatched<S>(MagmaRight, uplo, MagmaNoTrans, uplo == MagmaUpper ? MagmaNonUnit : MagmaUnit, n,
                                     n, lu.size_array(md, 0), lu.size_array(md, 1), S(1.0),
                                     lu.template pointer_array<S*>(md, 0), lu.size_array(md, 3),
                                     lu.template pointer_array<S*>(md, 2), lu.size_array(md, 5), 1, queue);
                }
                g.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
                g.gemm(md, MagmaNoTrans, MagmaTrans, 1.0, 1.0, queue);
                g.gemm(md, MagmaTrans, MagmaNoTrans, 1.0, 1.0, queue);
                check_cuda(cudaStreamSynchronize(stream), "warm-up");
                for (void* p : {static_cast<void*>(d_a), static_cast<void*>(d_b), static_cast<void*>(d_c),
                                static_cast<void*>(d_piv), static_cast<void*>(d_info), static_cast<void*>(d_piv_ptr)}) {
                    cudaFree(p);
                }
            }
            meta.clear();
            return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
    }

    // Before the first level of the process's first factorization: the
    // page-locked staging buffers the levels use (pinned_pool(), kept for
    // the whole run): the copiers' ring slots, the metadata images of the
    // waves and of their owner passes, the sketches' result buffer.  Page
    // locking costs about a second per GB and otherwise lands in the first
    // level that uses each buffer (the copier threads' ring slots held up
    // the level start by 0.1-0.2 s, and by 0.08 s again at the first level
    // with exchanges after a CA leaf).  Returns its time (0 when it already
    // ran), which the factor time leaves out.
    double warm_up_buffers() {
        static bool done = false;
        if (done || !color_gpu_enabled()) return 0.0;
        done = true;
        if constexpr (!gpu_data_type<DataType>) {
            return 0.0;
        } else {
            const auto t0 = std::chrono::steady_clock::now();
            Context::instance().activate();
            PinnedPool& pool = pinned_pool();
            pool.meta.reserve(1);
            pool.owner_meta.reserve(1);
            pool.result.reserve(1);
            const int rings = tree_->mpi_size > 1 ? 4 : 2;  // (the exchange copier's: multi-rank levels)
            for (int s = 0; s < rings; ++s) pool.ring[s].reserve(PinnedPool::kRingSlot);
            return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
    }

    // Before the first level: the application's evaluator, when a level may
    // take the device box path, which calls it (Evaluator::warm_up, once per
    // evaluator: an NVRTC compile, CuPy compiling its kernels, ...).  Returns
    // its time (0 when there is nothing to warm), which the factor time
    // leaves out.
    double warm_up_evaluator() {
        if constexpr (!gpu_data_type<DataType>) {
            return 0.0;
        } else {
            if (!color_gpu_enabled() || opt_.use_sketch != 2 || opt_.is_hermitian ||
                opt_.method != FactorizationMethod::LU) {
                return 0.0;
            }
            return warm_up_registered_evaluator<DataType>(tree_, kernel_->gpu_evaluator, opt_.is_symmetric);
        }
    }

    // Before the first level of the process's first factorization: a 64 KB
    // message with every rank this one exchanges with on a level (each
    // level's one-hop neighbour ranks), sent and received the way the levels
    // will (from and into the exchange arena with GPU-aware MPI).  A pair's
    // first device exchanges are slow (Laplace 192^3 on 64 ranks: the first
    // waves of the first exchanging level at 3-4 GB/s, ~12 GB/s from the
    // fourth on; ~0.1-0.2 s per factorization, more at larger sizes); this
    // takes that first contact out of the levels.  Then the same from and
    // into the main heap, where a level's larger messages land once the
    // arena is full: GPU-aware MPI meets a buffer's whole allocation on a
    // pair's first transfer from it, and the heap is one allocation of most
    // of the card (2026-10-12: after a CA0 leaf, which exchanges nothing,
    // the first four transports of the level above waited 3.4 s more than
    // after a Color leaf at 384^3 on 8 ranks, 6.5 s at 768^3 on 64 ranks on
    // 80 GB cards; a Color leaf paid the same in its own first waves).
    // Returns its time (0 when it already ran); collective over the tree's
    // ranks.
    double warm_up_exchange() {
        static bool done = false;
        if (done || !color_gpu_enabled() || tree_->mpi_size <= 1) return 0.0;
        done = true;
        const auto t0 = std::chrono::steady_clock::now();
        std::set<int> peer_set;
        for (int l = 2; l < tree_->num_levels; ++l) {
            for (int r : compute_one_hop_neighbor_ranks(tree_, tree_->levels[static_cast<size_t>(l)], l)) {
                if (r != tree_->mpi_rank) peer_set.insert(r);
            }
        }
        const std::vector<int> peers(peer_set.begin(), peer_set.end());
        constexpr size_t kBytes = size_t{64} << 10;
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        const size_t total = kBytes * std::max<size_t>(peers.size(), 1);
        char* send = device_exchange_enabled() ? arena.try_alloc(total) : nullptr;
        char* recv = send != nullptr ? arena.try_alloc(total) : nullptr;
        std::vector<char> host_send, host_recv;
        if (recv == nullptr) {  // (host messages: MPI is not GPU-aware, or no room)
            if (send != nullptr) arena.free(send);
            host_send.resize(total);
            host_recv.resize(total);
            send = host_send.data();
            recv = host_recv.data();
        }
        std::vector<MPI_Request> reqs;
        for (size_t i = 0; i < peers.size(); ++i) {
            reqs.emplace_back();
            MPI_Irecv(recv + i * kBytes, static_cast<int>(kBytes), MPI_BYTE, peers[i], 790, tree_->comm, &reqs.back());
        }
        for (size_t i = 0; i < peers.size(); ++i) {
            reqs.emplace_back();
            MPI_Isend(send + i * kBytes, static_cast<int>(kBytes), MPI_BYTE, peers[i], 790, tree_->comm, &reqs.back());
        }
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        if (host_send.empty()) {
            arena.free(send);
            arena.free(recv);
        }
        if (device_exchange_enabled()) {   // (the main heap: every rank takes part, a heap or not)
            DeviceHeap& heap = DeviceHeap::instance();
            heap.ensure_initialized();
            char* hs = heap.try_alloc(total);
            char* hr = hs != nullptr ? heap.try_alloc(total) : nullptr;
            int ok = hr != nullptr ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, tree_->comm);
            if (ok != 0) {
                reqs.clear();
                for (size_t i = 0; i < peers.size(); ++i) {
                    reqs.emplace_back();
                    MPI_Irecv(hr + i * kBytes, static_cast<int>(kBytes), MPI_BYTE, peers[i], 791, tree_->comm, &reqs.back());
                }
                for (size_t i = 0; i < peers.size(); ++i) {
                    reqs.emplace_back();
                    MPI_Isend(hs + i * kBytes, static_cast<int>(kBytes), MPI_BYTE, peers[i], 791, tree_->comm, &reqs.back());
                }
                MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
            }
            if (hs != nullptr) heap.free(hs);
            if (hr != nullptr) heap.free(hr);
        }
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }

    // ---- where the blocks of a level go next (tested where levels start)

    // The device box path runs level `lvl`.
    bool box_path_runs(int lvl) const {
        if (!color_gpu_enabled() || lvl <= 1 || tree_->level_uses_CA(lvl)) return false;
        if (!tree_->levels[static_cast<size_t>(lvl)].is_process_active) return false;
        // (a nonsymmetric matrix too: level_eliminator_would_run tests what it needs)
        const bool streamed = opt_.use_sketch == 2 && !opt_.is_hermitian;
        if (!streamed || opt_.lazy_schur == 0) return false;
        // (with the level's own lazy mode: after a CA level the runtime still
        // has that level's, capped at 1)
        return level_eliminator_would_run(tree_, lvl, kernel_, opt_.method, opt_.lazy_schur);
    }
    // Level 1 is not eliminated: its transition builds the root, which is
    // factored on the device when level 1 runs on one rank.
    bool root_path_runs() const {
        if (!color_gpu_enabled() || tree_->num_levels < 2) return false;
        const auto& l1 = tree_->levels[1];
        const auto& l0 = tree_->levels[0];
        if (!l1.is_process_active || l1.num_active_processes != 1) return false;
        if (!l0.is_process_active || l0.num_active_processes != 1) return false;
        if (opt_.method != FactorizationMethod::LU) return false;
        return level_eliminator_would_run(tree_, 1, kernel_, opt_.method);
    }
    bool keeps_blocks_of(int lvl) const {
        return lvl == 0 ? true : (lvl == 1 ? root_path_runs() : box_path_runs(lvl));
    }
    // Whether the parents' blocks can stay on the device across the process
    // reduction from level `lvl` into level `parent_lvl`: the parent level
    // would run on the device on its owners (the tests of keeps_blocks_of
    // without this rank's own activity there, which a sender lacks), every
    // region of the child level is whole parent boxes (device_reduction.hpp:
    // the equal cubes, or the Morton ranges of the unstructured path), and on
    // a CA parent level this rank's parents fit the device sketch.
    // Rank-independent but for the last test; the caller all-reduces it over
    // the transition's ranks.
    bool device_reduction_possible(int lvl, int parent_lvl, const std::vector<Box>& parents) const {
        if (!color_gpu_enabled() || parent_lvl < 1 || lvl < 2) return false;
        const Level& child = tree_->levels[static_cast<size_t>(lvl)];
        const Level& parent = tree_->levels[static_cast<size_t>(parent_lvl)];
        if (!regions_hold_whole_parents(child, tree_->dimension) || parent.region_start.size() < 2) return false;
        if (parent_lvl == 1) {
            // level 1 adopts the blocks on its rank; its transition builds the root there
            const Level& l0 = tree_->levels[0];
            return parent.num_active_processes == 1 && l0.num_active_processes == 1 &&
                   opt_.method == FactorizationMethod::LU && level_eliminator_would_run(tree_, 1, kernel_, opt_.method);
        }
        const bool streamed = opt_.use_sketch == 2 && !opt_.is_hermitian;
        if (!streamed) return false;
        if (tree_->level_uses_CA(parent_lvl)) {
            if (!ca_level_runs(parent_lvl) || !ca_level_supported(parent_lvl, nullptr)) return false;
            for (const Box& p : parents) {
                if (p.num_points > kEntryDestMask) return false;
            }
            return true;
        }
        if (opt_.lazy_schur == 0) return false;
        return level_eliminator_would_run(tree_, parent_lvl, kernel_, opt_.method, opt_.lazy_schur);
    }

    // ---- a Color level that eliminates

    // Level start (all ranks of level_comm): the device box path when it
    // covers the level (streamed sketches, lazy far fill, a symmetric
    // kernel).  Its sketches run on the host when a box of any rank is
    // wider than the device sketch's packed lists take.  `announce`: this
    // rank prints the decision.
    // Host staging of the solve factors.  A level's kept store is read by
    // nothing of the factorization once its transition is done (the kept-LU
    // logdet is summed here first), so at the start of a level the deeper
    // levels' stores go to pinned host memory, and come back at the first
    // solve.  What it buys is the store's bytes under every later level.
    // The bytes every level's host boxes still hold, by field, at a level's
    // end: the host side of 8x the boxes a rank (2026-10-10).
    void report_host_boxes(int at_level) {
        auto mat = [](const MatrixStorage<DataType>& m) { return static_cast<double>(m.data.size()) * sizeof(DataType); };
        auto vec = [](const auto& v) { return static_cast<double>(v.size()) * sizeof(v[0]); };
        struct Sum {
            double T = 0, xrr = 0, piv = 0, xrs = 0, xsr = 0, schur = 0, xrn = 0, xnr = 0, xrrf = 0, xrse = 0, t2 = 0,
                   orig = 0, pts = 0, lists = 0;
            int64_t n = 0;
            double total() const { return T + xrr + piv + xrs + xsr + schur + xrn + xnr + xrrf + xrse + t2 + orig + pts + lists; }
        };
        auto add = [&](Sum& s, const Box& b) {
            ++s.n;
            s.T += mat(b.interpolation_matrix); s.xrr += mat(b.X_RR); s.piv += vec(b.X_RR_pivots);
            s.xrs += mat(b.X_RS); s.xsr += mat(b.X_SR); s.schur += mat(b.schur_complement); s.xrn += mat(b.X_RN);
            s.xnr += mat(b.X_NR); s.xrrf += mat(b.X_RR_full); s.xrse += mat(b.X_RS_entry);
            s.t2 += vec(b.deferred_xnn_temp2); s.orig += vec(b.lazy_original_x_nr);
            s.pts += vec(b.point_coords) + vec(b.point_indices);
            s.lists += vec(b.redundant_indices) + vec(b.skeleton_indices) + vec(b.one_hop) + vec(b.two_hop) +
                       vec(b.use_full_set) + vec(b.deferred_xnn_neighbor_point_counts);
        };
        auto print = [&](int lvl, const char* what, const Sum& s) {
            if (s.n == 0 || s.total() < 1e6) return;
            std::printf("  [gpu] at level %d end, level %d host %s (%lld): %.2f GB = T %.2f, LU %.2f, piv %.2f, X_RS %.2f, "
                        "X_SR %.2f, Schur %.2f, X_RN %.2f, X_NR %.2f, X_RR_full %.2f, X_RS_entry %.2f, temp2 %.2f, "
                        "orig X_NR %.2f, points %.2f, lists %.2f\n", at_level, lvl, what, static_cast<long long>(s.n),
                        s.total() / 1e9, s.T / 1e9, s.xrr / 1e9, s.piv / 1e9, s.xrs / 1e9, s.xsr / 1e9, s.schur / 1e9,
                        s.xrn / 1e9, s.xnr / 1e9, s.xrrf / 1e9, s.xrse / 1e9, s.t2 / 1e9, s.orig / 1e9, s.pts / 1e9,
                        s.lists / 1e9);
        };
        for (size_t lvl = 0; lvl < tree_->levels.size(); ++lvl) {
            const Level& L = tree_->levels[lvl];
            Sum l, g, r;
            for (const Box& b : L.local_boxes) add(l, b);
            for (const Box& b : L.ghost_boxes) add(g, b);
            for (const Box& b : L.generator_boxes) add(r, b);
            print(static_cast<int>(lvl), "local boxes", l);
            print(static_cast<int>(lvl), "ghost boxes", g);
            print(static_cast<int>(lvl), "generator boxes", r);
        }
    }
    bool pressure_seen_ = false;   // a level packed or chunked for room: memory is tight from here on
    void stage_solve_factors(int lvl, bool announce, MPI_Comm level_comm) {
        // the host arenas' next chunks pin while the level sketches (H2 host staging (c2))
        pinned_chunk_pool().pin_ahead();
        DeviceSolveStore& store = device_solve_store();
        double pin_s = 0.0, copy_s = 0.0, bytes = 0.0;
        std::string which;
        // The decision is the level's, not the rank's (2026-10-11): a rank
        // that saw no pressure kept its upper levels' factors on the card
        // while its neighbours staged theirs, and at 768^3 on 64 ranks it
        // then met the level-3 transition and the 64 -> 8 reduction with 27
        // GB less room than they had.
        int want = pressure_seen_ || device_under_pressure() ? 1 : 0;
        const bool host_trace = lvl <= 7 && wave_trace_enabled();   // (8x diagnostics: every rank's host memory around the staging)
        auto host_line = [&](const char* when, double staged) {
            if (!host_trace) return;
            double rss = 0.0, hwm = 0.0;
            host_memory_gb(rss, hwm);
            int rank = 0;
            MPI_Comm_rank(tree_->comm, &rank);
            const PinnedChunkPool& pool = pinned_chunk_pool();
            std::printf("  [gpu] level %d start, rank %d, %s: staged %.2f GB, host RSS %.2f GB (high-water %.2f), pinned %.2f GB (%.2f idle)\n",
                        lvl, rank, when, staged / 1e9, rss, hwm, static_cast<double>(pool.pinned()) / 1e9,
                        static_cast<double>(pool.idle()) / 1e9);
            std::fflush(stdout);
        };
        host_line("before staging", 0.0);
        if (level_comm != MPI_COMM_NULL) MPI_Allreduce(MPI_IN_PLACE, &want, 1, MPI_INT, MPI_MAX, level_comm);
        for (auto& [level, kl] : store.kept) {
            if (level <= lvl || kl.on_host() || kl.blocks.empty()) continue;
            if (want == 0) continue;   // room enough on every rank: the factors stay
            double la = 0.0, ar = 0.0;
            device_logdet_kept_level(tree_, level, la, ar);  // (kept in the store before the LU leaves)
            stage_kept_level(kl, Context::instance().stream(), pin_s, copy_s);
            bytes += static_cast<double>(kl.staged_bytes);
            which += (which.empty() ? "" : ",") + std::to_string(level);
        }
        host_line("after staging", bytes);
        if (announce && bytes > 0.0) {
            std::printf("  [gpu] level %d start: the solve factors of level%s %s staged to the host, %.2f GB "
                        "(allocation %.2f s, copy %.2f s)\n", lvl, which.find(',') != std::string::npos ? "s" : "",
                        which.c_str(), bytes / 1e9, pin_s, copy_s);
        }
    }
    // (8x diagnostics) every rank's card at a level's start, the small levels
    void card_at_level_start(int lvl) {
        if (lvl > 4 || !wave_trace_enabled()) return;
        const DeviceHeap& h = DeviceHeap::instance();
        if (!h.initialized()) return;
        int rank = 0;
        MPI_Comm_rank(tree_->comm, &rank);
        std::string held;
        for (size_t k = 0; k < static_cast<size_t>(ResidentKind::Count); ++k) {
            const size_t bytes = h.resident_now(static_cast<ResidentKind>(k));
            if (bytes < (size_t{256} << 20)) continue;
            char b[96];
            std::snprintf(b, sizeof b, "%s%s %.2f", held.empty() ? "" : ", ", resident_kind_name(static_cast<ResidentKind>(k)),
                          static_cast<double>(bytes) / 1e9);
            held += b;
        }
        std::printf("  [gpu] level %d rank %d card at start: %.2f of %.2f GB in use; %s\n", lvl, rank,
                    static_cast<double>(h.used()) / 1e9, static_cast<double>(h.capacity()) / 1e9, held.c_str());
        std::fflush(stdout);
    }
    void start_level(int lvl, bool use_streamed_level, MPI_Comm level_comm, bool announce) {
        stage_solve_factors(lvl, announce, level_comm);
        card_at_level_start(lvl);
        {
            // H2 host memory: the pinned chunks the level above gave back, now
            // that the kept factors took what they could, stay pinned for the
            // levels below unless the job is short of host memory on its node
            // (2026-10-11).  Unpinning half of them at every level start made
            // the next level's stores pin again, and pinning runs 10x slower on
            // a node whose memory is fragmented (6 s against 62 s a rank at
            // 768^3, the level waiting on the slowest rank); short of room they
            // all go but one chunk's worth.
            PinnedChunkPool& pool = pinned_chunk_pool();
            double s = 0.0;
            const bool short_host = host_room_bytes() < kHostReserveBytes;
            const size_t bytes = short_host ? pool.trim(PinnedChunkPool::kIdleKeep, &s) : 0;
            if (announce && bytes > 0) {
                std::printf("  [gpu] level %d start: %.2f GB of idle pinned chunks unpinned (%.2f s), %.2f GB pinned now "
                            "(%.2f idle)\n", lvl, static_cast<double>(bytes) / 1e9, s,
                            static_cast<double>(pool.pinned()) / 1e9, static_cast<double>(pool.idle()) / 1e9);
            }
        }
        if (!(color_gpu_enabled() && use_streamed_level && lazy_far_field_mode() == LazyFarFieldMode::LAZY &&
              !opt_.is_hermitian)) {
            return;
        }
        if (wave_trace_enabled()) trace_.start_begin();
        int fits = 1;
        for (const Box& box : tree_->levels[static_cast<size_t>(lvl)].local_boxes) {
            if (box.num_points > kEntryDestMask) fits = 0;
        }
        MPI_Allreduce(MPI_IN_PLACE, &fits, 1, MPI_INT, MPI_MIN, level_comm);
        std::string reason;
        level_ = make_level_eliminator(tree_, lvl, kernel_, opt_.tolerance, opt_.method, &reason,
                                       std::move(blocks_), opt_.occupancy, fits != 0,
                                       /*ca=*/false, /*owner=*/false, opt_.compress_factors);
        level_index_ = lvl;
        decide_solve_keep(lvl, level_comm, false);
        if (wave_trace_enabled()) trace_.start_end();
        if (announce) {
            std::cout << "  GPU box path: "
                      << (level_ ? std::string(tensor_core_gemm() && std::is_same_v<DataType, double> ? "on (FP64 tensor-core GEMMs)" : "on")
                                 : "off (" + reason + ")")
                      << (level_ && !fits ? ", host sketches (boxes over " + std::to_string(kEntryDestMask) +
                                                " points)"
                                          : std::string())
                      << std::endl;
        }
    }
    // ---- a replicated CA level (owner component 0) that eliminates

    // The device runs CA level `lvl` (all ranks decide alike): streamed
    // sketches, lazy far fill, a symmetric kernel, the LU of X_RR.
    bool ca_level_runs(int lvl, int owner_component) const {
        if (!color_gpu_enabled() || lvl <= 1 || !tree_->level_uses_CA(lvl)) return false;
        // owner component 3 (ca_owner.hpp): its payloads move between the
        // devices, so it needs CUDA-aware MPI
        if (owner_component != 0 && !(owner_component == 3 && device_exchange_enabled())) return false;
        const bool streamed = opt_.use_sketch == 2 && opt_.is_symmetric && !opt_.is_hermitian;
        return streamed && opt_.lazy_schur > 0 && opt_.method == FactorizationMethod::LU && !opt_.occupancy;
    }
    // ---- a component-owner CA level (owner component 3; ca_owner.hpp)

    // The level's blocks a device transition kept (the owned components'
    // remote boxes' entry blocks are fetched into them before the start).
    bool has_device_blocks() const { return blocks_ != nullptr; }
    DeviceLevelBlocks<DataType>* device_blocks() { return blocks_.get(); }
    // Level start, after the host halo gather and the fetch (all ranks of
    // level_comm): the eliminator in owner mode with the level's plan.
    // Returns whether the device runs the level (every box fits the device
    // sketch's packed lists on every rank, and the eliminator supports it).
    bool start_owner_level(int lvl, MPI_Comm level_comm, const OwnerPlanView* view, bool announce) {
        stage_solve_factors(lvl, announce, level_comm);
        if (wave_trace_enabled()) trace_.start_begin();
        const Level& level = tree_->levels[static_cast<size_t>(lvl)];
        int fits = 1;
        for (const Box& box : level.local_boxes) {
            if (box.num_points > kEntryDestMask) fits = 0;
        }
        for (const Box& box : level.ghost_boxes) {
            if (box.num_points > kEntryDestMask) fits = 0;
        }
        std::string reason = "boxes over " + std::to_string(kEntryDestMask) + " points";
        if (fits && !ca_level_supported(lvl, &reason)) fits = 0;
        MPI_Allreduce(MPI_IN_PLACE, &fits, 1, MPI_INT, MPI_MIN, level_comm);
        if (blocks_ != nullptr && !fits) {
            throw std::runtime_error("start_owner_level: blocks kept on the device for a level the device declines");
        }
        if (fits) {
            level_ = make_level_eliminator(tree_, lvl, kernel_, opt_.tolerance, opt_.method, &reason,
                                           blocks_ != nullptr ? std::move(blocks_) : nullptr, opt_.occupancy, true,
                                           /*ca=*/true, /*owner=*/true, opt_.compress_factors);
            if (level_) level_->set_owner_plan(view);
        }
        level_index_ = lvl;
        decide_solve_keep(lvl, level_comm, false);  // (the factors kept here are this rank's boxes')
        if (wave_trace_enabled()) trace_.start_end();
        if (announce) {
            std::cout << "  GPU owner level: "
                      << (level_ ? std::string(tensor_core_gemm() && std::is_same_v<DataType, double> ? "on (FP64 tensor-core GEMMs)" : "on")
                                 : "off (" + reason + ")")
                      << std::endl;
        }
        return level_ != nullptr;
    }
    // After a colour's owned waves: its transport (ca_owner.hpp).
    void owner_exchange(int colour) {
        if (!level_) throw std::runtime_error("owner_exchange: no device level");
        level_->owner_exchange(colour);
    }
    bool ca_level_runs(int lvl) const { return ca_level_runs(lvl, opt_.ca_owner_component); }
    // The eliminator's CA mode supports level `lvl` (the same on every rank;
    // the box sizes are checked apart, collectively).
    bool ca_level_supported(int lvl, std::string* reason) const {
        const Level& level = tree_->levels[static_cast<size_t>(lvl)];
        // (the CA levels' lazy mode: at most 1)
        const int lazy = std::min(opt_.lazy_schur, 1);
        if (!level_eliminator_supported(level, kernel_, tree_->dimension, opt_.method, reason, true, lazy)) {
            return false;
        }
        if (!device_sketch_supported(tree_, *evaluator_of(kernel_->gpu_evaluator))) {
            if (reason) *reason = "CA levels need the device sketch";
            return false;
        }
        return true;
    }
    // Level start, after the host halo gather (all ranks of level_comm).
    // Returns whether the device runs the level: every box fits the device
    // sketch's packed lists on every rank, and the eliminator supports it.
    bool start_ca_level(int lvl, MPI_Comm level_comm, bool announce) {
        stage_solve_factors(lvl, announce, level_comm);
        if (wave_trace_enabled()) trace_.start_begin();
        const Level& level = tree_->levels[static_cast<size_t>(lvl)];
        int fits = 1;
        for (const Box& box : level.local_boxes) {
            if (box.num_points > kEntryDestMask) fits = 0;
        }
        for (const Box& box : level.ghost_boxes) {
            if (box.num_points > kEntryDestMask) fits = 0;
        }
        std::string reason = "boxes over " + std::to_string(kEntryDestMask) + " points";
        // decided on every rank before any builds its eliminator (which moves
        // the host blocks to the device)
        if (fits && !ca_level_supported(lvl, &reason)) fits = 0;
        MPI_Allreduce(MPI_IN_PLACE, &fits, 1, MPI_INT, MPI_MIN, level_comm);
        const bool device_halo = blocks_ != nullptr;  // the transition kept the blocks (M2)
        if (device_halo && !fits) {
            throw std::runtime_error("start_ca_level: blocks kept on the device for a level the device declines");
        }
        if (fits) {
            if (device_halo) {
                const CaHaloStats h = exchange_ca_ghost_blocks(tree_, lvl, *blocks_, level_comm);
                if (announce) {
                    std::printf("  [gpu] level %d device halo: %.2f s (plan %.2f, pack %.2f, MPI %.2f, unpack %.2f), "
                                "sent %.2f GB, received %.2f GB in %lld blocks%s\n",
                                lvl, h.total, h.plan, h.pack, h.mpi, h.unpack, h.bytes_sent / 1e9,
                                h.bytes_received / 1e9, static_cast<long long>(h.blocks_received),
                                h.staged ? " (through the host)" : "");
                    std::fflush(stdout);
                }
            }
            level_ = make_level_eliminator(tree_, lvl, kernel_, opt_.tolerance, opt_.method, &reason,
                                           device_halo ? std::move(blocks_) : nullptr, opt_.occupancy, true,
                                           /*ca=*/true, /*owner=*/false, opt_.compress_factors);
        }
        level_index_ = lvl;
        decide_solve_keep(lvl, level_comm, true);
        if (level_) {   // (every rank of the level: fits was decided together)
            int pack = level_->shared_packing_expected() ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &pack, 1, MPI_INT, MPI_MAX, level_comm);
            // H2_CA_pack_shared: auto (0) by the estimate above; always (1);
            // never (2) - and no rank packs them on its own when it runs
            // short either (that is how ghost copies came to differ), so a
            // card without the room for them raw runs out of heap instead
            if (opt_.ca_pack_shared == 1) pack = 1;
            if (opt_.ca_pack_shared == 2) pack = 0;
            level_->set_pack_shared(pack != 0, opt_.ca_pack_shared != 2);
            // H2_CA_solve home: the ghosts' factors are not kept (the solves run
            // the level from each box's home factors, in its elimination order)
            level_->set_home_solve(opt_.ca_solve == 1);
            if (announce && opt_.ca_solve == 1) {
                std::printf("  [gpu] level %d: solved from its home factors (H2_CA_solve home): the ghosts' factors "
                            "are not kept\n", lvl);
                std::fflush(stdout);
            }
            double est[2] = {level_->shared_packing_estimate(), -level_->shared_packing_room()};
            MPI_Allreduce(MPI_IN_PLACE, est, 2, MPI_DOUBLE, MPI_MAX, level_comm);
            if (announce) {
                static const char* const modes[3] = {"auto", "always", "never"};
                std::printf("  [gpu] level %d: the blocks on more than one rank %s through the ghosts' groups "
                            "(H2_CA_pack_shared %s; their near set raw: up to %.1f GB a rank against %.1f GB of room "
                            "under the hold limit on the tightest)\n", lvl,
                            pack != 0 ? "held packed" : (opt_.ca_pack_shared == 2 ? "never packed" : "left to each rank's pressure"),
                            modes[std::min(std::max(opt_.ca_pack_shared, 0), 2)], est[0] / 1e9, -est[1] / 1e9);
                std::fflush(stdout);
            }
        }
        if (wave_trace_enabled()) trace_.start_end();
        if (announce) {
            std::cout << "  GPU CA level: "
                      << (level_ ? std::string(tensor_core_gemm() && std::is_same_v<DataType, double> ? "on (FP64 tensor-core GEMMs)" : "on")
                                 : "off (" + reason + ")")
                      << std::endl;
        }
        return level_ != nullptr;
    }
    // Whether the level keeps its solve factors on the device (and skips
    // the host copies of their X_NR), before its waves and alike on every
    // rank of `level_comm`: their estimated bytes within the keep budget
    // (kSolveKeepFraction of the heap less what it holds now) on every rank.
    // The estimate is the upper bound (solve_keep_bound) times keep_ratio_:
    // the largest ratio of the factors' bytes to bound on the level below,
    // with a margin (kKeepMargin; the bound itself for the first level).  Too low an
    // estimate is caught by the waves, which keep factors only within the
    // budget: the level's keep then fails (its factors go to the host, as
    // without keeping) and later levels use the bound itself.  Collective
    // over level_comm.
    //
    // H2_compress_factors: X_NR is stored compressed, so the bound's X_NR
    // term is divided by kXnrBoundDivisor.  Without that the leaf never keeps
    // - its bound is tens of GB against a budget of what the heap has left
    // once the level's own working set is in it - and the leaf is where most
    // of X_NR lives, so its factors went to the host and came back
    // uncompressed at the first solve (GPU_CA_PLAN.md C4).  The divisor stays
    // far below the 18-19x the compression reaches, and a group of a wave
    // that does not fit is still caught by the waves' own budget check.
    void decide_solve_keep(int lvl, MPI_Comm level_comm, bool ca) {
        if (!color_gpu_enabled()) return;
        int ok = 1;
        double budget = 0.0;
        keep_bound_ = 0.0;
        if (level_) {
            const DeviceHeap& heap = DeviceHeap::instance();
            keep_bound_ = solve_keep_bound(tree_->levels[static_cast<size_t>(lvl)], ca,
                                           opt_.compress_factors >= 1 ? kXnrBoundDivisor : 1.0,
                                           !opt_.is_symmetric);
            budget = kSolveKeepFraction * static_cast<double>(heap.capacity()) - static_cast<double>(heap.used());
            if (compression_ <= 0.0) compression_ = opt_.compress_factors >= 1 ? kFirstCompression : 1.0;
            kept_expected_ = keep_ratio_ * keep_bound_ / compression_;   // what the store will hold
            // A nonsymmetric matrix's level also holds its fill sources raw
            // (temp2^T and temp4, nothing compressed yet), about the bound
            // itself (7.1 GB against a 9.2 GB bound at level 4 of 96^3 on 8
            // ranks): they must fit beside the kept factors.
            const double working = opt_.is_symmetric ? 0.0 : keep_bound_;
            ok = kept_expected_ + working <= budget ? 1 : 0;
        }
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, level_comm);
        if (level_) level_->set_solve_keep(ok != 0, budget, kept_expected_);
    }
    // After the host's post-elimination assisting gather: the assisting
    // boxes' skeletons, for the device transition.
    void refresh_remote_skeletons() {
        if (level_) level_->refresh_remote_skeletons();
    }
    // BPACK_CHECK=replica, after finish_level of a device CA level (all
    // ranks of level_comm): every ghost copy of a box against its owner's,
    // bitwise, by hashes of its factors.  Prints the copies checked and, per
    // rank, the first copies that differ and in which factors.
    void check_ca_replicas(int lvl, MPI_Comm level_comm, bool announce) const {
        if (!ca_replica_check_enabled() || opt_.ca_owner_component != 0) return;  // (an owner level has no replicas)
        const Level& level = tree_->levels[static_cast<size_t>(lvl)];
        constexpr int kParts = 8;
        static const char* const kPartNames[kParts] = {"skeleton", "T", "LU(X_RR)", "pivots",
                                                       "X_SR", "X_RS", "X_NR", "X_RR_full"};
        using Hashes = std::array<uint64_t, kParts>;
        auto hash_bytes = [](const void* p, size_t bytes) {
            uint64_t h = 0x9e3779b97f4a7c15ull ^ bytes;
            const unsigned char* c = static_cast<const unsigned char*>(p);
            size_t i = 0;
            for (; i + 8 <= bytes; i += 8) {
                uint64_t w;
                std::memcpy(&w, c + i, 8);
                h = (h ^ w) * 0xff51afd7ed558ccdull;
                h ^= h >> 32;
            }
            for (; i < bytes; ++i) h = (h ^ c[i]) * 0x100000001b3ull;
            return h;
        };
        auto hashes_of = [&](const Box& box) {
            auto vec = [&](const auto& v) { return hash_bytes(v.data(), v.size() * sizeof(v[0])); };
            return Hashes{vec(box.skeleton_indices), vec(box.interpolation_matrix.data), vec(box.X_RR.data),
                          vec(box.X_RR_pivots), vec(box.X_SR.data), vec(box.X_RS_entry.data),
                          vec(box.X_NR.data), vec(box.X_RR_full.data)};
        };
        auto eliminated = [&](const Box& box) {
            return level.eliminated_boxes.count(box.morton_index) != 0 && box.num_points > 0;
        };
        // the owners' hashes, gathered: (morton, hashes) per local box
        std::vector<uint64_t> mine;
        {
            const size_t nb = level.local_boxes.size();
            std::vector<Hashes> h(nb);
            #pragma omp parallel for schedule(dynamic, 8)
            for (int64_t b = 0; b < static_cast<int64_t>(nb); ++b) h[static_cast<size_t>(b)] = hashes_of(level.local_boxes[static_cast<size_t>(b)]);
            for (size_t b = 0; b < nb; ++b) {
                if (!eliminated(level.local_boxes[b])) continue;
                mine.push_back(static_cast<uint64_t>(level.local_boxes[b].morton_index));
                mine.insert(mine.end(), h[b].begin(), h[b].end());
                mine.push_back(static_cast<uint64_t>(level.local_boxes[b].X_NR.rows));
            }
        }
        int nranks = 1;
        MPI_Comm_size(level_comm, &nranks);
        int count = static_cast<int>(mine.size());
        std::vector<int> counts(static_cast<size_t>(nranks)), displs(static_cast<size_t>(nranks), 0);
        MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, level_comm);
        for (int r = 1; r < nranks; ++r) displs[static_cast<size_t>(r)] = displs[static_cast<size_t>(r - 1)] + counts[static_cast<size_t>(r - 1)];
        std::vector<uint64_t> all(static_cast<size_t>(displs.back() + counts.back()));
        MPI_Allgatherv(mine.data(), count, MPI_UINT64_T, all.data(), counts.data(), displs.data(), MPI_UINT64_T,
                       level_comm);
        std::unordered_map<int64_t, const uint64_t*> owner;
        for (size_t i = 0; i < all.size(); i += kParts + 2) owner[static_cast<int64_t>(all[i])] = &all[i + 1];
        // this rank's ghost copies against them
        const size_t ng = level.ghost_boxes.size();
        std::vector<Hashes> gh(ng);
        #pragma omp parallel for schedule(dynamic, 8)
        for (int64_t g = 0; g < static_cast<int64_t>(ng); ++g) gh[static_cast<size_t>(g)] = hashes_of(level.ghost_boxes[static_cast<size_t>(g)]);
        long long totals[2] = {0, 0};  // copies checked, copies that differ
        int rank = 0;
        MPI_Comm_rank(level_comm, &rank);
        for (size_t g = 0; g < ng; ++g) {
            const Box& box = level.ghost_boxes[g];
            if (!eliminated(box)) continue;
            auto it = owner.find(box.morton_index);
            if (it == owner.end()) {
                throw std::runtime_error("CA replica check: ghost " + std::to_string(box.morton_index) +
                                         " has no owner copy");
            }
            ++totals[0];
            std::string parts;
            for (int p = 0; p < kParts; ++p) {
                if (gh[g][static_cast<size_t>(p)] != it->second[p]) parts += std::string(parts.empty() ? "" : ", ") + kPartNames[p];
            }
            if (parts.empty()) continue;
            if (++totals[1] <= 3) {
                auto w = level.elimination_wave.find(box.morton_index);
                std::printf("  [gpu] level %d CA replica check, rank %d: ghost %lld (wave %d) differs from its owner's "
                            "copy in %s (X_NR rows %lld here, %lld on the owner)\n",
                            lvl, rank, static_cast<long long>(box.morton_index),
                            w == level.elimination_wave.end() ? -1 : static_cast<int>(w->second), parts.c_str(),
                            static_cast<long long>(box.X_NR.rows), static_cast<long long>(it->second[kParts]));
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, totals, 2, MPI_LONG_LONG, MPI_SUM, level_comm);
        if (announce) {
            std::printf("  [gpu] level %d CA replica check: %lld ghost copies, %lld differ from their owners' "
                        "(bitwise)\n",
                        lvl, totals[0], totals[1]);
        }
        std::fflush(stdout);
    }

    bool on() const { return level_ != nullptr; }
    // one active rank: the level has nothing to transport
    bool local() const {
        return level_ != nullptr && tree_->levels[static_cast<size_t>(level_index_)].num_active_processes == 1;
    }

    // A transport of a device level on several ranks: the device exchange,
    // or the generators of the last wave, the host transport
    // (host_transport(&installed) installs the remote generators and returns
    // its MPI time), and the device updates from the installed generators.
    // (Entries through a collective callback: the level's ranks quiesce their
    // rounds before the exchange and after the device updates, which fetch.)
    template<typename HostTransport>
    Duration transport(PendingFactorUpdates<DataType>& pending, HostTransport&& host_transport) {
        quiesce_entries();
        if (level_->device_exchange()) {
            const Duration duration = level_->exchange();
            quiesce_entries();
            return duration;
        }
        const auto t_emit = std::chrono::steady_clock::now();
        level_->emit_generators(pending);
        const auto t_exchange = std::chrono::steady_clock::now();
        transport_timers().presend += std::chrono::duration<double>(t_exchange - t_emit).count();
        std::vector<int64_t> installed;
        const Duration duration = host_transport(&installed);
        eliminator_stats().emit += std::chrono::duration<double>(t_exchange - t_emit).count();
        eliminator_stats().exchange +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_exchange).count();
        level_->receive_remote(installed);
        quiesce_entries();
        return duration;
    }

    // Box region, owner pass, mirror and finalize of a wave on the device;
    // returns its number of boundary boxes.  The loop records the wave's
    // boxes as eliminated.
    // A replicated CA level's interior group starts: no rank eliminates a
    // ghost from here on (see LevelEliminatorBase::replicas_done).
    void replicas_done() {
        if (level_) level_->replicas_done();
    }
    int eliminate_wave(const std::vector<int64_t>& wave, int counter) {
        if (!wave_trace_enabled()) return level_->eliminate_wave(wave, counter);
        trace_.wave_begin();
        const int boundary = level_->eliminate_wave(wave, counter);
        size_t ghosts = 0;   // (a replicated CA level's: boxes of other ranks eliminated here too)
        if (level_index_ >= 0) {
            const auto& L = tree_->levels[static_cast<size_t>(level_index_)];
            const int64_t lo = L.local_morton_start, hi = lo + static_cast<int64_t>(L.local_boxes.size());
            for (int64_t m : wave) ghosts += (m < lo || m >= hi) ? 1 : 0;
        }
        trace_.wave_end(counter, wave.size(), ghosts);
        return boundary;
    }

    // Level end, after the final transport (all ranks of level_comm): the
    // factors on the host, the blocks kept for the transition, and whether
    // the solve's factors stay on the device (only if every rank of the
    // level still holds all of its own).
    void finish_level(int lvl, const PendingFactorUpdates<DataType>& pending, MPI_Comm level_comm, bool announce) {
        if (level_) quiesce_entries();  // (before the collectives below)
        const auto t_finish = std::chrono::steady_clock::now();
        double finish_seconds = 0.0;  // this rank's level_->finish()
        if (level_) {
            // (every rank of level_comm has a device level here)
            const bool trace = wave_trace_enabled();
            if (trace) trace_.finish_begin();
            level_->finish();
            finish_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_finish).count();
            if (trace) {
                trace_.finish_end();
                print_wave_trace(lvl, level_comm, announce);
            }
            if (const Level& level = tree_->levels[static_cast<size_t>(lvl)]; lvl <= 3 && level.is_process_active && wave_trace_enabled()) {   // (8x diagnostics, BPACK_TRACE=wave) every rank's boxes at the small levels
                int rank = 0;
                MPI_Comm_rank(tree_->comm, &rank);
                std::string boxes;
                int64_t maxn = 0;
                int full = 0;
                for (const Box& b : level.local_boxes) {
                    maxn = std::max<int64_t>(maxn, b.num_points);
                    if (b.redundant_indices.empty() && b.num_points > 0) ++full;
                    if (boxes.size() < 400) {
                        boxes += (boxes.empty() ? "" : " ") + std::to_string(b.num_points) + "/" +
                                 std::to_string(b.skeleton_indices.size());
                    }
                }
                std::printf("  [gpu] level %d rank %d boxes (points/skeleton): %s; largest %lld points, %d kept whole\n", lvl, rank,
                            boxes.c_str(), static_cast<long long>(maxn), full);
                std::fflush(stdout);
            }
            if (local() && (!pending.replace_blocks.empty() || !pending.accumulated_deltas.empty() ||
                            !pending.generators.empty())) {
                throw std::runtime_error("device level left updates for other ranks");
            }
        }
        if (color_gpu_enabled()) {
            int kept = level_ && level_->solve_factors_kept() ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &kept, 1, MPI_INT, MPI_MIN, level_comm);
            // the next level's estimate (decide_solve_keep): the largest
            // ratio of the factors' bytes to bound, kept or not; the bound
            // itself after a failed keep
            const double bytes = level_ ? level_->level_factor_bytes() : 0.0;
            const int failure = level_ ? level_->solve_keep_failure() : 0;
            // (before the commit: it moves the kept state into the store)
            const double kept_now = level_ ? level_->level_kept_bytes() : 0.0;
            const double raw_now = bytes;
            // v[4]: minus the compression this rank measured, so the maximum
            // is the least compression of any rank (the next level's estimate)
            double v[5] = {keep_bound_ > 0.0 ? bytes / keep_bound_ : 0.0, bytes, failure == 1 ? 1.0 : 0.0,
                           failure == 2 ? 1.0 : 0.0, kept_now > 0.0 && raw_now > 0.0 ? -(raw_now / kept_now) : 0.0};
            MPI_Allreduce(MPI_IN_PLACE, v, 5, MPI_DOUBLE, MPI_MAX, level_comm);
            if (v[2] > 0.0 || v[3] > 0.0) {
                keep_ratio_ = 1.0;
                compression_ = 1.0;   // (the estimate was too small: the bound itself from here on)
            } else if (v[0] > 0.0) {
                keep_ratio_ = std::min(1.0, kKeepMargin * v[0]);
            }
            if (v[4] < 0.0 && !(v[2] > 0.0 || v[3] > 0.0)) compression_ = std::max(1.0, -v[4]);
            if (level_) level_->commit_solve_factors(kept != 0);
            // (the level's end: this rank's finish and the slowest rank's, the
            // collectives' wait included in the total)
            double fin[2] = {finish_seconds, finish_seconds};
            MPI_Allreduce(MPI_IN_PLACE, &fin[1], 1, MPI_DOUBLE, MPI_MAX, level_comm);
            if (announce) {
                std::printf("  [gpu] level %d end: finish %.2f s (slowest rank %.2f), total %.2f s\n", lvl, fin[0], fin[1],
                            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_finish).count());
                if (kept != 0) {
                    // stored vs the uncompressed size, so H2_compress_factors
                    // shows up (level_factor_bytes() is from the dimensions)
                    const double ratio = kept_now > 0.0 ? raw_now / kept_now : 1.0;
                    std::printf("  [gpu] level %d solve factors: kept on the device (up to %.2f GB per rank, at most "
                                "%.0f%% of their bound; rank 0 stored %.2f GB of %.2f GB, %.1fx)\n",
                                lvl, v[1] / 1e9, 100.0 * v[0], kept_now / 1e9, raw_now / 1e9, ratio);
                    const auto& es = eliminator_stats();
                    if (es.xnr_waves_compressed + es.xnr_waves_plain + es.xnr_not_requested > 0) {
                        if (es.xnr_raw_bytes > 0.0) {
                            std::printf("  [gpu] level %d X_NR is %.0f%% of the kept block; it compresses %.1fx "
                                        "(%.2f -> %.2f GB), T and X_SR %.1fx (%.2f -> %.2f GB), the LU and "
                                        "pivots are %.2f GB uncompressed\n", lvl,
                                        100.0 * es.xnr_raw_bytes / (es.xnr_raw_bytes + es.other_raw_bytes),
                                        es.xnr_stored_bytes > 0.0 ? es.xnr_raw_bytes / es.xnr_stored_bytes : 1.0,
                                        es.xnr_raw_bytes / 1e9, es.xnr_stored_bytes / 1e9,
                                        es.tx_stored_bytes > 0.0 ? es.tx_raw_bytes / es.tx_stored_bytes : 1.0,
                                        es.tx_raw_bytes / 1e9, es.tx_stored_bytes / 1e9,
                                        (es.other_raw_bytes - es.tx_raw_bytes) / 1e9);
                        }
                        std::printf("  [gpu] level %d X_NR compression: %lld batches compressed, %lld fell "
                                    "back, %lld waves not requested\n", lvl,
                                    static_cast<long long>(es.xnr_waves_compressed),
                                    static_cast<long long>(es.xnr_waves_plain),
                                    static_cast<long long>(es.xnr_not_requested));
                    }
                } else {
                    std::printf("  [gpu] level %d solve factors: to be uploaded at the first solve%s%s\n", lvl,
                                v[2] > 0.0 ? " (the keep budget ran out during the level)" : "",
                                v[3] > 0.0 ? " (no free range of the heap held a wave's factors)" : "");
                }
                std::fflush(stdout);
            }
        }
    }

    // BPACK_TRACE=wave (all ranks of level_comm): the spread of the
    // ranks' wave totals, and the wave tables of rank 0 and of the rank whose
    // waves took longest (a replicated CA level: the rank with most ghosts).
    void print_wave_trace(int lvl, MPI_Comm level_comm, bool announce) const {
        int rank = 0, size = 1;
        MPI_Comm_rank(level_comm, &rank);
        MPI_Comm_size(level_comm, &size);
        const auto mine = trace_.summary();
        constexpr int K = WaveTrace::kSummary;
        std::vector<double> all(static_cast<size_t>(size) * K);
        MPI_Gather(mine.data(), K, MPI_DOUBLE, all.data(), K, MPI_DOUBLE, 0, level_comm);
        struct { double v; int r; } local{mine[0], rank}, slowest{0.0, 0};
        MPI_Allreduce(&local, &slowest, 1, MPI_DOUBLE_INT, MPI_MAXLOC, level_comm);
        if (announce) {
            static const char* const names[K] = {"waves", "boxes", "sketch wait", "ID store", "plan",
                                                 "device sketch", "device elimination"};
            std::printf("  [gpu] level %d wave trace over %d ranks (min / median / max; rank of the max):", lvl, size);
            for (int k = 0; k < K; ++k) {
                std::vector<std::pair<double, int>> v;
                for (int r = 0; r < size; ++r) v.emplace_back(all[static_cast<size_t>(r) * K + k], r);
                std::sort(v.begin(), v.end());
                const double f = k == 1 ? 1.0 : 1e3;  // boxes, or ms
                std::printf(" %s %.0f / %.0f / %.0f (%d)%s", names[k], f * v.front().first,
                            f * v[v.size() / 2].first, f * v.back().first, v.back().second, k + 1 < K ? ";" : "\n");
            }
            trace_.print(lvl, rank);
        }
        MPI_Barrier(level_comm);
        if (rank == slowest.r && !announce) trace_.print(lvl, rank);
        MPI_Barrier(level_comm);
        trace_.print_rows(lvl, rank);   // (every rank, a line a wave)
    }

    // Level 1 (active, not eliminated) adopts the blocks of the level-2
    // transition; its own transition builds the root on the device.
    void adopt_level(int lvl) {
        if (!blocks_) return;
        std::string reason;
        level_ = make_level_eliminator(tree_, lvl, kernel_, opt_.tolerance, opt_.method, &reason,
                                       std::move(blocks_), opt_.occupancy);
        level_index_ = lvl;
        level_->adopt_without_elimination();
    }

    // ---- the level transition

    // This rank builds the parents of its boxes on the device (transition).
    bool builds_parents(const Level& level, const Level& parent_level) const {
        const bool reduction_ahead = parent_level.num_active_processes != level.num_active_processes;
        return level_ && level.is_process_active && level_->can_build_parent() &&
               (reduction_ahead || parent_level.is_process_active) && !opt_.is_hermitian;
    }

    // Device transition: this level's blocks never leave the device, and the
    // parent's go to the next level's eliminator (or to the host when that
    // level runs there).  Each rank builds the parents of its own boxes, with
    // its copies of the blocks shared with other ranks; at a process
    // reduction they move to their new owner's device when the parent level
    // runs there (device_reduction.hpp; every child rank decides alike over
    // `level_comm`), else to the host for the host reduction.
    // `build_structure()` returns the parent boxes without their blocks.
    // Returns whether the parents were built here; otherwise the level's
    // blocks are back in the host BoxData for the host transition.
    // `level_comm` (the level's active ranks, or MPI_COMM_NULL): a device CA
    // level next keeps the blocks on the device (its halo then moves between
    // the devices, M2) when every rank's parents fit the device sketch.
    // `transition_comm` (the ranks of the transition: the level's active
    // ranks and the parent level's, whose owners at a reduction may be ranks
    // idle at this level): the device reduction runs over it, the idle
    // owners receiving only.
    // `section_end` (every rank of the call, after the parents are built and
    // before any communication of the transition): entries through a
    // collective callback end the level's section of rounds there, so the
    // parent-only ranks serving them come to the transition (a rank that
    // builds its parents on the host does so before the call: builds_parents).
    template<typename StructureBuilder>
    bool transition(const Level& level, const Level& parent_level, std::vector<Box>& parents,
                    StructureBuilder&& build_structure, bool announce, MPI_Comm level_comm = MPI_COMM_NULL,
                    MPI_Comm transition_comm = MPI_COMM_NULL, const std::function<void()>& section_end = {}) {
        bool on_device = false;
        const bool reduction_ahead = parent_level.num_active_processes != level.num_active_processes;
        const bool builds = builds_parents(level, parent_level);
        if (builds) {
            parents = build_structure();
            {   // (8x diagnostics) a parent wider than the device sketch takes, with its children's sizes
                int rank = 0;
                MPI_Comm_rank(tree_->comm, &rank);
                for (const Box& p : parents) {
                    if (p.num_points <= kEntryDestMask) continue;
                    std::string kids;
                    for (int c = 0; c < 8; ++c) {
                        const int64_t cm = p.children_morton[c];
                        if (cm < 0) continue;
                        const int64_t ci = cm - level.local_morton_start;
                        const Box* cb = ci >= 0 && ci < static_cast<int64_t>(level.local_boxes.size()) ? &level.local_boxes[static_cast<size_t>(ci)] : nullptr;
                        kids += (kids.empty() ? "" : ", ") + std::to_string(cm) + ":" +
                                (cb ? std::to_string(cb->skeleton_indices.size()) + "/" + std::to_string(cb->num_points) : "?");
                    }
                    std::printf("  [gpu] level %d rank %d: parent %lld has %lld points (children skeleton/points: %s)\n",
                                level.level, rank, static_cast<long long>(p.morton_index),
                                static_cast<long long>(p.num_points), kids.c_str());
                }
                std::fflush(stdout);
            }
            // The parent blocks may be stored compressed only where the next
            // level's eliminator adopts them from this device and reads them
            // through begin()'s restore.  A process reduction sends them to
            // another rank and download_level_blocks() copies them to the
            // host, and both read the raw blocks: packing there left the
            // reduction with nothing to send ("adopted Schur block of box 1
            // is missing" at the level that follows).  `keep` below can only
            // be promoted from false to true (a CA parent level), so missing
            // that promotion costs a level's packing, not a result.
            level_->set_pack_parent_blocks(!reduction_ahead && keeps_blocks_of(level.level - 1));
            blocks_ = level_->build_parent(parents);
            level_.reset();
            on_device = true;
            // H2 8x (2026-10-11): a component-owner level's owned ghosts keep
            // host copies of their X_NR and X_RR_full only as fill sources the
            // transition may restore; their factors went home with the FULL
            // payload, and the solves run the level from its home boxes.  At
            // the 384^3 leaf on 8 ranks they were 6.2 GB a rank, held under
            // every later level, and the node's room for the parents' packed
            // copies ran out.
            if (Level& lv = tree_->levels[static_cast<size_t>(level.level)]; lv.ca_home_factors) {
                for (Box& g : lv.ghost_boxes) {
                    g.X_NR = MatrixStorage<DataType>{};
                    g.X_RR_full = MatrixStorage<DataType>{};
                }
            }
        }
        if (section_end) section_end();
        const int parent_lvl = level.level - 1;
        bool keep = false;
        bool reduced = false;
        ReductionStats red;
        if (reduction_ahead && transition_comm != MPI_COMM_NULL) {
            // every rank of the transition can take part: the builders with
            // their blocks, the owners as receivers
            const bool mine = level.is_process_active ? (builds && blocks_ != nullptr) : parent_level.is_process_active;
            int ok = mine && device_reduction_possible(level.level, parent_lvl, parents) ? 1 : 0;
            {
                // H2 8x (2026-10-11): the receivers' room.  A receiver holds the
                // merged parent blocks raw; when what its senders hold would not
                // fit beside what its card has, every rank takes the host
                // reduction instead (download_level_blocks below).  At 768^3 on
                // 64 ranks the 64 -> 8 reduction into level 2 exhausted the
                // receivers' cards.
                int rank = 0, tsize = 0, trank = 0;
                MPI_Comm_rank(tree_->comm, &rank);
                MPI_Comm_size(transition_comm, &tsize);
                MPI_Comm_rank(transition_comm, &trank);
                const int owner = level.parent_level_owner;
                const bool sender = level.is_process_active && owner != rank && blocks_ != nullptr;
                const bool receiver = parent_level.is_process_active;
                double my_bytes = 0.0;
                int my_blocks = 0, largest_r = 0, largest_c = 0;
                if (blocks_ != nullptr) {
                    auto take = [&](int rows, int cols) {
                        const double b = static_cast<double>(rows) * static_cast<double>(cols) * sizeof(DataType);
                        my_bytes += b;
                        ++my_blocks;
                        if (static_cast<double>(largest_r) * largest_c < static_cast<double>(rows) * cols) { largest_r = rows; largest_c = cols; }
                    };
                    for (const auto& kv : blocks_->schur) take(kv.second.rows, kv.second.cols);
                    for (const auto& kv : blocks_->edges) take(kv.second.rows, kv.second.cols);
                }
                std::vector<double> to_owner(static_cast<size_t>(tsize), 0.0);
                if (sender) {
                    MPI_Group wg, tg;
                    MPI_Comm_group(tree_->comm, &wg);
                    MPI_Comm_group(transition_comm, &tg);
                    int owner_t = MPI_UNDEFINED;
                    MPI_Group_translate_ranks(wg, 1, &owner, tg, &owner_t);
                    MPI_Group_free(&wg);
                    MPI_Group_free(&tg);
                    if (owner_t != MPI_UNDEFINED) to_owner[static_cast<size_t>(owner_t)] += my_bytes;
                }
                MPI_Allreduce(MPI_IN_PLACE, to_owner.data(), tsize, MPI_DOUBLE, MPI_SUM, transition_comm);
                int fits = 1;
                const DeviceHeap& heap = DeviceHeap::instance();
                const double incoming = to_owner[static_cast<size_t>(trank)];
                if (receiver && heap.initialized()) {
                    const double room = static_cast<double>(heap.capacity()) - static_cast<double>(heap.used());
                    const double margin = std::max(2.0e9, 0.1 * static_cast<double>(heap.capacity()));
                    if (incoming + margin > room) fits = 0;
                }
                if (sender && my_bytes > 1.0e9) {
                    std::printf("  [gpu] level %d rank %d reduction: sends %.2f GB in %d blocks to rank %d, largest %d x %d\n",
                                level.level, rank, my_bytes / 1e9, my_blocks, owner, largest_r, largest_c);
                }
                if (receiver) {
                    std::printf("  [gpu] level %d rank %d reduction: receives %.2f GB beside its own %.2f GB; card %.2f of %.2f GB in use; %s\n",
                                level.level, rank, incoming / 1e9, my_bytes / 1e9,
                                heap.initialized() ? static_cast<double>(heap.used()) / 1e9 : 0.0,
                                heap.initialized() ? static_cast<double>(heap.capacity()) / 1e9 : 0.0,
                                fits ? "on the device" : "through the host (no room)");
                }
                std::fflush(stdout);
                if (!fits) ok = 0;
            }
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, transition_comm);
            if (ok != 0) {
                red = reduce_parent_blocks(tree_, level, parent_level, blocks_, transition_comm);
                eliminator_stats().reduction += red.total;
                reduced = true;
                keep = true;
            }
        } else if (builds && !reduction_ahead) {
            keep = keeps_blocks_of(parent_lvl);
            if (!keep && level_comm != MPI_COMM_NULL && ca_level_runs(parent_lvl) && ca_level_supported(parent_lvl, nullptr)) {
                int fits = 1;
                for (const Box& p : parents) {
                    if (p.num_points > kEntryDestMask) fits = 0;
                }
                MPI_Allreduce(MPI_IN_PLACE, &fits, 1, MPI_INT, MPI_MIN, level_comm);
                keep = fits != 0;
            }
        }
        if (builds && !keep) {
            download_level_blocks(*blocks_, parents);
            blocks_.reset();
        }
        if (builds && announce) {
            const auto& e = eliminator_stats();
            std::printf("  [gpu] level %d device transition: %.2f s (%lld chunks: plan %.2f [restore %.2f, %.2f GB], "
                        "blocks %.2f, P %.2f, fill %.2f), %lld fill GEMMs (%.1f GFLOP), blocks %s\n",
                        level.level, e.transition, static_cast<long long>(e.tr_chunks), e.tr_plan,
                        e.tr_restore, e.tr_restore_bytes / 1e9, e.tr_blocks, e.tr_p, e.tr_fill,
                        static_cast<long long>(e.transition_fill_gemms), e.transition_flops / 1e9,
                        reduced ? (blocks_ ? "merged on this device (process reduction)"
                                           : "sent to their owner's device (process reduction)")
                                : (blocks_ ? "kept on the device" : "copied to the host"));
            if (e.ghost_lu_bytes > 0) {
                std::printf("  [gpu] level %d ghosts' LU held to the level's end (H2_CA_solve home): %.2f GB\n",
                            level.level, e.ghost_lu_bytes / 1e9);
            }
            if (e.rc_sources > 0 || e.tr_xrr_host > 0 || e.xr_restored > 0) {
                std::printf("  [gpu] level %d fill sources without host copies: %lld freed during the waves, X_RR_full of "
                            "%lld restored from the host at the transition, %lld rebuilt from the LU at a keep's failure\n",
                            level.level, static_cast<long long>(e.rc_sources), static_cast<long long>(e.tr_xrr_host),
                            static_cast<long long>(e.xr_restored));
            }
            if (e.nl_blocks > 0) {
                const double tot = e.nl_bin[0] + e.nl_bin[1] + e.nl_bin[2] + e.nl_bin[3] + e.nl_bin[4];
                std::printf("  [gpu] level %d compressed length over the level, %lld blocks: raw %.2f GB, "
                            "first %.2f, max %.2f, last %.2f GB; sum max/first %.3f, worst block %.2f\n",
                            level.level, static_cast<long long>(e.nl_blocks), e.nl_raw / 1e9,
                            e.nl_first / 1e9, e.nl_max / 1e9, e.nl_last / 1e9,
                            e.nl_first > 0.0 ? e.nl_max / e.nl_first : 0.0, e.nl_worst);
                std::printf("  [gpu] level %d ... by max/first: %.0f%% never grow, %.0f%% <=1.1, %.0f%% <=1.25, "
                            "%.0f%% <=1.5, %.0f%% more\n", level.level,
                            100.0 * e.nl_bin[0] / tot, 100.0 * e.nl_bin[1] / tot, 100.0 * e.nl_bin[2] / tot,
                            100.0 * e.nl_bin[3] / tot, 100.0 * e.nl_bin[4] / tot);
            }
            {
                DeviceHeap& hp = DeviceHeap::instance();
                static size_t calls0 = 0, scans0 = 0;
                std::printf("  [gpu] level %d heap: %lld free ranges now; %lld transient takes scanned %lld ranges from the top\n",
                            level.level, static_cast<long long>(hp.free_ranges()),
                            static_cast<long long>(hp.take_calls() - calls0), static_cast<long long>(hp.take_scans() - scans0));
                calls0 = hp.take_calls();
                scans0 = hp.take_scans();
            }
            if (e.np_pressure_groups > 0 || e.tr_pack_chunks_packed > 0 || e.np_waves_chunked > 0) pressure_seen_ = true;
            if (e.np_pressure_groups + e.np_raw_groups + e.np_waves > 0) {
                std::printf("  [gpu] level %d memory pressure (hold limit %.2f GB): %lld owner groups packed, %lld left raw; "
                            "%lld of %lld waves chunked; parents packed in %lld of %lld transition chunks; %lld fill sources "
                            "copied to the host on demand\n", level.level, g_hold_limit_bytes / 1e9,
                            static_cast<long long>(e.np_pressure_groups), static_cast<long long>(e.np_raw_groups),
                            static_cast<long long>(e.np_waves_chunked), static_cast<long long>(e.np_waves),
                            static_cast<long long>(e.tr_pack_chunks_packed), static_cast<long long>(e.tr_pack_chunks),
                            static_cast<long long>(e.fs_staged_late));
                if (e.np_shared_blocks > 0) {
                    std::printf("  [gpu] level %d near blocks on more than one rank, held packed always: %lld "
                                "writes, %.2f GB raw (replicated CA)\n", level.level,
                                static_cast<long long>(e.np_shared_blocks), e.np_shared_raw / 1e9);
                }
                if (e.fw_no_window > 0) {
                    std::printf("  [gpu] level %d gather window: %lld chunks got none, %lld blocks put back one by one\n",
                                level.level, static_cast<long long>(e.fw_no_window), static_cast<long long>(e.fw_fallbacks));
                }
            }
            {
                double rss = 0.0, hwm = 0.0;
                host_memory_gb(rss, hwm);
                std::printf("  [gpu] level %d host memory: RSS %.2f GB, high-water %.2f GB\n", level.level, rss, hwm);
                if (wave_trace_enabled()) report_host_boxes(level.level);
            }
            if (e.dead_pair_tasks + e.dead_schur_tasks + e.dead_beyond_tasks > 0) {
                std::printf("  [gpu] level %d owner pass, updates read nowhere on this rank and not made: pairs of two "
                            "eliminated ghosts %.0f GFLOP (%lld), eliminated ghosts' Schur blocks %.0f GFLOP (%lld), their "
                            "pairs with boxes beyond the halo %.0f GFLOP (%lld); made: %.0f GFLOP\n", level.level,
                            e.dead_pair_flops / 1e9, static_cast<long long>(e.dead_pair_tasks), e.dead_schur_flops / 1e9,
                            static_cast<long long>(e.dead_schur_tasks), e.dead_beyond_flops / 1e9,
                            static_cast<long long>(e.dead_beyond_tasks), e.owner_flops / 1e9);
            }
            if (e.og_mat + e.og_win + e.og_stage + e.og_launch + e.og_pack > 0.0) {
                std::printf("  [gpu] level %d owner pass host: materialize %.2f, window %.2f, stage %.2f, launch+sync %.2f, "
                            "pack %.2f s\n", level.level, e.og_mat, e.og_win, e.og_stage, e.og_launch, e.og_pack);
                double gemm_s = 0.0, gemm_flops = 0.0;
                for (int j = 0; j < EliminatorStats::kOwnerDepthBins; ++j) {
                    gemm_s += e.og_dev_gemm[j];
                    gemm_flops += e.og_gemm_flops[j];
                }
                std::printf("  [gpu] level %d owner pass device: kernel values %.2f s, GEMMs %.2f s (%.0f GF/s); by "
                            "source position j [launches, GEMMs, GFLOP, s, GF/s, (mean m n)^1/2, mean k]:", level.level, e.og_dev_eval,
                            gemm_s, gemm_s > 0.0 ? gemm_flops / gemm_s / 1e9 : 0.0);
                for (int j = 0; j < EliminatorStats::kOwnerDepthBins; ++j) {
                    if (e.og_gemm_launches[j] == 0) continue;
                    std::printf(" %d%s [%lld, %lld, %.0f, %.2f, %.0f, %.0f, %.0f]", j + 1,
                                j + 1 == EliminatorStats::kOwnerDepthBins ? "+" : "",
                                static_cast<long long>(e.og_gemm_launches[j]), static_cast<long long>(e.og_gemm_count[j]),
                                e.og_gemm_flops[j] / 1e9, e.og_dev_gemm[j],
                                e.og_dev_gemm[j] > 0.0 ? e.og_gemm_flops[j] / e.og_dev_gemm[j] / 1e9 : 0.0,
                                std::sqrt(e.og_gemm_mn[j] / static_cast<double>(std::max<int64_t>(e.og_gemm_count[j], 1))),
                                e.og_gemm_mn[j] > 0.0 ? e.og_gemm_flops[j] / (2.0 * e.og_gemm_mn[j]) : 0.0);
                }
                std::printf("; targets by sources:");
                for (int d = 1; d < EliminatorStats::kOwnerDepthBins; ++d) {
                    if (e.og_targets_by_depth[d] > 0) {
                        std::printf(" %d%s %lld", d, d + 1 == EliminatorStats::kOwnerDepthBins ? "+" : "",
                                    static_cast<long long>(e.og_targets_by_depth[d]));
                    }
                }
                std::printf("\n");
            }
            if (e.np_gen_blocks > 0 || e.np_gen_fresh_blocks > 0) {
                std::printf("  [gpu] level %d packed after the generator passes: %lld adopted blocks (%.2f GB raw), "
                            "%lld targets made from the kernel (%.2f GB raw)\n", level.level,
                            static_cast<long long>(e.np_gen_blocks), e.np_gen_raw / 1e9,
                            static_cast<long long>(e.np_gen_fresh_blocks), e.np_gen_fresh_raw / 1e9);
            }
            if (e.np_groups > 0 || e.np_groups_plain > 0) {
                std::printf("  [gpu] level %d near blocks held packed: %.2f -> %.2f GB (%.1fx) over %lld blocks "
                            "in %lld groups, %.2f s%s\n", level.level, e.np_raw / 1e9, e.np_stored / 1e9,
                            e.np_stored > 0.0 ? e.np_raw / e.np_stored : 0.0,
                            static_cast<long long>(e.np_blocks), static_cast<long long>(e.np_groups),
                            e.np_seconds,
                            e.np_groups_plain > 0 ? " (some groups stayed raw)" : "");
                std::printf("  [gpu] level %d near store: %lld chunks grew into a tail, %lld blocks found no "
                            "room and stayed uncompressed\n", level.level,
                            static_cast<long long>(e.np_tails), static_cast<long long>(e.np_no_room));
                std::printf("  [gpu] level %d packed groups reclaimed: %lld\n", level.level,
                            static_cast<long long>(e.np_reclaimed));
            }
            if (e.nb_comp_bytes > 0.0 && e.nb_comp_seconds > 0.0) {
                std::printf("  [gpu] level %d near-block compression: %.0f GB/s over %.0f MB (%.0f GB/s at the "
                            "tighter tolerance), against 262-365 GB/s decompressing\n", level.level,
                            e.nb_comp_bytes / e.nb_comp_seconds / 1e9, e.nb_comp_bytes / 1e6,
                            e.nb_tight_seconds > 0.0 ? e.nb_tight_bytes / e.nb_tight_seconds / 1e9 : 0.0);
            }
            if (e.nb_raw > 0.0 || e.ne_raw > 0.0) {
                std::printf("  [gpu] level %d compresses, at the factors' tolerance and four times tighter: "
                            "Schur %.1f MB of %lld -> %.1fx / %.1fx | edges %.1f MB of %lld -> %.1fx / %.1fx\n",
                            level.level,
                            e.nb_raw / 1e6, static_cast<long long>(e.nb_blocks),
                            e.nb_loose > 0.0 ? e.nb_raw / e.nb_loose : 0.0,
                            e.nb_tight > 0.0 ? e.nb_raw / e.nb_tight : 0.0,
                            e.ne_raw / 1e6, static_cast<long long>(e.ne_blocks),
                            e.ne_loose > 0.0 ? e.ne_raw / e.ne_loose : 0.0,
                            e.ne_tight > 0.0 ? e.ne_raw / e.ne_tight : 0.0);
            }
            if (e.tw_max > 0.0) {
                std::printf("  [gpu] level %d largest wave's transient %.2f GB: A_NS %.2f, X_NS %.2f (scratch, "
                            "%.0f%%) | X_NR %.2f, temp2 %.2f, X_NR orig %.2f (outputs, %.0f%%) | result %.2f, "
                            "T %.2f, small %.2f\n",
                            level.level, e.tw_max / 1e9, e.tw_ans / 1e9, e.tw_xns / 1e9,
                            100.0 * (e.tw_ans + e.tw_xns) / e.tw_max,
                            e.tw_xnr / 1e9, e.tw_temp2 / 1e9, e.tw_xnr_orig / 1e9,
                            100.0 * (e.tw_xnr + e.tw_temp2 + e.tw_xnr_orig) / e.tw_max,
                            e.tw_result / 1e9, e.tw_t / 1e9, e.tw_small / 1e9);
            }
            if (e.wr_final_q[0] + e.wr_final_q[1] + e.wr_final_q[2] + e.wr_final_q[3] > 0.0) {
                const double tot = e.wr_final_q[0] + e.wr_final_q[1] + e.wr_final_q[2] + e.wr_final_q[3];
                std::printf("  [gpu] level %d near blocks go final: %.0f%% / %.0f%% / %.0f%% / %.0f%% of %.2f GB "
                            "by the first, second, third and last quarter of the waves\n", level.level,
                            100.0 * e.wr_final_q[0] / tot, 100.0 * e.wr_final_q[1] / tot,
                            100.0 * e.wr_final_q[2] / tot, 100.0 * e.wr_final_q[3] / tot, tot / 1e9);
            }
            if (e.wr_windows > 0) {
                std::printf("  [gpu] level %d near-block writes: %.2f GB in the largest wave, %.2f GB over %lld "
                            "waves, %.2f GB at their largest over %lld blocks, %.1f writes each\n",
                            level.level, e.wr_max_window / 1e9, e.wr_sum / 1e9,
                            static_cast<long long>(e.wr_windows), e.wr_distinct / 1e9,
                            static_cast<long long>(e.wr_blocks),
                            e.wr_blocks > 0 ? static_cast<double>(e.wr_updates) / static_cast<double>(e.wr_blocks)
                                            : 0.0);
            }
            {   // What a kind holds at the peak, against everything ever
                // allocated to it.  The two answer different questions and
                // confusing them cost three wrong designs: 34.78 GB passed
                // through "adopted blocks" at the 192^3 leaf to leave 0.00
                // resident, while 6.16 GB sat resident at a level that
                // allocated none.  Churn is not occupancy.
                const DeviceHeap& h = DeviceHeap::instance();
                for (ResidentKind k : {ResidentKind::Adopted, ResidentKind::OwnerNear}) {
                    if (h.kind_alloc_bytes(k) == 0 && h.resident_now(k) == 0) continue;
                    std::printf("  [gpu] level %d kind %s: now %.2f GB, ever allocated %.2f GB in %zu, "
                                "freed %.2f GB in %zu\n", level.level, resident_kind_name(k),
                                static_cast<double>(h.resident_now(k)) / 1e9,
                                static_cast<double>(h.kind_alloc_bytes(k)) / 1e9, h.kind_alloc_count(k),
                                static_cast<double>(h.kind_free_bytes(k)) / 1e9, h.kind_free_count(k));
                }
            }
            if (e.np5_blocks > 0) {
                std::printf("  [gpu] level %d step 5's blocks packed as written: %.2f -> %.2f GB (%.1fx) over %lld "
                            "blocks in %lld groups, %.3f s\n", level.level, e.np5_raw / 1e9, e.np5_stored / 1e9,
                            e.np5_stored > 0.0 ? e.np5_raw / e.np5_stored : 1.0,
                            static_cast<long long>(e.np5_blocks), static_cast<long long>(e.np5_groups),
                            e.np5_seconds);
            }
            if (e.gw_products > 0) {
                std::printf("  [gpu] level %d remote generators: %lld partial products staged in %.2f GB, summed "
                            "into %lld targets of %.2f GB (%.1f partials a target, at most %lld)\n", level.level,
                            static_cast<long long>(e.gw_products), e.gw_bytes / 1e9,
                            static_cast<long long>(e.gw_targets), e.gw_target_bytes / 1e9,
                            e.gw_targets > 0 ? static_cast<double>(e.gw_products) / static_cast<double>(e.gw_targets)
                                             : 0.0,
                            static_cast<long long>(e.gw_max_parts));
            }
            if (e.ow_groups > 0) {
                std::printf("  [gpu] level %d owner source window: %.2f GB in the largest group, %.2f GB "
                            "decompressed over %lld groups in %lld slots\n", level.level, e.ow_max / 1e9,
                            e.ow_bytes / 1e9, static_cast<long long>(e.ow_groups),
                            static_cast<long long>(e.ow_slots));
            }
            if (e.xt_raw_bytes > 0.0) {
                std::printf("  [gpu] level %d X_NR^T store: %.2f -> %.2f GB (%.1fx) in %lld slots of %lld boxes"
                            "%s\n", level.level, e.xt_raw_bytes / 1e9, e.xt_stored_bytes / 1e9,
                            e.xt_stored_bytes > 0.0 ? e.xt_raw_bytes / e.xt_stored_bytes : 1.0,
                            static_cast<long long>(e.xt_slots), static_cast<long long>(e.xt_boxes),
                            e.xt_short > 0 ? " (some groups found no staging)" : "");
            }
            if (e.os_groups > 0) {
                std::printf("  [gpu] level %d owner source union: %.2f GB by the largest group, %.2f GB summed "
                            "over %lld groups in %lld slots (largest group with repeats %.2f GB)\n", level.level,
                            e.os_union_max / 1e9, e.os_union_sum / 1e9, static_cast<long long>(e.os_groups),
                            static_cast<long long>(e.os_slots), e.os_raw_max / 1e9);
            }
            if (e.tx_raw_bytes > 0.0) {
                std::printf("  [gpu] level %d kept T and X_SR: %.2f GB raw -> %.2f GB stored (%.1fx); X_NR %.2f -> "
                            "%.2f GB; the LU and pivots %.2f GB, uncompressed\n", level.level, e.tx_raw_bytes / 1e9,
                            e.tx_stored_bytes / 1e9, e.tx_stored_bytes > 0.0 ? e.tx_raw_bytes / e.tx_stored_bytes : 1.0,
                            e.xnr_raw_bytes / 1e9, e.xnr_stored_bytes / 1e9,
                            (e.other_raw_bytes - e.tx_raw_bytes) / 1e9);
            }
            if (e.kx_boxes > 0) {
                std::printf("  [gpu] level %d kept X_NR out of the fill source: %lld boxes in %lld groups, "
                            "%.2f GB raw -> %.2f GB stored, %.3f s%s\n", level.level,
                            static_cast<long long>(e.kx_boxes), static_cast<long long>(e.kx_groups),
                            e.kx_raw_bytes / 1e9, e.kx_stored_bytes / 1e9, e.kx_seconds,
                            e.kx_host > 0 ? " (some went round through the host)" : "");
            }
            if (e.mat_blocks > 0) {
                static const char* const kSites[6] = {"begin Schur", "owner join", "update group",
                                                      "host download", "generators", "build_parent"};
                std::printf("  [gpu] level %d decompressed by:", level.level);
                for (int i = 0; i < 6; ++i) {
                    if (e.mat_site_bytes[i] <= 0.0) continue;
                    std::printf(" %s %.2f GB in %lld", kSites[i], e.mat_site_bytes[i] / 1e9,
                                static_cast<long long>(e.mat_site_calls[i]));
                }
                std::printf("\n");
            }
            if (e.mat_blocks > 0 || e.mat_left > 0) {
                std::printf("  [gpu] level %d adopted blocks materialized: %.2f GB in %lld blocks over %lld batches, "
                            "%.3f s; %lld left packed; %lld group gathers prefetched\n",
                            level.level, e.mat_bytes / 1e9, static_cast<long long>(e.mat_blocks),
                            static_cast<long long>(e.mat_batches), e.mat_seconds,
                            static_cast<long long>(e.mat_left), static_cast<long long>(e.mat_prefetched));
            }
            if (e.pb_blocks > 0) {
                std::printf("  [gpu] level %d parent blocks compress: %.2f -> %.2f GB (%.1fx) over %lld blocks", level.level,
                        e.pb_raw_bytes / 1e9, e.pb_stored_bytes / 1e9,
                        e.pb_stored_bytes > 0.0 ? e.pb_raw_bytes / e.pb_stored_bytes : 1.0,
                        static_cast<long long>(e.pb_blocks));
                if (e.pb_packed_short > 0) std::printf(", %lld did not fit", static_cast<long long>(e.pb_packed_short));
                if (e.pb_host_short > 0) std::printf(", %lld groups left raw for want of host room", static_cast<long long>(e.pb_host_short));
            if (e.pb_check_tol > 0.0) {
                std::printf("; round trip %.3e against %.3e allowed%s", e.pb_check_err, e.pb_check_tol,
                            e.pb_check_err <= e.pb_check_tol ? " (ok)" : "  *** OVER TOLERANCE ***");
            }
                std::printf("\n");
            }
            if (e.tr_chunks > 0) {
                std::printf("  [gpu] level %d transition blocks: %.2f GB in all, %.2f GB written by the largest "
                            "chunk (%lld blocks), %.2f GB summed over %lld chunks, %lld revisits\n",
                            level.level,
                        static_cast<double>(DeviceHeap::instance().resident_now(ResidentKind::Transition)) / 1e9,
                        e.tr_block_bytes_max_chunk / 1e9,
                            static_cast<long long>(e.tr_blocks_max_chunk), e.tr_block_bytes_sum / 1e9,
                            static_cast<long long>(e.tr_chunks), static_cast<long long>(e.tr_block_revisits));
            }
            if (reduced) {
                std::printf("  [gpu] level %d device reduction: %.2f s (plan %.2f, pack %.2f, MPI %.2f, unpack %.2f; "
                            "%d rounds%s), sent %.2f GB in %lld blocks, received %.2f GB in %lld blocks, %lld own "
                            "copies dropped\n",
                            level.level, red.total, red.plan, red.pack, red.mpi, red.unpack, red.rounds,
                            red.staged ? ", through the host" : "", red.bytes_sent / 1e9,
                            static_cast<long long>(red.blocks_sent), red.bytes_received / 1e9,
                            static_cast<long long>(red.blocks_received), static_cast<long long>(red.blocks_dropped));
            }
            std::fflush(stdout);
        }
        if (level_) {
            level_->restore_host_copies();  // (X_NR kept without a host copy)
            level_->download_blocks();      // host transition
            level_.reset();
        }
        compact_device_memory(level.level, announce);
        return on_device;
    }

    // After the transition, the level's own blocks freed: the blocks that
    // outlive it, the solve factors kept on the device (every level's) and
    // the parents' blocks of a device transition, slide down the heap
    // (DeviceHeap::compact), so the free memory is one range again for the
    // next level instead of the pieces between them.  The pointers into them
    // (the kept boxes' factors, the parent matrices) follow.  Local to the
    // rank: no decision is shared.  The moves run on the stream after the
    // level's device work; an overlapping move goes through a staging block
    // from the top of the heap in pieces, or in pieces of its gap.
    void compact_device_memory(int lvl, bool announce) {
        if (!color_gpu_enabled()) return;
        DeviceHeap& heap = DeviceHeap::instance();
        if (!heap.initialized()) return;
        DeviceSolveStore& store = device_solve_store();
        std::vector<char*> blocks;
        for (auto& [level, kl] : store.kept) blocks.insert(blocks.end(), kl.blocks.begin(), kl.blocks.end());
        if (blocks_) {
            for (auto& kv : blocks_->schur) blocks.push_back(reinterpret_cast<char*>(kv.second.ptr));
            for (auto& kv : blocks_->edges) blocks.push_back(reinterpret_cast<char*>(kv.second.ptr));
        }
        if (blocks.empty()) return;
        const auto t0 = std::chrono::steady_clock::now();
        Context& ctx = Context::instance();
        ctx.activate();
        cudaStream_t stream = ctx.stream();
        constexpr size_t kStaging = size_t{256} << 20;
        char* staging = heap.try_alloc(kStaging);
        size_t moved = 0;
        auto copy = [&](char* dst, const char* src, size_t bytes) {
            check_cuda(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToDevice, stream), "compaction");
        };
        auto move = [&](char* dst, const char* src, size_t bytes) {
            moved += bytes;
            if (dst + bytes <= src) {
                copy(dst, src, bytes);
                return;
            }
            // (pieces in address order: a piece's destination is below what
            // the later pieces still read)
            const size_t gap = static_cast<size_t>(src - dst);
            if (staging != nullptr && gap < kStaging) {
                for (size_t off = 0; off < bytes; off += kStaging) {
                    const size_t n = std::min(kStaging, bytes - off);
                    copy(staging, src + off, n);
                    copy(dst + off, staging, n);
                }
                return;
            }
            for (size_t off = 0; off < bytes; off += gap) copy(dst + off, src + off, std::min(gap, bytes - off));
        };
        const std::vector<DeviceHeap::Move> moves = heap.compact(blocks, move);
        if (staging != nullptr) heap.free(staging);
        if (!moves.empty()) {
            for (auto& [level, kl] : store.kept) {
                // (a staged level's pointers are the addresses its blocks had,
                // kept for the upload's rebase; its old ranges may since have
                // been reused by blocks that just moved)
                if (kl.on_host()) continue;
                for (char*& b : kl.blocks) b = DeviceHeap::relocated(b, moves);
                for (auto& [morton, kb] : kl.boxes) {
                    kb.T = DeviceHeap::relocated(kb.T, moves);
                    kb.xsr = DeviceHeap::relocated(kb.xsr, moves);
                    kb.xnr = DeviceHeap::relocated(kb.xnr, moves);
                    kb.xrs_t = DeviceHeap::relocated(kb.xrs_t, moves);
                    kb.xrn_t = DeviceHeap::relocated(kb.xrn_t, moves);
                    kb.lu = DeviceHeap::relocated(kb.lu, moves);
                    kb.ipiv = DeviceHeap::relocated(kb.ipiv, moves);
                }
            }
            if (blocks_) {
                for (auto& kv : blocks_->schur) kv.second.ptr = DeviceHeap::relocated(kv.second.ptr, moves);
                for (auto& kv : blocks_->edges) kv.second.ptr = DeviceHeap::relocated(kv.second.ptr, moves);
            }
            check_cuda(cudaStreamSynchronize(stream), "compaction");
        }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        auto& e = eliminator_stats();
        e.compaction += seconds;
        e.compaction_bytes += static_cast<double>(moved);
        if (announce && !moves.empty()) {
            std::printf("  [gpu] level %d compaction: %zu blocks moved (%.2f GB) in %.3f s; %.2f GB free, "
                        "largest range %.2f GB\n",
                        lvl, moves.size(), static_cast<double>(moved) / 1e9, seconds,
                        static_cast<double>(heap.available()) / 1e9, static_cast<double>(heap.largest_free()) / 1e9);
            std::fflush(stdout);
        }
    }

    // ---- the root

    // The level-1 transition assembled the root block on the device.
    bool holds_root_block() const { return blocks_ != nullptr; }

    // The root's LU on the device: of the block the level-1 transition
    // assembled there, or of the host-assembled block (multi-rank runs)
    // when the heap has room.  Returns false when the host factors it.  A
    // symmetric problem factors the block's symmetric part (as the boxes'
    // X_RR, see compute_and_modify).
    bool factor_root(Box& root, bool announce) {
        const int64_t n = root.num_points;
        const bool symmetric = opt_.is_symmetric && !opt_.is_hermitian;
        if (blocks_) {
            if (announce) std::cout << "  Schur complement size: " << n << " × " << n << std::endl;
            factor_root_on_device(*blocks_, root, symmetric);
            blocks_.reset();
            if (announce) std::cout << "  ✓ Root LU factorization complete (GPU)" << std::endl;
            return true;
        }
        if (color_gpu_enabled() && opt_.method == FactorizationMethod::LU &&
            factor_host_root_on_device(root, symmetric)) {
            if (announce) {
                std::cout << "  Schur complement size: " << n << " × " << n << std::endl;
                std::cout << "  ✓ Root LU factorization complete (GPU)" << std::endl;
            }
            return true;
        }
        return false;
    }

    // Statistics of the level's device box path on this rank, if a wave ran
    // there (BPACK_TRACE=phase resets them per level); returns whether it
    // printed.
    static bool report_level(int lvl) {
        const auto& e = eliminator_stats();
        if (e.boxes == 0) return false;
        std::printf(
            "  [gpu] level %d box path (rank 0): begin %.2f s, sketch %.2f s, ID %.2f s, plan %.2f s, "
            "device %.2f s, download %.2f s, store %.2f s, finish %.2f s | up %.2f GB, "
            "down %.2f GB | owner %lld gemms in %lld batches, %lld new targets | heap peak %.2f GB "
            "(resident %.2f + transient %.2f), %lld reclaims\n",
            lvl, e.begin, e.sketch, e.id, e.plan, e.device, e.download, e.store, e.finish,
            e.bytes_up / 1e9, e.bytes_down / 1e9, static_cast<long long>(e.owner_gemms),
            static_cast<long long>(e.owner_batches), static_cast<long long>(e.new_targets),
            static_cast<double>(e.heap_peak) / 1e9,
            static_cast<double>(e.heap_peak_resident) / 1e9,
            static_cast<double>(e.heap_peak_transient) / 1e9,
            static_cast<long long>(e.heap_reclaims));
        {
            // What the level's resident peak was made of (the ceiling on
            // boxes per GPU; the transient part comes and goes within a wave)
            const DeviceHeap& heap = DeviceHeap::instance();
            std::printf("  [gpu] level %d resident at the peak: %.2f GB", lvl,
                        static_cast<double>(e.heap_peak_resident) / 1e9);
            size_t named = 0;
            for (size_t k = 1; k < static_cast<size_t>(ResidentKind::Count); ++k) {
                const auto kind = static_cast<ResidentKind>(k);
                const size_t bytes = heap.resident_at_peak(kind);
                if (bytes == 0) continue;
                named += bytes;
                std::printf(" | %s %.2f GB", resident_kind_name(kind), static_cast<double>(bytes) / 1e9);
            }
            const size_t other = heap.resident_at_peak(ResidentKind::Other);
            std::printf(" | other %.2f GB", static_cast<double>(other) / 1e9);
            const double accounted = static_cast<double>(named + other);
            std::printf(" (accounted %.2f GB)\n", accounted / 1e9);
            // And the transient half, which is now half the peak.  Without
            // this the only figure to hand was what a wave plans to allocate,
            // which is not what is live when the peak is reached.
            {
                std::printf("  [gpu] level %d transient at the peak: %.2f GB", lvl,
                            static_cast<double>(e.heap_peak_transient) / 1e9);
                size_t tnamed = 0;
                for (size_t k = 1; k < static_cast<size_t>(ResidentKind::Count); ++k) {
                    const auto kind = static_cast<ResidentKind>(k);
                    const size_t bytes = heap.transient_at_peak(kind);
                    if (bytes == 0) continue;
                    tnamed += bytes;
                    std::printf(" | %s %.2f GB", resident_kind_name(kind), static_cast<double>(bytes) / 1e9);
                }
                const size_t tother = heap.transient_at_peak(ResidentKind::Other);
                std::printf(" | other %.2f GB (accounted %.2f GB)\n", static_cast<double>(tother) / 1e9,
                            static_cast<double>(tnamed + tother) / 1e9);
            }
            if (e.persist_blocks > 0) {
                // The peak is live bytes and the other two are cumulative, so
                // subtracting them from it was meaningless - it printed a
                // negative temp2^T.  Each is reported as what it is.
                std::printf("  [gpu] level %d fill sources: %lld blocks, %.2f GB live at their peak; allocated "
                            "over the level: X_RR_full %.2f GB, generator boxes %.2f GB\n", lvl,
                            static_cast<long long>(e.persist_blocks), e.persist_peak / 1e9,
                            e.persist_xrr_bytes / 1e9, e.persist_gen_bytes / 1e9);
            }
            if (e.fs_host_bytes > 0.0 || e.ns_host_bytes > 0.0) {
                const PinnedChunkPool& pool = pinned_chunk_pool();
                std::printf("  [gpu] level %d on the host: fill store %.2f GB, near store %.2f GB at most (pinned chunks "
                            "%.2f GB, %.2f idle; pinning %.2f s so far, of which the levels waited %.2f s; %.2f GB "
                            "unpinned in %.2f s)\n", lvl, e.fs_host_bytes / 1e9, e.ns_host_bytes / 1e9,
                            static_cast<double>(pool.pinned()) / 1e9, static_cast<double>(pool.idle()) / 1e9,
                            pool.pin_seconds(), pool.pin_wait_seconds(), static_cast<double>(pool.unpinned()) / 1e9,
                            pool.unpin_seconds());
                if (e.ns_compact_held > 0.0) {
                    std::printf("  [gpu] level %d near store compacted at the level end: %.2f GB live of %.2f GB held "
                                "(%.2f s)\n", lvl, e.ns_compact_live / 1e9, e.ns_compact_held / 1e9, e.ns_compact_s);
                }
                if (wave_trace_enabled()) {  // which reader takes the host copies (the window must not)
                    std::printf("  [gpu] level %d fill store read as slots (device / host): window %lld / %lld, sketch %lld / "
                                "%lld, transition %lld / %lld\n", lvl, static_cast<long long>(e.fs_read_dev[0]),
                                static_cast<long long>(e.fs_read_host[0]), static_cast<long long>(e.fs_read_dev[1]),
                                static_cast<long long>(e.fs_read_host[1]), static_cast<long long>(e.fs_read_dev[2]),
                                static_cast<long long>(e.fs_read_host[2]));
                }
            }
        }
        if (e.fs_slots > 0 || e.fs_waves_plain > 0) {
            std::printf("  [gpu] level %d fill sources compressed: %.2f -> %.2f GB (%.1fx) in %lld slots of "
                        "%lld boxes, %lld groups fell back, %lld boxes unpacked", lvl, e.fs_raw_bytes / 1e9, e.fs_stored_bytes / 1e9,
                        e.fs_stored_bytes > 0.0 ? e.fs_raw_bytes / e.fs_stored_bytes : 1.0,
                        static_cast<long long>(e.fs_slots), static_cast<long long>(e.fs_boxes),
                        static_cast<long long>(e.fs_waves_plain),
                        static_cast<long long>(e.fs_packed_short));
            if (e.fs_check_tol > 0.0) {
                std::printf("; round trip worst error %.3e against %.3e allowed%s", e.fs_check_err, e.fs_check_tol,
                            e.fs_check_err <= e.fs_check_tol ? " (ok)" : "  *** OVER TOLERANCE ***");
            }
            std::printf("\n");
        }
        if (e.lsch_pairs > 0 || e.ledge_pairs > 0) {
            std::printf("  [gpu] level %d blocks a wave reads: Schur %.2f GB largest / %.2f GB summed over %lld "
                        "reads of %lld blocks; edges %.2f GB largest / %.2f GB summed over %lld reads of %lld\n",
                        lvl, e.lsch_union_max_wave / 1e9, e.lsch_union_sum / 1e9,
                        static_cast<long long>(e.lsch_pairs), static_cast<long long>(e.lsch_distinct),
                        e.ledge_union_max_wave / 1e9, e.ledge_union_sum / 1e9,
                        static_cast<long long>(e.ledge_pairs), static_cast<long long>(e.ledge_distinct));
        }
        if (e.fs_wave_slots > 0) {
            std::printf("  [gpu] level %d fill slots materialized: %.2f GB in the largest wave's scratch, "
                        "%.2f GB over the level in %lld slot chunks\n", lvl, e.fs_wave_scratch / 1e9,
                        e.fs_wave_bytes / 1e9, static_cast<long long>(e.fs_wave_slots));
        }
        if (e.sk_src_waves > 0) {
            std::printf("  [gpu] level %d fill sources read: %lld distinct over %lld waves, %.2f GB by the largest "
                        "wave, %.2f GB summed (%lld wave-source pairs, %.1f waves per source)\n", lvl,
                        static_cast<long long>(e.sk_src_distinct), static_cast<long long>(e.sk_src_waves),
                        e.sk_src_bytes_max_wave / 1e9, e.sk_src_bytes_sum / 1e9,
                        static_cast<long long>(e.sk_src_pairs),
                        e.sk_src_distinct > 0 ? static_cast<double>(e.sk_src_pairs) / e.sk_src_distinct : 0.0);
            std::printf("  [gpu] level %d fill source slices: %.2f GB by the largest wave, %.2f GB summed "
                        "(%.1f%% of the sources a wave reads)\n", lvl, e.sk_slice_bytes_max_wave / 1e9,
                        e.sk_slice_bytes_sum / 1e9,
                        e.sk_src_bytes_sum > 0.0 ? 100.0 * e.sk_slice_bytes_sum / e.sk_src_bytes_sum : 0.0);
            std::printf("  [gpu] level %d fill source runs (W side): %.2f GB by the largest wave, %.2f GB summed, "
                        "%lld runs; P+W largest wave %.2f GB, summed %.2f GB\n", lvl,
                        e.sk_run_bytes_max_wave / 1e9, e.sk_run_bytes_sum / 1e9,
                        static_cast<long long>(e.sk_runs),
                        (e.sk_slice_bytes_max_wave + e.sk_run_bytes_max_wave) / 1e9,
                        (e.sk_slice_bytes_sum + e.sk_run_bytes_sum) / 1e9);
            std::printf("  [gpu] level %d fill source union (distinct slots a wave reads): %.2f GB by the largest "
                        "wave, %.2f GB summed\n", lvl, e.sk_union_bytes_max_wave / 1e9,
                        e.sk_union_bytes_sum / 1e9);
        }
        std::printf("  [gpu] level %d detail: sketch plan %.2f s, sketch device %.2f s "
                    "(meta build %.2f, upload %.2f, rows %.2f, stored %.2f, P %.2f, fill %.2f, "
                    "ID %.2f, ranks down %.2f), background copies: busy %.2f s, level-end wait %.2f s\n",
                    lvl, e.sketch_plan, e.sketch_gpu, e.sk_meta, e.sk_upload, e.sk_rows,
                    e.sk_stored, e.sk_p, e.sk_fill, e.sk_id, e.sk_download, e.finish_store, e.finish_sources);
        if (e.adaptive_rounds > 0) {
            std::printf("  [gpu] level %d adaptive ID rows: %lld rounds over the waves, later rounds %.2f s, "
                        "%lld boxes sampled more than one node\n",
                        lvl, static_cast<long long>(e.adaptive_rounds), e.sk_adaptive,
                        static_cast<long long>(e.adaptive_multi));
        }
        std::printf("  [gpu] level %d background copier: busy %.2f s, of which waiting for device data "
                    "%.2f s\n", lvl, e.finish_store, e.copier_wait);
        std::printf("  [gpu] level %d elimination device: fills %.2f, X_RR/X_SR %.2f, LU %.2f, X_NR %.2f, "
                    "solves %.2f, Schur %.2f, near blocks+owner targets %.2f, owner GEMMs %.2f s\n",
                    lvl, e.el[0], e.el[1], e.el[2], e.el[3], e.el[4], e.el[5], e.el[6], e.el[7]);
        std::printf("  [gpu] level %d host: plan boxes %.2f s [buffers %.2f, items %.2f], near blocks %.2f s "
                    "(overlapped; allocations %.2f), owner pass %.2f s (overlapped) [candidates %.2f, pairs %.2f, join %.2f, "
                    "batches %.2f], launch %.2f s, sketch wait %.2f s, heap-reclaim wait %.2f s, exchange-buffer "
                    "wait %.2f s | GF/s: owner %.0f, solves %.0f\n",
                    lvl, e.plan_boxes, e.pb_buffers, e.pb_items, e.pb_near, e.pb_near_alloc, e.plan_owner, e.po_candidates,
                    e.po_pairs, e.po_order, e.po_batches, e.launch, e.sk_wait, e.reclaim_wait, e.exchange_wait,
                    e.owner_flops / std::max(e.el[7], 1e-9) / 1e9,
                    e.solve_flops / std::max(e.el[4], 1e-9) / 1e9);
        if (wave_trace_enabled()) {
            const DeviceHeap& hp = DeviceHeap::instance();
            std::printf("  [gpu] level %d heap work (the process so far): %zu allocations %.2f s, %zu frees %.2f s; "
                        "in the box region's buffers %.2f s (%lld calls); the box region planned while the stream "
                        "was idle in %lld of %lld chunks, %.2f s\n",
                        lvl, hp.alloc_calls(), hp.alloc_seconds(), hp.free_calls(), hp.free_seconds(), e.pb_heap,
                        static_cast<long long>(e.pb_heap_ops), static_cast<long long>(e.pb_idle_chunks),
                        static_cast<long long>(e.pb_chunks), e.pb_idle_plan);
            std::printf("  [gpu] level %d owner groups' raw blocks in one heap block a group: %lld groups, %lld blocks; "
                        "%lld groups found no range for it, %lld blocks moved out raw\n", lvl,
                        static_cast<long long>(e.og_arena_groups), static_cast<long long>(e.og_arena_blocks),
                        static_cast<long long>(e.og_arena_declined), static_cast<long long>(e.og_arena_moved));
        }
        if (e.remote_generators > 0 || e.exchange > 0.0) {
            const auto& tt = transport_timers();
            std::printf("  [gpu] level %d ranks (%s memory): emit %.2f s, exchange %.2f s [sizes %.2f, rest %.2f "
                        "(serialize %.2f), deserialize %.2f, assisting %.2f, install %.2f; sent %.2f GB, "
                        "received %.2f GB], receive %.2f s (%lld remote generators, %lld buffers "
                        "outside the exchange arena)\n",
                        lvl, e.device_exchange ? "device" : "host", e.emit, e.exchange,
                        tt.sizes, tt.payload, tt.serialize,
                        tt.deserialize, tt.assisting, tt.install, tt.bytes_sent / 1e9,
                        tt.bytes_received / 1e9, e.remote, static_cast<long long>(e.remote_generators),
                        static_cast<long long>(e.exchange_fallbacks));
            if (e.exchange_host > 0) {
                std::printf("  [gpu] level %d exchange buffers in pinned host memory (no device range fitted): %lld, %.2f GB\n",
                            lvl, static_cast<long long>(e.exchange_host), e.exchange_host_bytes / 1e9);
            }
        }
        std::fflush(stdout);
        return true;
    }

private:
    // Entries through a collective callback: rounds until no rank of the
    // level is fetching (collective_entries.hpp), before communication
    // between them.  Nothing otherwise.
    void quiesce_entries() {
        if (auto* entries = kernel_entries(kernel_, 0)) entries->quiesce();
    }

    Tree* tree_;
    KernelType* kernel_;
    Options opt_;
    int level_index_ = -1;
    std::unique_ptr<LevelEliminatorBase<CoordType, DataType>> level_;  // this level's device box path
    std::unique_ptr<DeviceLevelBlocks<DataType>> blocks_;              // the next level's blocks
    WaveTrace trace_;                                                  // BPACK_TRACE=wave
    // the keep estimate's scale of the bound (decide_solve_keep): the bound
    // itself for the first level (the leaf, whose large working set the
    // bound's slack leaves room for: its factors take about half of it),
    // then the level below's largest factors/bound times the margin
    static constexpr double kFirstKeepRatio = 1.0, kKeepMargin = 1.25;
    // H2_compress_factors: how much of the bound's X_NR term the keep expects
    // to store (decide_solve_keep).  The bound already overstates X_NR by
    // about n/r (levels come in at 53-68% of it) and the compression takes
    // another 18x, so anything up to about 8 would be safe on the ratio.
    //
    // It is 1 - the compression is not taken into account - because the
    // budget is the wrong quantity to relax while the leaf's working set is
    // what it is: decide_solve_keep is called at level start and measures the
    // heap as it is *then*, before the level has allocated the 26 GB it will
    // hold at its peak, so for the leaf the budget is structurally
    // over-optimistic whatever this divisor says.  Raising it does flip the
    // leaf's verdict, and the leaf then keeps about half its factors and
    // fails part way (GPU_CA_PLAN.md C4).  Raise it once the leaf's working
    // set leaves room - step 4 - or once the decision can estimate the
    // level's own working set.
    static constexpr double kXnrBoundDivisor = 4.0;
    double keep_ratio_ = kFirstKeepRatio;
    double keep_bound_ = 0.0;                                          // this level's (decide_solve_keep)
    // The kept store's raw / stored, measured on the last level (0: none
    // yet; kFirstCompression then, with the factors compressed).  The keep
    // estimate and the hold limit divide by it: without it the leaf at 8x
    // the boxes a rank declines on the full card (158 GB raw against a 25 GB
    // budget) and its X_NR goes to the host uncompressed (2026-10-10).
    static constexpr double kFirstCompression = 8.0;
    double compression_ = 0.0;
    double kept_expected_ = 0.0;                                       // this level's estimate, compressed
};

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
