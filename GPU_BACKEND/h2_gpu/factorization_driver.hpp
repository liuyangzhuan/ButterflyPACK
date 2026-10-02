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
#include "level_eliminator.hpp"
#include "wave_trace.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
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

    // Before the first level: the application's evaluator, when a level may
    // take the device box path, which calls it (Evaluator::warm_up, once per
    // evaluator: an NVRTC compile, CuPy compiling its kernels, ...).  Returns
    // its time (0 when there is nothing to warm), which the factor time
    // leaves out.
    double warm_up_evaluator() {
        if constexpr (!gpu_data_type<DataType>) {
            return 0.0;
        } else {
            if (!color_gpu_enabled() || opt_.use_sketch != 2 || !opt_.is_symmetric || opt_.is_hermitian ||
                opt_.method != FactorizationMethod::LU) {
                return 0.0;
            }
            return warm_up_registered_evaluator<DataType>(tree_, kernel_->gpu_evaluator);
        }
    }

    // Before the first level of the process's first factorization: a 64 KB
    // message with every rank this one exchanges with on a level (each
    // level's one-hop neighbour ranks), sent and received the way the levels
    // will (from and into the exchange arena with GPU-aware MPI).  A pair's
    // first device exchanges are slow (Laplace 192^3 on 64 ranks: the first
    // waves of the first exchanging level at 3-4 GB/s, ~12 GB/s from the
    // fourth on; ~0.1-0.2 s per factorization, more at larger sizes); this
    // takes that first contact out of the levels.  Returns its time (0 when
    // it already ran); collective over the tree's ranks.
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
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }

    // ---- where the blocks of a level go next (tested where levels start)

    // The device box path runs level `lvl`.
    bool box_path_runs(int lvl) const {
        if (!color_gpu_enabled() || lvl <= 1 || tree_->level_uses_CA(lvl)) return false;
        if (!tree_->levels[static_cast<size_t>(lvl)].is_process_active) return false;
        const bool streamed = opt_.use_sketch == 2 && opt_.is_symmetric && !opt_.is_hermitian;
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

    // ---- a Color level that eliminates

    // Level start (all ranks of level_comm): the device box path when it
    // covers the level (streamed sketches, lazy far fill, a symmetric
    // kernel).  Its sketches run on the host when a box of any rank is
    // wider than the device sketch's packed lists take.  `announce`: this
    // rank prints the decision.
    void start_level(int lvl, bool use_streamed_level, MPI_Comm level_comm, bool announce) {
        if (!(color_gpu_enabled() && use_streamed_level && lazy_far_field_mode() == LazyFarFieldMode::LAZY &&
              opt_.is_symmetric && !opt_.is_hermitian)) {
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
                                       std::move(blocks_), opt_.occupancy, fits != 0);
        level_index_ = lvl;
        decide_solve_keep(lvl, level_comm, false);
        if (wave_trace_enabled()) trace_.start_end();
        if (announce) {
            std::cout << "  GPU box path: "
                      << (level_ ? std::string(tensor_core_gemm() ? "on (FP64 tensor-core GEMMs)" : "on")
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
        if (!color_gpu_enabled() || lvl <= 1 || !tree_->level_uses_CA(lvl) || owner_component != 0) return false;
        const bool streamed = opt_.use_sketch == 2 && opt_.is_symmetric && !opt_.is_hermitian;
        return streamed && opt_.lazy_schur > 0 && opt_.method == FactorizationMethod::LU && !opt_.occupancy;
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
                                           /*ca=*/true);
        }
        level_index_ = lvl;
        decide_solve_keep(lvl, level_comm, true);
        if (wave_trace_enabled()) trace_.start_end();
        if (announce) {
            std::cout << "  GPU CA level: "
                      << (level_ ? std::string(tensor_core_gemm() ? "on (FP64 tensor-core GEMMs)" : "on")
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
    void decide_solve_keep(int lvl, MPI_Comm level_comm, bool ca) {
        if (!color_gpu_enabled()) return;
        int ok = 1;
        double budget = 0.0;
        keep_bound_ = 0.0;
        if (level_) {
            const DeviceHeap& heap = DeviceHeap::instance();
            keep_bound_ = solve_keep_bound(tree_->levels[static_cast<size_t>(lvl)], ca);
            budget = kSolveKeepFraction * static_cast<double>(heap.capacity()) - static_cast<double>(heap.used());
            ok = keep_ratio_ * keep_bound_ <= budget ? 1 : 0;
        }
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, level_comm);
        if (level_) level_->set_solve_keep(ok != 0, budget);
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
        if (!ca_replica_check_enabled()) return;
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
    template<typename HostTransport>
    Duration transport(PendingFactorUpdates<DataType>& pending, HostTransport&& host_transport) {
        if (level_->device_exchange()) return level_->exchange();
        const auto t_emit = std::chrono::steady_clock::now();
        level_->emit_generators(pending);
        const auto t_exchange = std::chrono::steady_clock::now();
        std::vector<int64_t> installed;
        const Duration duration = host_transport(&installed);
        eliminator_stats().emit += std::chrono::duration<double>(t_exchange - t_emit).count();
        eliminator_stats().exchange +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_exchange).count();
        level_->receive_remote(installed);
        return duration;
    }

    // Box region, owner pass, mirror and finalize of a wave on the device;
    // returns its number of boundary boxes.  The loop records the wave's
    // boxes as eliminated.
    int eliminate_wave(const std::vector<int64_t>& wave, int counter) {
        if (!wave_trace_enabled()) return level_->eliminate_wave(wave, counter);
        trace_.wave_begin();
        const int boundary = level_->eliminate_wave(wave, counter);
        trace_.wave_end(counter, wave.size());
        return boundary;
    }

    // Level end, after the final transport (all ranks of level_comm): the
    // factors on the host, the blocks kept for the transition, and whether
    // the solve's factors stay on the device (only if every rank of the
    // level still holds all of its own).
    void finish_level(int lvl, const PendingFactorUpdates<DataType>& pending, MPI_Comm level_comm, bool announce) {
        if (level_) {
            // (every rank of level_comm has a device level here)
            const bool trace = wave_trace_enabled();
            if (trace) trace_.finish_begin();
            level_->finish();
            if (trace) {
                trace_.finish_end();
                print_wave_trace(lvl, level_comm, announce);
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
            double v[4] = {keep_bound_ > 0.0 ? bytes / keep_bound_ : 0.0, bytes, failure == 1 ? 1.0 : 0.0,
                           failure == 2 ? 1.0 : 0.0};
            MPI_Allreduce(MPI_IN_PLACE, v, 4, MPI_DOUBLE, MPI_MAX, level_comm);
            if (v[2] > 0.0 || v[3] > 0.0) {
                keep_ratio_ = 1.0;
            } else if (v[0] > 0.0) {
                keep_ratio_ = std::min(1.0, kKeepMargin * v[0]);
            }
            if (level_) level_->commit_solve_factors(kept != 0);
            if (announce) {
                if (kept != 0) {
                    std::printf("  [gpu] level %d solve factors: kept on the device (up to %.2f GB per rank, at most "
                                "%.0f%% of their bound)\n",
                                lvl, v[1] / 1e9, 100.0 * v[0]);
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

    // Device transition: this level's blocks never leave the device, and the
    // parent's go to the next level's eliminator (or to the host when that
    // level runs there).  Each rank builds the parents of its own boxes, with
    // its copies of the blocks shared with other ranks; at a process
    // reduction they then go to the host, for their new owner.
    // `build_structure()` returns the parent boxes without their blocks.
    // Returns whether the parents were built here; otherwise the level's
    // blocks are back in the host BoxData for the host transition.
    // `level_comm` (the level's active ranks, or MPI_COMM_NULL): a device CA
    // level next keeps the blocks on the device (its halo then moves between
    // the devices, M2) when every rank's parents fit the device sketch.
    template<typename StructureBuilder>
    bool transition(const Level& level, const Level& parent_level, std::vector<Box>& parents,
                    StructureBuilder&& build_structure, bool announce, MPI_Comm level_comm = MPI_COMM_NULL) {
        bool on_device = false;
        const bool reduction_ahead = parent_level.num_active_processes != level.num_active_processes;
        if (level_ && level.is_process_active && level_->can_build_parent() &&
            (reduction_ahead || parent_level.is_process_active) && opt_.is_symmetric && !opt_.is_hermitian) {
            parents = build_structure();
            blocks_ = level_->build_parent(parents);
            level_.reset();
            const int parent_lvl = level_index_ - 1;
            bool keep = !reduction_ahead && keeps_blocks_of(parent_lvl);
            if (!keep && !reduction_ahead && level_comm != MPI_COMM_NULL && ca_level_runs(parent_lvl) && ca_level_supported(parent_lvl, nullptr)) {
                int fits = 1;
                for (const Box& p : parents) {
                    if (p.num_points > kEntryDestMask) fits = 0;
                }
                MPI_Allreduce(MPI_IN_PLACE, &fits, 1, MPI_INT, MPI_MIN, level_comm);
                keep = fits != 0;
            }
            if (!keep) {
                download_level_blocks(*blocks_, parents);
                blocks_.reset();
            }
            on_device = true;
            if (announce) {
                const auto& e = eliminator_stats();
                std::printf("  [gpu] level %d device transition: %.2f s (%lld chunks: plan %.2f [restore %.2f, %.2f GB], "
                            "blocks %.2f, P %.2f, fill %.2f), %lld fill GEMMs (%.1f GFLOP), blocks %s\n",
                            level_index_, e.transition, static_cast<long long>(e.tr_chunks), e.tr_plan,
                            e.tr_restore, e.tr_restore_bytes / 1e9, e.tr_blocks, e.tr_p, e.tr_fill,
                            static_cast<long long>(e.transition_fill_gemms), e.transition_flops / 1e9,
                            blocks_ ? "kept on the device" : "copied to the host");
                std::fflush(stdout);
            }
        }
        if (level_) {
            level_->restore_host_copies();  // (X_NR kept without a host copy)
            level_->download_blocks();      // host transition
            level_.reset();
        }
        return on_device;
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
            "down %.2f GB | owner %lld gemms in %lld batches, %lld new targets | heap peak %.2f GB, %lld reclaims\n",
            lvl, e.begin, e.sketch, e.id, e.plan, e.device, e.download, e.store, e.finish,
            e.bytes_up / 1e9, e.bytes_down / 1e9, static_cast<long long>(e.owner_gemms),
            static_cast<long long>(e.owner_batches), static_cast<long long>(e.new_targets),
            static_cast<double>(e.heap_peak) / 1e9, static_cast<long long>(e.heap_reclaims));
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
                    "solves %.2f, Schur+near %.2f, owner targets %.2f, owner GEMMs %.2f s\n",
                    lvl, e.el[0], e.el[1], e.el[2], e.el[3], e.el[4], e.el[5], e.el[6], e.el[7]);
        std::printf("  [gpu] level %d host: plan boxes %.2f s [buffers %.2f, items %.2f of which near-block "
                    "allocations %.2f], owner pass %.2f s (overlapped) [candidates %.2f, pairs %.2f, join %.2f, "
                    "batches %.2f], launch %.2f s, sketch wait %.2f s, heap-reclaim wait %.2f s, exchange-buffer "
                    "wait %.2f s | GF/s: owner %.0f, solves %.0f\n",
                    lvl, e.plan_boxes, e.pb_buffers, e.pb_items, e.pb_near_alloc, e.plan_owner, e.po_candidates,
                    e.po_pairs, e.po_order, e.po_batches, e.launch, e.sk_wait, e.reclaim_wait, e.exchange_wait,
                    e.owner_flops / std::max(e.el[7], 1e-9) / 1e9,
                    e.solve_flops / std::max(e.el[4], 1e-9) / 1e9);
        if (e.remote_generators > 0 || e.exchange > 0.0) {
            const auto& tt = transport_timers();
            std::printf("  [gpu] level %d ranks (%s memory): emit %.2f s, exchange %.2f s [sizes %.2f, payload %.2f "
                        "(serialize %.2f), deserialize %.2f, assisting %.2f, install %.2f; sent %.2f GB, "
                        "received %.2f GB], receive %.2f s (%lld remote generators, %lld buffers "
                        "outside the exchange arena)\n",
                        lvl, e.device_exchange ? "device" : "host", e.emit, e.exchange,
                        tt.sizes, tt.payload, tt.serialize,
                        tt.deserialize, tt.assisting, tt.install, tt.bytes_sent / 1e9,
                        tt.bytes_received / 1e9, e.remote, static_cast<long long>(e.remote_generators),
                        static_cast<long long>(e.exchange_fallbacks));
        }
        std::fflush(stdout);
        return true;
    }

private:
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
    double keep_ratio_ = kFirstKeepRatio;
    double keep_bound_ = 0.0;                                          // this level's (decide_solve_keep)
};

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
