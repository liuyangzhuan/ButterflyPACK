#pragma once

#include "gpu_common/bpack_env.hpp"
#include "core/butterfly_types.hpp"
#include "core/butterfly_ca_level.hpp"
#include "core/butterfly_communicators.hpp"
#include "core/butterfly_logdet.hpp"
#include "butterfly_solve.hpp"
#include "butterfly_verification.hpp"
#include "core/memory_diagnostics.hpp"
#include "h2_gpu/owner_pass_gpu.hpp"
#include "h2_gpu/factorization_driver.hpp"
#include "h2_gpu/ca_owner.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace butterfly {
using namespace fmm;

// ---------------------------------------------------------------------------
// BPACK_TRACE=skew (doc/environment_variables.md): where a level's waits for
// the neighbours come from.  Every rank records the elimination's wall time
// and the part of it spent waiting for or receiving neighbour data (the
// transports' `sizes` and `payload` timers, serialization excluded), per
// wave of the Color loop as well.  At the level's end the level's first rank
// gathers them and replays the Color waves with zero latency: a wave that
// begins with a transport starts on a rank when its own and its
// neighbours' previous waves are done (the data that transport carries),
// any other wave when its own previous one is, and takes the rank's
// measured busy time (wall minus wait); a wave's transport first does its
// send-side work (the last wave's device tail, the headers: the send-ready
// time), and the neighbours wait for that too, so a neighbour's data is
// ready at the end of its previous wave plus its send-ready time.  The model's level time, against the
// average busy time and against the measured level time, splits the level
// into compute, work skew (ranks with more work, and the chains of waits
// they cause) and latency (the messages, the transfers and the
// synchronisation itself, on the critical path).  The per-rank data go to
// skew_<job>.<step>_level<L>.csv in the working directory.
// ---------------------------------------------------------------------------
struct LevelSkewProbe {
    using clock = std::chrono::steady_clock;
    bool on = false;
    clock::time_point t_level, t_wave;
    double wait_level0 = 0.0, wait_wave0 = 0.0, bytes0 = 0.0, presend_wave0 = 0.0;
    double pre_ms = -1.0;              // level start to the first wave (-1: no wave loop)
    std::vector<double> wave_wall;     // ms, per iteration of the wave loop, then the final transport
    std::vector<double> wave_wait;     // ms, the part waiting for or receiving neighbour data
    std::vector<char> wave_dep;        // the iteration began with a transport
    std::vector<double> wave_presend;  // ms, its transport's send-ready time

    static double waited_ms() {
        const auto& t = transport_timers();
        return (t.sizes + t.payload - t.serialize) * 1e3;
    }
    static double ms_since(clock::time_point t) {
        return std::chrono::duration<double, std::milli>(clock::now() - t).count();
    }
    void level_begin(bool enabled) {
        on = enabled;
        pre_ms = -1.0;
        wave_wall.clear();
        wave_wait.clear();
        wave_dep.clear();
        wave_presend.clear();
        if (!on) return;
        t_level = clock::now();
        wait_level0 = waited_ms();
        bytes0 = transport_timers().bytes_received;
    }
    void loop_begin() {
        if (on) pre_ms = ms_since(t_level);
    }
    void wave_begin(bool dep) {
        if (!on) return;
        t_wave = clock::now();
        wait_wave0 = waited_ms();
        presend_wave0 = transport_timers().presend * 1e3;
        wave_dep.push_back(dep ? 1 : 0);
    }
    void wave_end() {
        if (!on) return;
        wave_wall.push_back(ms_since(t_wave));
        wave_wait.push_back(waited_ms() - wait_wave0);
        wave_presend.push_back(transport_timers().presend * 1e3 - presend_wave0);
    }
    // the wave loop's iterations end on every path (`continue` included)
    struct WaveScope {
        LevelSkewProbe& p;
        WaveScope(LevelSkewProbe& probe, bool dep) : p(probe) { p.wave_begin(dep); }
        ~WaveScope() { p.wave_end(); }
    };
};

// The level's report (collective over level_comm; `elim_ms` this rank's
// elimination wall time).
template<typename CoordType, typename DataType>
void level_skew_report(const LevelSkewProbe& p, fmm::ParallelTree<CoordType, DataType>* tree,
                       const fmm::TreeLevel<CoordType, DataType>& level, int lvl, MPI_Comm level_comm, double elim_ms) {
    int P = 1, me = 0;
    MPI_Comm_size(level_comm, &P);
    MPI_Comm_rank(level_comm, &me);
    std::vector<int> nbrs = compute_one_hop_neighbor_ranks(tree, level, lvl);
    nbrs.erase(std::remove(nbrs.begin(), nbrs.end(), tree->mpi_rank), nbrs.end());
    int dims[2] = {static_cast<int>(p.wave_wall.size()), static_cast<int>(nbrs.size())};
    int maxd[2] = {0, 0};
    MPI_Allreduce(dims, maxd, 2, MPI_INT, MPI_MAX, level_comm);
    const int NW = maxd[0], NN = maxd[1];
    // [tree rank, neighbours, elim, wait, MB received, pre, waves, (wall, wait, dep, send-ready) x NW, neighbour ranks x NN]
    const int H = 7, K = 4, REC = H + K * NW + NN;
    std::vector<double> rec(static_cast<size_t>(REC), -1.0);
    rec[0] = tree->mpi_rank;
    rec[1] = static_cast<double>(nbrs.size());
    rec[2] = elim_ms;
    rec[3] = LevelSkewProbe::waited_ms() - p.wait_level0;
    rec[4] = (transport_timers().bytes_received - p.bytes0) / 1e6;
    rec[5] = p.pre_ms;
    rec[6] = static_cast<double>(p.wave_wall.size());
    for (size_t w = 0; w < p.wave_wall.size(); ++w) {
        rec[static_cast<size_t>(H + K * w)] = p.wave_wall[w];
        rec[static_cast<size_t>(H + K * w + 1)] = p.wave_wait[w];
        rec[static_cast<size_t>(H + K * w + 2)] = p.wave_dep[w];
        rec[static_cast<size_t>(H + K * w + 3)] = p.wave_presend[w];
    }
    for (size_t i = 0; i < nbrs.size(); ++i) rec[static_cast<size_t>(H + K * NW) + i] = nbrs[i];
    std::vector<double> all(me == 0 ? static_cast<size_t>(REC) * static_cast<size_t>(P) : 0);
    MPI_Gather(rec.data(), REC, MPI_DOUBLE, all.data(), REC, MPI_DOUBLE, 0, level_comm);
    if (me != 0) return;
    auto at = [&](int r, int i) { return all[static_cast<size_t>(r) * static_cast<size_t>(REC) + static_cast<size_t>(i)]; };
    // per-rank totals
    double T = 0, c_sum = 0, c_min = 1e300, c_max = -1, w_sum = 0, w_of_busiest = 0, mb_sum = 0;
    int busiest = 0;
    for (int r = 0; r < P; ++r) {
        const double t = at(r, 2), w = at(r, 3), c = t - w;
        T = std::max(T, t);
        c_sum += c;
        w_sum += w;
        mb_sum += at(r, 4);
        c_min = std::min(c_min, c);
        if (c > c_max) {
            c_max = c;
            busiest = r;
            w_of_busiest = w;
        }
    }
    const double c_avg = c_sum / P;
    // by number of neighbour ranks
    std::map<int, std::array<double, 4>> cls;   // count, busy, wait, MB
    for (int r = 0; r < P; ++r) {
        auto& a = cls[static_cast<int>(at(r, 1))];
        a[0] += 1;
        a[1] += at(r, 2) - at(r, 3);
        a[2] += at(r, 3);
        a[3] += at(r, 4);
    }
    std::string by;
    for (const auto& [n, a] : cls) {
        char buf[160];
        std::snprintf(buf, sizeof(buf), " | %d neighbours: %d ranks, busy %.0f, wait %.0f ms, %.0f MB in", n,
                      static_cast<int>(a[0]), a[1] / a[0], a[2] / a[0], a[3] / a[0]);
        by += buf;
    }
    std::printf("  [skew] level %d (%d ranks): elimination %.0f ms (slowest rank); busy min/avg/max %.0f/%.0f/%.0f ms "
                "(the busiest: tree rank %d, %d neighbours, waits %.0f ms); waits avg %.0f ms\n",
                lvl, P, T, c_min, c_avg, c_max, static_cast<int>(at(busiest, 0)), static_cast<int>(at(busiest, 1)),
                w_of_busiest, w_sum / P);
    std::printf("  [skew] level %d by neighbour ranks:%s\n", lvl, by.c_str());
    // the zero-latency replay of the waves
    if (NW > 0 && at(0, 5) >= 0.0) {
        std::unordered_map<int, int> index_of;   // tree rank -> gathered index
        for (int r = 0; r < P; ++r) index_of[static_cast<int>(at(r, 0))] = r;
        std::vector<double> E(static_cast<size_t>(P)), next(static_cast<size_t>(P)), measured(static_cast<size_t>(P));
        for (int r = 0; r < P; ++r) {
            E[static_cast<size_t>(r)] = std::max(0.0, at(r, 5));
            measured[static_cast<size_t>(r)] = std::max(0.0, at(r, 5));
        }
        double busy_sum = 0;
        for (int r = 0; r < P; ++r) busy_sum += std::max(0.0, at(r, 5));
        for (int w = 0; w < NW; ++w) {
            for (int r = 0; r < P; ++r) {
                if (w >= static_cast<int>(at(r, 6))) {
                    next[static_cast<size_t>(r)] = E[static_cast<size_t>(r)];
                    continue;
                }
                const double wall = at(r, H + K * w), wait = at(r, H + K * w + 1);
                const bool dep = at(r, H + K * w + 2) > 0.5;
                const double ready = std::max(0.0, at(r, H + K * w + 3));   // send-ready time (busy)
                double start = E[static_cast<size_t>(r)] + ready;          // its own send side
                if (dep) {
                    for (int i = 0; i < static_cast<int>(at(r, 1)); ++i) {
                        const auto it = index_of.find(static_cast<int>(at(r, H + K * NW + i)));
                        if (it == index_of.end()) continue;
                        const int q = it->second;
                        const double q_ready = w < static_cast<int>(at(q, 6)) ? std::max(0.0, at(q, H + K * w + 3)) : 0.0;
                        start = std::max(start, E[static_cast<size_t>(q)] + q_ready);
                    }
                }
                next[static_cast<size_t>(r)] = start + std::max(0.0, wall - wait - ready);
                measured[static_cast<size_t>(r)] += wall;
                busy_sum += std::max(0.0, wall - wait);
            }
            E.swap(next);
        }
        const double model = *std::max_element(E.begin(), E.end());
        const double meas = *std::max_element(measured.begin(), measured.end());
        const double busy_avg = busy_sum / P;
        double end_avg = 0, end_max = 0;
        for (int r = 0; r < P; ++r) {
            const double end = at(r, 2) - measured[static_cast<size_t>(r)];
            end_avg += end / P;
            end_max = std::max(end_max, end);
        }
        std::printf("  [skew] level %d waves (%d, the last the final transport): measured %.0f ms to the final transport's end, "
                    "zero-latency replay %.0f ms, average busy %.0f ms => compute %.0f, work skew %.0f, latency %.0f ms "
                    "(then the level's end: avg %.0f, max %.0f ms)\n",
                    lvl, NW, meas, model, busy_avg, busy_avg, model - busy_avg, meas - model, end_avg, end_max);
    }
    // the raw data
    const char* job = std::getenv("SLURM_JOB_ID");
    const char* step = std::getenv("SLURM_STEP_ID");
    const std::string name = std::string("skew_") + (job ? job : "run") + "." + (step ? step : "0") + "_level" +
                             std::to_string(lvl) + ".csv";
    if (FILE* f = std::fopen(name.c_str(), "w")) {
        std::fprintf(f, "tree_rank,neighbours,elim_ms,wait_ms,mb_in,pre_ms,waves");
        for (int w = 0; w < NW; ++w) std::fprintf(f, ",wall%d,wait%d,dep%d,ready%d", w, w, w, w);
        for (int i = 0; i < NN; ++i) std::fprintf(f, ",nbr%d", i);
        std::fprintf(f, "\n");
        for (int r = 0; r < P; ++r) {
            for (int i = 0; i < REC; ++i) std::fprintf(f, i == 0 ? "%.6g" : ",%.6g", at(r, i));
            std::fprintf(f, "\n");
        }
        std::fclose(f);
    }
    std::fflush(stdout);
}


#ifdef H2_HAVE_GPU
// A replicated CA level on the device (h2_gpu/GPU_CA_PLAN.md, M1): the
// schedule of factorize_CA_level (its groups, each in Morton-parity waves,
// with the same wave numbers), every wave eliminated by the GPU driver.  The
// ghost boxes are eliminated along with the local ones; nothing is sent
// during the level.
template<typename CoordType, typename DataType, typename KernelType>
void factorize_CA_level_gpu(
    fmm::ParallelTree<CoordType, DataType>* tree,
    int current_level,
    gpu::ColorGpuDriver<CoordType, DataType, KernelType>& gpu_driver,
    int64_t& total_skeleton,
    int64_t& total_redundant,
    int64_t& local_max_skel,
    bool print_detail,
    int level_print_rank) {
    auto& level = tree->levels[current_level];
    const int dimension = tree->dimension;
    const int num_waves = 1 << dimension;

    std::vector<std::string> group_names;
    std::vector<std::vector<int64_t>> groups;
    std::vector<int64_t> interior;
    interior.reserve(level.interior_id.size());
    for (int64_t local_idx : level.interior_id) {
        interior.push_back(local_idx + level.local_morton_start);
    }
    if (level.num_active_processes == 1) {
        std::vector<int64_t> boundary;
        for (int64_t local_idx : level.boundary_id) {
            boundary.push_back(local_idx + level.local_morton_start);
        }
        group_names = {"boundary", "interior"};
        groups = {std::move(boundary), std::move(interior)};
    } else {
        std::vector<int64_t> blue_orange = level.blue;
        blue_orange.insert(blue_orange.end(), level.orange.begin(), level.orange.end());
        group_names.push_back("blue/orange");
        groups.push_back(std::move(blue_orange));
        if (dimension == 3) {
            group_names.push_back("purple");
            groups.push_back(level.purple);
        }
        group_names.push_back("green");
        groups.push_back(level.green);
        group_names.push_back("interior");
        groups.push_back(std::move(interior));
    }

    for (size_t group_idx = 0; group_idx < groups.size(); ++group_idx) {
        const auto& group = groups[group_idx];
        if (print_detail && tree->mpi_rank == level_print_rank) {
            std::cout << "  Processing CA " << group_names[group_idx] << " group ("
                      << group.size() << " boxes, device)..." << std::endl;
        }
        std::vector<std::vector<int64_t>> waves(static_cast<size_t>(num_waves));
        for (int64_t morton : group) {
            waves[static_cast<size_t>(morton & (num_waves - 1))].push_back(morton);
        }
        if (group_names[group_idx] == "interior") {
            clear_ghosts(level);
        }
        for (size_t wave_idx = 0; wave_idx < waves.size(); ++wave_idx) {
            const auto& wave = waves[wave_idx];
            if (wave.empty()) continue;
            const int32_t wave_sequence =
                static_cast<int32_t>(group_idx) * num_waves + static_cast<int32_t>(wave_idx);
            gpu_driver.eliminate_wave(wave, wave_sequence);
            for (int64_t morton : wave) {
                level.eliminated_boxes.insert(morton);
                level.elimination_wave[morton] = wave_sequence;
            }
        }
    }

    for (const auto& box : level.local_boxes) {
        total_skeleton += static_cast<int64_t>(box.skeleton_indices.size());
        total_redundant += static_cast<int64_t>(box.redundant_indices.size());
        local_max_skel = std::max<int64_t>(local_max_skel, box.skeleton_indices.size());
    }
}
#endif

/**
 * @brief Hierarchical factorization routine (parallel)
 * 
 * Performs bottom-up hierarchical factorization from leaf level to root
 * with MPI communication and process reduction.
 * 
 * Algorithm:
 * 1. For each level from leaf (N-1) down to 1:
 *    a) Gather ghost and assisting boxes (MPI communication)
 *    b) Eliminate boxes in colored order: blue → orange → purple (3D) → green → interior
 *    c) Gather assisting boxes post-elimination (MPI communication)
 *    d) Build parent level interactions
 *    e) Handle process reduction (send/receive parent boxes if needed)
 * 2. At level 1 → 0 transition:
 *    - No elimination at level 1 (only 4/8 boxes, no far-field)
 *    - Aggregate to level 0 (single box or few boxes)
 *    - Factorize level 0 schur complement for diagonal solve
 * 
 * @tparam CoordType Coordinate data type
 * @tparam DataType Matrix data type
 * @tparam KernelType Kernel evaluator type
 * @param tree Hierarchical tree structure
 * @param kernel Kernel evaluator
 * @param tolerance ID tolerance
 * @param is_symmetric Matrix symmetry flag
 * @param is_hermitian Matrix Hermitian flag
 * @param factorization_method Dense diagonal-block factorization method
 * @param unit_proxy_points Unit sphere proxy points
 * @param num_proxy Number of proxy points
 * @param proxy_radius Proxy sphere radius multiplier
 * @param verbosity -1: quiet, 0: summary, 1+: detailed progress
 */
template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_factorization_parallel(
    fmm::ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    int use_sketch,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
    int64_t* out_rankmax,
    size_t* memory_per_rank,
    int lazy_schur,
    int gemm_split,
    int ca_staged_halo,
    int ca_owner_component,
    int ca_owner_serial,
    int verbosity = 1) {

    // To Do: NEED TO FIX KERNEL!!!!!
    using clock = std::chrono::high_resolution_clock;
    
    int rank = tree->mpi_rank;
    int size = tree->mpi_size;
    const bool print_summary = verbosity >= 0;
    const bool print_detail = verbosity >= 1;
    const bool print_trace = verbosity >= 2;
    H2FactorizationMemoryDiagnostics memory_diagnostics;
#ifdef H2_HAVE_GPU
    // (a component-owner CA level needs a larger exchange arena: a fifth of
    // the GPU, 8 GB of 40; a colour's arrivals stay in it until its end)
    if (color_gpu_enabled()) gpu::begin_operator_build(tree, tree->comm, ca_owner_component == 3 ? 0.2 : 0.0);
    gpu::invalidate_device_solve();  // the device copies of the previous factors
#endif
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    FactorizationCommunicatorSet factorization_comms =
        make_factorization_communicators(tree);
    
    int dimension = tree->dimension;
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    const int num_children = morton::children_per_box(dimension);
    const int factorization_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);

    if (use_sketch < 0 || use_sketch > 2) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: use_sketch must be 0, 1, or 2");
    }
    if (lazy_schur < 0 || lazy_schur > 2) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: lazy_schur must be 0, 1, or 2");
    }
    if (lazy_schur != 0 && use_sketch != 2) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: lazy_schur requires use_sketch=2");
    }
    if (gemm_split < 0) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: gemm_split must be nonnegative");
    }
    if (ca_staged_halo != 0 && ca_staged_halo != 2) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: ca_staged_halo must be 0 or 2");
    }
    if (ca_owner_component != 0 && ca_owner_component != 3) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: ca_owner_component must be 0 or 3");
    }
    if (ca_owner_serial != 0 && ca_owner_serial != 1) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: ca_owner_serial must be 0 or 1");
    }
    if (ca_owner_component == 3 && lazy_schur == 0) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: ca_owner_component=3 requires lazy_schur=1 or 2");
    }

    struct FactorizationRuntimeReset {
        ~FactorizationRuntimeReset() {
            configure_color_factorization_runtime(false, 0, 0);
            configure_ca_factorization_runtime(0, 0, false, false);
        }
    } factorization_runtime_reset;
    owner_solve_records().clear();

    if (dynamic_threading.enabled &&
        !tree->levels[leaf_level].is_process_active) {
        park_inactive_rank_on_service_cpu(dynamic_threading);
    }
    
    if (print_summary && rank == factorization_header_rank) {
        const char* factorization_name = "Unknown";
        switch (factorization_method) {
            case FactorizationMethod::CHOLESKY:
                factorization_name = "Cholesky";
                break;
            case FactorizationMethod::LU:
                factorization_name = "LU";
                break;
            case FactorizationMethod::BUNCH_KAUFMAN:
                factorization_name = "Bunch-Kaufman";
                break;
            case FactorizationMethod::NONE:
                factorization_name = "None";
                break;
        }
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Factorization (Parallel)" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "MPI processes: " << size << std::endl;
        std::cout << "Levels: " << num_levels << std::endl;
        std::cout << "Leaf level: " << leaf_level << std::endl;
        std::cout << "Dimension: " << dimension << "D" << std::endl;
        std::cout << "Tolerance: " << tolerance << std::endl;
        std::cout << "Matrix property: " << (is_symmetric ? "Symmetric" : "Nonsymmetric") << std::endl;
        std::cout << "Factorization: " << factorization_name << std::endl;
        printf("max thread: %d\n", omp_get_max_threads());
        std::cout << "========================================\n" << std::endl;
    }
    
    // ===== Main factorization loop: leaf level down to level 1 =====
    auto total_time = clock::now();
    auto segment_start = clock::now();
    clock::duration total_data_exchange_time{};
    clock::duration total_reduction_time{};
    int64_t local_max_skel = 0;

    memory_diagnostics.record(tree, leaf_level, "setup", "factor_start");

#ifdef H2_HAVE_GPU
    // H2_use_gpu: the device box path of each level, and the blocks the
    // device transition keeps for the next one (h2_gpu/factorization_driver.hpp)
    typename gpu::ColorGpuDriver<CoordType, DataType, KernelType>::Options gpu_options;
    gpu_options.tolerance = tolerance;
    gpu_options.method = factorization_method;
    gpu_options.use_sketch = use_sketch;
    gpu_options.lazy_schur = lazy_schur;
    gpu_options.is_symmetric = is_symmetric;
    gpu_options.is_hermitian = is_hermitian;
    gpu_options.ca_owner_component = ca_owner_component;
    gpu::ColorGpuDriver<CoordType, DataType, KernelType> gpu_driver(tree, kernel, gpu_options);
    {
        const double kernels = gpu_driver.warm_up();
        const double buffers = gpu_driver.warm_up_buffers();
        const double exchanges = gpu_driver.warm_up_exchange();
        const double evaluator = gpu_driver.warm_up_evaluator();
        h2_warmup_seconds() = kernels + buffers + exchanges + evaluator;
        if (h2_warmup_seconds() > 0.0 && print_summary && rank == 0) {
            std::printf("  [gpu] warm-up: kernels %.2f s, buffers %.2f s, exchanges %.2f s, evaluator %.2f s, before "
                        "the levels (not in the factor time)\n",
                        kernels, buffers, exchanges, evaluator);
            std::fflush(stdout);
        }
    }
#endif

    for (int current_level = leaf_level; current_level >= 1; current_level--) {
        auto level_start = std::chrono::high_resolution_clock::now();
        clock::duration level_data_exchange{};
        clock::duration level_reduction{};
        // the level's host phases outside the elimination and transition
        // timers (ms): the owner schedule's setup, from the level's start to
        // the elimination, from the elimination to the transition, after it
        double host_phase_ms[4] = {0.0, 0.0, 0.0, 0.0};
        auto host_phase_since = [](std::chrono::high_resolution_clock::time_point t) {
            return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t).count();
        };
        auto elim_end_level = level_start;
        
        auto& level = tree->levels[current_level];
        auto& parent_level = tree->levels[current_level - 1];
        const bool participates_in_transition =
            level.is_process_active || parent_level.is_process_active;
        // BPACK_TRACE=phase: per-level elimination phase breakdown
        static const bool phase_report = fmm::env::trace("phase");
        if (phase_report) {
            FMM_PHASE_RESET();
#ifdef H2_HAVE_GPU
            gpu::owner_pass_stats() = gpu::OwnerPassStats{};
            gpu::eliminator_stats() = gpu::EliminatorStats{};
            transport_timers() = TransportTimers{};
#endif
        }
        if (!participates_in_transition) {
            continue;
        }

        MPI_Comm level_comm =
            factorization_comms.level[static_cast<size_t>(current_level)];
        MPI_Comm transition_comm =
            factorization_comms.transition[static_cast<size_t>(current_level)];
        const int level_print_rank = smallest_active_rank(level);
        LevelThreadPlan level_thread_plan;
        if (level.is_process_active) {
            level_thread_plan = configure_level_thread_plan(
                dynamic_threading, level_comm, level);
        }
        const bool CA_requested = tree->level_requests_CA(current_level);
        const bool use_CA_level = tree->level_uses_CA(current_level);
        const bool use_streamed_level =
            current_level > 1 && use_sketch == 2 &&
            is_symmetric && !is_hermitian;
        const int level_lazy_schur = use_streamed_level
            ? (use_CA_level ? std::min(lazy_schur, 1) : lazy_schur)
            : 0;
        configure_color_factorization_runtime(
            use_streamed_level, level_lazy_schur, gemm_split);
        configure_ca_factorization_runtime(
            use_CA_level ? ca_staged_halo : 0,
            use_CA_level ? ca_owner_component : 0,
            use_CA_level && ca_owner_serial != 0,
            use_CA_level && is_symmetric);

        OwnerScheduleState<CoordType, DataType> owner_schedule;
        StagedHaloState<CoordType, DataType> staged_state;
        bool staged_overlap_scheduling = false;
        // A replicated CA level the device may run (h2_gpu/GPU_CA_PLAN.md,
        // M1): its halo arrives whole before the elimination (no staged
        // overlap); start_ca_level confirms it after the gather.
        bool gpu_ca_candidate = false;
        bool ca_on_device = false;  // the device eliminates this CA level
        level.ca_home_factors = false;
#ifdef H2_HAVE_GPU
        gpu_ca_candidate = level.is_process_active && current_level > 1 &&
                           gpu_driver.ca_level_runs(current_level, ca_owner_component);
        gpu::OwnerPlanView owner_view;  // a component-owner level on the device (h2_gpu/ca_owner.hpp)
#endif
        if (level.is_process_active && use_CA_level &&
            current_level > 1 && level.num_active_processes > 1 &&
            ca_owner_component == 3) {
            if (!(is_symmetric && !is_hermitian)) {
                throw std::runtime_error(
                    "H2_CA_owner_component=3 requires a symmetric, non-Hermitian factorization");
            }
            if (!use_streamed_level ||
                lazy_far_field_mode() != LazyFarFieldMode::LAZY) {
                throw std::runtime_error(
                    "H2_CA_owner_component=3 requires streamed sketching and lazy Schur updates");
            }
            const auto t_setup = std::chrono::high_resolution_clock::now();
            owner_schedule_setup(
                tree, current_level, level_comm, owner_schedule,
                /*sync_colors=*/gpu_ca_candidate);  // (a device level: per-colour balancing)
            host_phase_ms[0] = host_phase_since(t_setup);
        }

        if (level.is_process_active && use_CA_level) {
            segment_start = clock::now();
            CAGatherTimes gather_times;
            FactorizationMemoryDiagnosticCallback gather_memory_diagnostic;
            if (memory_diagnostics.enabled()) {
                gather_memory_diagnostic =
                    [&](const char* phase, size_t pending, size_t communication) {
                        memory_diagnostics.record(
                            tree, current_level, "CA",
                            std::string("initial_gather_") + phase,
                            pending, 0, 0, communication);
                    };
            }
            if (ca_staged_halo == 2 && is_symmetric && !is_hermitian) {
                initiate_staged_halo_gather(
                    tree, current_level, staged_state, level_comm,
                    owner_schedule.active
                        ? &owner_schedule.block_ghosts
                        : nullptr);
                staged_halo_wait_stage(level, staged_state, 0);

                const bool overlap_ok =
                    current_level > 1 &&
                    lazy_far_field_mode() == LazyFarFieldMode::LAZY &&
                    !gpu_ca_candidate;
                if (overlap_ok) {
                    build_staged_stage_map(
                        level, level.staged_stage_map);
                    level.staged_unarrived.clear();
                    for (int stage = 1;
                         stage < STAGED_HALO_STAGES; ++stage) {
                        for (const auto& peer_mortons :
                             staged_state.recv_mortons[stage]) {
                            level.staged_unarrived.insert(
                                peer_mortons.begin(),
                                peer_mortons.end());
                        }
                    }
                    level.staged_pending = std::make_unique<
                        typename TreeLevel<CoordType, DataType>::
                            StagedOverlapPending>();
                    level.staged_overlap_on = true;
                    staged_overlap_scheduling = true;
                } else {
                    for (int stage = 1;
                         stage < STAGED_HALO_STAGES; ++stage) {
                        staged_halo_wait_stage(
                            level, staged_state, stage);
                    }
                    staged_halo_finish(staged_state);
                }
            } else {
                gather_CA_factorization_data(
                    tree, current_level, level_comm,
                    is_symmetric && !is_hermitian,
                    gather_memory_diagnostic, &gather_times);
            }
            const auto gather_duration = clock::now() - segment_start;
            level_data_exchange += gather_duration;
            // (the messages wait for the slowest rank's plan and packing)
            double gather_prepare_max = gather_times.plan + gather_times.serialize;
            if (print_detail && ca_staged_halo != 2) {
                MPI_Allreduce(MPI_IN_PLACE, &gather_prepare_max, 1, MPI_DOUBLE, MPI_MAX, level_comm);
            }
            if (print_detail && rank == level_print_rank) {
                std::cout << "  CA initial ghost/assisting gather time: "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(
                                 gather_duration).count()
                          << " ms";
                if (ca_staged_halo != 2) {
                    auto ms = [](double s) { return static_cast<long long>(s * 1e3 + 0.5); };
                    if (gather_times.geometric) {
                        std::cout << " (geometric plan " << ms(gather_times.plan) << ", serialize "
                                  << ms(gather_times.serialize) << " [both, slowest rank: "
                                  << ms(gather_prepare_max) << "], messages " << ms(gather_times.mpi)
                                  << ", unpack " << ms(gather_times.unpack) << ", symmetric halo "
                                  << ms(gather_times.reconstruct) << ")";
                    } else {
                        std::cout << " (plan " << ms(gather_times.plan) << " [first collective "
                                  << ms(gather_times.first_collective) << "], request sets "
                                  << ms(gather_times.request_sets) << ", payloads " << ms(gather_times.payloads)
                                  << ", symmetric halo " << ms(gather_times.reconstruct) << ", assisting "
                                  << ms(gather_times.assisting) << ")";
                    }
                }
                std::cout << std::endl;
            }
        }

        if (level.is_process_active) {
            kernel->register_level_coordinates(level);
        }

        memory_diagnostics.record(
            tree, current_level, use_CA_level ? "CA" : "color",
            level.is_process_active && use_CA_level ? "post_CA_gather" :
                                                       "level_start");

        // declare locks
        std::unordered_map<int64_t, omp_lock_t*> box_locks;

         // initialize locks
        if (level.is_process_active) {
            for (auto& box : level.local_boxes) {
                auto* lock = new omp_lock_t;
                omp_init_lock(lock);
                box_locks[box.morton_index] = lock;
            }
            if (use_CA_level) {
                for (auto& box : level.ghost_boxes) {
                    auto* lock = new omp_lock_t;
                    omp_init_lock(lock);
                    box_locks[box.morton_index] = lock;
                }
            }
            
            {
                auto* global_lock1 = new omp_lock_t;
                omp_init_lock(global_lock1);
                box_locks[-1] = global_lock1;
                auto* global_lock2 = new omp_lock_t;
                omp_init_lock(global_lock2);
                box_locks[-2] = global_lock2;
            }
            level.box_locks = std::move(box_locks);
        }
        
        if (print_summary && rank == level_print_rank) {
            std::cout << "\n===== Level " << current_level << " =====" << std::endl;
            std::cout << "Active processes: " << level.num_active_processes << " from rank: " << rank << std::endl;
            std::cout << "Boxes per process: " << level.num_boxes_local << std::endl;
        }
        print_level_thread_plan(
            dynamic_threading,
            current_level,
            "factorization",
            level_print_rank,
            rank,
            level_thread_plan,
            print_detail);
        if (print_detail && rank == level_print_rank && current_level > 1) {
            std::cout << "  Level algorithm: "
                      << (use_CA_level ? "CA" : "color");
            if (CA_requested && !use_CA_level) {
                std::cout << " (automatic CA fallback: "
                          << tree->boxes_per_active_process(current_level)
                          << " boxes per active process; distributed CA requires at least "
                          << tree->minimum_CA_boxes_per_active_process() << ")";
            }
            std::cout << std::endl;
        }
        
        // ===== Step 1: Gather ghost and assisting boxes =====
        

        
        // ===== Step 2: Eliminate boxes in colored order =====
        std::chrono::milliseconds elim_duration{0};
        int64_t total_skeleton = 0;
        int64_t total_redundant = 0;
        
        if (current_level > 1 && level.is_process_active) {
            // Regular levels: do ID and elimination in colored order

            auto elim_start = std::chrono::high_resolution_clock::now();
            host_phase_ms[1] = std::chrono::duration<double, std::milli>(elim_start - level_start).count();
            static const bool skew_trace = fmm::env::trace("skew");
            LevelSkewProbe skew_probe;
            skew_probe.level_begin(skew_trace);

#ifdef H2_HAVE_GPU
            if (use_CA_level && gpu_ca_candidate && ca_owner_component == 3 && owner_schedule.active &&
                gpu::start_owner_level_gpu(tree, current_level, gpu_driver, owner_schedule, owner_view, level_comm,
                                           print_detail && rank == level_print_rank)) {
                // a component-owner level on the device: its factors land at
                // their home ranks, so the solves run it as a Color level
                gpu::factorize_owner_level_gpu(
                    tree, current_level, gpu_driver, owner_schedule, owner_view,
                    total_skeleton, total_redundant, local_max_skel,
                    print_detail, level_print_rank);
                PendingFactorUpdates<DataType> no_updates;
                gpu_driver.finish_level(current_level, no_updates, level_comm,
                                        print_detail && rank == level_print_rank);
                level.ca_home_factors = true;
                ca_on_device = true;
            } else if (use_CA_level && gpu_ca_candidate && ca_owner_component == 0 &&
                       gpu_driver.start_ca_level(current_level, level_comm, print_detail && rank == level_print_rank)) {
                factorize_CA_level_gpu(
                    tree, current_level, gpu_driver,
                    total_skeleton, total_redundant, local_max_skel,
                    print_detail, level_print_rank);
                // the factors' host copies (the solve's, ghosts' too), as a Color level
                PendingFactorUpdates<DataType> no_updates;
                gpu_driver.finish_level(current_level, no_updates, level_comm,
                                        print_detail && rank == level_print_rank);
                gpu::device_solve_store().ca_levels.insert(current_level);  // the device solve runs it (M4b)
                ca_on_device = true;
            }
#endif
            if (use_CA_level && !ca_on_device) {
                factorize_CA_level(
                    tree, current_level, kernel, tolerance,
                    use_sketch,
                    is_symmetric, is_hermitian, factorization_method,
                    unit_proxy_points, num_proxy, proxy_radius,
                    total_skeleton, total_redundant, local_max_skel,
                    print_detail, level_print_rank, &memory_diagnostics,
                    &owner_schedule, &staged_state,
                    staged_overlap_scheduling);
            } else if (!use_CA_level) {
                const int num_colors = 1 << dimension;

            // ----------------------------------------------------------------
            // Build boundary color bins [0, num_colors-1]
            // Build interior sub-wave bins [num_colors, 2*num_colors-1]
            // Both use morton & (num_colors-1) coloring, guaranteeing that
            // boxes within the same wave are non-adjacent.
            // ----------------------------------------------------------------
            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx  = level.local_morton_start + local_idx;
                const auto&   box         = level.local_boxes[static_cast<size_t>(local_idx)];
                const int     color_id    = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }

            // Append interior sub-waves after boundary waves
            const int interior_start_loc = static_cast<int>(color_bins.size());  // == num_colors
            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }
            // color_bins layout:
            //   [0,           num_colors-1]: boundary waves
            //   [num_colors, 2*num_colors-1]: interior sub-waves

            PendingFactorUpdates<DataType> pending_updates;
            int boundary_count = 0;
            // After a wave, the assisting boxes (other ranks') of its color:
            // boundary ones in boundary waves, interior ones in interior
            // sub-waves.
            auto mark_assisting_eliminated = [&](int color, bool interior) {
                for (const auto& kv : level.assisting_box_points_for_kernel_evaluation) {
                    const int  assisting_color    = static_cast<int>(kv.first & (num_colors - 1));
                    const bool is_boundary_assist = level.assisting_boxes[kv.second].on_boundary;
                    if (assisting_color == color && is_boundary_assist != interior) {
                        level.eliminated_boxes.insert(kv.first);
                    }
                }
            };
            bool to_store = true;
            const bool store_interior_wave = true;

#ifdef H2_HAVE_GPU
            // H2_use_gpu=1: the device box path takes the whole level when it
            // covers it (see h2_gpu/level_eliminator.hpp); otherwise the
            // host box region runs and only the owner pass uses the GPU.
            gpu_driver.start_level(current_level, use_streamed_level, level_comm,
                                   print_detail && rank == level_print_rank);
            const bool gpu_level_local = gpu_driver.local();
#else
            const bool gpu_level_local = false;
#endif

            skew_probe.loop_begin();
            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {

                const int  color_id_mod    = counter % num_colors;
                const bool is_interior     = (counter >= interior_start_loc);
                const LevelSkewProbe::WaveScope skew_wave(
                    skew_probe, (!is_interior || counter == interior_start_loc) && !gpu_level_local);
                FMM_PHASE_LAP_BEGIN(wave_lap);  // master-thread wall laps of this wave

                // ----------------------------------------------------------------
                // Communication / transport step (single-threaded)
                // ----------------------------------------------------------------
                // (a single-rank device level has nothing to transport)
                if ((!is_interior || counter == interior_start_loc) && !gpu_level_local) {
                    if (memory_diagnostics.enabled()) {
                        memory_diagnostics.record(
                            tree, current_level, "color",
                            "wave" + std::to_string(counter) + "_pre_transport",
                            h2_diag_pending_bytes(pending_updates));
                    }
                    FactorizationMemoryDiagnosticCallback transport_memory_diagnostic;
                    if (memory_diagnostics.enabled()) {
                        transport_memory_diagnostic =
                            [&](const char* phase, size_t additional_pending,
                                size_t communication) {
                                memory_diagnostics.record(
                                    tree, current_level, "color",
                                    "wave" + std::to_string(counter) +
                                        "_transport_" + phase,
                                    h2_diag_pending_bytes(pending_updates) +
                                        additional_pending,
                                    0, 0, communication);
                            };
                    }
                    std::chrono::high_resolution_clock::duration comm_duration_raw{};
#ifdef H2_HAVE_GPU
                    if (gpu_driver.on()) {
                        // device level on several ranks
                        comm_duration_raw = gpu_driver.transport(
                            pending_updates, [&](std::vector<int64_t>* installed) {
                                return transport_and_apply_factor_updates_symmetric_onehop(
                                    tree, current_level, kernel, pending_updates, false,
                                    transport_memory_diagnostic, false, installed);
                            });
                    } else
#endif
                    {
                        comm_duration_raw =
                            transport_and_apply_factor_updates_symmetric_onehop(
                                tree, current_level, kernel, pending_updates, false,
                                transport_memory_diagnostic);
                    }
                    level_data_exchange += comm_duration_raw;
                    update_neighbor_slicing_for_level(level, is_symmetric);
                    if (memory_diagnostics.enabled()) {
                        memory_diagnostics.record(
                            tree, current_level, "color",
                            "wave" + std::to_string(counter) + "_post_transport",
                            h2_diag_pending_bytes(pending_updates));
                    }
                    auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_duration_raw);

                    if (print_detail && rank == level_print_rank) {
                        std::cout << "  Comm time: "
                                << comm_duration.count()
                                << " ms" << std::endl;
                    }
                }
                

              

                const auto& color_list = color_bins[static_cast<size_t>(counter)];

                if (print_detail && rank == level_print_rank) {
                    std::cout << "  Processing "
                            << (is_interior ? "interior sub-wave " : "boundary color ")
                            << color_id_mod
                            << " (" << color_list.size() << " boxes)..." << std::endl;
                }

                FMM_PHASE_LAP(wave_lap, WALL_STAGED);  // transport
                if (color_list.empty()) continue;
                FMM_PHASE_NOTE_WAVE(color_list.size());

                const bool enable_deferred_xnn =
                    to_store && (store_interior_wave || !is_interior);
                const bool use_owner_deferred_xnn =
                    enable_deferred_xnn && is_symmetric && !is_hermitian;

                std::unordered_set<int64_t> wave_box_set;
                wave_box_set.insert(color_list.begin(), color_list.end());

                for (int64_t morton_idx : color_list) {
                    BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
                    if (box_ptr == nullptr) {
                        throw std::runtime_error(
                            "Morton index " + std::to_string(morton_idx) +
                            " not found at level " + std::to_string(current_level) +
                            " (counter: " + std::to_string(counter) + ")");
                    }

                    for (int64_t neighbor_morton : box_ptr->one_hop) {
                        if (neighbor_morton != morton_idx &&
                            wave_box_set.count(neighbor_morton) != 0) {
                            std::ostringstream oss;
                            oss << "Wave construction error: one-hop conflict"
                                << " level=" << current_level
                                << " color_id_mod=" << color_id_mod
                                << " counter=" << counter
                                << " box=" << morton_idx
                                << " neighbor=" << neighbor_morton;
                            throw std::runtime_error(oss.str());
                        }
                    }
                }

#ifdef H2_HAVE_GPU
                if (gpu_driver.on()) {
                    // Box region, owner pass, mirror and finalize on the
                    // device.  One active rank: nothing is pending for other
                    // ranks and there are no assisting boxes to mark.
                    boundary_count += gpu_driver.eliminate_wave(color_list, counter);
                    for (int64_t morton_idx : color_list) {
                        level.eliminated_boxes.insert(morton_idx);
                        level.elimination_wave[morton_idx] =
                            static_cast<int32_t>(counter);
                    }
                    mark_assisting_eliminated(color_id_mod, is_interior);
                    FMM_PHASE_LAP(wave_lap, WALL_BOX);
                    FMM_PHASE_WAVE_BOX_DONE();
                    continue;
                }
#endif
                const int max_threads = std::max(1, omp_get_max_threads());
                std::vector<PendingFactorUpdates<DataType>> thread_pending(
                    static_cast<size_t>(max_threads));
                std::vector<std::vector<int64_t>> thread_xnn_candidate_boxes;
                std::vector<size_t> candidate_box_offsets;
                std::vector<int64_t> wave_xnn_candidate_boxes;
                std::vector<std::vector<DeferredXnnTargetKey>> wave_xnn_mirror_targets;
                std::vector<int> thread_boundary_counts(static_cast<size_t>(max_threads), 0);

                if (use_owner_deferred_xnn) {
                    thread_xnn_candidate_boxes.resize(static_cast<size_t>(max_threads));
                    candidate_box_offsets.resize(thread_xnn_candidate_boxes.size() + 1, 0);
                }

                std::exception_ptr wave_exception;
                std::mutex wave_exception_mutex;
                std::atomic<bool> wave_failed{false};
                size_t wave_scratch_bytes = 0;

                const int wave_team = std::max(1, omp_get_max_threads());
                const int wave_split = split_threads_for(
                    static_cast<int64_t>(color_list.size()), wave_team);

                #pragma omp parallel default(shared)
                {
                    const int tid = omp_get_thread_num();
                    FactorizationThreadScratch<CoordType, DataType> scratch;
                    scratch.split_threads = wave_split;

                    #pragma omp for schedule(dynamic)
                    for (int64_t bi = 0; bi < static_cast<int64_t>(color_list.size()); ++bi) {
                        if (wave_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }

                        try {
                            const int64_t morton_idx = color_list[static_cast<size_t>(bi)];
                            BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
                            if (box_ptr == nullptr) {
                                throw std::runtime_error(
                                    "Morton index " + std::to_string(morton_idx) +
                                    " not found at level " + std::to_string(current_level) +
                                    " (counter: " + std::to_string(counter) + ")");
                            }

                            auto& box = *box_ptr;
                            FMM_PHASE_BOX_SCOPE();
                            FMM_PHASE_LAP_BEGIN(sketch_lap);

                            if (use_streamed_level) {
                                gather_id_target_streamed(
                                    tree, &box, level, kernel, scratch,
                                    box.on_boundary, tolerance);
                            } else {
                                scratch.streamed_sketch_valid = false;
                                gather_id_workspace(
                                    tree,
                                    &box, level, kernel, tolerance,
                                    unit_proxy_points.data(), num_proxy,
                                    proxy_radius, is_symmetric,
                                    scratch.workspace, scratch.workspace_rows,
                                    scratch.workspace_cols, 0,
                                    box.on_boundary, false,
                                    &scratch.id_adaptive);
                            }

                            thread_boundary_counts[static_cast<size_t>(tid)] += box.on_boundary;
                            FMM_PHASE_LAP(sketch_lap, BOX_SKETCH);

                            compute_and_modify(dimension,
                                &box, level, kernel,
                                scratch,
                                tolerance, use_sketch, is_symmetric, is_hermitian,
                                &thread_pending[tid],
                                factorization_method, use_owner_deferred_xnn
                            );

                        } catch (...) {
                            if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(wave_exception_mutex);
                                wave_exception = std::current_exception();
                            }
                        }
                    }

                    if (memory_diagnostics.enabled()) {
                        const size_t thread_scratch_bytes =
                            h2_diag_factorization_scratch_bytes(scratch);
                        #pragma omp atomic update
                        wave_scratch_bytes += thread_scratch_bytes;
                        #pragma omp barrier
                        #pragma omp single
                        {
                            size_t all_pending_bytes =
                                h2_diag_pending_bytes(pending_updates);
                            for (const auto& thread_updates : thread_pending) {
                                all_pending_bytes +=
                                    h2_diag_pending_bytes(thread_updates);
                            }
                            memory_diagnostics.record(
                                tree, current_level, "color",
                                "wave" + std::to_string(counter) + "_compute",
                                all_pending_bytes, wave_scratch_bytes);
                        }
                    }
                }

                if (wave_exception) {
                    std::rethrow_exception(wave_exception);
                }
                FMM_PHASE_LAP(wave_lap, WALL_BOX);
                FMM_PHASE_WAVE_BOX_DONE();

                for (int local_boundary_count : thread_boundary_counts) {
                    boundary_count += local_boundary_count;
                }

                if (use_owner_deferred_xnn) {
                    for (int64_t morton_idx : color_list) {
                        level.eliminated_boxes.insert(morton_idx);
                        level.elimination_wave[morton_idx] =
                            static_cast<int32_t>(counter);
                    }

                    slice_far_field_blocks(level, is_symmetric, is_hermitian);

                    std::exception_ptr candidate_exception;
                    std::mutex candidate_exception_mutex;
                    std::atomic<bool> candidate_failed{false};

                    #pragma omp parallel default(shared) if (color_list.size() > 1)
                    {
                        const int tid = omp_get_thread_num();
                        auto& local_candidates =
                            thread_xnn_candidate_boxes[static_cast<size_t>(tid)];

                        #pragma omp for schedule(static)
                        for (int64_t bi = 0; bi < static_cast<int64_t>(color_list.size()); ++bi) {
                            if (candidate_failed.load(std::memory_order_relaxed)) {
                                continue;
                            }

                            try {
                                const int64_t morton_idx = color_list[static_cast<size_t>(bi)];
                                BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
                                if (box_ptr == nullptr) {
                                    continue;
                                }

                                collect_owner_deferred_xnn_candidates_for_source_box(
                                    box_ptr,
                                    level,
                                    wave_box_set,
                                    local_candidates);
                            } catch (...) {
                                if (!candidate_failed.exchange(true, std::memory_order_relaxed)) {
                                    std::lock_guard<std::mutex> lock(candidate_exception_mutex);
                                    candidate_exception = std::current_exception();
                                }
                            }
                        }
                    }

                    if (candidate_exception) {
                        std::rethrow_exception(candidate_exception);
                    }

                    for (size_t tid = 0; tid < thread_xnn_candidate_boxes.size(); ++tid) {
                        candidate_box_offsets[tid + 1] =
                            candidate_box_offsets[tid] + thread_xnn_candidate_boxes[tid].size();
                    }

                    wave_xnn_candidate_boxes.resize(candidate_box_offsets.back());
                    for (int tid_idx = 0; tid_idx < static_cast<int>(thread_xnn_candidate_boxes.size()); ++tid_idx) {
                        size_t out_idx = candidate_box_offsets[static_cast<size_t>(tid_idx)];
                        for (int64_t candidate_morton : thread_xnn_candidate_boxes[static_cast<size_t>(tid_idx)]) {
                            wave_xnn_candidate_boxes[out_idx++] = candidate_morton;
                        }
                    }

                    std::sort(
                        wave_xnn_candidate_boxes.begin(),
                        wave_xnn_candidate_boxes.end());
                    wave_xnn_candidate_boxes.erase(
                        std::unique(
                            wave_xnn_candidate_boxes.begin(),
                            wave_xnn_candidate_boxes.end()),
                        wave_xnn_candidate_boxes.end());
                    wave_xnn_mirror_targets.resize(wave_xnn_candidate_boxes.size());

                    std::exception_ptr owner_exception;
                    std::mutex owner_exception_mutex;
                    std::atomic<bool> owner_failed{false};
                    size_t owner_scratch_bytes = 0;

                    FMM_PHASE_LAP(wave_lap, WALL_CANDIDATES);
                    const int owner_team = std::max(1, omp_get_max_threads());
                    const int owner_split = split_threads_for(
                        static_cast<int64_t>(wave_xnn_candidate_boxes.size()),
                        owner_team);

                    bool owner_on_gpu = false;
#ifdef H2_HAVE_GPU
                    if (color_gpu_enabled()) {
                        gpu::run_owner_pass(
                            wave_xnn_candidate_boxes, level, kernel, wave_box_set,
                            wave_xnn_mirror_targets, thread_pending, owner_split);
                        owner_on_gpu = true;
                    }
#endif
                    if (!owner_on_gpu) {
                        #pragma omp parallel default(shared)
                        {
                            const int tid = omp_get_thread_num();
                            DeferredXnnOwnerScratch<DataType> scratch;
                            scratch.split_threads = owner_split;

                            #pragma omp for schedule(dynamic)
                            for (int64_t idx = 0; idx < static_cast<int64_t>(wave_xnn_candidate_boxes.size()); ++idx) {
                                if (owner_failed.load(std::memory_order_relaxed)) {
                                    continue;
                                }

                                try {
                                    const int64_t candidate_morton =
                                        wave_xnn_candidate_boxes[static_cast<size_t>(idx)];

                                    // Deferred store=true replay:
                                    //   local-local   -> local owner replay + local mirror
                                    //   local-remote  -> local replay + remote transport
                                    //   remote-remote -> remote transport only
                                    apply_owner_deferred_xnn_updates_for_candidate_box(
                                        candidate_morton,
                                        level,
                                        kernel,
                                        wave_box_set,
                                        scratch,
                                        wave_xnn_mirror_targets[static_cast<size_t>(idx)],
                                        &thread_pending[static_cast<size_t>(tid)]);
                                } catch (...) {
                                    if (!owner_failed.exchange(true, std::memory_order_relaxed)) {
                                        std::lock_guard<std::mutex> lock(owner_exception_mutex);
                                        owner_exception = std::current_exception();
                                    }
                                }
                            }

                            if (memory_diagnostics.enabled()) {
                                const size_t thread_scratch_bytes =
                                    h2_diag_deferred_owner_scratch_bytes(scratch);
                                #pragma omp atomic update
                                owner_scratch_bytes += thread_scratch_bytes;
                                #pragma omp barrier
                                #pragma omp single
                                {
                                    size_t all_pending_bytes =
                                        h2_diag_pending_bytes(pending_updates);
                                    for (const auto& thread_updates : thread_pending) {
                                        all_pending_bytes +=
                                            h2_diag_pending_bytes(thread_updates);
                                    }
                                    memory_diagnostics.record(
                                        tree, current_level, "color",
                                        "wave" + std::to_string(counter) + "_owner",
                                        all_pending_bytes, owner_scratch_bytes);
                                }
                            }
                        }
                    }

                    if (owner_exception) {
                        std::rethrow_exception(owner_exception);
                    }
                    FMM_PHASE_LAP(wave_lap, WALL_OWNER);

                    std::exception_ptr mirror_exception;
                    std::mutex mirror_exception_mutex;
                    std::atomic<bool> mirror_failed{false};

                    #pragma omp parallel default(shared) if (wave_xnn_candidate_boxes.size() > 1)
                    {
                        #pragma omp for schedule(static)
                        for (int64_t idx = 0; idx < static_cast<int64_t>(wave_xnn_candidate_boxes.size()); ++idx) {
                            if (mirror_failed.load(std::memory_order_relaxed)) {
                                continue;
                            }

                            try {
                                const int64_t candidate_morton =
                                    wave_xnn_candidate_boxes[static_cast<size_t>(idx)];
                                BoxData<CoordType, DataType>* candidate_box =
                                    level.find_local_box(candidate_morton);
                                if (candidate_box == nullptr) {
                                    candidate_box = level.find_ghost_box(candidate_morton);
                                }
                                if (candidate_box == nullptr) {
                                    continue;
                                }

                                apply_symmetric_owner_deferred_xnn_updates_for_candidate_box(
                                    candidate_box,
                                    wave_xnn_mirror_targets[static_cast<size_t>(idx)],
                                    level);
                            } catch (...) {
                                if (!mirror_failed.exchange(true, std::memory_order_relaxed)) {
                                    std::lock_guard<std::mutex> lock(mirror_exception_mutex);
                                    mirror_exception = std::current_exception();
                                }
                            }
                        }
                    }

                    if (mirror_exception) {
                        std::rethrow_exception(mirror_exception);
                    }
                    FMM_PHASE_LAP(wave_lap, WALL_MIRROR);

                    share_symmetric_level_edges(level);

                    // Package the generator while X_NR still contains the
                    // original one-hop coupling. Finalization replaces X_NR
                    // with temp2, which remains the representation retained
                    // by local and installed generator boxes.
                    if (lazy_far_field_mode() == LazyFarFieldMode::LAZY) {
                        emit_lazy_generators(
                            level, color_list, counter, pending_updates);
                    }

                    std::exception_ptr finalize_exception;
                    std::mutex finalize_exception_mutex;
                    std::atomic<bool> finalize_failed{false};

                    #pragma omp parallel default(shared) if (color_list.size() > 1)
                    {
                        #pragma omp for schedule(static)
                        for (int64_t bi = 0; bi < static_cast<int64_t>(color_list.size()); ++bi) {
                            if (finalize_failed.load(std::memory_order_relaxed)) {
                                continue;
                            }

                            try {
                                const int64_t source_morton = color_list[static_cast<size_t>(bi)];
                                BoxData<CoordType, DataType>* source_box = level.find_local_box(source_morton);
                                if (source_box == nullptr) {
                                    continue;
                                }

                                finalize_deferred_xnn_source_box(source_box);
                            } catch (...) {
                                if (!finalize_failed.exchange(true, std::memory_order_relaxed)) {
                                    std::lock_guard<std::mutex> lock(finalize_exception_mutex);
                                    finalize_exception = std::current_exception();
                                }
                            }
                        }
                    }

                    if (finalize_exception) {
                        std::rethrow_exception(finalize_exception);
                    }
                    FMM_PHASE_LAP(wave_lap, WALL_FINALIZE);

                } else {
                    slice_far_field_blocks(level, is_symmetric, is_hermitian);
                }

                for (int t = 0; t < max_threads; ++t) {
                    merge_pending(pending_updates, thread_pending[t]);
                    clear_pending_factor_updates_memory(thread_pending[t]);
                }

                if (!use_owner_deferred_xnn) {
                    for (size_t bi = 0; bi < color_list.size(); ++bi) {
                        const int64_t morton_idx = color_list[bi];
                        level.eliminated_boxes.insert(morton_idx);
                        level.elimination_wave[morton_idx] =
                            static_cast<int32_t>(counter);
                    }
                }

                mark_assisting_eliminated(color_id_mod, is_interior);

                if (memory_diagnostics.enabled()) {
                    memory_diagnostics.record(
                        tree, current_level, "color",
                        "wave" + std::to_string(counter) + "_complete",
                        h2_diag_pending_bytes(pending_updates));
                }

            }

            // Keep the final flush outside the wave body so empty trailing bins
            // cannot strand updates produced by the last nonempty wave.
            if (memory_diagnostics.enabled()) {
                memory_diagnostics.record(
                    tree, current_level, "color", "final_pre_transport",
                    h2_diag_pending_bytes(pending_updates));
            }
            FactorizationMemoryDiagnosticCallback final_transport_memory_diagnostic;
            if (memory_diagnostics.enabled()) {
                final_transport_memory_diagnostic =
                    [&](const char* phase, size_t additional_pending,
                        size_t communication) {
                        memory_diagnostics.record(
                            tree, current_level, "color",
                            std::string("final_transport_") + phase,
                            h2_diag_pending_bytes(pending_updates) +
                                additional_pending,
                            0, 0, communication);
                    };
            }
            if (!gpu_level_local) {
                const LevelSkewProbe::WaveScope skew_wave(skew_probe, true);
                std::chrono::high_resolution_clock::duration final_comm_duration{};
#ifdef H2_HAVE_GPU
                if (gpu_driver.on()) {
                    final_comm_duration = gpu_driver.transport(
                        pending_updates, [&](std::vector<int64_t>* installed) {
                            return transport_and_apply_factor_updates_symmetric_onehop(
                                tree, current_level, kernel, pending_updates, false,
                                final_transport_memory_diagnostic, false, installed);
                        });
                } else
#endif
                {
                    final_comm_duration =
                        transport_and_apply_factor_updates_symmetric_onehop(
                            tree, current_level, kernel, pending_updates, false,
                            final_transport_memory_diagnostic);
                }
                level_data_exchange += final_comm_duration;
                update_neighbor_slicing_for_level(level, is_symmetric);
            }
            if (memory_diagnostics.enabled()) {
                memory_diagnostics.record(
                    tree, current_level, "color", "final_post_transport",
                    h2_diag_pending_bytes(pending_updates));
            }
#ifdef H2_HAVE_GPU
            // after the final transport: its generators' host copies
            gpu_driver.finish_level(current_level, pending_updates, level_comm,
                                    print_detail && rank == level_print_rank);
#endif
            
            for (const auto& box : level.local_boxes) {
                total_skeleton += box.skeleton_indices.size();
                total_redundant += box.redundant_indices.size();
                local_max_skel = std::max<int64_t>(local_max_skel, box.skeleton_indices.size());
            }
            }

            const auto elim_end = std::chrono::high_resolution_clock::now();
            elim_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                elim_end - elim_start);
            elim_end_level = elim_end;
            if (skew_probe.on) {
                level_skew_report(skew_probe, tree, level, current_level, level_comm,
                                  std::chrono::duration<double, std::milli>(elim_end - elim_start).count());
            }
        } else {
            // Level 1: Skip elimination (only 4/8 boxes, no far-field)
            if (print_detail && rank == level_print_rank) {
                std::cout << "  Skipping elimination at level 1 (final coarsening step)" << std::endl;
            }
        }

        memory_diagnostics.record(
            tree, current_level, use_CA_level ? "CA" : "color",
            "post_elimination");

        if (current_level > 1 && level.is_process_active) {
            double min_elim_ms = 0.0;
            double max_elim_ms = 0.0;
            reduce_active_duration_bounds_ms(
                level_comm,
                0,
                true,
                elim_duration,
                min_elim_ms,
                max_elim_ms);
            
            const int64_t local_box_count =
                static_cast<int64_t>(level.local_boxes.size());
            const int64_t local_compression_counts[3] = {
                total_skeleton,
                total_redundant,
                local_box_count
            };
            int64_t global_compression_counts[3] = {0, 0, 0};
            MPI_Reduce(
                local_compression_counts,
                global_compression_counts,
                3,
                MPI_INT64_T,
                MPI_SUM,
                0,
                level_comm);

            if (print_detail && rank == level_print_rank) {
                std::cout << "  Elimination time: shortest=" << std::llround(min_elim_ms)
                          << " ms, longest=" << std::llround(max_elim_ms) << " ms" << std::endl;

                const int64_t global_skeleton = global_compression_counts[0];
                const int64_t global_redundant = global_compression_counts[1];
                const int64_t global_box_count = global_compression_counts[2];
                const int64_t global_original =
                    global_skeleton + global_redundant;
                const double compression = global_original > 0
                    ? static_cast<double>(global_skeleton) /
                        static_cast<double>(global_original)
                    : 0.0;
                std::cout << "  Compression ratio: " << compression
                          << " (" << (compression * 100) << "%)" << std::endl;
                std::cout << "  Average skeleton size: "
                          << (global_box_count > 0
                              ? static_cast<double>(global_skeleton) /
                                  static_cast<double>(global_box_count)
                              : 0.0)
                          << std::endl;
            }
            if (phase_report) {
                int level_rank = 0;
                MPI_Comm_rank(level_comm, &level_rank);
                FMM_PHASE_REPORT(level_comm, 0, true, level_rank, current_level,
                                 omp_get_max_threads(), print_detail);
#ifdef H2_HAVE_GPU
                if (color_gpu_enabled() && print_detail && level_rank == 0 &&
                    !gpu_driver.report_level(current_level)) {
                    // a host level: the GPU owner pass
                    const auto& g = gpu::owner_pass_stats();
                    std::printf(
                        "  [gpu] level %d owner pass (rank 0): record %.2f s, pack %.2f s, "
                        "upload %.2f s (%.2f GB), gemm %.2f s (%lld gemms, %lld launches, %lld chunks), "
                        "download %.2f s (%.2f GB), write-back %.2f s\n",
                        current_level, g.record, g.pack, g.transfer_in, g.bytes_in / 1e9, g.gemm,
                        static_cast<long long>(g.gemms), static_cast<long long>(g.launches),
                        static_cast<long long>(g.chunks), g.transfer_out, g.bytes_out / 1e9, g.write_back);
                    std::fflush(stdout);
                }
#endif
            }
        }

        if (current_level > 1) {
            // (a device owner level: its factors are at home, no owner solve)
            const bool owner_engine_level =
                use_CA_level && ca_owner_component == 3 &&
                level.num_active_processes > 1 && !ca_on_device;
            owner_solve_record_level(
                owner_schedule, current_level, owner_engine_level);
            if (level.is_process_active) {
                if (print_detail && rank == level_print_rank &&
                    owner_schedule.active) {
                    owner_schedule_print_timing(
                        owner_schedule, current_level);
                }
                level.pair_interest_filter = nullptr;
                level.sketch_eliminated_filter = nullptr;
            }
        }

        if (level.is_process_active && staged_overlap_scheduling) {
            if (level.staged_pending != nullptr &&
                (!level.staged_pending->deltas.empty() ||
                 !level.staged_pending->mirrors.empty())) {
                throw std::runtime_error(
                    "staged halo overlap: pending updates not drained at level end");
            }
            level.staged_overlap_on = false;
            level.staged_unarrived.clear();
            level.staged_stage_map.clear();
            level.staged_pending.reset();
        }
        if (level.is_process_active && staged_state.active) {
            if (print_detail && rank == level_print_rank) {
                staged_halo_print_timing(
                    staged_state, current_level);
            }
            staged_halo_finish(staged_state);
        }
        
        // ===== Step 3: Gather assisting boxes post-elimination =====
        if (current_level > 1 && level.is_process_active && use_CA_level) {
            segment_start = clock::now();
            FactorizationMemoryDiagnosticCallback assisting_memory_diagnostic;
            if (memory_diagnostics.enabled()) {
                assisting_memory_diagnostic =
                    [&](const char* phase, size_t pending, size_t communication) {
                        memory_diagnostics.record(
                            tree, current_level, "CA",
                            std::string("post_elimination_assisting_") + phase,
                            pending, 0, 0, communication);
                    };
            }
            // (only the skeletons are new: over the assisting gather's plan)
            if (!refresh_CA_assisting_skeletons(
                    tree, current_level, level_comm, 412)) {
                gather_CA_assisting_boxes_factorization(
                    tree, current_level, level_comm, 412,
                    assisting_memory_diagnostic);
            }
            const auto gather_duration = clock::now() - segment_start;
            level_data_exchange += gather_duration;
            if (!ca_on_device) {  // (the device transition slices its own blocks)
                update_neighbor_slicing_for_level(
                    level, is_symmetric,
                    /*use_received_assisting_skeletons=*/true);
            }
#ifdef H2_HAVE_GPU
            if (ca_on_device) {
                // the assisting boxes' skeletons, for the device transition,
                // and for the device solve (the level clears its copy)
                gpu_driver.refresh_remote_skeletons();
                auto& assisting = gpu::device_solve_store().ca_assisting[current_level];
                for (const auto& a : level.assisting_boxes) {
                    auto& rec = assisting[a.morton_index];
                    rec.n = static_cast<int>(a.indices.size());
                    rec.skeleton.assign(a.skel_indices.begin(), a.skel_indices.end());
                }
                gpu_driver.check_ca_replicas(current_level, level_comm,
                                             print_detail && rank == level_print_rank);
            }
#endif
            if (print_detail && rank == level_print_rank) {
                std::cout << "  CA post-elimination assisting gather time: "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(
                                 gather_duration).count()
                          << " ms" << std::endl;
            }
        } else {
            clear_ghosts(level);
        }

        memory_diagnostics.record(
            tree, current_level, use_CA_level ? "CA" : "color",
            "post_assisting_gather");
        
        

        
        // ===== Step 4: Build parent level interactions =====
        
        // Special handling for level 1: set all DOFs as skeleton (no elimination)
        if (current_level == 1 && level.is_process_active) {
            if (print_detail && rank == level_print_rank) {
                std::cout << "  Setting all DOFs as skeleton (no redundant DOFs at level 1)" << std::endl;
            }
            
            for (auto& box : level.local_boxes) {
                // All points are skeleton (no elimination at level 1)
                box.skeleton_indices.resize(box.num_points);
                for (int64_t i = 0; i < box.num_points; ++i) {
                    box.skeleton_indices[i] = i;
                }
                
                // No redundant DOFs
                box.redundant_indices.clear();
            }
        }
        
#ifdef H2_HAVE_GPU
        if (current_level == 1 && level.is_process_active) {
            // level 1 adopts the blocks of the level-2 transition; its own
            // transition builds the root on the device
            gpu_driver.adopt_level(current_level);
        }
#endif
        auto transition_start = std::chrono::high_resolution_clock::now();
        host_phase_ms[2] = std::chrono::duration<double, std::milli>(transition_start - elim_end_level).count();
        // if(rank == 1)
        // {
        //     build_parent_level_interactions<CoordType, DataType, KernelType>(
        //         level,          // child_level
        //         dimension,
        //         is_symmetric,
        //         is_hermitian,
        //         kernel,
        //         tree->global_bounds
        //     );
            
        // }
        // MPI_Barrier(tree->comm);
        // exit(0);
        std::vector<BoxData<CoordType, DataType>> parent_boxes;
        bool transition_on_device = false;
#ifdef H2_HAVE_GPU
        // Device transition (see ColorGpuDriver::transition), or this
        // level's blocks back in the host BoxData for the host transition.
        transition_on_device = gpu_driver.transition(
            level, parent_level, parent_boxes,
            [&] {
                return build_parent_level_structure(
                    level, tree->levels[current_level - 1], dimension, tree->global_bounds);
            },
            print_detail && rank == level_print_rank, level_comm, transition_comm);
#endif
        if (level.is_process_active && !transition_on_device) {
            parent_boxes = build_parent_level_interactions<CoordType, DataType, KernelType>(
                level,
                tree->levels[current_level - 1],
                dimension,
                is_symmetric,
                is_hermitian,
                kernel,
                tree->global_bounds
            );
        }

        if (memory_diagnostics.enabled()) {
            memory_diagnostics.record(
                tree, current_level, use_CA_level ? "CA" : "color",
                "post_parent_build", 0, 0,
                h2_diag_box_vector_bytes(parent_boxes));
        }
        
        
        
        
        auto transition_end = std::chrono::high_resolution_clock::now();
        auto transition_duration = std::chrono::duration_cast<std::chrono::milliseconds>(transition_end - transition_start);
        long long transition_max_ms = static_cast<long long>(transition_duration.count());
        if (print_detail) {
            MPI_Allreduce(MPI_IN_PLACE, &transition_max_ms, 1, MPI_LONG_LONG, MPI_MAX, transition_comm);
        }
        
        if (print_detail && rank == level_print_rank) {
            std::cout << "  Level transition time: " << transition_duration.count() << " ms (max over ranks "
                      << transition_max_ms << " ms)" << std::endl;
            std::cout << "  Parent boxes created: " << parent_boxes.size() << std::endl;
        }
        
        // ===== Step 5: Handle process reduction =====
        
        const bool reduction_occurred =
            (parent_level.num_active_processes != level.num_active_processes);
        const bool keep_local_parent_boxes =
            level.is_process_active && (level.parent_level_owner == rank);
        const bool sends_parent_boxes =
            reduction_occurred && level.is_process_active && level.parent_level_owner != rank;

        std::vector<char> send_buffer;
        int64_t send_buffer_size = 0;
        // wall-time split of the reduction (this rank): pack, barrier, receive,
        // unpack, edge sharing, send, clear
        double red_t[7] = {0, 0, 0, 0, 0, 0, 0};
        auto red_since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
        auto red_mark = clock::now();
        if (sends_parent_boxes) {
            // Pre-pack parent boxes before the synchronization point so the
            // communication timer excludes sender-side serialization work.
            send_buffer = serialize_boxes(parent_boxes);
            send_buffer_size = static_cast<int64_t>(send_buffer.size());
        }
        red_t[0] = red_since(red_mark);

        if (memory_diagnostics.enabled()) {
            memory_diagnostics.record(
                tree, current_level, use_CA_level ? "CA" : "color",
                "post_parent_pack", 0, 0,
                h2_diag_box_vector_bytes(parent_boxes),
                h2_diag_vector_bytes(send_buffer));
        }

        if (reduction_occurred) {
            // Synchronize after parent-box construction/packing so the reduction
            // communication timer does not count time spent waiting for slower
            // ranks that are still preparing their payloads.
            red_mark = clock::now();
            MPI_Barrier(transition_comm);
            red_t[1] = red_since(red_mark);
        }

        if (!reduction_occurred) {
            if (parent_level.is_process_active) {
                parent_level.local_boxes = std::move(parent_boxes);
            }
        } else if (parent_level.is_process_active) {
            std::vector<BoxData<CoordType, DataType>> all_parent_boxes;

            for (int child_rank : parent_level.children_senders) {
                if (keep_local_parent_boxes && child_rank == rank) {
                    all_parent_boxes.insert(
                        all_parent_boxes.end(),
                        std::make_move_iterator(parent_boxes.begin()),
                        std::make_move_iterator(parent_boxes.end())
                    );
                    continue;
                }
                
                std::vector<char> recv_buffer;
                int64_t buffer_size = 0;
                MPI_Status status;
                segment_start = clock::now();
                MPI_Recv(&buffer_size, 1, MPI_INT64_T, child_rank, 0, tree->comm, &status);
                
                recv_buffer.resize(buffer_size);
                MPI_Recv_large(recv_buffer.data(), buffer_size, MPI_CHAR, child_rank, 1, tree->comm, &status);
                level_reduction += (clock::now() - segment_start);
                red_t[2] += red_since(segment_start);
                red_mark = clock::now();

                if (memory_diagnostics.enabled()) {
                    memory_diagnostics.record(
                        tree, current_level, use_CA_level ? "CA" : "color",
                        "parent_receive_buffer", 0, 0,
                        h2_diag_box_vector_bytes(parent_boxes) +
                            h2_diag_box_vector_bytes(all_parent_boxes),
                        h2_diag_vector_bytes(recv_buffer));
                }
                
                std::vector<BoxData<CoordType, DataType>> child_parent_boxes =
                    deserialize_boxes<CoordType, DataType>(recv_buffer);
                red_t[3] += red_since(red_mark);

                if (memory_diagnostics.enabled()) {
                    memory_diagnostics.record(
                        tree, current_level, use_CA_level ? "CA" : "color",
                        "parent_deserialized", 0, 0,
                        h2_diag_box_vector_bytes(parent_boxes) +
                            h2_diag_box_vector_bytes(all_parent_boxes) +
                            h2_diag_box_vector_bytes(child_parent_boxes),
                        h2_diag_vector_bytes(recv_buffer));
                }
                
                all_parent_boxes.insert(
                    all_parent_boxes.end(),
                    std::make_move_iterator(child_parent_boxes.begin()),
                    std::make_move_iterator(child_parent_boxes.end())
                );
            }

            parent_level.local_boxes = std::move(all_parent_boxes);
        }

        red_mark = clock::now();
        if (parent_level.is_process_active &&
            is_symmetric && !is_hermitian) {
            share_symmetric_level_edges(parent_level);
        }
        red_t[4] = red_since(red_mark);

        red_mark = clock::now();
        if (sends_parent_boxes) {
            // This process does NOT own the parent - send to parent_level_owner
            segment_start = clock::now();
            MPI_Send(&send_buffer_size, 1, MPI_INT64_T, level.parent_level_owner, 0, tree->comm);

            // Send data
            MPI_Send_large(send_buffer.data(), send_buffer_size, MPI_CHAR, level.parent_level_owner, 1, tree->comm);
            level_reduction += (clock::now() - segment_start);

            // This process will become inactive at parent level
            parent_level.local_boxes.clear();

            if (print_trace) {
                std::cout << "  Process " << rank << " sending " << parent_boxes.size() 
                          << " parent boxes to process " << level.parent_level_owner << std::endl;
            }
        }

        if (memory_diagnostics.enabled()) {
            memory_diagnostics.record(
                tree, current_level, use_CA_level ? "CA" : "color",
                "post_parent_redistribution", 0, 0,
                h2_diag_box_vector_bytes(parent_boxes),
                h2_diag_vector_bytes(send_buffer));
        }
        
        red_t[5] = red_since(red_mark);
        // Clear modified interaction matrices to free memory
        red_mark = clock::now();
        if (level.is_process_active) {
            clear_modified_interaction_matrices(level, use_CA_level);
        }
        red_t[6] = red_since(red_mark);
        if (reduction_occurred && print_detail) {
            double red_max[7];
            MPI_Reduce(red_t, red_max, 7, MPI_DOUBLE, MPI_MAX, 0, transition_comm);
            int transition_rank = 0;
            MPI_Comm_rank(transition_comm, &transition_rank);
            if (transition_rank == 0) {
                std::printf("  process reduction (max over ranks): pack %.2f s, barrier %.2f s, receive %.2f s, "
                            "unpack %.2f s, share edges %.2f s, send %.2f s, clear %.2f s\n",
                            red_max[0], red_max[1], red_max[2], red_max[3], red_max[4], red_max[5], red_max[6]);
                std::fflush(stdout);
            }
        }

        memory_diagnostics.record(
            tree, current_level, use_CA_level ? "CA" : "color",
            "post_level_clear");

        // Teardown
        if (level.is_process_active) {
            for (auto& [morton, lock] : level.box_locks) {
                omp_destroy_lock(lock);
                delete lock;
            }
        }

        {
            double local_deltas[2] = {
                std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(level_data_exchange).count(),
                std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(level_reduction).count()
            };
            double max_deltas[2] = {0.0, 0.0};
            int transition_rank = -1;
            MPI_Comm_rank(transition_comm, &transition_rank);
            MPI_Reduce(
                local_deltas,
                max_deltas,
                2,
                MPI_DOUBLE,
                MPI_MAX,
                0,
                transition_comm);
            if (transition_rank == 0) {
                total_data_exchange_time += std::chrono::duration_cast<clock::duration>(
                    std::chrono::duration<double, std::milli>(max_deltas[0]));
                total_reduction_time += std::chrono::duration_cast<clock::duration>(
                    std::chrono::duration<double, std::milli>(max_deltas[1]));
            }
        }

        auto level_end = std::chrono::high_resolution_clock::now();
        auto level_duration = std::chrono::duration_cast<std::chrono::milliseconds>(level_end - level_start);
        host_phase_ms[3] = std::chrono::duration<double, std::milli>(level_end - transition_end).count();
        
        if (print_detail && rank == level_print_rank) {
            std::cout << "  Total level time: " << level_duration.count() << " ms" << std::endl;
        }
        if (print_detail && current_level > 1 && level.is_process_active) {
            double host_phase_max[4];
            MPI_Allreduce(host_phase_ms, host_phase_max, 4, MPI_DOUBLE, MPI_MAX, level_comm);
            if (rank == level_print_rank) {
                std::printf("  level host phases (ms, this rank / slowest): schedule setup %.0f / %.0f, start to elimination %.0f / %.0f, "
                            "elimination to transition %.0f / %.0f, after the transition %.0f / %.0f\n",
                            host_phase_ms[0], host_phase_max[0], host_phase_ms[1], host_phase_max[1], host_phase_ms[2], host_phase_max[2],
                            host_phase_ms[3], host_phase_max[3]);
                std::fflush(stdout);
            }
        }

        if (dynamic_threading.enabled &&
            level.is_process_active && !parent_level.is_process_active) {
            park_inactive_rank_on_service_cpu(dynamic_threading);
        }
        
    }


    if (out_rankmax) *out_rankmax = local_max_skel;

    
    // ===== Special handling for level 0 (root) =====
    
    const int root_print_rank = smallest_active_rank(tree->levels[0]);
    if (print_detail && rank == root_print_rank) {
        std::cout << "\n===== Level 0 (Root) =====" << std::endl;
    }
    
    auto& root_level = tree->levels[0];
    
    if (!root_level.is_process_active || root_level.local_boxes.empty()) {
        if (print_trace) {
            std::cout << "  Process " << rank << " has no boxes at root level" << std::endl;
        }
    } else {
        if (root_level.local_boxes.size() != 1) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Expected exactly 1 box at root level, got " +
                std::to_string(root_level.local_boxes.size()));
        }
        
        auto& root_box = root_level.local_boxes[0];
        
        if (print_detail && rank == root_print_rank) {
            std::cout << "  Root box points: " << root_box.num_points << std::endl;
        }
        
        int64_t n = root_box.num_points;
        bool root_on_device = false;
#ifdef H2_HAVE_GPU
        // root block assembled on the device by the level-1 transition, or
        // on the host (multi-rank runs)
        root_on_device = gpu_driver.factor_root(root_box, print_detail && rank == root_print_rank);
#endif
        if (!root_on_device) {
            // At level 0, the assembled matrix is just the schur complement
            if (!root_box.schur_complement.is_allocated()) {
                throw std::runtime_error(
                    "hierarchical_factorization_parallel: Root box schur complement not allocated");
            }
        
            n = root_box.schur_complement.rows;
        
            if (print_detail && rank == root_print_rank) {
                std::cout << "  Schur complement size: " << n << " × " << n << std::endl;
            }
        
            // Factorize the root schur complement for diagonal solve
            root_box.X_RR.allocate(n, n, MatrixStorage<DataType>::FULL);
        
            if (factorization_method == FactorizationMethod::CHOLESKY) {
                root_box.X_RR_pivots.clear();
                // Copy schur complement to X_RR
                root_box.X_RR.data = root_box.schur_complement.data;
            
                // Perform Cholesky factorization in-place
                char uplo = 'L';
                int nn = n;
                int info = 0;
            
                if constexpr (std::is_same_v<DataType, double>) {
                    dpotrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn, &info);
                } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                    zsychol_(&uplo, &nn, root_box.X_RR.data.data(), &nn, &info);
                }
            
                if (info != 0) {
                    throw std::runtime_error(
                        "hierarchical_factorization_parallel: Cholesky factorization of root failed at pivot " +
                        std::to_string(info));
                }
            
                root_box.X_RR.format = MatrixStorage<DataType>::CHOLESKY_L;
            
                if (print_detail && rank == root_print_rank) {
                    std::cout << "  ✓ Root Cholesky factorization complete" << std::endl;
                }
            
            } else if (factorization_method == FactorizationMethod::LU) {
                root_box.X_RR.data = root_box.schur_complement.data;
                root_box.X_RR_pivots.resize(static_cast<size_t>(n));
                // LU reads both triangles: a symmetric problem factors the root's
                // symmetric part, as its boxes' X_RR (see compute_and_modify)
                if (is_symmetric || is_hermitian) {
                    symmetrize_in_place(root_box.X_RR, n, is_hermitian && !is_symmetric);
                }

                int nn = n;
                int info = 0;

                if constexpr (std::is_same_v<DataType, double>) {
                    dgetrf_(&nn, &nn, root_box.X_RR.data.data(), &nn,
                            root_box.X_RR_pivots.data(), &info);
                } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                    zgetrf_(&nn, &nn, root_box.X_RR.data.data(), &nn,
                            root_box.X_RR_pivots.data(), &info);
                }

                if (info != 0) {
                    throw std::runtime_error(
                        "hierarchical_factorization_parallel: LU factorization of root failed at pivot " +
                        std::to_string(info));
                }

                root_box.X_RR.format = MatrixStorage<DataType>::LU_FACTORED;

                if (print_detail && rank == root_print_rank) {
                    std::cout << "  ✓ Root LU factorization complete" << std::endl;
                }

            } else if (factorization_method == FactorizationMethod::BUNCH_KAUFMAN) {
                root_box.X_RR.data = root_box.schur_complement.data;
                root_box.X_RR_pivots.resize(static_cast<size_t>(n));

                char uplo = 'L';
                int nn = n;
                int lwork = -1;
                int info = 0;
                std::vector<DataType> work(1);
                sytrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn,
                       root_box.X_RR_pivots.data(), work.data(), &lwork, &info);
                if (info != 0) {
                    throw std::runtime_error(
                        "hierarchical_factorization_parallel: Bunch-Kaufman root workspace query failed with INFO = " +
                        std::to_string(info));
                }
                lwork = std::max(1, static_cast<int>(std::real(work[0])));
                work.resize(static_cast<size_t>(lwork));
                sytrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn,
                       root_box.X_RR_pivots.data(), work.data(), &lwork, &info);

                if (info != 0) {
                    throw std::runtime_error(
                        "hierarchical_factorization_parallel: Bunch-Kaufman factorization of root failed at pivot " +
                        std::to_string(info));
                }

                root_box.X_RR.format = MatrixStorage<DataType>::BUNCH_KAUFMAN;

                if (print_detail && rank == root_print_rank) {
                    std::cout << "  ✓ Root Bunch-Kaufman factorization complete" << std::endl;
                }

            } else {
                // No factorization: just copy schur complement to X_RR
                root_box.X_RR.data = root_box.schur_complement.data;
                root_box.X_RR.format = MatrixStorage<DataType>::FULL;
                root_box.X_RR_pivots.clear();
            
                if (print_detail && rank == root_print_rank) {
                    std::cout << "  ✓ Root matrix copied (no factorization)" << std::endl;
                }
            }
        
        }

        // Mark root as skeleton only (no redundant DOFs at this level)
        root_box.skeleton_indices.resize(n);
        for (int64_t i = 0; i < n; ++i) {
            root_box.skeleton_indices[i] = i;
        }
        root_box.redundant_indices.clear();
    }

    memory_diagnostics.record(tree, 0, "root", "post_root_factorization");
    
    passive_mpi_barrier(tree->comm);

    const double local_timing_ms[2] = {
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            total_data_exchange_time).count(),
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            total_reduction_time).count()
    };
    double global_timing_ms[2] = {0.0, 0.0};
    MPI_Reduce(
        local_timing_ms,
        global_timing_ms,
        2,
        MPI_DOUBLE,
        MPI_SUM,
        root_print_rank,
        tree->comm);

    if (print_summary && rank == root_print_rank) {
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - total_time);
        std::cout << "\n========================================" << std::endl;
        std::cout << "✓ Hierarchical Factorization Complete" << std::endl;
        std::cout << "  total time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  data exchange communication time: " << std::llround(global_timing_ms[0]) << " ms" << std::endl;
        std::cout << "  process reduction communication time: " << std::llround(global_timing_ms[1]) << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    destroy_factorization_communicators(factorization_comms);
    restore_base_process_affinity();
    clear_runtime_fmm_thread_count();
    destroy_dynamic_threading_context(dynamic_threading);


    size_t local_memory_usage = 0;
    size_t CA_halo_memory_usage = 0;
    for (int current_level = leaf_level; current_level >= 1; current_level--) {
    
        auto& level = tree->levels[current_level];
        for (auto& box : level.local_boxes) {
            local_memory_usage += calculate_box_data_size(box);
        }
        if (tree->level_uses_CA(current_level)) {
            for (const auto& box : level.ghost_boxes) {
                CA_halo_memory_usage += calculate_box_data_size(box);
            }
        }
    }
    if (print_detail) {
        printf("factorization memory usage on rank %d: %.10f GB local + %.10f GB CA halo\n",
               rank,
               local_memory_usage / (1024.0 * 1024.0 * 1024.0),
               CA_halo_memory_usage / (1024.0 * 1024.0 * 1024.0));
        fflush(stdout);
    }

    *memory_per_rank = local_memory_usage + CA_halo_memory_usage;

    memory_diagnostics.record(tree, 0, "retained", "factorization_complete");

    const auto verification_start = clock::now();
    if (print_summary) {
        double logabsdet;
        DataType phase;
        hierarchical_logdet_parallel(tree, &logabsdet, &phase);
        if (rank == 0) {
            std::cout.precision(17);
            std::cout << "logdet: " << phase << " " << logabsdet << std::endl;
        }

        (void)butterfly::h2_quick_verification(tree, kernel);
    }
    h2_verification_seconds() = std::chrono::duration<double>(clock::now() - verification_start).count();

    memory_diagnostics.record(tree, 0, "verification", "post_quick_verification");
    memory_diagnostics.print(tree->comm, rank, size);
}

template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_factorization_parallel_if_supported(
    fmm::ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    int use_sketch,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
    int64_t* out_rankmax = nullptr,
    size_t* memory_per_rank = nullptr,
    int lazy_schur = 0,
    int gemm_split = 16,
    int ca_staged_halo = 0,
    int ca_owner_component = 0,
    int ca_owner_serial = 0,
    int verbosity = 1) {
    if constexpr (std::is_same_v<DataType, double> ||
                  std::is_same_v<DataType, std::complex<double>>) {
        hierarchical_factorization_parallel<CoordType, DataType, KernelType>(
            tree, 
            kernel, 
            tolerance, 
            use_sketch,
            is_symmetric, 
            is_hermitian, 
            factorization_method, 
            unit_proxy_points, 
            num_proxy, 
            proxy_radius, 
            out_rankmax,
            memory_per_rank,
            lazy_schur,
            gemm_split,
            ca_staged_halo,
            ca_owner_component,
            ca_owner_serial,
            verbosity);   // instantiated ONLY for double types
    } else {
        throw std::runtime_error("H2/FMM only supports double / std::complex<double>");
    }
}

template<typename CoordType, typename DataType>
void butterfly_factorization_parallel(H2<CoordType,DataType>* solver, double* factorization_time, double* entryeval_time) {

  if (solver->build_state == H2BuildState::RS_FACTORIZED) {
    if (factorization_time) *factorization_time = 0.0;
    if (entryeval_time) *entryeval_time = 0.0;
    return;
  }
  if (solver->build_state == H2BuildState::H2_COMPRESSED) {
    throw std::runtime_error(
      "butterfly_factorization_parallel: solver already contains a compression-only H2 representation");
  }

  int rank = 0;
  MPI_Comm_rank(solver->comm, &rank);

  const auto factorization_method =
    h2_xrr_factorization_method(solver->options);
  configure_color_gpu(solver->options.use_gpu != 0);
  fmm::env::report_environment(solver->comm);
#ifdef H2_HAVE_GPU
  fmm::gpu::tensor_core_gemm() = solver->options.use_gpu == 2;
#endif
  fmm::HierarchicalFactorization<CoordType, DataType, butterfly::H2Kernel<CoordType, DataType>> factorizer(
    solver->options.N,
    fmm::MatrixProperty::SYMMETRIC,
    &solver->kernel,
    solver->options.dimension,
    factorization_method,
    solver->options.num_proxy);

  const auto& unit_proxy = factorizer.get_unit_proxy_points();
  const int num_proxy = factorizer.get_num_proxy_points();

  const bool is_symmetric = true;
  const bool is_hermitian = false;
  
  solver->kernel.entryeval_time_per_thread.assign(omp_get_max_threads(), 0.0);
  auto& ev = solver->kernel.entryeval_time_per_thread;

//   auto total_start = std::chrono::high_resolution_clock::now();
  double t0 = MPI_Wtime();
  butterfly::h2_verification_seconds() = 0.0;
  butterfly::h2_warmup_seconds() = 0.0;
  butterfly::hierarchical_factorization_parallel_if_supported(
    solver->tree.get(),
    &solver->kernel,
    solver->options.tolerance,
    solver->options.use_sketch,
    is_symmetric,
    is_hermitian,
    factorization_method,
    unit_proxy,
    num_proxy,
    0.0, // proxy_radius = 0
    &solver->last_factor_rankmax,
    &solver->factorization_memory,
    solver->options.lazy_schur,
    solver->options.gemm_split,
    solver->options.ca_staged_halo,
    solver->options.ca_owner_component,
    solver->options.ca_owner_serial,
    solver->options.verbosity);
  // the logdet and quick verification that end the call, and the GPU
  // kernel warm-up before its first level, are not factorization
  double tf = MPI_Wtime() - t0 - butterfly::h2_verification_seconds() - butterfly::h2_warmup_seconds();
  MPI_Allreduce(MPI_IN_PLACE, &tf, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
  *factorization_time = tf;
  double tv = butterfly::h2_verification_seconds();
  MPI_Allreduce(MPI_IN_PLACE, &tv, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
  if (solver->tree->mpi_rank == 0 && solver->options.verbosity >= 0) {
    std::printf("  logdet and quick verification: %.2f s (not in the factor time)\n", tv);
    std::fflush(stdout);
  }
  
  double t_entry = 0.0;
  if (!ev.empty()) {
    t_entry = *std::max_element(ev.begin(), ev.end());
  }
  MPI_Allreduce(MPI_IN_PLACE, &t_entry, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
  *entryeval_time = t_entry;
  solver->factorized = true;
  solver->build_state = H2BuildState::RS_FACTORIZED;

  MPI_Allreduce(MPI_IN_PLACE, &solver->last_factor_rankmax, 1, MPI_INT64_T, MPI_MAX, solver->comm);

//   auto total_end = std::chrono::high_resolution_clock::now();
//   auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
//     total_end - total_start);


//   if (rank == 0) {
//     std::cout << "\n========================================" << std::endl;
//     std::cout << "Total factorization time: " << total_duration.count() << " ms" << std::endl;
//     std::cout << "========================================\n" << std::endl;
//   }

}

/**
 * @brief Hierarchical solve routine (parallel MPI version)
 */

} // namespace butterfly
