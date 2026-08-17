#pragma once

// Runtime threading / CPU-affinity support helpers for the parallel
// hierarchical factorization.
//
// These were originally defined inline inside distributed_routine_all.cpp
// (namespace fmm). They are extracted here so that headers which host a copy
// of hierarchical_factorization_parallel (e.g. butterfly_h2/butterfly_integration.hpp)
// can see the declarations without including a .cpp translation unit.
//
// NOTE: distributed_routine_all.cpp still contains its own copy of these
// definitions, but that file is not part of the ButterflyPACK build, so there
// is no ODR conflict. If it is ever revived, replace its copy with an include
// of this header.

#include <atomic>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <thread>
#include <vector>

#include <mpi.h>
#include <omp.h>
#include <sched.h>

#include "factorization.hpp"

namespace fmm {

template<typename CoordType, typename DataType>
int smallest_active_rank(const TreeLevel<CoordType, DataType>& level) {
    if (level.active_process_ranks.empty()) {
        return 0;
    }
    return *std::min_element(
        level.active_process_ranks.begin(),
        level.active_process_ranks.end());
}

template<typename DataType>
MPI_Datatype mpi_datatype_for() {
    if constexpr (std::is_same_v<DataType, double>) {
        return MPI_DOUBLE;
    } else if constexpr (std::is_same_v<DataType, float>) {
        return MPI_FLOAT;
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        return MPI_CXX_DOUBLE_COMPLEX;
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        return MPI_CXX_FLOAT_COMPLEX;
    } else {
        throw std::runtime_error("Unsupported DataType for MPI");
    }
}

inline int parse_positive_thread_count(const char* value) {
    if (value == nullptr || *value == '\0') {
        return 0;
    }

    const long parsed = std::strtol(value, nullptr, 10);
    if (parsed <= 0) {
        return 0;
    }
    if (parsed > static_cast<long>(std::numeric_limits<int>::max())) {
        return std::numeric_limits<int>::max();
    }
    return static_cast<int>(parsed);
}

inline bool query_current_affinity_mask(cpu_set_t& mask) {
    CPU_ZERO(&mask);
    return sched_getaffinity(0, sizeof(mask), &mask) == 0;
}

inline std::vector<int> cpu_list_from_mask(const cpu_set_t& mask) {
    std::vector<int> cpus;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        if (CPU_ISSET(cpu, &mask)) {
            cpus.push_back(cpu);
        }
    }
    return cpus;
}

inline std::atomic<int>& runtime_thread_override_storage() {
    static std::atomic<int> override_count{0};
    return override_count;
}

inline const std::vector<int>& base_process_cpu_list() {
    static std::once_flag init_flag;
    static std::vector<int> cpus;
    std::call_once(init_flag, []() {
        cpu_set_t mask;
        if (query_current_affinity_mask(mask)) {
            cpus = cpu_list_from_mask(mask);
        }
    });
    return cpus;
}

inline int visible_process_cpu_count() {
    const auto& cpus = base_process_cpu_list();
    if (!cpus.empty()) {
        return static_cast<int>(cpus.size());
    }
    const unsigned int hardware_threads = std::thread::hardware_concurrency();
    return hardware_threads == 0 ? 1 : static_cast<int>(hardware_threads);
}

inline bool build_cpu_affinity_mask(const std::vector<int>& cpus, cpu_set_t& mask) {
    CPU_ZERO(&mask);
    bool any_cpu = false;
    for (int cpu : cpus) {
        if (cpu >= 0 && cpu < CPU_SETSIZE) {
            CPU_SET(cpu, &mask);
            any_cpu = true;
        }
    }
    return any_cpu;
}

inline bool set_process_cpu_subset(const std::vector<int>& cpus) {
    if (cpus.empty()) {
        return false;
    }
    cpu_set_t mask;
    if (!build_cpu_affinity_mask(cpus, mask)) {
        return false;
    }
    return sched_setaffinity(0, sizeof(mask), &mask) == 0;
}

inline bool set_runtime_cpu_subset(const std::vector<int>& cpus, int thread_count) {
    if (!set_process_cpu_subset(cpus)) {
        return false;
    }

    const int team_size = std::max(1, std::min<int>(thread_count, static_cast<int>(cpus.size())));
    std::atomic<bool> pin_failed{false};

    // Pin each OpenMP worker to a specific CPU inside the rank-local slice so
    // the team does not float around within the whole subset.
    #pragma omp parallel num_threads(team_size) default(shared)
    {
        const int tid = omp_get_thread_num();
        const int cpu = cpus[static_cast<std::size_t>(std::min(tid, team_size - 1))];
        cpu_set_t thread_mask;
        CPU_ZERO(&thread_mask);
        CPU_SET(cpu, &thread_mask);
        if (sched_setaffinity(0, sizeof(thread_mask), &thread_mask) != 0) {
            pin_failed.store(true, std::memory_order_relaxed);
        }
    }

    return !pin_failed.load(std::memory_order_relaxed);
}

inline bool restore_base_process_affinity() {
    return set_process_cpu_subset(base_process_cpu_list());
}

inline void set_runtime_fmm_thread_count(int thread_count) {
    thread_count = std::max(1, thread_count);
    runtime_thread_override_storage().store(thread_count, std::memory_order_relaxed);
    omp_set_num_threads(thread_count);
}

inline void clear_runtime_fmm_thread_count() {
    runtime_thread_override_storage().store(0, std::memory_order_relaxed);
}

inline void passive_mpi_barrier(MPI_Comm comm) {
    if (comm == MPI_COMM_NULL) {
        return;
    }

    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Ibarrier(comm, &request);
    int complete = 0;
    while (!complete) {
        MPI_Test(&request, &complete, MPI_STATUS_IGNORE);
        if (!complete) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }
}

inline int configured_fmm_thread_count() {
    if (const int override_count =
            runtime_thread_override_storage().load(std::memory_order_relaxed);
        override_count > 0) {
        return override_count;
    }
    if (const char* omp_threads = std::getenv("OMP_NUM_THREADS")) {
        const int parsed = parse_positive_thread_count(omp_threads);
        if (parsed > 0) {
            return parsed;
        }
    }
    return std::max(1, omp_get_max_threads());
}

struct DynamicThreadingContext {
    bool enabled = false;
    MPI_Comm shared_comm = MPI_COMM_NULL;
    int shared_rank = 0;
    int shared_size = 1;
    int cpu_cap_per_node = 0;
    int original_thread_count = 1;
};

struct LevelThreadPlan {
    bool active = false;
    int threads = 1;
    int active_on_node = 0;
    int max_active_on_any_node = 1;
    int active_slot = -1;
    int cpu_begin = -1;
    int cpu_end = -1;
};

inline DynamicThreadingContext make_dynamic_threading_context(MPI_Comm comm) {
    DynamicThreadingContext context;
    context.original_thread_count = std::max(1, omp_get_max_threads());
    const int requested_cpus_per_node =
        parse_positive_thread_count(std::getenv("FMM_MAX_CPUS_PER_NODE"));
    if (requested_cpus_per_node <= 0) {
        return context;
    }

    context.enabled = true;
    MPI_Comm_split_type(
        comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &context.shared_comm);
    MPI_Comm_rank(context.shared_comm, &context.shared_rank);
    MPI_Comm_size(context.shared_comm, &context.shared_size);

    int visible_cpus = visible_process_cpu_count();
    int min_visible_cpus = visible_cpus;
    MPI_Allreduce(&visible_cpus, &min_visible_cpus, 1, MPI_INT, MPI_MIN, comm);
    context.cpu_cap_per_node = std::max(
        1,
        std::min(requested_cpus_per_node, min_visible_cpus));
    return context;
}

inline void destroy_dynamic_threading_context(DynamicThreadingContext& context) {
    if (context.enabled) {
        omp_set_num_threads(std::max(1, context.original_thread_count));
    }
    if (context.shared_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&context.shared_comm);
        context.shared_comm = MPI_COMM_NULL;
    }
}

inline LevelThreadPlan configure_static_process_thread_plan(
    const DynamicThreadingContext& context) {
    LevelThreadPlan plan;
    plan.active = true;

    if (!context.enabled) {
        return plan;
    }

    const auto& base_cpus = base_process_cpu_list();
    const int usable_cpu_count = std::min<int>(
        context.cpu_cap_per_node,
        static_cast<int>(base_cpus.size()));
    const int threads = std::max(1, usable_cpu_count / std::max(1, context.shared_size));

    plan.threads = threads;
    plan.active_on_node = context.shared_size;
    plan.max_active_on_any_node = context.shared_size;
    plan.active_slot = context.shared_rank;

    if (usable_cpu_count > 0) {
        int begin_offset = context.shared_rank * threads;
        if (begin_offset >= usable_cpu_count) {
            begin_offset = usable_cpu_count - 1;
        }
        const int end_offset = std::min(begin_offset + threads, usable_cpu_count);
        std::vector<int> assigned_cpus(
            base_cpus.begin() + begin_offset,
            base_cpus.begin() + end_offset);
        if (!assigned_cpus.empty()) {
            plan.cpu_begin = assigned_cpus.front();
            plan.cpu_end = assigned_cpus.back();
            set_runtime_fmm_thread_count(plan.threads);
            set_runtime_cpu_subset(assigned_cpus, plan.threads);
        } else {
            set_runtime_fmm_thread_count(plan.threads);
        }
    } else {
        set_runtime_fmm_thread_count(plan.threads);
    }

    return plan;
}

template<typename CoordType, typename DataType>
LevelThreadPlan configure_level_thread_plan(
    const DynamicThreadingContext& context,
    MPI_Comm level_comm,
    const TreeLevel<CoordType, DataType>& level) {
    LevelThreadPlan plan;
    plan.active = level.is_process_active;

    if (!context.enabled || !plan.active) {
        return plan;
    }
    if (level_comm == MPI_COMM_NULL) {
        throw std::runtime_error(
            "configure_level_thread_plan: active rank has no level communicator");
    }

    MPI_Comm active_shared_comm = MPI_COMM_NULL;
    MPI_Comm_split_type(
        level_comm,
        MPI_COMM_TYPE_SHARED,
        0,
        MPI_INFO_NULL,
        &active_shared_comm);
    MPI_Comm_rank(active_shared_comm, &plan.active_slot);
    MPI_Comm_size(active_shared_comm, &plan.active_on_node);
    MPI_Allreduce(
        &plan.active_on_node,
        &plan.max_active_on_any_node,
        1,
        MPI_INT,
        MPI_MAX,
        level_comm);
    plan.max_active_on_any_node = std::max(1, plan.max_active_on_any_node);

    const auto& base_cpus = base_process_cpu_list();
    const int usable_cpu_count = std::min<int>(
        context.cpu_cap_per_node,
        static_cast<int>(base_cpus.size()));
    plan.threads = std::max(
        1, usable_cpu_count / std::max(1, plan.active_on_node));

    if (usable_cpu_count > 0) {
        int begin_offset = plan.active_slot * plan.threads;
        if (begin_offset >= usable_cpu_count) {
            begin_offset = usable_cpu_count - 1;
        }
        const int end_offset = std::min(
            begin_offset + plan.threads, usable_cpu_count);
        std::vector<int> assigned_cpus(
            base_cpus.begin() + begin_offset,
            base_cpus.begin() + end_offset);
        if (!assigned_cpus.empty()) {
            plan.cpu_begin = assigned_cpus.front();
            plan.cpu_end = assigned_cpus.back();
            set_runtime_fmm_thread_count(plan.threads);
            set_runtime_cpu_subset(assigned_cpus, plan.threads);
        } else {
            set_runtime_fmm_thread_count(plan.threads);
        }
    } else {
        set_runtime_fmm_thread_count(plan.threads);
    }

    MPI_Comm_free(&active_shared_comm);

    return plan;
}

inline void print_level_thread_plan(
    const DynamicThreadingContext& context,
    int level_index,
    const char* phase_name,
    int level_print_rank,
    int rank,
    const LevelThreadPlan& local_plan,
    bool verbose) {
    if (!context.enabled || !verbose || rank != level_print_rank) {
        return;
    }

    std::cout << "  Thread plan [" << phase_name << "] level " << level_index
              << ": cpu_cap_per_node=" << context.cpu_cap_per_node
              << ", max_active_ranks_on_any_node=" << local_plan.max_active_on_any_node
              << ", active_threads=" << std::max(1, local_plan.threads)
              << ", node_active=" << local_plan.active_on_node
              << std::endl;
    if (local_plan.cpu_begin >= 0 && local_plan.cpu_end >= 0) {
        std::cout << "    Rank " << rank
                  << ": active=" << (local_plan.active ? "yes" : "no")
                  << ", threads=" << local_plan.threads
                  << ", cpus=[" << local_plan.cpu_begin << ", " << local_plan.cpu_end << "]"
                  << std::endl;
    }
}

template<typename DurationType>
inline void reduce_active_duration_bounds_ms(
    MPI_Comm comm,
    int root_rank,
    bool active,
    DurationType duration,
    double& min_ms,
    double& max_ms) {
    const double local_ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(duration).count();
    const double local_min_ms =
        active ? local_ms : std::numeric_limits<double>::infinity();
    const double local_max_ms = active ? local_ms : 0.0;
    MPI_Reduce(&local_min_ms, &min_ms, 1, MPI_DOUBLE, MPI_MIN, root_rank, comm);
    MPI_Reduce(&local_max_ms, &max_ms, 1, MPI_DOUBLE, MPI_MAX, root_rank, comm);
}

}  // namespace fmm
