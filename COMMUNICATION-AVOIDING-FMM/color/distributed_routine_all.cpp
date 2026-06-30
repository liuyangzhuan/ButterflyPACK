#include "factorization.hpp"
#include "solver.hpp"
#include "apply_mul.hpp"
#include "tree_impl.hpp"
#include "kernel.hpp"
#include "id_decomposition.hpp"
#include <iostream>
#include <iomanip>
#include <array>
#include <chrono>
#include <cmath>
#include <atomic>
#include <mutex>
#include <omp.h>
#include <random>
#include <complex>
#include <sched.h>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace fmm {

enum class KernelKind {
    LAPLACE,
    HELMHOLTZ,
    MATERN52,
    YUKAWA
};

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
    MPI_Comm comm,
    const TreeLevel<CoordType, DataType>& level) {
    LevelThreadPlan plan;
    plan.active = level.is_process_active;

    if (!context.enabled) {
        return plan;
    }

    std::vector<int> node_active_flags(static_cast<size_t>(context.shared_size), 0);
    const int local_active = plan.active ? 1 : 0;
    MPI_Allgather(
        &local_active,
        1,
        MPI_INT,
        node_active_flags.data(),
        1,
        MPI_INT,
        context.shared_comm);

    for (int active_flag : node_active_flags) {
        plan.active_on_node += active_flag;
    }
    const int inactive_on_node = context.shared_size - plan.active_on_node;
    MPI_Allreduce(
        &plan.active_on_node,
        &plan.max_active_on_any_node,
        1,
        MPI_INT,
        MPI_MAX,
        comm);
    plan.max_active_on_any_node = std::max(1, plan.max_active_on_any_node);

    int max_inactive_on_any_node = 0;
    MPI_Allreduce(
        &inactive_on_node,
        &max_inactive_on_any_node,
        1,
        MPI_INT,
        MPI_MAX,
        comm);

    const int reserved_idle_cpus =
        std::max(0, std::min(context.cpu_cap_per_node - 1, max_inactive_on_any_node));
    const int active_cpu_budget =
        std::max(1, context.cpu_cap_per_node - reserved_idle_cpus);
    const int active_threads =
        std::max(1, active_cpu_budget / plan.max_active_on_any_node);
    plan.threads = plan.active ? active_threads : 1;

    const auto& base_cpus = base_process_cpu_list();
    const int usable_cpu_count = std::min<int>(
        context.cpu_cap_per_node,
        static_cast<int>(base_cpus.size()));
    const int reserved_cpu_count = std::min(
        usable_cpu_count,
        plan.max_active_on_any_node * active_threads);

    if (plan.active) {
        plan.active_slot = 0;
        for (int local_rank = 0; local_rank < context.shared_rank; ++local_rank) {
            if (node_active_flags[static_cast<size_t>(local_rank)] != 0) {
                ++plan.active_slot;
            }
        }

        if (usable_cpu_count > 0) {
            int begin_offset = plan.active_slot * active_threads;
            if (begin_offset >= usable_cpu_count) {
                begin_offset = usable_cpu_count - 1;
            }
            const int end_offset = std::min(begin_offset + active_threads, usable_cpu_count);
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
    } else {
        if (usable_cpu_count > 0) {
            int inactive_slot = 0;
            for (int local_rank = 0; local_rank < context.shared_rank; ++local_rank) {
                if (node_active_flags[static_cast<size_t>(local_rank)] == 0) {
                    ++inactive_slot;
                }
            }

            int chosen_offset = 0;
            const int idle_span = std::max(1, usable_cpu_count - reserved_cpu_count);
            chosen_offset =
                reserved_cpu_count + std::min(inactive_slot, idle_span - 1);
            std::vector<int> idle_cpu{base_cpus[static_cast<size_t>(chosen_offset)]};
            plan.cpu_begin = idle_cpu.front();
            plan.cpu_end = idle_cpu.front();
            set_runtime_fmm_thread_count(1);
            set_runtime_cpu_subset(idle_cpu, 1);
        } else {
            set_runtime_fmm_thread_count(1);
        }
    }

    return plan;
}

inline void print_level_thread_plan(
    const DynamicThreadingContext& context,
    MPI_Comm comm,
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

inline uint64_t splitmix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

template<typename RealType>
RealType centered_uniform_from_hash(uint64_t key) {
    constexpr long double scale =
        1.0L / static_cast<long double>(std::numeric_limits<uint64_t>::max());
    const long double u =
        static_cast<long double>(splitmix64(key)) * scale;
    return static_cast<RealType>(u - 0.5L);
}

template<typename DataType>
DataType make_verification_entry(int64_t global_index, uint64_t seed) {
    const uint64_t base =
        splitmix64(seed ^ static_cast<uint64_t>(global_index));
    if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        return DataType(
            centered_uniform_from_hash<double>(base),
            centered_uniform_from_hash<double>(base ^ 0x6a09e667f3bcc909ULL));
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        return DataType(
            centered_uniform_from_hash<float>(base),
            centered_uniform_from_hash<float>(base ^ 0x6a09e667f3bcc909ULL));
    } else {
        return centered_uniform_from_hash<DataType>(base);
    }
}

template<typename CoordType, typename DataType>
std::vector<DataType> build_local_verification_vector(
    const ParallelTree<CoordType, DataType>* tree,
    uint64_t seed) {
    const int leaf_level = tree->num_levels - 1;
    const auto& level = tree->levels[leaf_level];

    int64_t total_points = 0;
    for (const auto& box : level.local_boxes) {
        total_points += box.num_points;
    }

    std::vector<DataType> values;
    values.reserve(static_cast<size_t>(total_points));
    for (const auto& box : level.local_boxes) {
        for (int64_t i = 0; i < box.num_points; ++i) {
            values.push_back(
                make_verification_entry<DataType>(box.point_indices[i], seed));
        }
    }
    return values;
}

template<typename DataType>
std::vector<DataType> build_global_verification_vector(
    int64_t num_points,
    uint64_t seed) {
    std::vector<DataType> values(static_cast<size_t>(num_points));
    for (int64_t i = 0; i < num_points; ++i) {
        values[static_cast<size_t>(i)] =
            make_verification_entry<DataType>(i, seed);
    }
    return values;
}

template<typename CoordType, typename DataType>
std::vector<DataType> extract_local_leaf_vector_from_global(
    const ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& global_vector) {
    const int leaf_level = tree->num_levels - 1;
    const auto& level = tree->levels[leaf_level];

    int64_t total_points = 0;
    for (const auto& box : level.local_boxes) {
        total_points += box.num_points;
    }

    std::vector<DataType> local_values;
    local_values.reserve(static_cast<size_t>(total_points));
    for (const auto& box : level.local_boxes) {
        for (int64_t i = 0; i < box.num_points; ++i) {
            const int64_t global_idx = box.point_indices[i];
            local_values.push_back(global_vector[static_cast<size_t>(global_idx)]);
        }
    }
    return local_values;
}

template<typename CoordType, typename DataType>
size_t compute_fft_matvec(
    const std::vector<DataType>& input,
    std::vector<DataType>& output,
    const std::vector<CoordType>& grid_points,
    int64_t grid_size,
    KernelKind kernel_kind,
    int dimension,
    double wave_divisor = 32.0,
    double length_scale = 0.1,
    double nugget = 1e-6,
    double kappa = 10.0) {
    const int64_t N = static_cast<int64_t>(input.size());
    output.assign(static_cast<size_t>(N), DataType{0.0});

    if (kernel_kind == KernelKind::LAPLACE) {
        if constexpr (std::is_same_v<DataType, double>) {
            if (dimension == 2) {
                kernel::LaplaceKernel2D_FFT<double> fft_kernel(
                    grid_points.data(), grid_size, N);
                fft_kernel.matvec(input.data(), output.data());
                return fft_kernel.memory_usage();
            }
            if (dimension == 3) {
                kernel::LaplaceKernel3D_FFT<double> fft_kernel(
                    grid_points.data(), grid_size, N);
                fft_kernel.matvec(input.data(), output.data());
                return fft_kernel.memory_usage();
            }
            throw std::invalid_argument(
                "compute_fft_matvec: Laplace FFT only supports dimension 2 or 3.");
        }
        throw std::runtime_error(
            "compute_fft_matvec: Laplace FFT currently supports only real double precision.");
    }

    if (kernel_kind == KernelKind::HELMHOLTZ) {
        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            if (dimension == 2) {
                kernel::HelmholtzKernel2D_FFT<CoordType, double> fft_kernel(
                    grid_points.data(), grid_size, N, wave_divisor);
                fft_kernel.matvec(input.data(), output.data());
                return fft_kernel.memory_usage();
            }
            if (dimension == 3) {
                kernel::HelmholtzKernel3D_FFT<CoordType, double> fft_kernel(
                    grid_points.data(), grid_size, N, wave_divisor);
                fft_kernel.matvec(input.data(), output.data());
                return fft_kernel.memory_usage();
            }
            throw std::invalid_argument(
                "compute_fft_matvec: Helmholtz FFT only supports dimension 2 or 3.");
        }
        throw std::runtime_error(
            "compute_fft_matvec: Helmholtz FFT currently supports only complex<double>.");
    }

    if (kernel_kind == KernelKind::MATERN52) {
        if constexpr (std::is_same_v<DataType, double>) {
            if (dimension == 2) {
                kernel::Matern52Kernel2D_FFT<double> fft_kernel(
                    grid_points.data(), grid_size, N, length_scale, nugget);
                fft_kernel.matvec(input.data(), output.data());
                return fft_kernel.memory_usage();
            }
            if (dimension == 3) {
                kernel::Matern52Kernel3D_FFT<double> fft_kernel(
                    grid_points.data(), grid_size, N, length_scale, nugget);
                fft_kernel.matvec(input.data(), output.data());
                return fft_kernel.memory_usage();
            }
            throw std::invalid_argument(
                "compute_fft_matvec: Matern52 FFT only supports dimension 2 or 3.");
        }
        throw std::runtime_error(
            "compute_fft_matvec: Matern52 FFT currently supports only real double precision.");
    }

    if (kernel_kind == KernelKind::YUKAWA) {
        if constexpr (std::is_same_v<DataType, double>) {
            if (dimension == 2) {
                kernel::YukawaKernel2D_FFT<double> fft_kernel(
                    grid_points.data(), grid_size, N, kappa);
                fft_kernel.matvec(input.data(), output.data());
                return fft_kernel.memory_usage();
            }
            if (dimension == 3) {
                kernel::YukawaKernel3D_FFT<double> fft_kernel(
                    grid_points.data(), grid_size, N, kappa);
                fft_kernel.matvec(input.data(), output.data());
                return fft_kernel.memory_usage();
            }
            throw std::invalid_argument(
                "compute_fft_matvec: Yukawa FFT only supports dimension 2 or 3.");
        }
        throw std::runtime_error(
            "compute_fft_matvec: Yukawa FFT currently supports only real double precision.");
    }

    throw std::invalid_argument("compute_fft_matvec: Unknown kernel kind.");
}


/**
* @brief Print vector values for debugging
*/
template<typename DataType>
void print_vector(
    const std::vector<DataType>& vec,
    const std::string& name,
    int64_t max_print = 20) {
    
    std::cout << "\n" << name << " (size " << vec.size() << "):" << std::endl;
    
    // Print first few entries
    int64_t n_print = std::min(max_print, static_cast<int64_t>(vec.size()));
    
    std::cout << "  First " << n_print << " entries:" << std::endl;
    for (int64_t i = 0; i < n_print; ++i) {
        std::cout << "    [" << std::setw(4) << i << "] = " 
                  << std::scientific << std::setprecision(6) << vec[i] << std::endl;
    }
    
    // // Print last few entries
    // if (vec.size() > max_print * 2) {
    //     std::cout << "  ..." << std::endl;
    //     std::cout << "  Last " << n_print << " entries:" << std::endl;
    //     for (int64_t i = vec.size() - n_print; i < vec.size(); ++i) {
    //         std::cout << "    [" << std::setw(4) << i << "] = " 
    //                   << std::scientific << std::setprecision(6) << vec[i] << std::endl;
    //     }
    // }
    
    // Compute statistics
    DataType min_val = *std::min_element(vec.begin(), vec.end());
    DataType max_val = *std::max_element(vec.begin(), vec.end());
    DataType sum = std::accumulate(vec.begin(), vec.end(), DataType{0.0});
    DataType mean = sum / vec.size();
    
    DataType variance = 0.0;
    for (const auto& val : vec) {
        variance += (val - mean) * (val - mean);
    }
    variance /= vec.size();
    DataType std_dev = std::sqrt(variance);
    
    std::cout << "\n  Statistics:" << std::endl;
    std::cout << "    Min:    " << std::scientific << min_val << std::endl;
    std::cout << "    Max:    " << max_val << std::endl;
    std::cout << "    Mean:   " << mean << std::endl;
    std::cout << "    Std:    " << std_dev << std::endl;
    std::cout << "    L2norm: " << std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0)) << std::endl;
}


/**
 * @brief Verify solution using direct BLAS matrix-vector product
 * 
 * Builds full dense matrix A and computes A*x using BLAS.
 * WARNING: O(N²) memory and O(N²) time - only for small problems!
 * 
 * @tparam CoordType Coordinate data type
 * @tparam DataType Matrix data type
 * @tparam KernelType Kernel evaluator type
 * @param kernel Kernel evaluator
 * @param rhs Original right-hand side vector
 * @param solution Computed solution vector
 * @param grid_points Grid point coordinates
 * @param N Number of points
 * @param verbose Print detailed output
 * @return Relative residual norm
 */
template<typename CoordType, typename DataType, typename KernelType>
double verify_solution_direct(
    KernelType* kernel,
    const std::vector<DataType>& rhs,
    const std::vector<DataType>& solution,
    const std::vector<CoordType>& grid_points,
    int64_t N,
    int dimension = 2,
    bool verbose = true) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank != 0) {
        return 0.0;  // Only verify on rank 0
    }
    
    if (verbose) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Solution Verification (Direct BLAS)" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Total DOFs: " << N << std::endl;
        
        double matrix_memory_mb = (static_cast<double>(N) * N * sizeof(DataType)) / (1024.0 * 1024.0);
        std::cout << "Matrix memory: " << matrix_memory_mb << " MB" << std::endl;
        
        if (N > 10000) {
            std::cout << "⚠ WARNING: N = " << N << " is large for direct method!" << std::endl;
        }
    }
    
    // ===== Step 1: Build full matrix A =====
    
    auto build_start = std::chrono::high_resolution_clock::now();
    
    std::vector<DataType> A(N * N);
    
    #pragma omp parallel for collapse(2) if(N > 1000)
    for (int64_t i = 0; i < N; ++i) {
        for (int64_t j = 0; j < N; ++j) {
            if (i == j) {
                // Diagonal entry
                A[i * N + j] = kernel->evaluate_diagonal(N);
            } else {
                // Off-diagonal entry
                CoordType xi[3], xj[3];
                for (int d = 0; d < dimension; ++d) {
                    xi[d] = grid_points[i * dimension + d];
                    xj[d] = grid_points[j * dimension + d];
                }
                
                A[i * N + j] = kernel->evaluate(xi, xj) / static_cast<DataType>(N);
            }
        }
    }
    
    auto build_end = std::chrono::high_resolution_clock::now();
    auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start);
    
    if (verbose) {
        std::cout << "Matrix assembly time: " << build_duration.count() << " ms" << std::endl;
    }
    
    // ===== Step 2: Compute A*x using BLAS =====
    
    auto matvec_start = std::chrono::high_resolution_clock::now();
    
    std::vector<DataType> Ax(N, DataType{0.0});
    
    if constexpr (std::is_same_v<DataType, double>) {
        // DGEMV: y = alpha*A*x + beta*y
        char trans = 'N';
        int m = static_cast<int>(N);
        int n = static_cast<int>(N);
        double alpha = 1.0;
        double beta = 0.0;
        int lda = static_cast<int>(N);
        int incx = 1;
        int incy = 1;

        dgemv_(&trans, &m, &n, &alpha, A.data(), &lda,
            solution.data(), &incx, &beta, Ax.data(), &incy);

    } else if constexpr (std::is_same_v<DataType, float>) {
        // SGEMV
        char trans = 'N';
        int m = static_cast<int>(N);
        int n = static_cast<int>(N);
        float alpha = 1.0f;
        float beta = 0.0f;
        int lda = static_cast<int>(N);
        int incx = 1;
        int incy = 1;

        sgemv_(&trans, &m, &n, &alpha, A.data(), &lda,
            solution.data(), &incx, &beta, Ax.data(), &incy);

    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        // ZGEMV
        char trans = 'N';
        int m = static_cast<int>(N);
        int n = static_cast<int>(N);
        std::complex<double> alpha(1.0, 0.0);
        std::complex<double> beta(0.0, 0.0);
        int lda = static_cast<int>(N);
        int incx = 1;
        int incy = 1;

        zgemv_(&trans, &m, &n, &alpha, A.data(), &lda,
            solution.data(), &incx, &beta, Ax.data(), &incy);

    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        // CGEMV
        char trans = 'N';
        int m = static_cast<int>(N);
        int n = static_cast<int>(N);
        std::complex<float> alpha(1.0f, 0.0f);
        std::complex<float> beta(0.0f, 0.0f);
        int lda = static_cast<int>(N);
        int incx = 1;
        int incy = 1;

        cgemv_(&trans, &m, &n, &alpha, A.data(), &lda,
            solution.data(), &incx, &beta, Ax.data(), &incy);

    } else {
        throw std::runtime_error("Unsupported DataType for BLAS");
    }
    
    auto matvec_end = std::chrono::high_resolution_clock::now();
    auto matvec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(matvec_end - matvec_start);
    
    if (verbose) {
        std::cout << "BLAS matvec time: " << matvec_duration.count() << " ms" << std::endl;
    }
    
    // ===== Step 3: Compute residual =====
    
    // Deduce the underlying real type
    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<
            std::is_same_v<DataType, std::complex<float>>, float,
            DataType
        >
    >;

    std::vector<DataType> residual(N);
    for (int64_t i = 0; i < N; ++i) {
        residual[i] = Ax[i] - rhs[i];
    }

    // Compute norms using BLAS
    RealType residual_norm = 0.0;
    RealType rhs_norm = 0.0;

    if constexpr (std::is_same_v<DataType, double>) {
        int n = static_cast<int>(N);
        int inc = 1;
        residual_norm = dnrm2_(&n, residual.data(), &inc);
        rhs_norm = dnrm2_(&n, rhs.data(), &inc);
    } else if constexpr (std::is_same_v<DataType, float>) {
        int n = static_cast<int>(N);
        int inc = 1;
        residual_norm = snrm2_(&n, residual.data(), &inc);
        rhs_norm = snrm2_(&n, rhs.data(), &inc);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        int n = static_cast<int>(N);
        int inc = 1;
        residual_norm = dznrm2_(&n, residual.data(), &inc);
        rhs_norm = dznrm2_(&n, rhs.data(), &inc);
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        int n = static_cast<int>(N);
        int inc = 1;
        residual_norm = scnrm2_(&n, residual.data(), &inc);
        rhs_norm = scnrm2_(&n, rhs.data(), &inc);
    }

    RealType relative_error = residual_norm / rhs_norm;

    if (verbose) {
        std::cout << "\nBackward Error Analysis:" << std::endl;
        std::cout << "  ||Ax - b||₂ = " << std::scientific << std::setprecision(6)
                  << residual_norm << std::endl;
        std::cout << "  ||b||₂      = " << rhs_norm << std::endl;
        std::cout << "  Relative residual = " << relative_error << std::endl;

        // Additional statistics
        RealType max_residual = 0.0;
        int64_t max_idx = 0;
        for (int64_t i = 0; i < N; ++i) {
            RealType abs_res = std::abs(residual[i]);
            if (abs_res > max_residual) {
                max_residual = abs_res;
                max_idx = i;
            }
        }

        std::cout << "  Max residual    = " << max_residual
                  << " (at index " << max_idx << ")" << std::endl;

        // Solution quality assessment
        if (relative_error < static_cast<RealType>(1e-10)) {
            std::cout << "\n  ✓ EXCELLENT: Relative error < 1e-10" << std::endl;
        } else if (relative_error < static_cast<RealType>(1e-6)) {
            std::cout << "\n  ✓ VERY GOOD: Relative error < 1e-6" << std::endl;
        } else if (relative_error < static_cast<RealType>(1e-3)) {
            std::cout << "\n  ✓ GOOD: Relative error < 1e-3" << std::endl;
        } else if (relative_error < static_cast<RealType>(1e-1)) {
            std::cout << "\n  ⚠ WARNING: Relative error > 1e-3" << std::endl;
        } else {
            std::cout << "\n  ✗ ERROR: Relative error > 1e-1 (solution likely incorrect)" << std::endl;
        }

        std::cout << "========================================\n" << std::endl;
    }
    
    return relative_error;
}


/**
 * @brief Verify solution using FFT-based matrix-vector product
 * 
 * Computes relative residual ||Ax - b|| / ||b|| where A*x is computed via FFT.
 * 
 * @tparam CoordType Coordinate data type
 * @tparam DataType Matrix data type
 * @param tree Hierarchical tree with solution
 * @param rhs Original right-hand side vector
 * @param solution Computed solution vector
 * @param grid_points Grid point coordinates (for FFT kernel)
 * @param grid_size Grid dimension n (n² = N)
 * @param verbose Print detailed output
 * @return Struct with relative_residual and backward_error
 */
struct SolveErrorMetrics {
    double relative_residual;  // ||Ax - b|| / ||b||
    double backward_error;     // ||Ax - b|| / (||A|| * ||x||), with ||A|| estimated as ||Ax||/||x||
};

template<typename CoordType, typename DataType>
SolveErrorMetrics verify_solution_fft(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& rhs,
    const std::vector<DataType>& solution,
    const std::vector<CoordType>& grid_points,
    int64_t grid_size,
    KernelKind kernel_kind,
    int dimension,
    double wave_divisor,
    double length_scale = 0.1,
    double nugget = 1e-6,
    double kappa = 10.0,
    bool verbose = true) {

    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<
            std::is_same_v<DataType, std::complex<float>>, float,
            DataType
        >
    >;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank != 0) {
        return {0.0, 0.0};  // Only verify on rank 0
    }

    int64_t N = solution.size();

    if (verbose) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Solution Verification (FFT)" << std::endl;
        std::cout << "========================================" << std::endl;
        if(dimension == 2) {
            std::cout << "Grid size: " << grid_size << " × " << grid_size << std::endl;
        } else if (dimension == 3) {
            std::cout << "Grid size: " << grid_size << " × " << grid_size << " × " << grid_size << std::endl;
        }
        std::cout << "Total DOFs: " << N << std::endl;
    }

    std::vector<DataType> Ax(N);
    size_t fft_memory_usage = 0;

    auto matvec_start = std::chrono::high_resolution_clock::now();
    fft_memory_usage = compute_fft_matvec(
        solution, Ax, grid_points, grid_size, kernel_kind, dimension,
        wave_divisor, length_scale, nugget, kappa);

    auto matvec_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> matvec_duration = matvec_end - matvec_start;

    if (verbose) {
        std::cout << "FFT kernel memory: "
                  << std::scientific << std::setprecision(6)
                  << (fft_memory_usage / 1024.0 / 1024.0)
                  << " MB" << std::endl;
        std::cout << "FFT matvec time: "
                  << std::defaultfloat << std::setprecision(6)
                  << matvec_duration.count() << " ms" << std::endl;
    }

    std::vector<DataType> residual(N);
    for (int64_t i = 0; i < N; ++i) {
        residual[i] = Ax[i] - rhs[i];
    }

    auto squared_magnitude = [](const DataType& value) -> RealType {
        if constexpr (std::is_same_v<DataType, std::complex<double>> ||
                      std::is_same_v<DataType, std::complex<float>>) {
            return static_cast<RealType>(std::norm(value));
        } else {
            return static_cast<RealType>(value * value);
        }
    };

    RealType residual_norm_sq = 0.0;
    RealType rhs_norm_sq = 0.0;
    RealType solution_norm_sq = 0.0;
    RealType Ax_norm_sq = 0.0;

    for (int64_t i = 0; i < N; ++i) {
        residual_norm_sq += squared_magnitude(residual[i]);
        rhs_norm_sq += squared_magnitude(rhs[i]);
        solution_norm_sq += squared_magnitude(solution[i]);
        Ax_norm_sq += squared_magnitude(Ax[i]);
    }

    RealType residual_norm = std::sqrt(residual_norm_sq);
    RealType rhs_norm = std::sqrt(rhs_norm_sq);
    RealType solution_norm = std::sqrt(solution_norm_sq);
    RealType Ax_norm = std::sqrt(Ax_norm_sq);
    double relative_error = std::numeric_limits<double>::infinity();
    if (rhs_norm > static_cast<RealType>(0.0)) {
        relative_error = static_cast<double>(residual_norm / rhs_norm);
    }
    double backward_error = std::numeric_limits<double>::infinity();
    if (Ax_norm > static_cast<RealType>(0.0)) {
        backward_error = static_cast<double>(residual_norm / Ax_norm);
    }

    if (verbose) {
        std::cout << "\nBackward Error Analysis:" << std::endl;
        std::cout << "  ||Ax - b||₂ = " << std::scientific << std::setprecision(6)
                  << residual_norm << std::endl;
        std::cout << "  ||b||₂      = " << rhs_norm << std::endl;
        std::cout << "  ||x||₂      = " << solution_norm << std::endl;
        std::cout << "  ||Ax||₂     = " << Ax_norm << std::endl;
        std::cout << "  Relative residual = " << relative_error << std::endl;
        std::cout << "  Backward error    = " << backward_error
                  << "  (using ||A|| ≈ ||Ax||/||x||)" << std::endl;

        RealType max_residual = 0.0;
        int64_t max_idx = 0;
        for (int64_t i = 0; i < N; ++i) {
            RealType abs_res = std::abs(residual[i]);
            if (abs_res > max_residual) {
                max_residual = abs_res;
                max_idx = i;
            }
        }

        std::cout << "  Max residual    = " << max_residual
                  << " (at index " << max_idx << ")" << std::endl;

        if (relative_error < 1e-10) {
            std::cout << "\n  ✓ EXCELLENT: Relative error < 1e-10" << std::endl;
        } else if (relative_error < 1e-6) {
            std::cout << "\n  ✓ EXCELLENT: Relative error < 1e-6" << std::endl;
        } else if (relative_error < 1e-3) {
            std::cout << "\n  ✓ GOOD: Relative error < 1e-3" << std::endl;
        } else if (relative_error < 1e-1) {
            std::cout << "\n  ⚠ WARNING: Relative error > 1e-3" << std::endl;
        } else {
            std::cout << "\n  ✗ ERROR: Relative error > 1e-1 (solution likely incorrect)" << std::endl;
        }

        std::cout << "========================================\n" << std::endl;
    }

    return {relative_error, backward_error};
}



/**
 * @brief Build grid points for FFT verification
 * 
 * Extracts grid coordinates from leaf-level boxes across all processes.
 * Gathers to rank 0. Assumes uniform grid structure.
 */
template<typename CoordType, typename DataType>
std::vector<CoordType> build_grid_points_from_tree(
    ParallelTree<CoordType, DataType>* tree,
    int64_t grid_size) {
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int leaf_level = tree->num_levels - 1;
    auto& leaf_level_ref = tree->levels[leaf_level];
    int dimension = tree->dimension;
    
    // ===== Step 1: Collect local indices and coordinates =====
    
    std::vector<int64_t> local_indices;
    std::vector<CoordType> local_coords;
    
    for (int64_t box_idx = 0; box_idx < leaf_level_ref.num_boxes_local; ++box_idx) {
        auto& box = leaf_level_ref.local_boxes[box_idx];
        
        for (int64_t i = 0; i < box.num_points; ++i) {
            local_indices.push_back(box.point_indices[i]);
            
            // Extract coordinates (dimension values per point)
            for (int d = 0; d < dimension; ++d) {
                local_coords.push_back(box.point_coords[i * dimension + d]);
            }
        }
    }
    
    // ===== Step 2: Gather counts to root =====
    
    int local_count = static_cast<int>(local_indices.size());
    std::vector<int> recv_counts(rank == 0 ? size : 0);
    std::vector<int> recv_displs(rank == 0 ? size : 0);
    
    MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // ===== Step 3: Calculate displacements on rank 0 =====
    
    int total_count = 0;
    if (rank == 0) {
        recv_displs[0] = 0;
        for (int i = 0; i < size; ++i) {
            if (i > 0) {
                recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
            }
            total_count += recv_counts[i];
        }
    }
    
    // ===== Step 4: Gather indices =====
    
    std::vector<int64_t> all_indices(rank == 0 ? total_count : 0);
    MPI_Gatherv(
        local_indices.data(), local_count, MPI_INT64_T,
        all_indices.data(), recv_counts.data(), recv_displs.data(), MPI_INT64_T,
        0, MPI_COMM_WORLD
    );
    
    // ===== Step 5: Gather coordinates =====
    
    MPI_Datatype mpi_coord_type;
    if constexpr (std::is_same_v<CoordType, double>) {
        mpi_coord_type = MPI_DOUBLE;
    } else if constexpr (std::is_same_v<CoordType, float>) {
        mpi_coord_type = MPI_FLOAT;
    } else {
        throw std::runtime_error("Unsupported CoordType for MPI");
    }
    
    // Adjust counts/displs for coordinate data (dimension values per point)
    std::vector<int> coord_recv_counts(rank == 0 ? size : 0);
    std::vector<int> coord_recv_displs(rank == 0 ? size : 0);
    
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            coord_recv_counts[i] = recv_counts[i] * dimension;
            coord_recv_displs[i] = recv_displs[i] * dimension;
        }
    }
    
    int local_coord_count = local_count * dimension;
    std::vector<CoordType> all_coords(rank == 0 ? total_count * dimension : 0);
    
    MPI_Gatherv(
        local_coords.data(), local_coord_count, mpi_coord_type,
        all_coords.data(), coord_recv_counts.data(), coord_recv_displs.data(), mpi_coord_type,
        0, MPI_COMM_WORLD
    );
    
    // ===== Step 6: Assemble grid points on rank 0 =====
    
    std::vector<CoordType> grid_points;
    
    if (rank == 0) {
        int64_t N = (dimension == 2) ? (grid_size * grid_size) : (grid_size * grid_size * grid_size);
        grid_points.resize(N * dimension, CoordType{0.0});
        
        for (int i = 0; i < total_count; ++i) {
            int64_t global_idx = all_indices[i];
            
            if (global_idx < 0 || global_idx >= N) {
                throw std::runtime_error(
                    "build_grid_points_from_tree: Invalid global index " + 
                    std::to_string(global_idx));
            }
            
            // Copy coordinates (interleaved: [x0, y0, x1, y1, ...] or [x0, y0, z0, x1, y1, z1, ...])
            for (int d = 0; d < dimension; ++d) {
                grid_points[global_idx * dimension + d] = all_coords[i * dimension + d];
            }
        }
        // print_vector(grid_points, "Grid Points", std::min<int64_t>(10, grid_points.size()));
    }
    
    
    return grid_points;
}


// /**
//  * @brief Serialize vector of BoxData into buffer for MPI communication
//  * @param boxes Vector of boxes to serialize
//  * @return Buffer containing serialized data
//  */
// template<typename CoordType, typename DataType>
// std::vector<char> serialize_boxes(const std::vector<BoxData<CoordType, DataType>>& boxes) {
//     // First pass: calculate total size needed
//     int64_t num_boxes = boxes.size();
//     int64_t total_size = sizeof(int64_t);  // Space for num_boxes
    
//     // Calculate size for each box (we need to do this by actually serializing, 
//     // or having a size calculation function - for now, use a temporary buffer approach)
//     std::vector<int64_t> box_sizes(num_boxes);
    
//     for (int64_t i = 0; i < num_boxes; ++i) {
//         // Estimate size conservatively (can be refined with actual size calculation)
//         // For now, serialize to temporary buffer to get exact size
//         std::vector<char> temp_buffer(1024 * 1024);  // 1MB temp buffer per box
//         char* end_ptr = serialize(boxes[i], temp_buffer.data());
//         box_sizes[i] = end_ptr - temp_buffer.data();
//         total_size += box_sizes[i];
//     }
    
//     // Allocate buffer
//     std::vector<char> buffer(total_size);
//     char* ptr = buffer.data();
    
//     // Write number of boxes
//     std::memcpy(ptr, &num_boxes, sizeof(int64_t));
//     ptr += sizeof(int64_t);
    
//     // Serialize each box
//     for (int64_t i = 0; i < num_boxes; ++i) {
//         ptr = serialize(boxes[i], ptr);
//     }
    
//     // Resize to actual used size
//     buffer.resize(ptr - buffer.data());
    
//     return buffer;
// }

/**
 * @brief Serialize vector of BoxData into buffer for MPI communication
 * @param boxes Vector of boxes to serialize
 * @return Buffer containing serialized data
 */
template<typename CoordType, typename DataType>
std::vector<char> serialize_boxes(const std::vector<BoxData<CoordType, DataType>>& boxes) {
    // First pass: calculate total size needed using get_serialized_size
    int64_t num_boxes = boxes.size();
    size_t total_size = sizeof(int64_t);  // Space for num_boxes
    
    std::vector<size_t> box_sizes(num_boxes);
    
    for (int64_t i = 0; i < num_boxes; ++i) {
        box_sizes[i] = get_serialized_size(boxes[i]);
        total_size += box_sizes[i];
    }
    
    // Allocate exact buffer size
    std::vector<char> buffer(total_size);
    char* ptr = buffer.data();
    
    // Write number of boxes
    std::memcpy(ptr, &num_boxes, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    // Serialize each box
    for (int64_t i = 0; i < num_boxes; ++i) {
        char* ptr_before = ptr;
        ptr = serialize(boxes[i], ptr);
        
        // Verify we wrote the expected size
        size_t bytes_written = ptr - ptr_before;
        if (bytes_written != box_sizes[i]) {
            throw std::runtime_error("Box " + std::to_string(i) + 
                                   " size mismatch: expected " + 
                                   std::to_string(box_sizes[i]) + 
                                   " but wrote " + std::to_string(bytes_written));
        }
    }
    
    // Verify total size
    size_t actual_size = ptr - buffer.data();
    if (actual_size != total_size) {
        throw std::runtime_error("Total size mismatch: expected " + 
                               std::to_string(total_size) + 
                               " but wrote " + std::to_string(actual_size));
    }
    
    return buffer;
}

/**
 * @brief Deserialize vector of BoxData from buffer
 * @param buffer Buffer containing serialized data
 * @return Vector of deserialized boxes
 */
template<typename CoordType, typename DataType>
std::vector<BoxData<CoordType, DataType>> deserialize_boxes(const std::vector<char>& buffer) {
    const char* ptr = buffer.data();
    
    // Read number of boxes
    int64_t num_boxes = 0;
    std::memcpy(&num_boxes, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    // Deserialize each box
    std::vector<BoxData<CoordType, DataType>> boxes;
    boxes.reserve(num_boxes);
    
    for (int64_t i = 0; i < num_boxes; ++i) {
        BoxData<CoordType, DataType> box;
        ptr = deserialize(box, ptr);
        boxes.push_back(std::move(box));
    }
    
    return boxes;
}




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
 * @param factorization_method Cholesky or None
 * @param unit_proxy_points Unit sphere proxy points
 * @param num_proxy Number of proxy points
 * @param proxy_radius Proxy sphere radius multiplier
 * @param verbose Print progress messages
 */
template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_factorization_parallel(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
    bool verbose = true) {
    using clock = std::chrono::high_resolution_clock;
    
    int rank = tree->mpi_rank;
    int size = tree->mpi_size;
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    
    int dimension = tree->dimension;
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    int num_children = (dimension == 2) ? 4 : 8;
    const int factorization_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);
    
    if (verbose && rank == factorization_header_rank) {
        const char* factorization_name =
            factorization_method == FactorizationMethod::CHOLESKY ? "Cholesky" :
            factorization_method == FactorizationMethod::LU ? "LU" : "None";
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

    for (int current_level = leaf_level; current_level >= 1; current_level--) {
        auto level_start = std::chrono::high_resolution_clock::now();
        clock::duration level_data_exchange{};
        clock::duration level_reduction{};
        
        auto& level = tree->levels[current_level];
        const int level_print_rank = smallest_active_rank(level);
        const LevelThreadPlan level_thread_plan =
            configure_level_thread_plan(dynamic_threading, tree->comm, level);
        
        if (!level.is_process_active) {
            if (verbose && rank == level_print_rank) {
                std::cout << "\n===== Level " << current_level << " (Process " << rank << " inactive) =====" << std::endl;
            }
        }

        // declare locks
        std::unordered_map<int64_t, omp_lock_t*> box_locks;

         // initialize locks
        if (level.is_process_active) {
            for (auto& box : level.local_boxes) {
                auto* lock = new omp_lock_t;
                omp_init_lock(lock);
                box_locks[box.morton_index] = lock;
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
        
        if (verbose && rank == level_print_rank) {
            std::cout << "\n===== Level " << current_level << " =====" << std::endl;
            std::cout << "Active processes: " << level.num_active_processes << " from rank: " << rank << std::endl;
            std::cout << "Boxes per process: " << level.num_boxes_local << std::endl;
        }
        print_level_thread_plan(
            dynamic_threading,
            tree->comm,
            current_level,
            "factorization",
            level_print_rank,
            rank,
            level_thread_plan,
            verbose);
        
        // ===== Step 1: Gather ghost and assisting boxes =====
        

        
        // ===== Step 2: Eliminate boxes in colored order =====
        std::chrono::milliseconds elim_duration{0};
        int64_t total_skeleton = 0;
        int64_t total_redundant = 0;
        
        if (current_level > 1 && level.is_process_active) {
            // Regular levels: do ID and elimination in colored order

            auto elim_start = std::chrono::high_resolution_clock::now();

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
            bool to_store = true;
            const bool store_interior_wave = true;

            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {

                const int  color_id_mod    = counter % num_colors;
                const bool is_interior     = (counter >= interior_start_loc);

                // ----------------------------------------------------------------
                // Communication / transport step (single-threaded)
                // ----------------------------------------------------------------
                if(!is_interior){
                    const auto comm_duration_raw = transport_and_apply_factor_updates_symmetric_onehop(
                    tree, current_level, kernel, pending_updates, false);
                    level_data_exchange += comm_duration_raw;
                    update_neighbor_slicing_for_level(level, is_symmetric);
                    auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_duration_raw);

                    if (verbose && rank == level_print_rank) {
                        std::cout << "  Comm time: "
                                << comm_duration.count()
                                << " ms" << std::endl;
                    }
                }else if(counter == interior_start_loc){
                    const std::vector<int> neighbor_ranks = compute_one_hop_neighbor_ranks(tree, level, current_level);
                    std::vector<int64_t> need_assist;
                    for (const auto& kv : level.assisting_box_points_for_kernel_evaluation) {
                        if(level.eliminated_boxes.count(kv.first) != 0){
                            need_assist.push_back(kv.first);
                        }
                    }
                    if (!need_assist.empty()) {
                        const auto comm_duration_raw = exchange_assisting_for_mortons_onehop(
                            tree, level, current_level, neighbor_ranks, need_assist);
                        level_data_exchange += comm_duration_raw;
                        auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_duration_raw);
                        if (verbose && rank == level_print_rank) {
                            std::cout << "  Comm time: "
                                    << comm_duration.count()
                                    << " ms" << std::endl;
                        }
                    }
                }
                

              

                const auto& color_list = color_bins[static_cast<size_t>(counter)];

                if (verbose && rank == level_print_rank) {
                    std::cout << "  Processing "
                            << (is_interior ? "interior sub-wave " : "boundary color ")
                            << color_id_mod
                            << " (" << color_list.size() << " boxes)..." << std::endl;
                }

                if (color_list.empty()) continue;

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

                #pragma omp parallel default(shared) if (color_list.size() > 1)
                {
                    const int tid = omp_get_thread_num();
                    FactorizationThreadScratch<CoordType, DataType> scratch;

                    #pragma omp for schedule(static)
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

                            gather_id_workspace(
                                &box, level, kernel,
                                unit_proxy_points.data(), num_proxy,
                                proxy_radius, is_symmetric,
                                scratch.workspace, scratch.workspace_rows, scratch.workspace_cols,
                                0,
                                box.on_boundary
                            );

                            thread_boundary_counts[static_cast<size_t>(tid)] += box.on_boundary;

                            compute_and_modify(dimension,
                                &box, level, kernel,
                                scratch,
                                tolerance, is_symmetric, is_hermitian,
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
                }

                if (wave_exception) {
                    std::rethrow_exception(wave_exception);
                }

                for (int local_boundary_count : thread_boundary_counts) {
                    boundary_count += local_boundary_count;
                }

                if (use_owner_deferred_xnn) {
                    for (int64_t morton_idx : color_list) {
                        level.eliminated_boxes.insert(morton_idx);
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

                    #pragma omp parallel default(shared) if (wave_xnn_candidate_boxes.size() > 1)
                    {
                        const int tid = omp_get_thread_num();
                        DeferredXnnOwnerScratch<DataType> scratch;

                        #pragma omp for schedule(static)
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
                    }

                    if (owner_exception) {
                        std::rethrow_exception(owner_exception);
                    }

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
                    }
                }

                // ----------------------------------------------------------------
                // Mark assisting boxes as eliminated after this wave.
                // Boundary assisting boxes are marked by their boundary color wave.
                // Interior assisting boxes are marked by their interior sub-wave.
                // ----------------------------------------------------------------
                for (const auto& kv : level.assisting_box_points_for_kernel_evaluation) {
                    const int  assisting_color    = static_cast<int>(kv.first & (num_colors - 1));
                    const bool is_boundary_assist = level.assisting_boxes[kv.second].on_boundary;

                    const bool mark =
                        (!is_interior && is_boundary_assist  && assisting_color == color_id_mod) ||
                        ( is_interior && !is_boundary_assist && assisting_color == color_id_mod);

                    if (mark) {
                        level.eliminated_boxes.insert(kv.first);
                    }
                }

                // Final transport after the very last wave
                if (counter == static_cast<int>(color_bins.size()) - 1) {
                    const auto final_comm_duration = transport_and_apply_factor_updates_symmetric_onehop(
                        tree, current_level, kernel, pending_updates, false);
                    level_data_exchange += final_comm_duration;
                    update_neighbor_slicing_for_level(level, is_symmetric);
                }
            }
            // // Regular levels: do ID and elimination in colored order
            
            // auto elim_start = std::chrono::high_resolution_clock::now();
            
            // std::vector<double> workspace;
            // std::vector<double> X_NN_full;
            // std::vector<double> sketch_storage;
            // PendingFactorUpdates<double> pending_updates;
            // int64_t ws_rows, ws_cols;
            
            // // Define color order: blue → orange → purple (3D) → green → interior
            // // New coloring scheme: 2^d colors (4 in 2D, 8 in 3D)
            // // color_id = morton_idx mod 2^d  (use bitmask since 2^d is power of 2)
            // const int num_colors = 1 << dimension;

            // // Bin ALL local boxes by color
            // // Boundary-only colored bins
            // std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));

            // // Interior boxes (local, in increasing Morton order by construction)
            // std::vector<int64_t> interior_boxes;
            // interior_boxes.reserve(level.local_boxes.size());

            // for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(level.local_boxes.size()); ++local_idx) {
            //     const int64_t morton_idx = level.local_morton_start + local_idx;
            //     const auto& box = level.local_boxes[static_cast<size_t>(local_idx)];

            //     if (box.on_boundary) {
            //         // 4-color in 2D / 8-color in 3D, based on morton % num_colors
            //         const int color_id = static_cast<int>(morton_idx & (num_colors - 1));
            //         color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
            //     } else {
            //         // Interior: keep a simple list (already in morton order)
            //         interior_boxes.push_back(morton_idx);
            //     }
            // }
            // color_bins.push_back(std::move(interior_boxes));
            // int64_t interior_loc = color_bins.size() - 1;

            // // for(auto &each_color : color_bins){
            // //     each_color.clear();
            // // }

            // // for (size_t c = 0; c < num_colors; ++c) {
            // //     int64_t start = (static_cast<int64_t>(c)     * level.local_boxes.size()) / static_cast<int64_t>(num_colors);
            // //     int64_t end   = (static_cast<int64_t>(c + 1) * level.local_boxes.size()) / static_cast<int64_t>(num_colors);

            // //     color_bins[c].reserve(static_cast<size_t>(end - start));
            // //     for (int64_t x = start; x < end; ++x) {
            // //         color_bins[c].push_back(x);
            // //     }
            // // }

            // // Process each color in order
            // int boundary_count = 0;
            // for (int counter = 0; counter < color_bins.size(); ++counter) {
            //     int color_id = counter;
                
                
            //     // If your algorithm requires communication/synchronization between colors,
            //     // this is the natural place to do it (per-color).
            //     // e.g. gather_ghost_and_assisting_boxes_factorization(...);
            //     // e.g. transport_and_apply_factor_updates_symmetric_onehop(...);
            //     auto comm_start = std::chrono::high_resolution_clock::now();
            //     transport_and_apply_factor_updates_symmetric_onehop(tree, current_level, kernel, pending_updates, false);
            //     update_neighbor_slicing_for_level(level, is_symmetric);
            //     auto comm_end = std::chrono::high_resolution_clock::now();
            //     auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_end - comm_start);
            //     // if(counter == 1 && rank == 0){
            //     //     fflush(stdout);
            //     //     assert(false);
            //     // }
            //     if (verbose && rank == 0) {
            //         std::cout << "  Ghost/assisting gather time: " << comm_duration.count() << " ms" << std::endl;
            //     }
            //     const auto& color_list = color_bins[static_cast<size_t>(color_id)];

            //     if (verbose && rank == 0) {
            //         std::cout << "  Processing color " << color_id
            //                 << " (" << color_list.size() << " boxes)..." << std::endl;
            //     }

            //     for (int64_t morton_idx : color_list) {
            //         // For this loop morton_idx should always be local, but keep the ghost fallback for safety
            //         BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
            //         if (box_ptr == nullptr) {
            //             box_ptr = level.find_ghost_box(morton_idx);
            //         }
            //         if (box_ptr == nullptr) {
            //             throw std::runtime_error(
            //                 "hierarchical_factorization_parallel: Morton index " +
            //                 std::to_string(morton_idx) + " not found in local_boxes or ghost_boxes at level " +
            //                 std::to_string(current_level) + " (color: " + std::to_string(color_id) + ")");
            //         }

            //         auto& box = *box_ptr;

            //         gather_id_workspace(
            //             &box, level, kernel,
            //             unit_proxy_points.data(), num_proxy,
            //             proxy_radius, is_symmetric,
            //             workspace, ws_rows, ws_cols,
            //             0,  // DEBUG
            //             box.on_boundary
            //         );

            //         boundary_count += box.on_boundary;

            //         compute_and_modify(dimension,
            //             &box, level, kernel,
            //             workspace, ws_rows, ws_cols,
            //             tolerance, is_symmetric, is_hermitian, X_NN_full, sketch_storage, &pending_updates,
            //             factorization_method
            //         );
                    
                    
            //         // if(box.morton_index == 65){
            //         //     fflush(stdout);
            //         //     assert(false);
            //         // }
                    
            //     }

            //     // update eliminated_boxes to include the incoming assisting box ids, otherwise it's going to think that they aren't eliminated
            //     for(const auto& kv : level.assisting_box_points_for_kernel_evaluation){
            //         if(counter == interior_loc){
            //             level.eliminated_boxes.insert(kv.first);
            //         }else{
            //             if ((kv.first % (1ULL << tree->dimension)) == color_id && level.assisting_boxes[kv.second].on_boundary){
            //                 level.eliminated_boxes.insert(kv.first);
            //             }
            //         }
                    
            //     }


            //     auto assist_start = std::chrono::high_resolution_clock::now();
                
                
            //     // verify_modified_interaction_symmetry(tree->levels[current_level], 1e-16, true);
            //     if(counter == color_bins.size() - 1) {
            //         transport_and_apply_factor_updates_symmetric_onehop(tree, current_level, kernel, pending_updates, false);
            //         update_neighbor_slicing_for_level(level, is_symmetric);
            //     }
                
            //     auto assist_end = std::chrono::high_resolution_clock::now();
            //     auto assist_duration = std::chrono::duration_cast<std::chrono::milliseconds>(assist_end - assist_start);
            //     if (verbose && rank == 0) {
            //         std::cout << "  Assisting boxes gather time: " << assist_duration.count() << " ms" << std::endl;
            //     }
                
                
            // }
            
            // int world_rank=0, world_size=0;
            // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

            // for (int counter = 0; counter < num_colors; ++counter) {
            //     int color_id = perm[counter];
            //     if (current_level == 4) color_id = counter;

            //     const auto& color_list = color_bins[(size_t)color_id];

            //     // (Optional) flush anything already pending before starting this color
            //     transport_and_apply_factor_updates_symmetric_onehop(tree, current_level, kernel,
            //                                                         pending_updates, /*is_hermitian=*/false);
                

            //     for (int turn = 0; turn < level.num_active_processes; ++turn) {
                    
            //         // ---------------- compute phase: ONLY rank==turn ----------------
            //         if (world_rank == turn) {
            //             if (verbose) {
            //                 std::cout << "[rank " << world_rank << "] color " << color_id
            //                         << " turn " << turn << " computing...\n";
            //             }
                        
            //             for (int64_t morton_idx : color_list) {
            //                 BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);

            //                 // IMPORTANT: do NOT fall back to ghost here
            //                 if (!box_ptr) continue; // not owned by this rank

            //                 auto& box = *box_ptr;
            //                 printf("rank %d processing box id: %d, color: %d, num_points: %d\n", world_rank, box.morton_index, color_id, static_cast<int>(box.num_points));
            //                 fflush(stdout);
                           
            //                 gather_id_workspace(&box, level, kernel,
            //                                     unit_proxy_points.data(), num_proxy,
            //                                     proxy_radius, is_symmetric,
            //                                     workspace, ws_rows, ws_cols,
            //                                     0, box.on_boundary);

            //                 compute_and_modify(dimension, &box, level, kernel,
            //                                 workspace, ws_rows, ws_cols,
            //                                 tolerance, is_symmetric, is_hermitian,
            //                                 X_NN_full, sketch_storage,
            //                                 &pending_updates,
            //                                 factorization_method);
            //                 // if(current_level == 4 && color_id == 2 && box.morton_index == 22){
            //                 //     fflush(stdout);
            //                 //     assert(false);
            //                 // }
            //                 if(box.morton_index == 65){
            //                     print_workspace_stats(level.local_boxes[1].far_field_modified_interactions[level.local_boxes[1].far_field_interaction_map[21]].A_NS.data, "A_NSaaaaaaaaaa (after gather_id_workspace)");
            //                     fflush(stdout);
            //                     // assert(false);
            //                 }
                            
            //             }
                        

            //             // Optional: inspect what I'm about to send
            //             // print_pending_factor_updates(pending_updates, std::cout,
            //             //     ("rank " + std::to_string(world_rank) + " outgoing: ").c_str());
            //         }
            //         if(world_rank == 1 && turn != 0){
            //                     print_workspace_stats(level.local_boxes[1].far_field_modified_interactions[level.local_boxes[1].far_field_interaction_map[21]].A_NS.data, "A_NSaaaaaaaaaa (after gather_id_workspace)");
            //                     fflush(stdout);
                                
            //             }
                    
            //         // update eliminated_boxes to include the incoming assisting box ids, otherwise it's going to think that they aren't eliminated
            //         // for(const auto& kv : level.assisting_box_points_for_kernel_evaluation){
            //         //     if ((kv.first % (1ULL << tree->dimension)) == color_id && (kv.first < ((turn + 1) * level.num_boxes_local))){
            //         //         level.eliminated_boxes.insert(kv.first);
            //         //     }
            //         // }
            //         if(current_level == 4){
            //             MPI_Barrier(tree->comm); // synchronize after each turn (optional, depends on algorithm requirements)
            //         }
                    

            //     }
                
            //     for(const auto& kv : level.assisting_box_points_for_kernel_evaluation){
            //         if ((kv.first % (1ULL << tree->dimension)) == color_id){
            //             level.eliminated_boxes.insert(kv.first);
            //         }
            //     }
 
                
            //     // ---------------- comm/apply phase: ALL ranks participate ----------------
            //     transport_and_apply_factor_updates_symmetric_onehop(tree, current_level, kernel,
            //                                                         pending_updates, /*is_hermitian=*/false);
                
            //     // If your slicing update depends on newly-arrived data, do it here (per turn)
            //     update_neighbor_slicing_for_level(level, is_symmetric);

                

            // }
            // printf("level: %d, on boundary sum: %d from rank: %d\n", current_level, boundary_count, rank);
            
            auto elim_end = std::chrono::high_resolution_clock::now();
            elim_duration = std::chrono::duration_cast<std::chrono::milliseconds>(elim_end - elim_start);

            for (const auto& box : level.local_boxes) {
                total_skeleton += box.skeleton_indices.size();
                total_redundant += box.redundant_indices.size();
            }
        } else {
            // Level 1: Skip elimination (only 4/8 boxes, no far-field)
            if (verbose && rank == level_print_rank) {
                std::cout << "  Skipping elimination at level 1 (final coarsening step)" << std::endl;
            }
        }

        if (current_level > 1) {
            double min_elim_ms = 0.0;
            double max_elim_ms = 0.0;
            reduce_active_duration_bounds_ms(
                tree->comm,
                level_print_rank,
                level.is_process_active,
                elim_duration,
                min_elim_ms,
                max_elim_ms);
            
            if (verbose && rank == level_print_rank) {
                std::cout << "  Elimination time: shortest=" << std::llround(min_elim_ms)
                          << " ms, longest=" << std::llround(max_elim_ms) << " ms" << std::endl;
                
                const double compression =
                    static_cast<double>(total_skeleton) / (total_skeleton + total_redundant);
                std::cout << "  Compression ratio: " << compression
                          << " (" << (compression * 100) << "%)" << std::endl;
                std::cout << "  Average skeleton size: "
                          << static_cast<double>(total_skeleton) / level.local_boxes.size()
                          << std::endl;
            }
        }
        
        // ===== Step 3: Gather assisting boxes post-elimination =====
        clear_ghosts(level);
        
        

        
        // ===== Step 4: Build parent level interactions =====
        
        // Special handling for level 1: set all DOFs as skeleton (no elimination)
        if (current_level == 1 && level.is_process_active) {
            if (verbose && rank == level_print_rank) {
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
        
        auto transition_start = std::chrono::high_resolution_clock::now();
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
        if (level.is_process_active) {
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
        
        
        
        
        auto transition_end = std::chrono::high_resolution_clock::now();
        auto transition_duration = std::chrono::duration_cast<std::chrono::milliseconds>(transition_end - transition_start);
        
        if (verbose && rank == level_print_rank) {
            std::cout << "  Level transition time: " << transition_duration.count() << " ms" << std::endl;
            std::cout << "  Parent boxes created: " << parent_boxes.size() << std::endl;
        }
        
        // ===== Step 5: Handle process reduction =====
        
        auto& parent_level = tree->levels[current_level - 1];
        
        // Error check: Level 2 → 1 should NOT trigger reduction
        if (current_level == 2 && level.parent_level_owner != rank && level.is_process_active) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Level 2 → 1 triggered process reduction! "
                "This is not allowed. parent_level_owner = " + std::to_string(level.parent_level_owner) +
                ", current rank = " + std::to_string(rank));
        }


        const bool reduction_occurred =
            (parent_level.num_active_processes != level.num_active_processes);
        const bool keep_local_parent_boxes =
            level.is_process_active && (level.parent_level_owner == rank);
        const bool sends_parent_boxes =
            reduction_occurred && level.is_process_active && level.parent_level_owner != rank;

        std::vector<char> send_buffer;
        int64_t send_buffer_size = 0;
        if (sends_parent_boxes) {
            // Pre-pack parent boxes before the synchronization point so the
            // communication timer excludes sender-side serialization work.
            send_buffer = serialize_boxes(parent_boxes);
            send_buffer_size = static_cast<int64_t>(send_buffer.size());
        }

        if (reduction_occurred) {
            // Synchronize after parent-box construction/packing so the reduction
            // communication timer does not count time spent waiting for slower
            // ranks that are still preparing their payloads.
            MPI_Barrier(tree->comm);
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
                
                std::vector<BoxData<CoordType, DataType>> child_parent_boxes =
                    deserialize_boxes<CoordType, DataType>(recv_buffer);
                
                all_parent_boxes.insert(
                    all_parent_boxes.end(),
                    std::make_move_iterator(child_parent_boxes.begin()),
                    std::make_move_iterator(child_parent_boxes.end())
                );
            }

            parent_level.local_boxes = std::move(all_parent_boxes);
        }

        if (sends_parent_boxes) {
            // This process does NOT own the parent - send to parent_level_owner
            segment_start = clock::now();
            MPI_Send(&send_buffer_size, 1, MPI_INT64_T, level.parent_level_owner, 0, tree->comm);

            // Send data
            MPI_Send_large(send_buffer.data(), send_buffer_size, MPI_CHAR, level.parent_level_owner, 1, tree->comm);
            level_reduction += (clock::now() - segment_start);

            // This process will become inactive at parent level
            parent_level.local_boxes.clear();

            if (verbose) {
                std::cout << "  Process " << rank << " sending " << parent_boxes.size() 
                          << " parent boxes to process " << level.parent_level_owner << std::endl;
            }
        }
        
        // Clear modified interaction matrices to free memory
        if (level.is_process_active) {
            clear_modified_interaction_matrices(level);
        }

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
            MPI_Allreduce(local_deltas, max_deltas, 2, MPI_DOUBLE, MPI_MAX, tree->comm);
            total_data_exchange_time += std::chrono::duration_cast<clock::duration>(
                std::chrono::duration<double, std::milli>(max_deltas[0]));
            total_reduction_time += std::chrono::duration_cast<clock::duration>(
                std::chrono::duration<double, std::milli>(max_deltas[1]));
        }

        auto level_end = std::chrono::high_resolution_clock::now();
        auto level_duration = std::chrono::duration_cast<std::chrono::milliseconds>(level_end - level_start);
        
        if (verbose && rank == level_print_rank) {
            std::cout << "  Total level time: " << level_duration.count() << " ms" << std::endl;
        }
        
    }
    
    // ===== Special handling for level 0 (root) =====
    
    const int root_print_rank = smallest_active_rank(tree->levels[0]);
    if (verbose && rank == root_print_rank) {
        std::cout << "\n===== Level 0 (Root) =====" << std::endl;
    }
    
    auto& root_level = tree->levels[0];
    
    if (!root_level.is_process_active || root_level.local_boxes.empty()) {
        if (verbose) {
            std::cout << "  Process " << rank << " has no boxes at root level" << std::endl;
        }
    } else {
        if (root_level.local_boxes.size() != 1) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Expected exactly 1 box at root level, got " +
                std::to_string(root_level.local_boxes.size()));
        }
        
        auto& root_box = root_level.local_boxes[0];
        
        if (verbose && rank == root_print_rank) {
            std::cout << "  Root box points: " << root_box.num_points << std::endl;
        }
        
        // At level 0, the assembled matrix is just the schur complement
        if (!root_box.schur_complement.is_allocated()) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Root box schur complement not allocated");
        }
        
        int64_t n = root_box.schur_complement.rows;
        
        if (verbose && rank == root_print_rank) {
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
            
            if (verbose && rank == root_print_rank) {
                std::cout << "  ✓ Root Cholesky factorization complete" << std::endl;
            }
            
        } else if (factorization_method == FactorizationMethod::LU) {
            root_box.X_RR.data = root_box.schur_complement.data;
            root_box.X_RR_pivots.resize(static_cast<size_t>(n));

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

            if (verbose && rank == root_print_rank) {
                std::cout << "  ✓ Root LU factorization complete" << std::endl;
            }

        } else if (factorization_method == FactorizationMethod::COMPLEX_SYM) {
            if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                root_box.X_RR.data = root_box.schur_complement.data;
                root_box.X_RR_pivots.resize(static_cast<size_t>(n));

                char uplo = 'L';
                int nn = n;
                int lwork = -1;
                int info = 0;
                std::vector<DataType> work(1);
                zsytrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn,
                        root_box.X_RR_pivots.data(), work.data(), &lwork, &info);
                lwork = static_cast<int>(work[0].real());
                work.resize(static_cast<size_t>(lwork));
                zsytrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn,
                        root_box.X_RR_pivots.data(), work.data(), &lwork, &info);

                if (info != 0) {
                    throw std::runtime_error(
                        "hierarchical_factorization_parallel: Bunch-Kaufman factorization of root failed at pivot " +
                        std::to_string(info));
                }

                root_box.X_RR.format = MatrixStorage<DataType>::BUNCH_KAUFMAN;

                if (verbose && rank == root_print_rank) {
                    std::cout << "  ✓ Root Bunch-Kaufman factorization complete" << std::endl;
                }
            } else {
                throw std::runtime_error("COMPLEX_SYM factorization only supported for complex<double>");
            }

        } else {
            // No factorization: just copy schur complement to X_RR
            root_box.X_RR.data = root_box.schur_complement.data;
            root_box.X_RR.format = MatrixStorage<DataType>::FULL;
            root_box.X_RR_pivots.clear();
            
            if (verbose && rank == root_print_rank) {
                std::cout << "  ✓ Root matrix copied (no factorization)" << std::endl;
            }
        }
        
        // Mark root as skeleton only (no redundant DOFs at this level)
        root_box.skeleton_indices.resize(n);
        for (int64_t i = 0; i < n; ++i) {
            root_box.skeleton_indices[i] = i;
        }
        root_box.redundant_indices.clear();
    }
    
    MPI_Barrier(tree->comm);

    const double data_exchange_ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            total_data_exchange_time).count();
    const double reduction_ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            total_reduction_time).count();

    if (verbose && rank == root_print_rank) {
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - total_time);
        std::cout << "\n========================================" << std::endl;
        std::cout << "✓ Hierarchical Factorization Complete" << std::endl;
        std::cout << "  total time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  data exchange communication time: " << std::llround(data_exchange_ms) << " ms" << std::endl;
        std::cout << "  process reduction communication time: " << std::llround(reduction_ms) << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    restore_base_process_affinity();
    clear_runtime_fmm_thread_count();
    destroy_dynamic_threading_context(dynamic_threading);


    size_t local_memory_usage = 0;
    for (int current_level = leaf_level; current_level >= 1; current_level--) {
    
        auto& level = tree->levels[current_level];
        for (auto& box : level.local_boxes) {
            local_memory_usage += calculate_box_data_size(box);
        }
    }
    printf("factorization memory usage on rank %d: %.2f GB\n", rank, local_memory_usage / (1024.0 * 1024.0 * 1024.0));
    fflush(stdout);
    
}


/**
 * @brief Hierarchical solve routine (parallel MPI version)
 */
template<typename CoordType, typename DataType>
void hierarchical_solve_parallel(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& rhs,
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> &solve_data,
    bool verbose = true) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    using clock = std::chrono::high_resolution_clock;
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    auto get_data_start = clock::now();
    clock::duration communication_total_forward{};
    clock::duration communication_total_backward{};
    
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    const int solve_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);
    const int solve_root_rank =
        smallest_active_rank(tree->levels[0]);
    const LevelThreadPlan solve_thread_plan =
        configure_static_process_thread_plan(dynamic_threading);
    
    if (verbose && rank == solve_header_rank) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Solve (Parallel MPI)" << std::endl;
        std::cout << "========================================" << std::endl;
        if (dynamic_threading.enabled) {
            std::cout << "Solve thread plan: threads=" << solve_thread_plan.threads;
            if (solve_thread_plan.cpu_begin >= 0 && solve_thread_plan.cpu_end >= 0) {
                std::cout << ", cpus=[" << solve_thread_plan.cpu_begin
                          << ", " << solve_thread_plan.cpu_end << "]";
            }
            std::cout << std::endl;
        }
    }
    
    // ===== Initialize solve data structures =====
    
    
    for (int level = 0; level < num_levels; ++level) {
        auto& tree_level = tree->levels[level];
        
        if (!tree_level.is_process_active) {
            continue;
        }
        
        solve_data[level].resize(tree_level.num_boxes_local);
        int64_t global_idx = 0;
        for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
            auto& box = tree_level.local_boxes[box_idx];
            auto& solve_box = solve_data[level][box_idx];
            
            // solve_box.morton_index = box.morton_index;
            // solve_box.source_rank = rank;
            // solve_box.right_side.resize(box.num_points, DataType{0.0});
            // solve_box.left_side.resize(box.num_points, DataType{0.0});
            solve_box.initialize(box.morton_index, rank, box.num_points);
            
            solve_box.skeleton_indices = box.skeleton_indices;
            solve_box.redundant_indices = box.redundant_indices;

            
            
            
            // Set RHS from global vector (only at leaf level)
            if (level == leaf_level) {
                for (int64_t i = 0; i < box.num_points; ++i) {
                    // int64_t global_idx = box.point_indices[i];
                    // solve_box.right_side[i] = rhs[global_idx];
                    solve_box.right_side[i] = rhs[global_idx];
                    global_idx++;
                }
                solve_box.left_side = solve_box.right_side;
            }
        }
    }

    
    if (verbose && rank == solve_header_rank) {
        std::cout << "Initialized solve data structures" << std::endl;
    }
    
    // ===== Phase 1: Forward Sweep (V^{-1}) with level transitions =====
    
    if (verbose && rank == solve_header_rank) {
        std::cout << "\n===== Phase 1: Forward Sweep (V^{-1}) =====" << std::endl;
    }
    
    auto forward_start = std::chrono::high_resolution_clock::now();
    PendingSolveUpdates<DataType> pending_solve;

    for (int level = leaf_level; level >= 1; level--) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);
        // Error check: level 2→1 should not trigger reduction
        if (level == 2) {
            auto& parent_level = tree->levels[1];
            if (tree_level.num_active_processes != parent_level.num_active_processes) {
                throw std::runtime_error("Reduction between level 2 and level 1 is not allowed");
            }
        }
        if (verbose && rank == level_print_rank) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local << " boxes" << " on rank: " << rank << std::endl;
        }
       
        
        // // 1. Apply V^{-1} by color (skip level 1 - no far-field elimination)
        // if (level >= 2) {
        //     // 2^d color-wave forward sweep (local boxes only)
        //     const int num_colors = 1 << tree->dimension;

        //     // Bin ALL local boxes by color
        //     // Boundary-only colored bins
        //     std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));

        //     // Interior boxes (local, in increasing Morton order by construction)
        //     std::vector<int64_t> interior_boxes;
        //     interior_boxes.reserve(tree_level.local_boxes.size());

        //     for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
        //         const int64_t morton_idx = tree_level.local_morton_start + local_idx;
        //         const auto& box = tree_level.local_boxes[static_cast<size_t>(local_idx)];

        //         if (box.on_boundary) {
        //             // 4-color in 2D / 8-color in 3D, based on morton % num_colors
        //             const int color_id = static_cast<int>(morton_idx & (num_colors - 1));
        //             color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
        //         } else {
        //             // Interior: keep a simple list (already in morton order)
        //             interior_boxes.push_back(morton_idx);
        //         }
        //     }
        //     color_bins.push_back(std::move(interior_boxes));

        //     // for(auto &each_color : color_bins){
        //     //     each_color.clear();
        //     // }

        //     // for (size_t c = 0; c < num_colors; ++c) {
        //     //     int64_t start = (static_cast<int64_t>(c)     * tree_level.local_boxes.size()) / static_cast<int64_t>(num_colors);
        //     //     int64_t end   = (static_cast<int64_t>(c + 1) * tree_level.local_boxes.size()) / static_cast<int64_t>(num_colors);

        //     //     color_bins[c].reserve(static_cast<size_t>(end - start));
        //     //     for (int64_t x = start; x < end; ++x) {
        //     //         color_bins[c].push_back(x);
        //     //     }
                
        //     // }

        //     // pending updates generated by forward elimination to nonlocal neighbors
            

        //     // Process waves color-by-color
        //     for (int color_id = 0; color_id < color_bins.size(); ++color_id) {

        //         // Apply elimination for all local boxes of this color
        //         for (int64_t morton_idx : color_bins[(size_t)color_id]) {
        //             const int64_t local_idx = morton_idx - tree_level.local_morton_start;

        //             apply_forward_elimination(
        //                 tree_level,
        //                 solve_data[level][(size_t)local_idx],
        //                 solve_data[level],
        //                 fmm::MatrixProperty::SYMMETRIC,
        //                 pending_solve,
        //                 /*is_ghost=*/false
        //             );
        //         }

        //         // After each color wave, ship accumulated nonlocal neighbor updates and apply received ones.
        //         // (This matches the "wave" semantics you used for factor updates.)
        //         transport_and_apply_solve_updates_onehop(
        //             tree,
        //             tree_level,
        //             level,
        //             pending_solve,
        //             solve_data
        //         );
                
        //     }
            
           
        // }

        if (level >= 2 && tree_level.is_process_active) {
            const int num_colors = 1 << tree->dimension;
            const int max_forward_threads = std::max(1, omp_get_max_threads());

            // ----------------------------------------------------------------
            // Build boundary color bins [0, num_colors-1]
            // Build interior sub-wave bins [num_colors, 2*num_colors-1]
            // Both use morton & (num_colors-1), matching factorization order.
            // ----------------------------------------------------------------
            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx = tree_level.local_morton_start + local_idx;
                const auto&   box        = tree_level.local_boxes[static_cast<size_t>(local_idx)];
                const int     color_id   = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }

            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }
            // color_bins layout:
            //   [0,           num_colors-1]: boundary waves
            //   [num_colors, 2*num_colors-1]: interior sub-waves

            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {
                const auto& color_list = color_bins[static_cast<size_t>(counter)];

                if (color_list.empty()) {
                    // Still need to transport even if no local boxes in this wave
                    communication_total_forward += transport_and_apply_solve_updates_onehop(
                        tree, tree_level, level, pending_solve, solve_data);
                    continue;
                }

                std::vector<PendingSolveUpdates<DataType>> thread_pending(
                    static_cast<size_t>(max_forward_threads));
                std::exception_ptr forward_wave_exception;
                std::mutex forward_wave_exception_mutex;
                std::atomic<bool> forward_wave_failed{false};

                #pragma omp parallel default(shared) if (color_list.size() > 1)
                {
                    const int tid = omp_get_thread_num();
                    auto& local_pending = thread_pending[static_cast<size_t>(tid)];

                    #pragma omp for schedule(static)
                    for (int64_t idx = 0; idx < static_cast<int64_t>(color_list.size()); ++idx) {
                        if (forward_wave_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }

                        try {
                            const int64_t morton_idx = color_list[static_cast<size_t>(idx)];
                            const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                            apply_forward_elimination(
                                tree_level,
                                solve_data[level][static_cast<size_t>(local_idx)],
                                solve_data[level],
                                fmm::MatrixProperty::SYMMETRIC,
                                local_pending,
                                /*is_ghost=*/false,
                                /*defer_local_updates=*/true
                            );
                        } catch (...) {
                            if (!forward_wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(forward_wave_exception_mutex);
                                forward_wave_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (forward_wave_exception) {
                    std::rethrow_exception(forward_wave_exception);
                }

                for (int t = 0; t < max_forward_threads; ++t) {
                    merge_pending_solve(pending_solve, thread_pending[static_cast<size_t>(t)]);
                    clear_pending_solve_updates_memory(thread_pending[static_cast<size_t>(t)]);
                }

                communication_total_forward += transport_and_apply_solve_updates_onehop(
                    tree, tree_level, level, pending_solve, solve_data);
            }
        }
        
        
        // Gather skeleton to parent level
        auto& parent_level = tree->levels[level - 1];
        get_data_start = clock::now();
        gather_skeleton_to_parent(
            tree_level,
            parent_level,
            solve_data[level],
            solve_data[level - 1],
            tree->dimension
        );
        communication_total_forward += (clock::now() - get_data_start);

       
   
    }
    

    auto forward_end = std::chrono::high_resolution_clock::now();
    auto forward_duration = std::chrono::duration_cast<std::chrono::milliseconds>(forward_end - forward_start - communication_total_forward);
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "  Forward sweep time: " << forward_duration.count() << " ms" << std::endl;
    }
    
    // ===== Phase 2: Diagonal Solve (D^{-1}) =====
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "\n===== Phase 2: Diagonal Solve (D^{-1}) =====" << std::endl;
    }
    
    auto diagonal_start = std::chrono::high_resolution_clock::now();
    
    int64_t num_diagonal_solves = 0;
    
    // Solve all X_RR blocks (levels N-1 down to 2)
    for (int level = leaf_level; level >= 2; level--) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);
        // std::cout << "  Diagonal solves before from rank: " << rank << std::endl;
        if (!tree_level.is_process_active) {
            continue;
        }

        std::exception_ptr diagonal_exception;
        std::mutex diagonal_exception_mutex;
        std::atomic<bool> diagonal_failed{false};

        #pragma omp parallel default(shared) if (tree_level.num_boxes_local > 1)
        {
            int64_t local_diagonal_solves = 0;

            #pragma omp for schedule(static)
            for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
                if (diagonal_failed.load(std::memory_order_relaxed)) {
                    continue;
                }

                try {
                    auto& box = tree_level.local_boxes[static_cast<size_t>(box_idx)];
                    auto& solve_box = solve_data[level][static_cast<size_t>(box_idx)];

                    if (box.redundant_indices.empty()) {
                        continue;
                    }

                    int64_t r = static_cast<int64_t>(box.redundant_indices.size());
                    std::vector<DataType> b_R(static_cast<size_t>(r));
                    for (int64_t i = 0; i < r; ++i) {
                        b_R[static_cast<size_t>(i)] =
                            solve_box.left_side[box.redundant_indices[static_cast<size_t>(i)]];
                    }

                    if (box.X_RR.format == MatrixStorage<DataType>::CHOLESKY_L) {
                        char uplo = 'L';
                        int n = static_cast<int>(r), nrhs = 1;
                        int lda = static_cast<int>(r), ldb = static_cast<int>(r), info = 0;

                        if constexpr (std::is_same_v<DataType, double>) {
                            dpotrs_(&uplo, &n, &nrhs,
                                    box.X_RR.data.data(), &lda,
                                    b_R.data(), &ldb, &info);
                        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            zsychol_solve_(&uplo, &n, &nrhs,
                                           box.X_RR.data.data(), &lda,
                                           b_R.data(), &ldb, &info);
                        }

                        if (info != 0) {
                            throw std::runtime_error("Diagonal solve failed for X_RR");
                        }
                    } else if (box.X_RR.format == MatrixStorage<DataType>::LU_FACTORED) {
                        if (box.X_RR_pivots.size() < static_cast<size_t>(r)) {
                            throw std::runtime_error("Diagonal solve missing LU pivots for X_RR");
                        }

                        char trans = 'N';
                        int n = static_cast<int>(r), nrhs = 1;
                        int lda = static_cast<int>(r), ldb = static_cast<int>(r), info = 0;

                        if constexpr (std::is_same_v<DataType, double>) {
                            dgetrs_(&trans, &n, &nrhs,
                                    box.X_RR.data.data(), &lda,
                                    box.X_RR_pivots.data(),
                                    b_R.data(), &ldb, &info);
                        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            zgetrs_(&trans, &n, &nrhs,
                                    box.X_RR.data.data(), &lda,
                                    box.X_RR_pivots.data(),
                                    b_R.data(), &ldb, &info);
                        }

                        if (info != 0) {
                            throw std::runtime_error("LU diagonal solve failed for X_RR");
                        }
                    } else if (box.X_RR.format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
                        if (box.X_RR_pivots.size() < static_cast<size_t>(r)) {
                            throw std::runtime_error("Diagonal solve missing Bunch-Kaufman pivots for X_RR");
                        }

                        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            char uplo = 'L';
                            int n = static_cast<int>(r), nrhs = 1;
                            int lda = static_cast<int>(r), ldb = static_cast<int>(r), info = 0;
                            zsytrs_(&uplo, &n, &nrhs,
                                    box.X_RR.data.data(), &lda,
                                    box.X_RR_pivots.data(),
                                    b_R.data(), &ldb, &info);
                            if (info != 0) {
                                throw std::runtime_error("Bunch-Kaufman diagonal solve failed for X_RR");
                            }
                        } else {
                            throw std::runtime_error("BUNCH_KAUFMAN format only supported for complex<double>");
                        }
                    } else {
                        throw std::runtime_error("Unsupported X_RR format in diagonal solve");
                    }

                    for (int64_t i = 0; i < r; ++i) {
                        solve_box.left_side[box.redundant_indices[static_cast<size_t>(i)]] =
                            b_R[static_cast<size_t>(i)];
                    }

                    local_diagonal_solves++;
                } catch (...) {
                    if (!diagonal_failed.exchange(true, std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lock(diagonal_exception_mutex);
                        diagonal_exception = std::current_exception();
                    }
                }
            }

            #pragma omp atomic
            num_diagonal_solves += local_diagonal_solves;
        }

        if (diagonal_exception) {
            std::rethrow_exception(diagonal_exception);
        }
    }
    
    // Solve root X_RR (level 0)
    auto& root_level = tree->levels[0];
    if (root_level.is_process_active && !root_level.local_boxes.empty()) {
        apply_diagonal_solve(root_level, solve_data[0][0], false);
        num_diagonal_solves++;
    }
    
    auto diagonal_end = std::chrono::high_resolution_clock::now();
    auto diagonal_duration = std::chrono::duration_cast<std::chrono::milliseconds>(diagonal_end - diagonal_start);
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "  Diagonal solves: " << num_diagonal_solves  << " from rank: " << rank << std::endl;
        std::cout << "  Diagonal solve time: " << diagonal_duration.count() << " ms" << std::endl;
    }
    
    
    // ===== Phase 3: Backward Sweep (W^{-1}) with level transitions =====
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "\n===== Phase 3: Backward Sweep (W^{-1}) =====" << std::endl;
    }
    
    auto backward_start = std::chrono::high_resolution_clock::now();
    
    for (int level = 1; level <= leaf_level; level++) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);
        
        if (verbose && rank == level_print_rank) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local << " boxes" << " from rank: " << rank << std::endl;
        }

        // for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
        //     auto& box = tree_level.local_boxes[box_idx];
        //     auto& solve_box = solve_data[level][box_idx];
            
        //     for(int gg = 0; gg < box.skeleton_indices.size(); gg++)
        //     {
        //         assert(solve_box.skeleton_indices[gg] == box.skeleton_indices[gg]);
        //     }
        // }
        
        // 1. Scatter solution from parent
        auto& parent_level = tree->levels[level - 1];
        get_data_start = clock::now();
        if (tree_level.is_process_active || parent_level.is_process_active) {
            scatter_solution_to_children(
                tree_level,
                parent_level,
                solve_data[level],
                solve_data[level - 1],
                tree->dimension
            );
        }
        communication_total_backward += (clock::now() - get_data_start);

        // for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
        //     auto& box = tree_level.local_boxes[box_idx];
        //     auto& solve_box = solve_data[level][box_idx];
            
        //     for(int gg = 0; gg < box.skeleton_indices.size(); gg++)
        //     {
        //         assert(solve_box.skeleton_indices[gg] == box.skeleton_indices[gg]);
        //     }
        // }
        
        if (verbose && rank == level_print_rank) {
            std::cout << "    ← Scattered from level " << (level - 1) << std::endl;
        }
        
        if (level >= 2 && tree_level.is_process_active) {
            const int num_colors = 1 << tree->dimension;

            // ----------------------------------------------------------------
            // Build boundary color bins [0, num_colors-1]
            // Build interior sub-wave bins [num_colors, 2*num_colors-1]
            // Matches factorization binning exactly - reversed on iteration.
            // ----------------------------------------------------------------
            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx = tree_level.local_morton_start + local_idx;
                const auto&   box        = tree_level.local_boxes[static_cast<size_t>(local_idx)];
                const int     color_id   = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }

            // Append interior sub-waves after boundary waves, matching forward order
            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }
            // color_bins layout (same as forward):
            //   [0,           num_colors-1]: boundary waves
            //   [num_colors, 2*num_colors-1]: interior sub-waves
            //
            // Backward sweep: iterate in reverse, interior sub-waves first,
            // then boundary waves. Within each wave, reverse Morton order.

            for (int counter = static_cast<int>(color_bins.size()) - 1; counter >= 0; --counter) {
                // Gather new updates before each wave (same position as forward transport)
                get_data_start = clock::now();
                gather_boxes_solve(tree, level, solve_data[level], true);
                communication_total_backward += (clock::now() - get_data_start);

                auto& bins = color_bins[static_cast<size_t>(counter)];

                if (bins.empty()) continue;

                std::exception_ptr backward_wave_exception;
                std::mutex backward_wave_exception_mutex;
                std::atomic<bool> backward_wave_failed{false};

                #pragma omp parallel default(shared) if (bins.size() > 1)
                {
                    #pragma omp for schedule(static)
                    for (int64_t idx = 0; idx < static_cast<int64_t>(bins.size()); ++idx) {
                        if (backward_wave_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }

                        try {
                            const int64_t morton_idx =
                                bins[bins.size() - 1 - static_cast<size_t>(idx)];
                            const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                            apply_backward_substitution(
                                tree_level,
                                solve_data[level][static_cast<size_t>(local_idx)],
                                solve_data[level],
                                fmm::MatrixProperty::SYMMETRIC,
                                /*is_ghost=*/false
                            );
                        } catch (...) {
                            if (!backward_wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(backward_wave_exception_mutex);
                                backward_wave_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (backward_wave_exception) {
                    std::rethrow_exception(backward_wave_exception);
                }
            }
        }
        // if (level >= 2) {
            
        //    // 2^d color-wave backward sweep (reverse order vs forward):
        //     // colors: (2^d - 1) ... 0
        //     // within each color: reverse Morton order
            
        //     const int num_colors = 1 << tree->dimension;

        //     // Bin ALL local boxes by color
        //     // Boundary-only colored bins
        //     std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));

        //     // Interior boxes (local, in increasing Morton order by construction)
        //     std::vector<int64_t> interior_boxes;
        //     interior_boxes.reserve(tree_level.local_boxes.size());

        //     for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
        //         const int64_t morton_idx = tree_level.local_morton_start + local_idx;
        //         const auto& box = tree_level.local_boxes[static_cast<size_t>(local_idx)];

        //         if (box.on_boundary) {
        //             // 4-color in 2D / 8-color in 3D, based on morton % num_colors
        //             const int color_id = static_cast<int>(morton_idx & (num_colors - 1));
        //             color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
        //         } else {
        //             // Interior: keep a simple list (already in morton order)
        //             interior_boxes.push_back(morton_idx);
        //         }
        //     }
        //     color_bins.push_back(std::move(interior_boxes));

        //     // for(auto &each_color : color_bins){
        //     //     each_color.clear();
        //     // }

        //     // for (size_t c = 0; c < num_colors; ++c) {
        //     //     int64_t start = (static_cast<int64_t>(c)     * tree_level.local_boxes.size()) / static_cast<int64_t>(num_colors);
        //     //     int64_t end   = (static_cast<int64_t>(c + 1) * tree_level.local_boxes.size()) / static_cast<int64_t>(num_colors);

        //     //     color_bins[c].reserve(static_cast<size_t>(end - start));
        //     //     for (int64_t x = start; x < end; ++x) {
        //     //         color_bins[c].push_back(x);
        //     //     }
        //     // }

        //     // Backward: reverse color order, and reverse within-color order
        //     for (int color_id = color_bins.size() - 1; color_id >= 0; --color_id) {
        //         // gather new updates after every wave
        //         gather_boxes_solve(tree, level, solve_data[level], true);
                
        //         auto& bins = color_bins[(size_t)color_id];
        //         for (auto it = bins.rbegin(); it != bins.rend(); ++it) {
        //             const int64_t morton_idx = *it;
        //             const int64_t local_idx  = morton_idx - tree_level.local_morton_start;

        //             apply_backward_substitution(
        //                 tree_level,
        //                 solve_data[level][(size_t)local_idx],
        //                 solve_data[level],
        //                 fmm::MatrixProperty::SYMMETRIC,
        //                 /*is_ghost=*/false
        //             );
        //         }

                
        //     }
            
           
        // }
    }
    
    auto backward_end = std::chrono::high_resolution_clock::now();
    auto backward_duration = std::chrono::duration_cast<std::chrono::milliseconds>(backward_end - backward_start - communication_total_backward);
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "  Backward sweep time: " << backward_duration.count() << " ms" << std::endl;
    }
    
    
    if (verbose && rank == solve_root_rank) {
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - forward_start);
        auto communication_duration_forward = std::chrono::duration_cast<std::chrono::milliseconds>(communication_total_forward);
        auto communication_duration_backward = std::chrono::duration_cast<std::chrono::milliseconds>(communication_total_backward);
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "Total solve time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  Forward:  " << forward_duration.count() << " ms" << std::endl;
        std::cout << "  Diagonal: " << diagonal_duration.count() << " ms" << std::endl;
        std::cout << "  Backward: " << backward_duration.count() << " ms" << std::endl;
        std::cout << "  Communication (forward):  " << communication_duration_forward.count() << " ms" << std::endl;
        std::cout << "  Communication (backward):  " << communication_duration_backward.count() << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    restore_base_process_affinity();
    clear_runtime_fmm_thread_count();
    destroy_dynamic_threading_context(dynamic_threading);
}


/**
 * @brief Hierarchical multiply: y = F * x, with F the approximate factorization.
 *
 * F = (prod_i V_i) * P_l * D * P_l^* * (prod_i W_i)
 *
 * Phase 1 (forward, leaf -> root): apply W operators. W reads from neighbors
 *   but does not write to them -> use the same gather-before-wave pattern as
 *   solve's backward sweep (W^{-1}).
 *
 * Phase 2 (diagonal): multiply X_RR blocks.
 *
 * Phase 3 (backward, root -> leaf): apply V operators. V writes to neighbors
 *   -> use the same pending-update / transport pattern as solve's forward
 *   sweep (V^{-1}). Waves are traversed in reverse order.
 */
template<typename CoordType, typename DataType>
void hierarchical_mul_parallel(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& input_vec,
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> &solve_data,
    bool verbose = true) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    using clock = std::chrono::high_resolution_clock;
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    auto get_data_start = clock::now();
    clock::duration communication_total_forward{};
    clock::duration communication_total_backward{};

    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    const int mul_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);
    const int mul_root_rank =
        smallest_active_rank(tree->levels[0]);
    const LevelThreadPlan mul_thread_plan =
        configure_static_process_thread_plan(dynamic_threading);

    if (verbose && rank == mul_header_rank) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Multiply (Parallel MPI)" << std::endl;
        std::cout << "========================================" << std::endl;
        if (dynamic_threading.enabled) {
            std::cout << "Multiply thread plan: threads=" << mul_thread_plan.threads;
            if (mul_thread_plan.cpu_begin >= 0 && mul_thread_plan.cpu_end >= 0) {
                std::cout << ", cpus=[" << mul_thread_plan.cpu_begin
                          << ", " << mul_thread_plan.cpu_end << "]";
            }
            std::cout << std::endl;
        }
    }

    // ===== Initialize solve data structures =====
    for (int level = 0; level < num_levels; ++level) {
        auto& tree_level = tree->levels[level];
        if (!tree_level.is_process_active) continue;

        solve_data[level].resize(tree_level.num_boxes_local);
        int64_t global_idx = 0;
        for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
            auto& box = tree_level.local_boxes[box_idx];
            auto& solve_box = solve_data[level][box_idx];

            solve_box.initialize(box.morton_index, rank, box.num_points);
            solve_box.skeleton_indices = box.skeleton_indices;
            solve_box.redundant_indices = box.redundant_indices;

            if (level == leaf_level) {
                for (int64_t i = 0; i < box.num_points; ++i) {
                    solve_box.right_side[i] = input_vec[global_idx];
                    global_idx++;
                }
                solve_box.left_side = solve_box.right_side;
            }
        }
    }

    if (verbose && rank == mul_header_rank) {
        std::cout << "Initialized multiply data structures" << std::endl;
    }

    // ===== Phase 1: Forward Sweep (W) =====
    if (verbose && rank == mul_header_rank) {
        std::cout << "\n===== Phase 1: Forward Sweep (W) =====" << std::endl;
    }

    auto forward_start = std::chrono::high_resolution_clock::now();

    for (int level = leaf_level; level >= 1; level--) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);

        if (level == 2) {
            auto& parent_level = tree->levels[1];
            if (tree_level.num_active_processes != parent_level.num_active_processes) {
                throw std::runtime_error("Reduction between level 2 and level 1 is not allowed");
            }
        }
        if (verbose && rank == level_print_rank) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local
                      << " boxes on rank: " << rank << std::endl;
        }

        if (level >= 2 && tree_level.is_process_active) {
            const int num_colors = 1 << tree->dimension;

            // Build boundary-then-interior sub-wave bins (same as solve forward).
            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx = tree_level.local_morton_start + local_idx;
                const auto& box = tree_level.local_boxes[static_cast<size_t>(local_idx)];
                const int color_id = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }
            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }

            // Apply W in forward color order (reversed of solve backward which
            // applies W^{-1}). Within each wave, gather fresh neighbor data
            // BEFORE applying (W is a gather, not scatter).
            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {
                // Fresh neighbor data for this wave.
                get_data_start = clock::now();
                gather_boxes_solve(tree, level, solve_data[level], true);
                communication_total_forward += (clock::now() - get_data_start);

                auto& bins = color_bins[static_cast<size_t>(counter)];
                if (bins.empty()) continue;

                std::exception_ptr wave_exception;
                std::mutex wave_exception_mutex;
                std::atomic<bool> wave_failed{false};

                #pragma omp parallel default(shared) if (bins.size() > 1)
                {
                    #pragma omp for schedule(static)
                    for (int64_t idx = 0; idx < static_cast<int64_t>(bins.size()); ++idx) {
                        if (wave_failed.load(std::memory_order_relaxed)) continue;

                        try {
                            const int64_t morton_idx = bins[static_cast<size_t>(idx)];
                            const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                            fmm::apply_mul_forward_W(
                                tree_level,
                                solve_data[level][static_cast<size_t>(local_idx)],
                                solve_data[level],
                                fmm::MatrixProperty::SYMMETRIC,
                                /*is_ghost=*/false);
                        } catch (...) {
                            if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(wave_exception_mutex);
                                wave_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (wave_exception) {
                    std::rethrow_exception(wave_exception);
                }
            }
        }

        // Gather skeleton values up to parent level (same as solve forward).
        auto& parent_level = tree->levels[level - 1];
        get_data_start = clock::now();
        gather_skeleton_to_parent(
            tree_level,
            parent_level,
            solve_data[level],
            solve_data[level - 1],
            tree->dimension);
        communication_total_forward += (clock::now() - get_data_start);
    }

    auto forward_end = std::chrono::high_resolution_clock::now();
    auto forward_duration = std::chrono::duration_cast<std::chrono::milliseconds>(forward_end - forward_start - communication_total_forward);

    if (verbose && rank == mul_root_rank) {
        std::cout << "  Forward sweep time: " << forward_duration.count() << " ms" << std::endl;
    }

    // ===== Phase 2: Diagonal Multiply (D) =====
    if (verbose && rank == mul_root_rank) {
        std::cout << "\n===== Phase 2: Diagonal Multiply (D) =====" << std::endl;
    }

    auto diagonal_start = std::chrono::high_resolution_clock::now();

    int64_t num_diagonal_muls = 0;

    for (int level = leaf_level; level >= 2; level--) {
        auto& tree_level = tree->levels[level];
        if (!tree_level.is_process_active) continue;

        std::exception_ptr diagonal_exception;
        std::mutex diagonal_exception_mutex;
        std::atomic<bool> diagonal_failed{false};

        #pragma omp parallel default(shared) if (tree_level.num_boxes_local > 1)
        {
            int64_t local_diagonal_muls = 0;

            #pragma omp for schedule(static)
            for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
                if (diagonal_failed.load(std::memory_order_relaxed)) continue;

                try {
                    auto& box = tree_level.local_boxes[static_cast<size_t>(box_idx)];
                    auto& solve_box = solve_data[level][static_cast<size_t>(box_idx)];

                    if (box.redundant_indices.empty()) continue;

                    int64_t r = static_cast<int64_t>(box.redundant_indices.size());
                    std::vector<DataType> x_R(static_cast<size_t>(r));
                    for (int64_t i = 0; i < r; ++i) {
                        x_R[static_cast<size_t>(i)] =
                            solve_box.left_side[box.redundant_indices[static_cast<size_t>(i)]];
                    }

                    int n = static_cast<int>(r);
                    int incx = 1;

                    if (box.X_RR.format == MatrixStorage<DataType>::CHOLESKY_L) {
                        char uplo = 'L', diag_N = 'N';
                        int lda = static_cast<int>(r);

                        if constexpr (std::is_same_v<DataType, double>) {
                            char trans_T = 'T';
                            dtrmv_(&uplo, &trans_T, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            char trans_N = 'N';
                            dtrmv_(&uplo, &trans_N, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            char trans_T = 'T';
                            ztrmv_(&uplo, &trans_T, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            char trans_N = 'N';
                            ztrmv_(&uplo, &trans_N, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                        }
                    } else if (box.X_RR.format == MatrixStorage<DataType>::LU_FACTORED) {
                        if (box.X_RR_pivots.size() < static_cast<size_t>(r)) {
                            throw std::runtime_error("Diagonal multiply missing LU pivots");
                        }
                        int lda = static_cast<int>(r);

                        if constexpr (std::is_same_v<DataType, double>) {
                            char uplo_U = 'U', trans_N = 'N', diag_N = 'N';
                            dtrmv_(&uplo_U, &trans_N, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            char uplo_L = 'L', diag_U = 'U';
                            dtrmv_(&uplo_L, &trans_N, &diag_U, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            // P * z (not P^T * z) — use INCX=-1.
                            int nrhs = 1, k1 = 1, k2 = n, inc_rev = -1;
                            dlaswp_(&nrhs, x_R.data(), &n,
                                    &k1, &k2, box.X_RR_pivots.data(), &inc_rev);
                        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            char uplo_U = 'U', trans_N = 'N', diag_N = 'N';
                            ztrmv_(&uplo_U, &trans_N, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            char uplo_L = 'L', diag_U = 'U';
                            ztrmv_(&uplo_L, &trans_N, &diag_U, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            int nrhs = 1, k1 = 1, k2 = n, inc_rev = -1;
                            zlaswp_(&nrhs, x_R.data(), &n,
                                    &k1, &k2, box.X_RR_pivots.data(), &inc_rev);
                        }
                    } else if (box.X_RR.format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
                        if (box.X_RR_pivots.size() < static_cast<size_t>(r)) {
                            throw std::runtime_error("Diagonal multiply missing Bunch-Kaufman pivots");
                        }
                        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            fmm::bunch_kaufman_multiply(n, box.X_RR.data.data(), n,
                                                       box.X_RR_pivots.data(), x_R.data());
                        } else {
                            throw std::runtime_error("BUNCH_KAUFMAN format only supported for complex<double>");
                        }
                    } else {
                        throw std::runtime_error("Diagonal multiply: unsupported X_RR format");
                    }

                    for (int64_t i = 0; i < r; ++i) {
                        solve_box.left_side[box.redundant_indices[static_cast<size_t>(i)]] =
                            x_R[static_cast<size_t>(i)];
                    }

                    local_diagonal_muls++;
                } catch (...) {
                    if (!diagonal_failed.exchange(true, std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lock(diagonal_exception_mutex);
                        diagonal_exception = std::current_exception();
                    }
                }
            }

            #pragma omp atomic
            num_diagonal_muls += local_diagonal_muls;
        }

        if (diagonal_exception) {
            std::rethrow_exception(diagonal_exception);
        }
    }

    // Root X_RR multiply
    auto& root_level = tree->levels[0];
    if (root_level.is_process_active && !root_level.local_boxes.empty()) {
        fmm::apply_diagonal_multiply(root_level, solve_data[0][0], false);
        num_diagonal_muls++;
    }

    auto diagonal_end = std::chrono::high_resolution_clock::now();
    auto diagonal_duration = std::chrono::duration_cast<std::chrono::milliseconds>(diagonal_end - diagonal_start);

    if (verbose && rank == mul_root_rank) {
        std::cout << "  Diagonal multiplies: " << num_diagonal_muls << " from rank: " << rank << std::endl;
        std::cout << "  Diagonal multiply time: " << diagonal_duration.count() << " ms" << std::endl;
    }

    // ===== Phase 3: Backward Sweep (V) =====
    if (verbose && rank == mul_root_rank) {
        std::cout << "\n===== Phase 3: Backward Sweep (V) =====" << std::endl;
    }

    auto backward_start = std::chrono::high_resolution_clock::now();

    for (int level = 1; level <= leaf_level; level++) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);

        if (verbose && rank == level_print_rank) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local
                      << " boxes from rank: " << rank << std::endl;
        }

        // Scatter solution from parent (same as solve backward).
        auto& parent_level = tree->levels[level - 1];
        get_data_start = clock::now();
        if (tree_level.is_process_active || parent_level.is_process_active) {
            scatter_solution_to_children(
                tree_level,
                parent_level,
                solve_data[level],
                solve_data[level - 1],
                tree->dimension);
        }
        communication_total_backward += (clock::now() - get_data_start);

        if (verbose && rank == level_print_rank) {
            std::cout << "    <- Scattered from level " << (level - 1) << std::endl;
        }

        if (level >= 2 && tree_level.is_process_active) {
            const int num_colors = 1 << tree->dimension;

            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx = tree_level.local_morton_start + local_idx;
                const auto& box = tree_level.local_boxes[static_cast<size_t>(local_idx)];
                const int color_id = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }
            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }

            // Apply V in REVERSE wave order (mirror of solve forward's V^{-1}).
            // V is a scatter: generate pending updates, then transport.
            const int max_backward_threads = std::max(1, omp_get_max_threads());
            PendingSolveUpdates<DataType> pending_mul;

            for (int counter = static_cast<int>(color_bins.size()) - 1; counter >= 0; --counter) {
                auto& bins = color_bins[static_cast<size_t>(counter)];

                if (bins.empty()) {
                    // Still participate in transport so neighbors don't deadlock.
                    communication_total_backward += transport_and_apply_solve_updates_onehop(
                        tree, tree_level, level, pending_mul, solve_data);
                    continue;
                }

                std::vector<PendingSolveUpdates<DataType>> thread_pending(
                    static_cast<size_t>(max_backward_threads));
                std::exception_ptr wave_exception;
                std::mutex wave_exception_mutex;
                std::atomic<bool> wave_failed{false};

                #pragma omp parallel default(shared) if (bins.size() > 1)
                {
                    const int tid = omp_get_thread_num();
                    auto& local_pending = thread_pending[static_cast<size_t>(tid)];

                    #pragma omp for schedule(static)
                    for (int64_t idx = 0; idx < static_cast<int64_t>(bins.size()); ++idx) {
                        if (wave_failed.load(std::memory_order_relaxed)) continue;

                        try {
                            // Mirror solve backward's reverse-within-wave order.
                            const int64_t morton_idx =
                                bins[bins.size() - 1 - static_cast<size_t>(idx)];
                            const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                            fmm::apply_mul_backward_V_with_pending(
                                tree_level,
                                solve_data[level][static_cast<size_t>(local_idx)],
                                solve_data[level],
                                fmm::MatrixProperty::SYMMETRIC,
                                local_pending,
                                /*is_ghost=*/false);
                        } catch (...) {
                            if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(wave_exception_mutex);
                                wave_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (wave_exception) {
                    std::rethrow_exception(wave_exception);
                }

                for (int t = 0; t < max_backward_threads; ++t) {
                    merge_pending_solve(pending_mul, thread_pending[static_cast<size_t>(t)]);
                    clear_pending_solve_updates_memory(thread_pending[static_cast<size_t>(t)]);
                }

                communication_total_backward += transport_and_apply_solve_updates_onehop(
                    tree, tree_level, level, pending_mul, solve_data);
            }
        }
    }

    auto backward_end = std::chrono::high_resolution_clock::now();
    auto backward_duration = std::chrono::duration_cast<std::chrono::milliseconds>(backward_end - backward_start - communication_total_backward);

    if (verbose && rank == mul_root_rank) {
        std::cout << "  Backward sweep time: " << backward_duration.count() << " ms" << std::endl;
    }

    if (verbose && rank == mul_root_rank) {
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - forward_start);
        auto communication_duration_forward = std::chrono::duration_cast<std::chrono::milliseconds>(communication_total_forward);
        auto communication_duration_backward = std::chrono::duration_cast<std::chrono::milliseconds>(communication_total_backward);
        std::cout << "\n========================================" << std::endl;
        std::cout << "Total multiply time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  Forward:  " << forward_duration.count() << " ms" << std::endl;
        std::cout << "  Diagonal: " << diagonal_duration.count() << " ms" << std::endl;
        std::cout << "  Backward: " << backward_duration.count() << " ms" << std::endl;
        std::cout << "  Communication (forward):  " << communication_duration_forward.count() << " ms" << std::endl;
        std::cout << "  Communication (backward):  " << communication_duration_backward.count() << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    restore_base_process_affinity();
    clear_runtime_fmm_thread_count();
    destroy_dynamic_threading_context(dynamic_threading);
}


/**
 * @brief Gather solution and RHS from all processes to process 0
 * 
 * @tparam CoordType Coordinate type (float or double)
 * @tparam DataType Data type for matrix entries
 * 
 * @param tree Parallel tree structure
 * @param solve_data Solve data at all levels (only leaf level used)
 * @param solution Output solution vector (only populated on rank 0)
 * @param aggregated_rhs Output RHS vector (only populated on rank 0)
 * 
 * Gathers point_indices, left_side, and right_side from all leaf boxes 
 * across all processes and assembles into global vectors on rank 0.
 */
template<typename CoordType, typename DataType>
void gather_solution_to_root(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& solve_data,
    std::vector<DataType>& solution,
    std::vector<DataType>& aggregated_rhs) {
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int leaf_level = tree->num_levels - 1;
    auto& leaf_level_ref = tree->levels[leaf_level];

    const int64_t local_solve_levels = static_cast<int64_t>(solve_data.size());
    if (leaf_level < 0 || leaf_level >= static_cast<int>(solve_data.size())) {
        std::fprintf(stderr,
                     "gather_solution_to_root: rank=%d invalid leaf level %d for solve_data size %lld\n",
                     rank,
                     leaf_level,
                     static_cast<long long>(local_solve_levels));
        std::fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const int64_t expected_leaf_boxes = leaf_level_ref.num_boxes_local;
    const int64_t actual_leaf_boxes = static_cast<int64_t>(solve_data[leaf_level].size());
    if (actual_leaf_boxes != expected_leaf_boxes) {
        std::fprintf(stderr,
                     "gather_solution_to_root: rank=%d solve_data leaf mismatch: expected %lld boxes, got %lld entries\n",
                     rank,
                     static_cast<long long>(expected_leaf_boxes),
                     static_cast<long long>(actual_leaf_boxes));
        std::fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // ===== Step 1: Collect local indices, solution values, and RHS values =====
    
    std::vector<int64_t> local_indices;
    std::vector<DataType> local_solution;
    std::vector<DataType> local_rhs;
    
    for (int64_t box_idx = 0; box_idx < leaf_level_ref.num_boxes_local; ++box_idx) {
        auto& box = leaf_level_ref.local_boxes[box_idx];
        auto& solve_box = solve_data[leaf_level][box_idx];
        
        for (int64_t i = 0; i < box.num_points; ++i) {
            local_indices.push_back(box.point_indices[i]);
            local_solution.push_back(solve_box.left_side[i]);
            local_rhs.push_back(solve_box.right_side[i]);
        }
    }
    
    // ===== Step 2: Gather counts on root =====

    const size_t local_count = local_indices.size();
    if (local_count > static_cast<size_t>(INT_MAX)) {
        throw std::runtime_error(
            "gather_solution_to_root: local verification gather count exceeds INT_MAX on rank " +
            std::to_string(rank) + ": " + std::to_string(local_count));
    }

    const int local_count_int = static_cast<int>(local_count);
    std::vector<int> recv_counts(rank == 0 ? size : 0, 0);
    int err = MPI_Gather(&local_count_int,
                         1,
                         MPI_INT,
                         rank == 0 ? recv_counts.data() : nullptr,
                         1,
                         MPI_INT,
                         0,
                         MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            "gather_solution_to_root: MPI_Gather(counts) failed with error " +
            std::string(errbuf, errlen));
    }
    
    // ===== Step 3: Calculate displacements on rank 0 =====
    
    int total_count = 0;
    std::vector<int> recv_displs(rank == 0 ? size : 0, 0);
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            if (i > 0) {
                recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
            }
            if (recv_counts[i] < 0 || total_count > INT_MAX - recv_counts[i]) {
                throw std::runtime_error(
                    "gather_solution_to_root: total verification gather count exceeds INT_MAX");
            }
            total_count += recv_counts[i];
        }
    }
    
    // ===== Step 4: Gather indices, solution, and RHS on root =====
    
    const MPI_Datatype mpi_datatype = mpi_datatype_for<DataType>();

    std::vector<int64_t> all_indices(rank == 0 ? static_cast<size_t>(total_count) : 0);
    std::vector<DataType> all_solution(rank == 0 ? static_cast<size_t>(total_count) : 0);
    std::vector<DataType> all_rhs(rank == 0 ? static_cast<size_t>(total_count) : 0);

    err = MPI_Gatherv(local_indices.empty() ? nullptr : local_indices.data(),
                      local_count_int,
                      MPI_INT64_T,
                      rank == 0 ? all_indices.data() : nullptr,
                      rank == 0 ? recv_counts.data() : nullptr,
                      rank == 0 ? recv_displs.data() : nullptr,
                      MPI_INT64_T,
                      0,
                      MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            "gather_solution_to_root: MPI_Gatherv(indices) failed with error " +
            std::string(errbuf, errlen));
    }

    err = MPI_Gatherv(local_solution.empty() ? nullptr : local_solution.data(),
                      local_count_int,
                      mpi_datatype,
                      rank == 0 ? all_solution.data() : nullptr,
                      rank == 0 ? recv_counts.data() : nullptr,
                      rank == 0 ? recv_displs.data() : nullptr,
                      mpi_datatype,
                      0,
                      MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            "gather_solution_to_root: MPI_Gatherv(solution) failed with error " +
            std::string(errbuf, errlen));
    }

    err = MPI_Gatherv(local_rhs.empty() ? nullptr : local_rhs.data(),
                      local_count_int,
                      mpi_datatype,
                      rank == 0 ? all_rhs.data() : nullptr,
                      rank == 0 ? recv_counts.data() : nullptr,
                      rank == 0 ? recv_displs.data() : nullptr,
                      mpi_datatype,
                      0,
                      MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            "gather_solution_to_root: MPI_Gatherv(rhs) failed with error " +
            std::string(errbuf, errlen));
    }
    
    // ===== Step 5: Assemble solution and RHS on rank 0 =====
    
    if (rank == 0) {
        solution.resize(tree->num_points, DataType{0.0});
        aggregated_rhs.resize(tree->num_points, DataType{0.0});
        
        for (int i = 0; i < total_count; ++i) {
            int64_t global_idx = all_indices[i];
            
            if (global_idx < 0 || global_idx >= tree->num_points) {
                throw std::runtime_error(
                    "gather_solution_to_root: Invalid global index " + 
                    std::to_string(global_idx));
            }
            
            solution[global_idx] = all_solution[i];
            aggregated_rhs[global_idx] = all_rhs[i];
        }
    }
}

/**
 * @brief Verify forward (multiply) error: ||F*x - A*x|| / ||A*x||
 *
 * Generates random x, computes y_true = A*x via FFT, computes y_hat = F*x
 * via hierarchical_mul_parallel, and measures the relative error.
 * This directly measures factorization quality (epsilon) without
 * condition number amplification.
 */
template<typename CoordType, typename DataType>
double verify_forward_error_fft(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<CoordType>& grid_points,
    int64_t grid_size,
    KernelKind kernel_kind,
    int dimension,
    double wave_divisor,
    double length_scale = 0.1,
    double nugget = 1e-6,
    double kappa = 10.0,
    bool verbose = true) {
    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<
            std::is_same_v<DataType, std::complex<float>>, float,
            DataType
        >
    >;

    auto squared_magnitude = [](const DataType& value) -> RealType {
        if constexpr (std::is_same_v<DataType, std::complex<double>> ||
                      std::is_same_v<DataType, std::complex<float>>) {
            return static_cast<RealType>(std::norm(value));
        } else {
            return static_cast<RealType>(value * value);
        }
    };

    constexpr uint64_t forward_seed = 0xa1b2c3d4e5f60718ULL;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0 && verbose) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Forward Error Verification (FFT)" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Using random test vector";
        if constexpr (std::is_same_v<DataType, std::complex<double>> ||
                      std::is_same_v<DataType, std::complex<float>>) {
            std::cout << " with independent real/imaginary parts";
        }
        std::cout << std::endl;
    }

    std::vector<DataType> x_test;
    std::vector<DataType> y_true;
    size_t fft_memory_usage = 0;
    double fft_matvec_ms = 0.0;

    if (rank == 0) {
        x_test = build_global_verification_vector<DataType>(tree->num_points, forward_seed);
        const auto fft_start = std::chrono::high_resolution_clock::now();
        fft_memory_usage = compute_fft_matvec(
            x_test, y_true, grid_points, grid_size, kernel_kind, dimension,
            wave_divisor, length_scale, nugget, kappa);
        const auto fft_end = std::chrono::high_resolution_clock::now();
        fft_matvec_ms =
            std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
                fft_end - fft_start).count();
    } else {
        x_test.resize(static_cast<size_t>(tree->num_points));
    }

    if (rank == 0 && verbose) {
        std::cout << "FFT kernel memory: "
                  << std::scientific << std::setprecision(6)
                  << (fft_memory_usage / 1024.0 / 1024.0)
                  << " MB" << std::endl;
        std::cout << "FFT matvec time: "
                  << std::defaultfloat << std::setprecision(6)
                  << fft_matvec_ms << " ms" << std::endl;
    }

    MPI_Bcast(x_test.data(),
              static_cast<int>(tree->num_points),
              mpi_datatype_for<DataType>(),
              0,
              MPI_COMM_WORLD);

    const std::vector<DataType> local_x =
        extract_local_leaf_vector_from_global(tree, x_test);

    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> mul_data(
        tree->num_levels);
    hierarchical_mul_parallel(tree, local_x, mul_data, verbose);

    std::vector<DataType> y_hat;
    std::vector<DataType> aggregated_input;
    gather_solution_to_root(tree, mul_data, y_hat, aggregated_input);

    if (rank != 0) {
        return 0.0;
    }

    if (y_hat.size() != y_true.size()) {
        throw std::runtime_error("verify_forward_error_fft: gathered result size mismatch");
    }

    RealType diff_norm_sq = 0.0;
    RealType y_true_norm_sq = 0.0;
    RealType max_abs_diff = 0.0;
    int64_t max_idx = 0;
    for (int64_t i = 0; i < tree->num_points; ++i) {
        const DataType diff =
            y_hat[static_cast<size_t>(i)] - y_true[static_cast<size_t>(i)];
        diff_norm_sq += squared_magnitude(diff);
        y_true_norm_sq += squared_magnitude(y_true[static_cast<size_t>(i)]);
        const RealType abs_diff = static_cast<RealType>(std::abs(diff));
        if (abs_diff > max_abs_diff) {
            max_abs_diff = abs_diff;
            max_idx = i;
        }
    }

    const RealType diff_norm = std::sqrt(diff_norm_sq);
    const RealType y_true_norm = std::sqrt(y_true_norm_sq);
    double relative_error = std::numeric_limits<double>::infinity();
    if (y_true_norm > static_cast<RealType>(0.0)) {
        relative_error = static_cast<double>(diff_norm / y_true_norm);
    }

    if (verbose) {
        std::cout << "\nForward Error Analysis:" << std::endl;
        std::cout << "  ||F*x - A*x||₂ = " << std::scientific
                  << std::setprecision(6) << diff_norm << std::endl;
        std::cout << "  ||A*x||₂       = " << y_true_norm << std::endl;
        std::cout << "  Relative forward error = " << relative_error << std::endl;
        std::cout << "  Max |F*x - A*x| = " << max_abs_diff
                  << " (at index " << max_idx << ")" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    return relative_error;
}

/**
 * @brief Estimate condition number via power iteration.
 *
 * Uses power iteration on A (via FFT matvec) to estimate ||A||,
 * and inverse power iteration via F^{-1} (hierarchical solve) to estimate ||A^{-1}||.
 */
template<typename CoordType, typename DataType>
double estimate_condition_number(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<CoordType>& grid_points,
    int64_t grid_size,
    KernelKind kernel_kind,
    int dimension,
    int num_samples,
    double wave_divisor,
    double length_scale = 0.1,
    double nugget = 1e-6,
    double kappa = 10.0) {

    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<
            std::is_same_v<DataType, std::complex<float>>, float,
            DataType
        >
    >;

    constexpr uint64_t cond_seed = 0xc04d1710e5714a7eULL;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Condition Number Estimation" << std::endl;
        std::cout << "  Power iteration samples: " << num_samples << std::endl;
        std::cout << "========================================" << std::endl;
    }

    const int64_t N = tree->num_points;

    std::vector<DataType> v_global;
    if (rank == 0) {
        v_global = build_global_verification_vector<DataType>(N, cond_seed);
        RealType norm_sq = 0;
        for (int64_t i = 0; i < N; ++i)
            norm_sq += static_cast<RealType>(std::norm(v_global[static_cast<size_t>(i)]));
        RealType norm = std::sqrt(norm_sq);
        for (int64_t i = 0; i < N; ++i)
            v_global[static_cast<size_t>(i)] /= static_cast<DataType>(norm);
    } else {
        v_global.resize(static_cast<size_t>(N));
    }

    double sigma_max = 0.0;
    for (int iter = 0; iter < num_samples; ++iter) {
        MPI_Bcast(v_global.data(), static_cast<int>(N),
                  mpi_datatype_for<DataType>(), 0, MPI_COMM_WORLD);

        std::vector<DataType> Av;
        if (rank == 0) {
            compute_fft_matvec(v_global, Av, grid_points, grid_size,
                               kernel_kind, dimension, wave_divisor,
                               length_scale, nugget, kappa);
            RealType Av_norm_sq = 0;
            for (int64_t i = 0; i < N; ++i)
                Av_norm_sq += static_cast<RealType>(std::norm(Av[static_cast<size_t>(i)]));
            sigma_max = std::sqrt(static_cast<double>(Av_norm_sq));

            for (int64_t i = 0; i < N; ++i)
                v_global[static_cast<size_t>(i)] = Av[static_cast<size_t>(i)] / static_cast<DataType>(sigma_max);

            std::cout << "  ||A|| iteration " << (iter + 1) << ": " << std::scientific
                      << std::setprecision(6) << sigma_max << std::endl;
        }
    }

    if (rank == 0) {
        v_global = build_global_verification_vector<DataType>(N, cond_seed + 1);
        RealType norm_sq = 0;
        for (int64_t i = 0; i < N; ++i)
            norm_sq += static_cast<RealType>(std::norm(v_global[static_cast<size_t>(i)]));
        RealType norm = std::sqrt(norm_sq);
        for (int64_t i = 0; i < N; ++i)
            v_global[static_cast<size_t>(i)] /= static_cast<DataType>(norm);
    }

    double sigma_min_inv = 0.0;
    for (int iter = 0; iter < num_samples; ++iter) {
        MPI_Bcast(v_global.data(), static_cast<int>(N),
                  mpi_datatype_for<DataType>(), 0, MPI_COMM_WORLD);

        const std::vector<DataType> local_v =
            extract_local_leaf_vector_from_global(tree, v_global);

        std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> solve_data(
            tree->num_levels);
        hierarchical_solve_parallel(tree, local_v, solve_data, false);

        std::vector<DataType> result_global;
        std::vector<DataType> agg_input;
        gather_solution_to_root(tree, solve_data, result_global, agg_input);

        if (rank == 0) {
            RealType result_norm_sq = 0;
            for (int64_t i = 0; i < N; ++i)
                result_norm_sq += static_cast<RealType>(std::norm(result_global[static_cast<size_t>(i)]));
            sigma_min_inv = std::sqrt(static_cast<double>(result_norm_sq));

            for (int64_t i = 0; i < N; ++i)
                v_global[static_cast<size_t>(i)] = result_global[static_cast<size_t>(i)] / static_cast<DataType>(sigma_min_inv);

            std::cout << "  ||A^{-1}|| iteration " << (iter + 1) << ": " << std::scientific
                      << std::setprecision(6) << sigma_min_inv << std::endl;
        }
    }

    double kappa_est = sigma_max * sigma_min_inv;

    if (rank == 0) {
        std::cout << "\n  ||A||      ≈ " << std::scientific << std::setprecision(6) << sigma_max << std::endl;
        std::cout << "  ||A^{-1}|| ≈ " << sigma_min_inv << std::endl;
        std::cout << "  κ(A)       ≈ " << kappa_est << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    return kappa_est;
}

} // namespace fmm


namespace {

enum class NumberKind {
    REAL,
    COMPLEX
};

struct ProgramOptions {
    int num_levels = 0;
    int64_t N = 0;
    int64_t grid_size = 0;
    double tolerance = 0.0;
    fmm::KernelKind kernel_kind = fmm::KernelKind::LAPLACE;
    NumberKind number_kind = NumberKind::REAL;
    int dimension = 3;
    int64_t reduction_threshold = 0;
    int num_proxy = -1;
    double wave_divisor = 32.0;
    double length_scale = 0.1;   // Matérn length scale ℓ
    double nugget = 1e-6;        // Matérn diagonal nugget σ_n²
    double kappa = 10.0;         // Yukawa screening parameter κ
    int cond_samples = 0;        // Power iteration samples for condition number estimate (0 = skip)
};

std::string to_lower_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

int parse_int_arg(const std::string& text, const std::string& name) {
    size_t idx = 0;
    int value = 0;
    try {
        value = std::stoi(text, &idx);
    } catch (const std::exception&) {
        throw std::invalid_argument("Invalid integer for " + name + ": '" + text + "'");
    }
    if (idx != text.size()) {
        throw std::invalid_argument("Invalid integer for " + name + ": '" + text + "'");
    }
    return value;
}

int64_t parse_int64_arg(const std::string& text, const std::string& name) {
    size_t idx = 0;
    int64_t value = 0;
    try {
        value = std::stoll(text, &idx);
    } catch (const std::exception&) {
        throw std::invalid_argument("Invalid integer for " + name + ": '" + text + "'");
    }
    if (idx != text.size()) {
        throw std::invalid_argument("Invalid integer for " + name + ": '" + text + "'");
    }
    return value;
}

double parse_double_arg(const std::string& text, const std::string& name) {
    size_t idx = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &idx);
    } catch (const std::exception&) {
        throw std::invalid_argument("Invalid floating-point value for " + name + ": '" + text + "'");
    }
    if (idx != text.size()) {
        throw std::invalid_argument("Invalid floating-point value for " + name + ": '" + text + "'");
    }
    return value;
}

int64_t default_reduction_threshold_for_dimension(int dimension) {
    return (dimension == 2) ? 256 : 4096;
}

int64_t min_reduction_threshold_for_dimension(int dimension) {
    return (dimension == 2) ? 64 : 512;
}

int default_num_proxy_for_dimension(int dimension) {
    return (dimension == 2) ? 32 : 256;
}

bool is_power_of_base(int64_t value, int64_t base) {
    if (value < 1 || base < 2) {
        return false;
    }
    while (value % base == 0) {
        value /= base;
    }
    return value == 1;
}

bool is_valid_reduction_threshold(int64_t reduction_threshold, int dimension) {
    const int64_t base = (dimension == 2) ? 4 : 8;
    return is_power_of_base(reduction_threshold, base);
}

std::string reduction_threshold_pattern(int dimension) {
    if (dimension == 2) {
        return "4^k (1, 4, 16, 64, 256, ...)";
    }
    return "8^k (1, 8, 64, 512, 4096, ...)";
}

fmm::KernelKind parse_kernel_kind(const std::string& text) {
    const std::string value = to_lower_copy(text);
    if (value == "laplace") {
        return fmm::KernelKind::LAPLACE;
    }
    if (value == "helmholtz") {
        return fmm::KernelKind::HELMHOLTZ;
    }
    if (value == "matern52" || value == "matern5/2" || value == "matern_52") {
        return fmm::KernelKind::MATERN52;
    }
    if (value == "yukawa" || value == "screened_coulomb" || value == "screened") {
        return fmm::KernelKind::YUKAWA;
    }
    throw std::invalid_argument("Unsupported kernel type: '" + text + "'");
}

NumberKind parse_number_kind(const std::string& text) {
    const std::string value = to_lower_copy(text);
    if (value == "real" || value == "double") {
        return NumberKind::REAL;
    }
    if (value == "complex" || value == "complex_double" || value == "complex128") {
        return NumberKind::COMPLEX;
    }
    throw std::invalid_argument("Unsupported number type: '" + text + "'");
}

const char* kernel_kind_to_string(fmm::KernelKind kind) {
    switch (kind) {
        case fmm::KernelKind::LAPLACE:
            return "laplace";
        case fmm::KernelKind::HELMHOLTZ:
            return "helmholtz";
        case fmm::KernelKind::MATERN52:
            return "matern52";
        case fmm::KernelKind::YUKAWA:
            return "yukawa";
        default:
            return "unknown";
    }
}

const char* number_kind_to_string(NumberKind kind) {
    switch (kind) {
        case NumberKind::REAL:
            return "real";
        case NumberKind::COMPLEX:
            return "complex";
        default:
            return "unknown";
    }
}

void print_usage(const char* executable_name) {
    std::cerr << "Usage:\n";
    std::cerr << "  " << executable_name
              << " <num_levels> <N> <grid_size> <tolerance>\n";
    std::cerr << "  " << executable_name
              << " <num_levels> <N> <grid_size> <tolerance> [<kernel> <number_type> <dimension>]\n";
    std::cerr << "  " << executable_name
              << " <num_levels> <N> <grid_size> <tolerance>"
              << " [--kernel <laplace|helmholtz|matern52|yukawa>]"
              << " [--number-type <real|complex>]"
              << " [--dimension <2|3>]"
              << " [--reduction-threshold <count>]"
              << " [--num-proxy <count>]"
              << " [--wave-divisor <value>]"
              << " [--length-scale <value>]"
              << " [--nugget <value>]"
              << " [--kappa <value>]"
              << " [--cond-samples <n>]\n\n";
    std::cerr << "Defaults preserve the previous behavior:\n";
    std::cerr << "  kernel=laplace, number_type=real, dimension=3\n";
    std::cerr << "  reduction_threshold=256 (2D) or 4096 (3D)\n";
    std::cerr << "  num_proxy=-1 uses default 32 (2D) or 256 (3D)\n";
    std::cerr << "  num_proxy=0 disables proxy rows\n";
    std::cerr << "  wave_divisor=32 (Helmholtz only)\n";
    std::cerr << "  length_scale=0.1 (Matern52 only)\n";
    std::cerr << "  nugget=1e-6 (Matern52 only)\n";
    std::cerr << "  kappa=10 (Yukawa only, screening parameter)\n";
    std::cerr << "  reduction_threshold must be >= 64 in 2D or >= 512 in 3D\n";
    std::cerr << "  reduction_threshold must be "
              << "4^k in 2D or 8^k in 3D\n";
}

ProgramOptions parse_program_options(int argc, char* argv[]) {
    if (argc < 5) {
        throw std::invalid_argument("Expected at least 4 positional arguments.");
    }

    ProgramOptions options;
    options.num_levels = parse_int_arg(argv[1], "num_levels");
    options.N = parse_int64_arg(argv[2], "N");
    options.grid_size = parse_int64_arg(argv[3], "grid_size");
    options.tolerance = parse_double_arg(argv[4], "tolerance");

    int extra_positionals = 0;
    for (int i = 5; i < argc; ++i) {
        const std::string token = argv[i];
        if (token == "--kernel" || token == "-k") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --kernel.");
            }
            options.kernel_kind = parse_kernel_kind(argv[++i]);
            continue;
        }
        if (token == "--number-type" || token == "--number" || token == "-n") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --number-type.");
            }
            options.number_kind = parse_number_kind(argv[++i]);
            continue;
        }
        if (token == "--dimension" || token == "--dim" || token == "-d") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --dimension.");
            }
            options.dimension = parse_int_arg(argv[++i], "dimension");
            continue;
        }
        if (token == "--reduction-threshold") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --reduction-threshold.");
            }
            options.reduction_threshold =
                parse_int64_arg(argv[++i], "reduction_threshold");
            continue;
        }
        if (token == "--num-proxy") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --num-proxy.");
            }
            options.num_proxy = parse_int_arg(argv[++i], "num_proxy");
            continue;
        }
        if (token == "--wave-divisor") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --wave-divisor.");
            }
            options.wave_divisor = parse_double_arg(argv[++i], "wave_divisor");
            continue;
        }
        if (token == "--length-scale" || token == "--ls") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --length-scale.");
            }
            options.length_scale = parse_double_arg(argv[++i], "length_scale");
            continue;
        }
        if (token == "--nugget") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --nugget.");
            }
            options.nugget = parse_double_arg(argv[++i], "nugget");
            continue;
        }
        if (token == "--kappa" || token == "--screening") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --kappa.");
            }
            options.kappa = parse_double_arg(argv[++i], "kappa");
            continue;
        }

        if (token == "--cond-samples") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value after --cond-samples.");
            }
            options.cond_samples = parse_int_arg(argv[++i], "cond-samples");
            continue;
        }

        if (extra_positionals == 0) {
            options.kernel_kind = parse_kernel_kind(token);
        } else if (extra_positionals == 1) {
            options.number_kind = parse_number_kind(token);
        } else if (extra_positionals == 2) {
            options.dimension = parse_int_arg(token, "dimension");
        } else {
            throw std::invalid_argument("Unexpected extra argument: '" + token + "'");
        }
        ++extra_positionals;
    }

    if (options.num_levels <= 0) {
        throw std::invalid_argument("num_levels must be positive.");
    }
    if (options.N <= 0) {
        throw std::invalid_argument("N must be positive.");
    }
    if (options.grid_size <= 0) {
        throw std::invalid_argument("grid_size must be positive.");
    }
    if (!(options.tolerance > 0.0)) {
        throw std::invalid_argument("tolerance must be positive.");
    }
    if (options.dimension != 2 && options.dimension != 3) {
        throw std::invalid_argument("dimension must be either 2 or 3.");
    }

    if (options.reduction_threshold == 0) {
        options.reduction_threshold =
            default_reduction_threshold_for_dimension(options.dimension);
    }
    if (options.num_proxy == -1) {
        options.num_proxy = default_num_proxy_for_dimension(options.dimension);
    }

    if (options.reduction_threshold <= 0) {
        throw std::invalid_argument("reduction_threshold must be positive.");
    }
    if (options.reduction_threshold < min_reduction_threshold_for_dimension(options.dimension)) {
        throw std::invalid_argument(
            "reduction_threshold must be at least " +
            std::to_string(min_reduction_threshold_for_dimension(options.dimension)) +
            " for dimension " + std::to_string(options.dimension) + ".");
    }
    if (!is_valid_reduction_threshold(options.reduction_threshold, options.dimension)) {
        throw std::invalid_argument(
            "reduction_threshold must be " +
            reduction_threshold_pattern(options.dimension) +
            " for dimension " + std::to_string(options.dimension) + ".");
    }
    if (options.num_proxy < -1) {
        throw std::invalid_argument(
            "num_proxy must be -1 (default), 0, or a positive integer.");
    }
    if (options.dimension == 3 && options.num_proxy == 1) {
        throw std::invalid_argument(
            "num_proxy=1 is invalid in 3D. Use 0 to disable proxies or >= 2.");
    }
    if (!(options.wave_divisor > 0.0)) {
        throw std::invalid_argument("wave_divisor must be positive.");
    }

    if (options.kernel_kind == fmm::KernelKind::LAPLACE &&
        options.number_kind != NumberKind::REAL) {
        throw std::invalid_argument(
            "Unsupported combination: Laplace currently supports only real number type.");
    }

    if (options.kernel_kind == fmm::KernelKind::HELMHOLTZ &&
        options.number_kind != NumberKind::COMPLEX) {
        throw std::invalid_argument(
            "Unsupported combination: Helmholtz currently supports only complex number type.");
    }

    if (options.kernel_kind == fmm::KernelKind::MATERN52) {
        if (options.number_kind != NumberKind::REAL) {
            throw std::invalid_argument(
                "Unsupported combination: Matern52 currently supports only real number type.");
        }
        if (!(options.length_scale > 0.0)) {
            throw std::invalid_argument("length_scale must be positive.");
        }
        if (!(options.nugget >= 0.0)) {
            throw std::invalid_argument("nugget must be non-negative.");
        }
    }

    if (options.kernel_kind == fmm::KernelKind::YUKAWA) {
        if (options.number_kind != NumberKind::REAL) {
            throw std::invalid_argument(
                "Unsupported combination: Yukawa currently supports only real number type.");
        }
        if (!(options.kappa > 0.0)) {
            throw std::invalid_argument("kappa must be positive.");
        }
    }

    int64_t expected_points = options.grid_size;
    for (int d = 1; d < options.dimension; ++d) {
        expected_points *= options.grid_size;
    }
    if (expected_points != options.N) {
        throw std::invalid_argument(
            "N must equal grid_size^dimension for FFT verification.");
    }

    return options;
}

template<typename CoordType, typename DataType, typename KernelType>
int run_with_configuration(const ProgramOptions& options, KernelType& kernel,
                           int rank, int size) {
    if (rank == 0) {
        std::cout << "=== Hierarchical Factorization Test ("
                  << options.dimension << "D "
                  << kernel_kind_to_string(options.kernel_kind) << " kernel, "
                  << number_kind_to_string(options.number_kind) << ") ===" << std::endl;
    }

    double bounds[6] = {0.0, 1.0, 0.0, 1.0, 0.0, (options.dimension == 2) ? 0.0 : 1.0};

    auto tree = std::unique_ptr<fmm::ParallelTree<CoordType, DataType>>(
        fmm::create_uniform_tree<CoordType, DataType>(
            nullptr,
            options.N,
            options.num_levels,
            bounds,
            options.dimension,
            MPI_COMM_WORLD,
            options.reduction_threshold));

    const int leaf_level = options.num_levels - 1;

    if (rank == 0) {
        std::cout << "\nTree Setup:" << std::endl;
        std::cout << "Levels: " << options.num_levels << std::endl;
        std::cout << "Leaf level: " << leaf_level << std::endl;
        std::cout << "Total points: " << options.N << std::endl;
        std::cout << "MPI processes: " << size << std::endl;
        std::cout << "Reduction threshold: " << options.reduction_threshold << std::endl;
        std::cout << "Num proxy points: " << options.num_proxy << std::endl;
        if (options.kernel_kind == fmm::KernelKind::HELMHOLTZ) {
            std::cout << "Wave divisor: " << options.wave_divisor << std::endl;
        }
        if (options.kernel_kind == fmm::KernelKind::MATERN52) {
            std::cout << "Length scale: " << options.length_scale << std::endl;
            std::cout << "Nugget: " << options.nugget << std::endl;
        }
        if (options.kernel_kind == fmm::KernelKind::YUKAWA) {
            std::cout << "Kappa (screening): " << options.kappa << std::endl;
        }
        if (const int dynamic_cpu_cap =
                fmm::parse_positive_thread_count(std::getenv("FMM_MAX_CPUS_PER_NODE"));
            dynamic_cpu_cap > 0) {
            std::cout << "Dynamic thread cpu cap per node: " << dynamic_cpu_cap << std::endl;
        }
    }

    const auto factorization_method =
        (options.number_kind == NumberKind::COMPLEX)
            ? fmm::FactorizationMethod::COMPLEX_SYM
            : fmm::FactorizationMethod::LU;
    fmm::HierarchicalFactorization<CoordType, DataType, KernelType> factorizer(
        options.N,
        fmm::MatrixProperty::SYMMETRIC,
        &kernel,
        options.dimension,
        factorization_method,
        options.num_proxy);

    const auto& unit_proxy = factorizer.get_unit_proxy_points();
    const int num_proxy = factorizer.get_num_proxy_points();

    const bool is_symmetric = true;
    const bool is_hermitian = false;

    MPI_Barrier(MPI_COMM_WORLD);

    auto total_start = std::chrono::high_resolution_clock::now();

    fmm::hierarchical_factorization_parallel(
        tree.get(),
        &kernel,
        options.tolerance,
        is_symmetric,
        is_hermitian,
        factorization_method,
        unit_proxy,
        num_proxy,
        2.5,
        true);

    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        total_end - total_start);

    if (rank == 0) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Total factorization time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\n=== Factorization Complete ===" << std::endl;
    }

    constexpr uint64_t backward_seed = 0x123456789abcdef0ULL;
    std::vector<DataType> rhs =
        fmm::build_local_verification_vector(tree.get(), backward_seed);

    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>> solve_data(
        options.num_levels);
    fmm::hierarchical_solve_parallel(tree.get(), rhs, solve_data, true);

    std::vector<DataType> solution;
    std::vector<DataType> aggregated_rhs;
    const auto gather_verify_start = std::chrono::high_resolution_clock::now();
    fmm::gather_solution_to_root(tree.get(), solve_data, solution, aggregated_rhs);
    const auto gather_verify_end = std::chrono::high_resolution_clock::now();
    const double gather_verify_ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            gather_verify_end - gather_verify_start).count();
    double gather_verify_max_ms = 0.0;
    MPI_Reduce(&gather_verify_ms, &gather_verify_max_ms, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Gather communication time for solve verification: "
                  << gather_verify_max_ms << " ms" << std::endl;
    }

    auto grid_points = fmm::build_grid_points_from_tree(tree.get(), options.grid_size);
    int can_run_fft_verification = 0;
    if (rank == 0 && !grid_points.empty()) {
        can_run_fft_verification = 1;
    }
    MPI_Bcast(&can_run_fft_verification, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (can_run_fft_verification) {
        const auto verify_fft_start = std::chrono::high_resolution_clock::now();
        const auto solve_metrics = fmm::verify_solution_fft(
            tree.get(),
            aggregated_rhs,
            solution,
            grid_points,
            options.grid_size,
            options.kernel_kind,
            options.dimension,
            options.wave_divisor,
            options.length_scale,
            options.nugget,
            options.kappa,
            true);
        const auto verify_fft_end = std::chrono::high_resolution_clock::now();
        const auto verify_fft_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            verify_fft_end - verify_fft_start);

        if (rank == 0) {
            std::cout << "FFT verification time: " << verify_fft_duration.count() << " ms" << std::endl;
            fflush(stdout);
        }

        // Forward error via multiply: ||F*x - A*x|| / ||A*x|| (no condition number amplification)
        const double forward_relative_error = fmm::verify_forward_error_fft(
            tree.get(),
            grid_points,
            options.grid_size,
            options.kernel_kind,
            options.dimension,
            options.wave_divisor,
            options.length_scale,
            options.nugget,
            options.kappa,
            true);
        if (rank == 0) {
            const std::streamsize fwd_old_precision = std::cout.precision();
            const std::ios::fmtflags fwd_old_flags = std::cout.flags();
            std::cout << "Forward relative error: " << std::scientific
                      << std::setprecision(std::numeric_limits<double>::max_digits10)
                      << forward_relative_error
                      << std::setprecision(fwd_old_precision);
            std::cout.flags(fwd_old_flags);
            std::cout << std::endl;
            fflush(stdout);
        }

        if (options.N <= 4096) {
            const double direct_error = fmm::verify_solution_direct(
                &kernel,
                aggregated_rhs,
                solution,
                grid_points,
                tree->num_points,
                options.dimension,
                true);
            if (rank == 0) {
                std::cout << "verifying with direct matrix vector multiply since N <= 4096" << std::endl;
                std::cout << "Direct error: " << direct_error << std::endl;
            }
        }

        if (rank == 0) {
            const std::streamsize old_precision = std::cout.precision();
            const std::ios::fmtflags old_flags = std::cout.flags();
            std::cout << "\nBackward residual: " << std::scientific
                      << std::setprecision(std::numeric_limits<double>::max_digits10)
                      << solve_metrics.relative_residual
                      << std::setprecision(old_precision);
            std::cout.flags(old_flags);
            std::cout << std::endl;
        }

        if (options.cond_samples > 0) {
            const double kappa_est = fmm::estimate_condition_number(
                tree.get(),
                grid_points,
                options.grid_size,
                options.kernel_kind,
                options.dimension,
                options.cond_samples,
                options.wave_divisor,
                options.length_scale,
                options.nugget,
                options.kappa);
            if (rank == 0) {
                const std::streamsize be_old_precision = std::cout.precision();
                const std::ios::fmtflags be_old_flags = std::cout.flags();
                std::cout << "Backward relative error (rough estimate, uses ||A|| approx ||Ax||/||x||): "
                          << std::scientific
                          << std::setprecision(std::numeric_limits<double>::max_digits10)
                          << solve_metrics.backward_error
                          << std::setprecision(be_old_precision);
                std::cout.flags(be_old_flags);
                std::cout << std::endl;
                std::cout << "Estimated condition number: " << std::scientific
                          << std::setprecision(6) << kappa_est << std::endl;
            }
        }
    }

    return 0;
}

int dispatch_and_run(const ProgramOptions& options, int rank, int size) {
    using CoordType = double;
    if (options.kernel_kind == fmm::KernelKind::LAPLACE) {
        if (options.dimension == 2) {
            kernel::LaplaceKernel2D<double> kernel(options.N);
            return run_with_configuration<CoordType, double>(options, kernel, rank, size);
        }
        kernel::LaplaceKernel3D<double> kernel(options.N);
        return run_with_configuration<CoordType, double>(options, kernel, rank, size);
    }

    if (options.kernel_kind == fmm::KernelKind::HELMHOLTZ) {
        if (options.dimension == 2) {
            kernel::HelmholtzKernel2D<double, double> kernel(options.N, options.wave_divisor);
            return run_with_configuration<CoordType, std::complex<double>>(options, kernel, rank, size);
        }
        kernel::HelmholtzKernel3D<double, double> kernel(options.N, options.wave_divisor);
        return run_with_configuration<CoordType, std::complex<double>>(options, kernel, rank, size);
    }

    if (options.kernel_kind == fmm::KernelKind::MATERN52) {
        if (options.dimension == 2) {
            kernel::Matern52Kernel2D<double> kernel(options.N, options.length_scale, options.nugget);
            return run_with_configuration<CoordType, double>(options, kernel, rank, size);
        }
        kernel::Matern52Kernel3D<double> kernel(options.N, options.length_scale, options.nugget);
        return run_with_configuration<CoordType, double>(options, kernel, rank, size);
    }

    if (options.kernel_kind == fmm::KernelKind::YUKAWA) {
        if (options.dimension == 2) {
            kernel::YukawaKernel2D<double> kernel(options.N, options.kappa);
            return run_with_configuration<CoordType, double>(options, kernel, rank, size);
        }
        kernel::YukawaKernel3D<double> kernel(options.N, options.kappa);
        return run_with_configuration<CoordType, double>(options, kernel, rank, size);
    }

    throw std::invalid_argument("Unsupported kernel selection.");
}

}  // namespace


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        std::cout << "Thread runtime: fmm_threads=" << fmm::configured_fmm_thread_count()
                  << ", omp_max_threads=" << omp_get_max_threads()
                  << ", openblas_threads=" << openblas_get_num_threads()
                  << ", visible_cpus=" << fmm::visible_process_cpu_count();
        if (const int dynamic_cpu_cap =
                fmm::parse_positive_thread_count(std::getenv("FMM_MAX_CPUS_PER_NODE"));
            dynamic_cpu_cap > 0) {
            std::cout << ", dynamic_cpu_cap_per_node=" << dynamic_cpu_cap;
        }
        std::cout << std::endl;
    }

    {
        const auto& base_cpus = fmm::base_process_cpu_list();
        const int cap = fmm::parse_positive_thread_count(std::getenv("FMM_MAX_CPUS_PER_NODE"));
        MPI_Comm shared_comm;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                            MPI_INFO_NULL, &shared_comm);
        int shared_rank = 0, shared_size = 1;
        MPI_Comm_rank(shared_comm, &shared_rank);
        MPI_Comm_size(shared_comm, &shared_size);
        MPI_Comm_free(&shared_comm);

        int slice_first = -1, slice_last = -1, slice_count = 0;
        if (!base_cpus.empty()) {
            const int usable = std::min<int>(cap > 0 ? cap : (int)base_cpus.size(),
                                              (int)base_cpus.size());
            const int threads = std::max(1, usable / std::max(1, shared_size));
            int slice_begin_idx = shared_rank * threads;
            if (slice_begin_idx >= usable) slice_begin_idx = std::max(0, usable - 1);
            const int slice_end_idx = std::min(slice_begin_idx + threads, usable);
            if (slice_end_idx > slice_begin_idx) {
                slice_first = base_cpus[slice_begin_idx];
                slice_last = base_cpus[slice_end_idx - 1];
                slice_count = slice_end_idx - slice_begin_idx;
            }
        }

        constexpr int kSocketCores = 64;
        constexpr int kCcdCores = 8;
        const bool socket_ok = (slice_first >= 0) &&
            ((slice_count >= kSocketCores && slice_first % kSocketCores == 0 && slice_count % kSocketCores == 0) ||
             (slice_count < kSocketCores && (slice_first / kSocketCores == slice_last / kSocketCores)));
        const bool ccd_ok = (slice_first >= 0) &&
            ((slice_count >= kCcdCores && slice_first % kCcdCores == 0 && slice_count % kCcdCores == 0) ||
             (slice_count < kCcdCores && (slice_first / kCcdCores == slice_last / kCcdCores)));

        struct AffinityInfo { int rank; int shared_rank; int first; int last; int count; int socket_ok; int ccd_ok; };
        AffinityInfo local{rank, shared_rank, slice_first, slice_last, slice_count,
                            socket_ok ? 1 : 0, ccd_ok ? 1 : 0};
        std::vector<AffinityInfo> all(size);
        MPI_Gather(&local, sizeof(AffinityInfo), MPI_BYTE,
                   all.data(), sizeof(AffinityInfo), MPI_BYTE,
                   0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "Projected per-rank CPU slice (Perlmutter: 64 cores/socket, 8 cores/CCD):" << std::endl;
            bool any_socket_bad = false, any_ccd_bad = false;
            for (const auto& a : all) {
                std::cout << "  rank " << a.rank
                          << " (node-local " << a.shared_rank << ")"
                          << ": cpus=[" << a.first << "," << a.last << "]"
                          << " count=" << a.count
                          << (a.socket_ok ? " socket-ok" : " SOCKET-CROSSES")
                          << (a.ccd_ok ? " ccd-ok" : " ccd-crosses")
                          << std::endl;
                if (!a.socket_ok) any_socket_bad = true;
                if (!a.ccd_ok) any_ccd_bad = true;
            }
            if (any_socket_bad) {
                std::cout << "  WARNING: at least one rank's slice crosses a socket boundary"
                          << " — NUMA traffic expected." << std::endl;
            } else if (any_ccd_bad) {
                std::cout << "  NOTE: slices cross CCD boundaries (cross-L3 within socket is OK but not optimal)."
                          << std::endl;
            }
        }
    }

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            if (rank == 0) {
                print_usage(argv[0]);
            }
            MPI_Finalize();
            return 0;
        }
    }

    ProgramOptions options;
    try {
        options = parse_program_options(argc, argv);
    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Argument error: " << e.what() << std::endl;
            print_usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    try {
        if (rank == 0) {
            std::cout << "Run configuration: kernel=" << kernel_kind_to_string(options.kernel_kind)
                      << ", number_type=" << number_kind_to_string(options.number_kind)
                      << ", dimension=" << options.dimension
                      << ", reduction_threshold=" << options.reduction_threshold
                      << ", num_proxy=" << options.num_proxy;
            if (options.kernel_kind == fmm::KernelKind::MATERN52) {
                std::cout << ", length_scale=" << options.length_scale
                          << ", nugget=" << options.nugget;
            }
            if (options.kernel_kind == fmm::KernelKind::YUKAWA) {
                std::cout << ", kappa=" << options.kappa;
            }
            std::cout << std::endl;
        }
        dispatch_and_run(options, rank, size);
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();
    return 0;
}
