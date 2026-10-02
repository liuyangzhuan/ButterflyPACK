#pragma once
// Environment variables of ButterflyPACK's H2 code and of its GPU backends
// (H2 and HODLR).  doc/environment_variables.md lists all of them and
// explains every value of BPACK_CHECK and BPACK_TRACE; keep it in step with
// the value tables below.
//
//   BPACK_GPU_HEAP_FRACTION, BPACK_GPU_EXCHANGE_MB   h2_gpu/device_heap.hpp
//   BPACK_GPU_AWARE_MPI                              h2_gpu/gpu_runtime.hpp
//   BPACK_MAX_CPUS_PER_NODE                          core/runtime_thread_support.hpp
//   BPACK_CHECK, BPACK_TRACE                         here: comma-separated
//                                                    lists of values, read once
//
// report_environment() names, once per process, the variables earlier
// versions read (they are ignored now) and the values of the two lists that
// are not known.

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace fmm {
namespace env {

namespace detail {

struct Value {
    const char* name;
    bool argument;  // given as name:argument
};

// The values of BPACK_CHECK and BPACK_TRACE (doc/environment_variables.md).
inline constexpr Value kCheckValues[] = {
    {"solve", false},  {"matvec", false}, {"kernel", false},          {"replica", false},
    {"parity", false}, {"hodlr", false},  {"hodlr-transpose", false}, {"hodlr-exact", false},
};
inline constexpr Value kTraceValues[] = {
    {"wave", false}, {"exchange", false}, {"phase", false},    {"bk", false},
    {"memory", false}, {"id", true},      {"hodlr-qr", false}, {"sync", false},
};

struct List {
    std::vector<std::pair<std::string, std::string>> values;  // name, argument
    std::vector<std::string> unknown;
};

inline std::string trim(std::string s) {
    const size_t a = s.find_first_not_of(" \t");
    if (a == std::string::npos) return {};
    return s.substr(a, s.find_last_not_of(" \t") - a + 1);
}

// (no MPI here: the first query may come from any thread)
template<size_t N>
List parse(const char* variable, const Value (&known)[N]) {
    List list;
    const char* v = std::getenv(variable);
    if (v == nullptr) return list;
    const std::string s(v);
    for (size_t start = 0; start <= s.size();) {
        size_t end = s.find(',', start);
        if (end == std::string::npos) end = s.size();
        const std::string item = trim(s.substr(start, end - start));
        start = end + 1;
        if (item.empty()) continue;
        const size_t colon = item.find(':');
        const std::string name = trim(item.substr(0, colon));
        const std::string argument = colon == std::string::npos ? std::string() : trim(item.substr(colon + 1));
        bool valid = false;
        for (const Value& k : known) valid = valid || (name == k.name && k.argument == !argument.empty());
        if (valid) {
            list.values.emplace_back(name, argument);
        } else {
            list.unknown.push_back(item);
        }
    }
    return list;
}

inline const List& check_list() {
    static const List list = parse("BPACK_CHECK", kCheckValues);
    return list;
}
inline const List& trace_list() {
    static const List list = parse("BPACK_TRACE", kTraceValues);
    return list;
}

inline bool listed(const List& list, std::string_view name) {
    for (const auto& [n, a] : list.values) {
        if (n == name) return true;
    }
    return false;
}

}  // namespace detail

// Whether BPACK_CHECK lists `value` (e.g. check("solve")).
inline bool check(std::string_view value) { return detail::listed(detail::check_list(), value); }

// Whether BPACK_TRACE lists `value`; trace_argument: the argument of a value
// given as value:argument (id:<file prefix>), empty if not listed.
inline bool trace(std::string_view value) { return detail::listed(detail::trace_list(), value); }
inline std::string trace_argument(std::string_view value) {
    for (const auto& [n, a] : detail::trace_list().values) {
        if (n == value) return a;
    }
    return {};
}

// Once per process, collective over comm (its rank 0 prints): the checks and
// traces that are on, the unknown values of the two lists, and the
// variables of earlier versions that are set (no longer read).
inline void report_environment(MPI_Comm comm) {
    static bool done = false;
    if (done) return;
    done = true;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank != 0) return;
    struct Old {
        const char* name;
        const char* instead;
    };
    static constexpr Old kOld[] = {
        {"H2_GPU_HEAP_FRACTION", "use BPACK_GPU_HEAP_FRACTION"},
        {"H2_GPU_HEAP_GB", "ranks that share a GPU now split its memory by themselves"},
        {"H2_GPU_EXCHANGE_MB", "use BPACK_GPU_EXCHANGE_MB"},
        {"HODLR_GPU_EXCHANGE_MB", "use BPACK_GPU_EXCHANGE_MB"},
        {"H2_GPU_AWARE_MPI", "use BPACK_GPU_AWARE_MPI"},
        {"H2_GPU_DEVICE_EXCHANGE", "use BPACK_GPU_AWARE_MPI=0"},
        {"FMM_MAX_CPUS_PER_NODE", "use BPACK_MAX_CPUS_PER_NODE"},
        {"H2_GPU_SOLVE_CHECK", "use BPACK_CHECK=solve"},
        {"H2_GPU_MATVEC_CHECK", "use BPACK_CHECK=matvec"},
        {"H2_GPU_KERNEL_CHECK", "use BPACK_CHECK=kernel"},
        {"H2_CA_REPLICA_CHECK", "use BPACK_CHECK=replica"},
        {"H2_OWNER_PARITY_WAVES", "use BPACK_CHECK=parity"},
        {"HODLR_GPU_CHECK", "use BPACK_CHECK=hodlr (2: hodlr-transpose)"},
        {"HODLR_GPU_SPLIT", "use BPACK_CHECK=hodlr-exact for HODLR_GPU_SPLIT=0"},
        {"H2_GPU_WAVE_TRACE", "use BPACK_TRACE=wave"},
        {"H2_GPU_EXCHANGE_TRACE", "use BPACK_TRACE=exchange"},
        {"H2_PHASE_REPORT", "use BPACK_TRACE=phase"},
        {"H2_BK_DIAGNOSTICS", "use BPACK_TRACE=bk"},
        {"FMM_MEMORY_DIAGNOSTICS", "use BPACK_TRACE=memory"},
        {"H2_ID_TRACE", "use BPACK_TRACE=id:<file prefix>"},
        {"HODLR_GPU_DEBUG", "use BPACK_TRACE=hodlr-qr (2: hodlr-qr,sync)"},
        {"H2_GPU_SOLVE", "always on"},
        {"H2_GPU_SKETCH", "always on"},
        {"H2_GPU_SOLVE_DIRECT", "always on with GPU-aware MPI"},
        {"H2_GPU_CA_SOLVE_HALO", "always on"},
        {"H2_GPU_CA_DEVICE_HALO", "always on"},
        {"H2_GPU_MATVEC", "always on"},
        {"H2_GPU_MATVEC_DIRECT", "always on with GPU-aware MPI"},
        {"H2_GPU_MATVEC_SYMMETRIC", "always on"},
        {"H2_GPU_SOLVE_KEEP", "always on"},
        {"H2_GPU_SOLVE_KEEP_FRACTION", "fixed at 0.8 of the device pool"},
        {"H2_GPU_MATVEC_KEEP_FRACTION", "fixed at 0.9 of the device pool"},
        {"H2_GPU_KEEP_OPERATORS", "always on"},
        {"H2_GPU_COPY_THREADS", "fixed"},
        {"H2_GPU_WARMUP", "always on, not in the factor time"},
        {"HODLR_GPU_KEEP_FACTORS", "always on"},
        {"HODLR_GPU_DEFER_HOST", "always on (off with BPACK_CHECK=hodlr or ErrFillFull=1)"},
    };
    auto print_list = [](const char* variable, const detail::List& list) {
        if (!list.values.empty()) {
            std::string names;
            for (const auto& [n, a] : list.values) names += (names.empty() ? "" : ", ") + n + (a.empty() ? "" : ":" + a);
            std::printf("BPACK: %s: %s\n", variable, names.c_str());
        }
        for (const std::string& u : list.unknown) {
            std::printf("BPACK: %s: unknown value '%s' ignored (doc/environment_variables.md)\n", variable, u.c_str());
        }
    };
    print_list("BPACK_CHECK", detail::check_list());
    print_list("BPACK_TRACE", detail::trace_list());
    for (const Old& o : kOld) {
        if (std::getenv(o.name) != nullptr) {
            std::printf("BPACK: %s is no longer read: %s (doc/environment_variables.md)\n", o.name, o.instead);
        }
    }
    std::fflush(stdout);
}

}  // namespace env
}  // namespace fmm
