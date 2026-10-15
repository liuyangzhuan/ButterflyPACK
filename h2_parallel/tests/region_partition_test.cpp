// Self-test of region_partition.hpp (the regions of any number of ranks on the
// unstructured Color path; GPU_BACKEND/h2_gpu/GPU_CA_PLAN.md, "Any number of
// ranks"): for 2D and 3D trees of every depth and 1 to 300 ranks, the layout
// of plan_regions against the invariants the tree code relies on -
//   - every level covered by its regions, none empty;
//   - where the count does not change, a region's range one level finer is
//     its range times 2^d (a box's children on its rank);
//   - where it does, each region's child regions tile its children, and
//     (counts other than (2^d)^k) every child region is whole parent boxes;
//   - (2^d)^k ranks: the caller's threshold, the counts of the old reduction
//     rule and the equal ranges (the bricks), so these runs are unchanged;
// then the worst leaf imbalance, and the layouts of a few counts.
//
//   make h2_region_partition_test && ./h2_region_partition_test     (in a build directory), or
//   g++ -std=c++17 -O2 -I../core region_partition_test.cpp -o region_partition_test && ./region_partition_test
//
// (no MPI calls; well under a second)

#include "region_partition.hpp"

#include <cstdio>
#include <exception>
#include <tuple>
#include <vector>

using namespace fmm;

static int failures = 0;

static void fail(const char* what, int dim, int levels, int procs, int level) {
    ++failures;
    std::printf("FAIL %s: dim %d levels %d procs %d level %d\n", what, dim, levels, procs, level);
}

// The counts of the reduction rule before any count was allowed ((2^d)^k only).
static std::vector<int32_t> cubic_counts(int dim, int levels, int procs, int64_t threshold) {
    std::vector<int32_t> count(static_cast<size_t>(levels), 1);
    count.back() = procs;
    const int32_t f = int32_t{1} << dim;
    for (int level = levels - 2; level >= 0; --level) {
        const int64_t boxes = int64_t{1} << (dim * level);
        const int32_t prev = count[static_cast<size_t>(level + 1)];
        count[static_cast<size_t>(level)] = boxes / prev < threshold ? std::max(1, prev / f) : prev;
    }
    return count;
}

int main() {
    const int64_t threshold = 8;
    long checked = 0, raised = 0, over = 0;
    double worst = 0.0;
    int worst_dim = 0, worst_levels = 0, worst_procs = 0;
    for (int dim = 2; dim <= 3; ++dim) {
        const int64_t f = int64_t{1} << dim;
        for (int levels = 2; levels <= (dim == 3 ? 8 : 11); ++levels) {
            const int64_t leaf_boxes = int64_t{1} << (dim * (levels - 1));
            for (int procs = 1; procs <= 300; ++procs) {
                if (leaf_boxes < procs) continue;
                RegionPlan p;
                try {
                    p = plan_regions(dim, levels, procs, threshold);
                } catch (const std::exception&) {
                    continue;  // no layout (counted by the caller as an error)
                }
                if (p.count[1] != 1) continue;  // the tree rejects these (level 1 has one rank)
                ++checked;
                const auto& c = p.count;
                const auto& s = p.start;
                const bool cubic = region_count_is_cubic(dim, procs);
                for (int l = 0; l < levels; ++l) {
                    const auto& t = s[static_cast<size_t>(l)];
                    const int32_t n = c[static_cast<size_t>(l)];
                    if (static_cast<int32_t>(t.size()) != n + 1 || t[0] != 0 || t.back() != (int64_t{1} << (dim * l))) {
                        fail("cover", dim, levels, procs, l);
                    }
                    for (int32_t i = 0; i < n; ++i) {
                        if (t[static_cast<size_t>(i) + 1] <= t[static_cast<size_t>(i)]) {
                            fail("empty region", dim, levels, procs, l);
                            break;
                        }
                    }
                    if (l + 1 == levels) continue;
                    const auto& u = s[static_cast<size_t>(l) + 1];
                    const int32_t below = c[static_cast<size_t>(l) + 1];
                    for (int32_t i = 0; i <= n; ++i) {
                        const int32_t child = n == below ? i : region_first_child(i, n, below);
                        if (u[static_cast<size_t>(child)] != t[static_cast<size_t>(i)] * f) {
                            fail(n == below ? "aligned children" : "grouped children", dim, levels, procs, l);
                            break;
                        }
                    }
                    if (n != below && !cubic) {
                        for (int32_t ch = 0; ch <= below; ++ch) {
                            if (u[static_cast<size_t>(ch)] % f != 0) {
                                fail("whole parent boxes", dim, levels, procs, l);
                                break;
                            }
                        }
                    }
                }
                if (cubic) {
                    if (cubic_counts(dim, levels, procs, threshold) != c || p.threshold != threshold) {
                        fail("(2^d)^k counts", dim, levels, procs, -1);
                    }
                    for (int l = 0; l < levels; ++l) {
                        const int64_t boxes = int64_t{1} << (dim * l);
                        const int32_t n = c[static_cast<size_t>(l)];
                        for (int32_t i = 0; i <= n; ++i) {
                            if (s[static_cast<size_t>(l)][static_cast<size_t>(i)] != i * (boxes / n)) {
                                fail("(2^d)^k ranges", dim, levels, procs, l);
                                break;
                            }
                        }
                    }
                } else {
                    if (p.threshold != threshold) ++raised;
                    if (p.leaf_imbalance > kRegionImbalance) ++over;
                    // the worst among layouts with 64 leaf boxes a rank or more (real runs)
                    if (leaf_boxes >= 64 * procs && p.leaf_imbalance > worst) {
                        worst = p.leaf_imbalance;
                        worst_dim = dim;
                        worst_levels = levels;
                        worst_procs = procs;
                    }
                }
            }
        }
    }
    std::printf("checked %ld layouts; other than (2^d)^k: threshold raised in %ld, leaf imbalance over %.2f in %ld\n",
                checked, raised, kRegionImbalance, over);
    std::printf("worst leaf imbalance with 64 leaf boxes a rank or more: %.3f (%dD, %d levels, %d ranks)\n", worst,
                worst_dim, worst_levels, worst_procs);
    std::printf("dim levels ranks: threshold, leaf imbalance, ranks by level\n");
    for (const auto& [dim, levels, procs] : std::vector<std::tuple<int, int, int>>{
             {3, 6, 8}, {3, 6, 12}, {3, 6, 20}, {3, 6, 24}, {3, 6, 50}, {3, 7, 64}, {3, 7, 100}, {3, 8, 192},
             {2, 8, 50}, {2, 10, 100}}) {
        const RegionPlan p = plan_regions(dim, levels, procs, threshold);
        std::printf("%dD %d %4d: %3ld %.3f  ", dim, levels, procs, static_cast<long>(p.threshold), p.leaf_imbalance);
        for (const int32_t x : p.count) std::printf("%d ", x);
        std::printf("\n");
    }
    if (failures) {
        std::printf("%d failures\n", failures);
        return 1;
    }
    std::printf("all checks passed\n");
    return 0;
}
