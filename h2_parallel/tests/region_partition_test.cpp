// Self-test of region_partition.hpp (the regions of any number of ranks on the
// unstructured Color path; GPU_BACKEND/h2_gpu/GPU_CA_PLAN.md, "Any number of
// ranks"): for 2D and 3D trees of every depth and 1 to 300 ranks, the layout
// of plan_regions against the invariants the tree code relies on -
//   - every level covered by its regions, none empty;
//   - where the count does not change, a region's range one level finer is
//     its range times 2^d (a box's children on its rank);
//   - where it does, each region's child regions tile its children, and
//     (unless the layout is the bricks) every child region is whole parent
//     boxes;
//   - (2^d)^k ranks without weights: the caller's threshold, the counts of
//     the old reduction rule and the equal ranges (the bricks), so these runs
//     are unchanged; with the same weight in every box, the same layouts as
//     without weights, for every count;
// and the same invariants with the boxes weighted by points (the unstructured
// path) - a thin shell (a sphere's surface: most boxes empty, the rest
// clipped by different amounts), a dense cluster and random weights with
// empty boxes - with the leaf regions' imbalance by points and the bricks
// flag; and with the counts keeping more ranks on the coarse levels (the
// unstructured path): the same invariants, level 1 on one rank, (2^d)^k in
// 3D at threshold 8 still the bricks, and how many coarse ranks are kept;
// the reduction groups (region_parent against region_first_child) for every
// pair of counts; then the worst imbalances, and the layouts of a few
// counts.
//
//   make h2_region_partition_test && ./h2_region_partition_test     (in a build directory), or
//   g++ -std=c++17 -O2 -I../core region_partition_test.cpp -o region_partition_test && ./region_partition_test
//
// (no MPI calls; a few seconds)

#include "morton.hpp"
#include "region_partition.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <random>
#include <tuple>
#include <vector>

using namespace fmm;

static int failures = 0;

static void fail(const char* what, const char* weights, int dim, int levels, int procs, int level) {
    ++failures;
    std::printf("FAIL %s (%s): dim %d levels %d procs %d level %d\n", what, weights, dim, levels, procs, level);
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

// The invariants of a layout (see the top).
static void check_layout(const RegionPlan& p, const char* weights, int dim, int levels, int procs) {
    const int64_t f = int64_t{1} << dim;
    const auto& c = p.count;
    const auto& s = p.start;
    bool bricks = region_count_is_cubic(dim, procs);
    for (int l = 0; l < levels; ++l) {
        const auto& t = s[static_cast<size_t>(l)];
        const int32_t n = c[static_cast<size_t>(l)];
        if (static_cast<int32_t>(t.size()) != n + 1 || t[0] != 0 || t.back() != (int64_t{1} << (dim * l))) {
            fail("cover", weights, dim, levels, procs, l);
            return;
        }
        for (int32_t i = 0; i < n; ++i) {
            if (t[static_cast<size_t>(i) + 1] <= t[static_cast<size_t>(i)]) {
                fail("empty range", weights, dim, levels, procs, l);
                break;
            }
        }
        bricks = bricks && region_starts_are_bricks(dim, l, t);
        if (l + 1 == levels) continue;
        const auto& u = s[static_cast<size_t>(l) + 1];
        const int32_t below = c[static_cast<size_t>(l) + 1];
        for (int32_t i = 0; i <= n; ++i) {
            const int32_t child = n == below ? i : region_first_child(i, n, below);
            if (u[static_cast<size_t>(child)] != t[static_cast<size_t>(i)] * f) {
                fail(n == below ? "aligned children" : "grouped children", weights, dim, levels, procs, l);
                break;
            }
        }
    }
    if (bricks != p.bricks) fail("bricks flag", weights, dim, levels, procs, -1);
    if (!p.bricks) {
        for (int l = 0; l + 1 < levels; ++l) {
            if (c[static_cast<size_t>(l)] == c[static_cast<size_t>(l) + 1]) continue;
            for (const int64_t x : s[static_cast<size_t>(l) + 1]) {
                if (x % f != 0) {
                    fail("whole parent boxes", weights, dim, levels, procs, l);
                    break;
                }
            }
        }
    }
}

// Points per box of `level` for three densities on the unit cube.
static std::vector<int64_t> weights_of(int kind, int dim, int level, std::mt19937_64& rng) {
    const int64_t boxes = int64_t{1} << (dim * level);
    const double width = 1.0 / static_cast<double>(1u << level);
    std::vector<int64_t> w(static_cast<size_t>(boxes), 0);
    std::uniform_int_distribution<int> pick(0, 99);
    for (int64_t m = 0; m < boxes; ++m) {
        uint32_t x = 0, y = 0, z = 0;
        morton::decode_nd(dim, static_cast<uint64_t>(m), x, y, z);
        const double cx = (x + 0.5) * width - 0.5, cy = (y + 0.5) * width - 0.5, cz = dim == 3 ? (z + 0.5) * width - 0.5 : 0.0;
        const double r = std::sqrt(cx * cx + cy * cy + cz * cz);
        double v = 0.0;
        if (kind == 0) {
            // a shell of radius 0.4: the boxes it crosses, more points the
            // nearer their centre is to it
            v = 1000.0 * std::max(0.0, 1.0 - std::abs(r - 0.4) / width);
        } else if (kind == 1) {
            // a cluster about (-0.2, -0.2, -0.2) over a sparse background
            const double dx = cx + 0.2, dy = cy + 0.2, dz = dim == 3 ? cz + 0.2 : 0.0;
            v = 1.0e4 * std::exp(-(dx * dx + dy * dy + dz * dz) / (2 * 0.1 * 0.1)) + 1.0;
        } else {
            // random, a third of the boxes empty
            v = pick(rng) < 33 ? 0.0 : static_cast<double>(pick(rng) + 1);
        }
        w[static_cast<size_t>(m)] = static_cast<int64_t>(std::llround(v));
    }
    return w;
}

// Keeping ranks against the rule, both with these weights: the levels with
// more and with fewer ranks, and the leaf imbalance.
struct KeepStats {
    long layouts = 0, more = 0, fewer = 0, over = 0;
    double worst = 0.0;
};
static void check_keeping(const RegionPlan& rule, int dim, int levels, int procs, int64_t threshold,
                          const RegionWeights& w, const char* weights, KeepStats& st) {
    RegionPlan k;
    try {
        k = plan_regions(dim, levels, procs, threshold, w, true);
    } catch (const std::exception&) {
        fail("no layout keeping ranks", weights, dim, levels, procs, -1);
        return;
    }
    ++st.layouts;
    check_layout(k, weights, dim, levels, procs);
    if (k.count[1] != 1) fail("keeping ranks: level 1 on one rank", weights, dim, levels, procs, 1);
    if (dim == 3 && threshold == 8 && region_count_is_cubic(dim, procs) && rule.bricks &&
        (k.count != rule.count || !k.bricks)) {
        fail("keeping ranks: (2^d)^k the bricks", weights, dim, levels, procs, -1);
    }
    bool more = false, fewer = false;
    for (int l = 0; l < levels; ++l) {
        more = more || k.count[static_cast<size_t>(l)] > rule.count[static_cast<size_t>(l)];
        fewer = fewer || k.count[static_cast<size_t>(l)] < rule.count[static_cast<size_t>(l)];
    }
    st.more += more;
    st.fewer += fewer;
    if (k.leaf_imbalance > kRegionImbalance && k.leaf_imbalance > rule.leaf_imbalance) ++st.over;
    if (!k.empty_region) st.worst = std::max(st.worst, k.leaf_imbalance);
}

int main() {
    const int64_t threshold = 8;
    // the reduction groups: region c below is in the group of its parent
    for (int32_t below = 1; below <= 300; ++below) {
        for (int32_t n = 1; n <= below; ++n) {
            for (int32_t c = 0; c < below; ++c) {
                const int32_t q = region_parent(c, n, below);
                if (q < 0 || q >= n || region_first_child(q, n, below) > c || c >= region_first_child(q + 1, n, below)) {
                    fail("region_parent", "none", 0, 0, below, n);
                    break;
                }
            }
        }
    }
    KeepStats keep_none, keep_shell;
    const char* const kind_name[3] = {"shell", "cluster", "random"};
    long checked = 0, raised = 0, over = 0, weighted = 0, w_over[3] = {0, 0, 0}, w_empty[3] = {0, 0, 0},
         left_bricks[3] = {0, 0, 0};
    double worst = 0.0, w_worst[3] = {0.0, 0.0, 0.0};
    int worst_dim = 0, worst_levels = 0, worst_procs = 0;
    std::mt19937_64 rng(20261015);
    for (int dim = 2; dim <= 3; ++dim) {
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
                const bool cubic = region_count_is_cubic(dim, procs);
                check_layout(p, "none", dim, levels, procs);
                if (cubic) {
                    if (!p.bricks || cubic_counts(dim, levels, procs, threshold) != p.count || p.threshold != threshold) {
                        fail("(2^d)^k: the bricks", "none", dim, levels, procs, -1);
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
                // the same weight everywhere: the same layout
                const int h = region_weight_level(dim, levels, procs);
                const size_t hb = size_t{1} << (dim * h);
                if (static_cast<int64_t>(hb) > kRegionWeightBoxes * procs && h + 1 < levels) {
                    fail("weight level too fine", "uniform", dim, levels, procs, h);
                }
                const RegionPlan u = plan_regions(dim, levels, procs, threshold,
                                                  RegionWeights::of_counts(h, std::vector<int64_t>(hb, 5)));
                if (u.count != p.count || u.start != p.start || u.threshold != p.threshold || u.bricks != p.bricks) {
                    fail("uniform weights change the layout", "uniform", dim, levels, procs, -1);
                }
                check_keeping(p, dim, levels, procs, threshold, RegionWeights{}, "none, keeping ranks", keep_none);
                for (int kind = 0; kind < 3; ++kind) {
                    const RegionWeights w = RegionWeights::of_counts(h, weights_of(kind, dim, h, rng));
                    if (w.prefix.back() == 0) continue;
                    RegionPlan q;
                    try {
                        q = plan_regions(dim, levels, procs, threshold, w);
                    } catch (const std::exception&) {
                        fail("no weighted layout", kind_name[kind], dim, levels, procs, -1);
                        continue;
                    }
                    ++weighted;
                    check_layout(q, kind_name[kind], dim, levels, procs);
                    if (kind == 0) check_keeping(q, dim, levels, procs, threshold, w, "shell, keeping ranks", keep_shell);
                    if (cubic && !q.bricks) ++left_bricks[kind];
                    if (q.empty_region) ++w_empty[kind];
                    if (q.leaf_imbalance > kRegionImbalance) ++w_over[kind];
                    if (leaf_boxes >= 512 * procs && !q.empty_region) {
                        w_worst[kind] = std::max(w_worst[kind], q.leaf_imbalance);
                    }
                }
            }
        }
    }
    std::printf("checked %ld layouts; other than (2^d)^k: threshold raised in %ld, leaf imbalance over %.2f in %ld\n",
                checked, raised, kRegionImbalance, over);
    std::printf("worst leaf imbalance with 64 leaf boxes a rank or more: %.3f (%dD, %d levels, %d ranks)\n", worst,
                worst_dim, worst_levels, worst_procs);
    std::printf("weighted by points: %ld layouts\n", weighted);
    for (int kind = 0; kind < 3; ++kind) {
        std::printf("  %-7s: over %.2f in %ld, a leaf region without points in %ld, (2^d)^k off the bricks in %ld; "
                    "worst with 512 leaf boxes a rank or more %.3f\n",
                    kind_name[kind], kRegionImbalance, w_over[kind], w_empty[kind], left_bricks[kind], w_worst[kind]);
    }
    for (const auto& [name, st] : {std::pair<const char*, const KeepStats&>{"none", keep_none},
                                   std::pair<const char*, const KeepStats&>{"shell", keep_shell}}) {
        std::printf("keeping ranks (%s): %ld layouts, more ranks on some level in %ld, fewer on some level in %ld, "
                    "leaf imbalance over %.2f and over the rule's in %ld, worst %.3f\n",
                    name, st.layouts, st.more, st.fewer, kRegionImbalance, st.over, st.worst);
    }
    std::printf("dim levels ranks: threshold, leaf imbalance, ranks by level (none | shell | none, keeping ranks)\n");
    for (const auto& [dim, levels, procs] : std::vector<std::tuple<int, int, int>>{
             {3, 6, 8}, {3, 6, 12}, {3, 6, 16}, {3, 6, 20}, {3, 6, 24}, {3, 6, 50}, {3, 7, 64}, {3, 7, 100}, {3, 8, 192},
             {2, 8, 50}, {2, 10, 100}}) {
        const RegionPlan p = plan_regions(dim, levels, procs, threshold);
        const int h = region_weight_level(dim, levels, procs);
        const RegionPlan q =
            plan_regions(dim, levels, procs, threshold, RegionWeights::of_counts(h, weights_of(0, dim, h, rng)));
        std::printf("%dD %d %4d: %3ld %.3f  ", dim, levels, procs, static_cast<long>(p.threshold), p.leaf_imbalance);
        for (const int32_t x : p.count) std::printf("%d ", x);
        std::printf("| %3ld %.3f%s  ", static_cast<long>(q.threshold), q.leaf_imbalance, q.bricks ? " bricks" : "");
        for (const int32_t x : q.count) std::printf("%d ", x);
        const RegionPlan k = plan_regions(dim, levels, procs, threshold, RegionWeights{}, true);
        std::printf("| %3ld %.3f%s  ", static_cast<long>(k.threshold), k.leaf_imbalance, k.bricks ? " bricks" : "");
        for (const int32_t x : k.count) std::printf("%d ", x);
        std::printf("\n");
    }
    if (failures) {
        std::printf("%d failures\n", failures);
        return 1;
    }
    std::printf("all checks passed\n");
    return 0;
}
