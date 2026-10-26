#pragma once
// The regions of a level for any number of active ranks (the unstructured
// Color path; the structured and CA levels keep (2^dimension)^k).
//
// A region is the boxes of one active rank on a level: a contiguous range of
// the level's Morton order, as before, but no longer one of equal cubes.
// The counts follow the reduction rule of compute_active_processes - all
// ranks at the leaf, divided by 2^dimension (rounded up) where a level has
// fewer than the threshold's boxes per rank - and the regions of the level
// below are grouped under the reduced level's as evenly as possible
// (region_first_child: 50 under 7 is 7 groups of 7 or 8, not six of 8 and
// one of 2, whose parent's share of a coarse level would be a box or two).
//
// The ranges are split top down.  Where the count does not change, a
// region's range one level finer is its range times 2^d, so a box's children
// are on its rank.  Where it grows (going down; a reduction going up) the
// parent region's boxes are divided among its child regions in whole boxes
// of the parent's level, in proportion to the leaf regions under each and at
// least one each: every child rank builds whole parent boxes from its own
// children, as the transitions (CPU and GPU), the H2 construction and the
// solve's gathers expect.  With a count of (2^d)^k every group is 2^d
// regions and every split even: the ranges are the equal ones of
// distribute_boxes, the bricks the rest of the code knew.
//
// Whole parent boxes are coarse where a reduction happens (fewer than the
// threshold's boxes per rank of that level by construction): with 50 ranks
// in 3D the leaf regions are whole level-2 boxes, 64 for 50 ranks, and the
// largest leaf region is 2.3x the mean.  plan_regions therefore takes a
// count's first reduction a level finer (the threshold times 2^d) while the
// largest leaf region is more than kRegionImbalance of the mean: 50 ranks
// then reduce at level 3 (512 boxes) and the leaf regions are within 1.07.
//
// On the unstructured path the counts keep more ranks on the coarse levels
// (region_counts_keeping_ranks): a reduction divides the count by the
// smallest k <= 2^d that keeps the leaf regions even, not always by 2^d.
//
// The splits weigh the boxes by a load (RegionWeights): on the unstructured
// path the points in each box, so that empty boxes count for nothing and
// every leaf region gets about the same number of points; without weights
// every box counts one, as above.  (2^d)^k ranks keep the bricks unless
// their points are more uneven than kRegionImbalance (or a brick has none).
//
// Storage: per level, count + 1 range starts - the size of the process count,
// like morton_to_rank; the weights are the points per box of one level with at
// most kRegionWeightBoxes boxes a rank; nothing grows with the number of
// leaf boxes.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace fmm {

// Active ranks per level (index = level, leaf last) for `num_procs` ranks:
// the leaf has all of them; a level with fewer than `threshold` boxes per
// rank of the level below has that count divided by 2^dimension, rounded up.
inline std::vector<int32_t> region_counts(int dimension, int num_levels, int num_procs, int64_t threshold) {
    if (dimension < 1 || dimension > 3) throw std::invalid_argument("region_counts: dimension must be 1, 2 or 3");
    if (num_levels <= 0 || num_procs <= 0) throw std::invalid_argument("region_counts: no levels or no ranks");
    const int32_t f = int32_t{1} << dimension;
    const int leaf = num_levels - 1;
    std::vector<int32_t> count(static_cast<size_t>(num_levels), 1);
    count[static_cast<size_t>(leaf)] = num_procs;
    const int64_t leaf_boxes = int64_t{1} << (dimension * leaf);
    if (leaf_boxes < num_procs) {
        throw std::invalid_argument("H2 tree: " + std::to_string(num_procs) + " ranks but only " +
                                    std::to_string(leaf_boxes) + " leaf boxes");
    }
    for (int level = leaf - 1; level >= 0; --level) {
        const int64_t boxes = int64_t{1} << (dimension * level);
        const int32_t prev = count[static_cast<size_t>(level + 1)];
        const int64_t per_proc = boxes / prev;
        count[static_cast<size_t>(level)] = per_proc < threshold ? std::max(1, (prev + f - 1) / f) : prev;
    }
    return count;
}

// The child regions of region `i` of a level with `n` regions over `below`
// regions one level finer: [first, last).  Equal groups when n == below or
// below = 2^d n; else groups of floor or ceil(below / n).
inline int32_t region_first_child(int32_t i, int32_t n, int32_t below) {
    if (n == below) return i;
    return static_cast<int32_t>(static_cast<int64_t>(i) * below / n);
}
// The region of the reduced level over region `c` of the level below.
inline int32_t region_parent(int32_t c, int32_t n, int32_t below) {
    if (n == below) return c;
    return static_cast<int32_t>((static_cast<int64_t>(c + 1) * n + below - 1) / below - 1);
}

// The load of the boxes: the points in each box of one level (`level`, as
// prefix sums), spread evenly over the boxes of a finer level; without a
// level (empty()) every box of every level counts one.
struct RegionWeights {
    int level = -1;
    std::vector<int64_t> prefix;  // boxes of `level` + 1
    bool empty() const { return level < 0; }
    static RegionWeights of_counts(int level, const std::vector<int64_t>& per_box) {
        RegionWeights w;
        w.level = level;
        w.prefix.assign(per_box.size() + 1, 0);
        for (size_t i = 0; i < per_box.size(); ++i) {
            if (per_box[i] < 0) throw std::invalid_argument("RegionWeights: a negative weight");
            w.prefix[i + 1] = w.prefix[i] + per_box[i];
        }
        return w;
    }
    // The load of boxes [a, b) of level `lvl`: the points (in units of
    // 1 / 2^(dimension (lvl - level)) of a point on a level finer than `level`;
    // loads are compared within one level only).
    __int128 load(int dimension, int lvl, int64_t a, int64_t b) const {
        if (level < 0) return b - a;
        if (lvl <= level) {
            const int s = dimension * (level - lvl);
            return prefix[static_cast<size_t>(b << s)] - prefix[static_cast<size_t>(a << s)];
        }
        const int s = dimension * (lvl - level);
        const int64_t g = int64_t{1} << s;
        auto upto = [&](int64_t x) {
            const size_t k = static_cast<size_t>(x >> s);
            const int64_t r = x & (g - 1);
            __int128 v = static_cast<__int128>(prefix[k]) * g;
            if (r != 0) v += static_cast<__int128>(prefix[k + 1] - prefix[k]) * r;
            return v;
        };
        return upto(b) - upto(a);
    }
};

// The boxes a rank the weights' level may have at most: the finest level with
// no more than this many boxes a rank (the leaf if it has fewer) carries them.
inline constexpr int64_t kRegionWeightBoxes = 512;
inline int region_weight_level(int dimension, int num_levels, int num_procs) {
    int level = 0;
    while (level + 1 < num_levels && (int64_t{1} << (dimension * (level + 1))) <= kRegionWeightBoxes * num_procs) ++level;
    return level;
}

// The region starts of every level: start[level] has count[level] + 1
// entries, region i holding boxes [start[i], start[i + 1]) of the level.
// Empty when a reduction's parent region has fewer boxes than child regions
// (no split in whole parent boxes).
inline std::vector<std::vector<int64_t>> region_starts(int dimension, const std::vector<int32_t>& count,
                                                       const RegionWeights& weights = RegionWeights{}) {
    const int num_levels = static_cast<int>(count.size());
    const int64_t f = int64_t{1} << dimension;
    if (num_levels == 0 || count[0] != 1) throw std::invalid_argument("region_starts: the root level needs one region");
    // the leaf regions under each region (weights of the splits)
    std::vector<std::vector<int64_t>> weight(static_cast<size_t>(num_levels));
    weight[static_cast<size_t>(num_levels - 1)].assign(static_cast<size_t>(count.back()), 1);
    for (int level = num_levels - 2; level >= 0; --level) {
        const int32_t n = count[static_cast<size_t>(level)], below = count[static_cast<size_t>(level + 1)];
        auto& w = weight[static_cast<size_t>(level)];
        if (n == below) {
            w = weight[static_cast<size_t>(level + 1)];
            continue;
        }
        if (n < 1 || n > below) {
            throw std::invalid_argument("region_starts: level " + std::to_string(level) + " has " + std::to_string(n) +
                                        " regions for " + std::to_string(below) + " below");
        }
        w.assign(static_cast<size_t>(n), 0);
        for (int32_t c = 0; c < below; ++c) {
            w[static_cast<size_t>(region_parent(c, n, below))] += weight[static_cast<size_t>(level + 1)][static_cast<size_t>(c)];
        }
    }
    std::vector<std::vector<int64_t>> start(static_cast<size_t>(num_levels));
    start[0] = {0, 1};
    for (int level = 0; level + 1 < num_levels; ++level) {
        const auto& s = start[static_cast<size_t>(level)];
        const int32_t n = count[static_cast<size_t>(level)], below = count[static_cast<size_t>(level + 1)];
        auto& t = start[static_cast<size_t>(level + 1)];
        t.assign(static_cast<size_t>(below) + 1, 0);
        t[static_cast<size_t>(below)] = s[static_cast<size_t>(n)] * f;
        for (int32_t i = 0; i < n; ++i) {
            const int64_t first = s[static_cast<size_t>(i)];                                // (boxes of this level)
            const int64_t boxes = s[static_cast<size_t>(i) + 1] - s[static_cast<size_t>(i)];
            if (n == below) {
                t[static_cast<size_t>(i)] = first * f;
                continue;
            }
            // its children: whole boxes of this level, each boundary where the
            // load before it is nearest to the children's share (by the leaf
            // regions under each) of the region's load, one box at least each
            const int32_t c0 = region_first_child(i, n, below);
            const int32_t c1 = region_first_child(i + 1, n, below);
            const int64_t m = c1 - c0;
            if (m < 1 || boxes < m) return {};
            int64_t total = 0;
            for (int32_t c = c0; c < c1; ++c) total += weight[static_cast<size_t>(level + 1)][static_cast<size_t>(c)];
            const __int128 region_load = weights.load(dimension, level, first, first + boxes);
            int64_t cum = 0, prev = -1;
            for (int32_t c = c0; c < c1; ++c) {
                // the load of the first `at` boxes against the share, times total
                auto miss = [&](int64_t at) {
                    return weights.load(dimension, level, first, first + at) * total - region_load * cum;
                };
                // the first box past the share, and the one before; the nearer,
                // the later on a tie (without weights: boxes cum / total rounded)
                int64_t lo = 0, hi = boxes + 1;
                while (lo < hi) {
                    const int64_t mid = lo + (hi - lo) / 2;
                    if (mid <= boxes && miss(mid) <= 0) lo = mid + 1; else hi = mid;
                }
                int64_t at = lo - 1;
                if (lo <= boxes && miss(lo) <= -miss(at)) at = lo;
                if (c == c0) at = 0;
                at = std::max(at, prev + 1);                    // a box each so far ...
                at = std::min(at, boxes - (c1 - c));            // ... and one left for each to come
                t[static_cast<size_t>(c)] = (first + at) * f;
                prev = at;
                cum += weight[static_cast<size_t>(level + 1)][static_cast<size_t>(c)];
            }
        }
    }
    return start;
}

// Whether a count is (2^dimension)^k, the layout of equal cubes (bricks) the
// structured and CA code require.
inline bool region_count_is_cubic(int dimension, int64_t count) {
    if (count <= 0 || (count & (count - 1)) != 0) return false;
    return __builtin_ctzll(static_cast<unsigned long long>(count)) % dimension == 0;
}

// The largest leaf region's load over the mean; and whether one has none.
inline double region_leaf_imbalance(int dimension, const std::vector<std::vector<int64_t>>& start,
                                    const RegionWeights& weights = RegionWeights{}, bool* empty_region = nullptr) {
    if (start.empty()) return 0.0;
    const int leaf = static_cast<int>(start.size()) - 1;
    const auto& t = start.back();
    const int64_t n = static_cast<int64_t>(t.size()) - 1;
    __int128 largest = 0, total = 0;
    bool empty = false;
    for (int64_t i = 0; i < n; ++i) {
        const __int128 w = weights.load(dimension, leaf, t[static_cast<size_t>(i)], t[static_cast<size_t>(i) + 1]);
        largest = std::max(largest, w);
        total += w;
        empty = empty || w == 0;
    }
    if (empty_region != nullptr) *empty_region = empty;
    if (total == 0) return 0.0;
    return static_cast<double>(largest) * static_cast<double>(n) / static_cast<double>(total);
}

// Whether a level's regions are the equal cubes: a count of (2^d)^k and the
// equal ranges (the layout the structured, CA and device-reduction code need).
inline bool region_starts_are_bricks(int dimension, int level, const std::vector<int64_t>& start) {
    const int64_t n = static_cast<int64_t>(start.size()) - 1;
    if (!region_count_is_cubic(dimension, n)) return false;
    const int64_t boxes = int64_t{1} << (dimension * level);
    if (boxes % n != 0) return false;
    for (int64_t i = 0; i <= n; ++i) {
        if (start[static_cast<size_t>(i)] != i * (boxes / n)) return false;
    }
    return true;
}

// The largest leaf region a layout may have against the mean before
// plan_regions takes the first reduction a level finer.
inline constexpr double kRegionImbalance = 1.15;

// The limit in use: kRegionImbalance, or the environment's H2_REGION_IMBALANCE
// (>= 1).  A larger limit keeps more ranks on the coarse levels of uneven
// geometries (whole parent boxes as leaf regions), at the price of a less even
// leaf level.  Every rank must see the same value.
inline double region_imbalance_limit() {
    static const double limit = [] {
        if (const char* v = std::getenv("H2_REGION_IMBALANCE")) {
            char* end = nullptr;
            const double x = std::strtod(v, &end);
            if (end != v && x >= 1.0) return x;
        }
        return kRegionImbalance;
    }();
    return limit;
}

// The counts of the levels above `level` by the reduction rule of
// region_counts, from count[level].
inline void region_counts_above(int dimension, int level, int64_t threshold, std::vector<int32_t>& count) {
    const int32_t f = int32_t{1} << dimension;
    for (int l = level - 1; l >= 0; --l) {
        const int64_t boxes = int64_t{1} << (dimension * l);
        const int32_t prev = count[static_cast<size_t>(l + 1)];
        count[static_cast<size_t>(l)] = boxes / prev < threshold ? std::max(1, (prev + f - 1) / f) : prev;
    }
}

// Active ranks per level keeping more of them on the coarse levels (the
// unstructured path).  Where region_counts divides a count by 2^d, take
// instead the largest ceil(count / k), k = 2 .. 2^d, that leaves a rank of
// the level `threshold` boxes or more (and at most 2^d (level - 1) ranks, so
// that level 1 is one) and whose layout - the levels above by the rule -
// has its leaf regions within kRegionImbalance of the mean or no further
// off than with 2^d (and no region without load that 2^d avoids); else
// 2^d.  Whole parent boxes make this a matter of factors: 12 ranks keep 4 on
// level 2 of a 3D tree (groups of 3 under 2 level-1 boxes each), not 8
// (groups of 1 and 2 under one box each: one leaf region 1.5x the mean) nor
// the rule's 2.  Leaf first, each level's choice fixed before the next.
inline std::vector<int32_t> region_counts_keeping_ranks(int dimension, int num_levels, int num_procs,
                                                        int64_t threshold, const RegionWeights& weights) {
    const int32_t f = int32_t{1} << dimension;
    std::vector<int32_t> count = region_counts(dimension, num_levels, num_procs, threshold);
    for (int level = num_levels - 2; level >= 2; --level) {
        const int64_t boxes = int64_t{1} << (dimension * level);
        const int32_t prev = count[static_cast<size_t>(level + 1)];
        if (boxes / prev >= threshold) continue;  // no reduction here (count[level] == prev)
        auto layout = [&](int32_t n, double& imbalance, bool& empty) {
            std::vector<int32_t> c = count;
            c[static_cast<size_t>(level)] = n;
            region_counts_above(dimension, level, threshold, c);
            const auto start = region_starts(dimension, c, weights);
            imbalance = start.empty() ? -1.0 : region_leaf_imbalance(dimension, start, weights, &empty);
            return c;
        };
        const int32_t plain = std::max(1, (prev + f - 1) / f);
        double plain_imbalance = 0.0;
        bool plain_empty = false;
        std::vector<int32_t> chosen = layout(plain, plain_imbalance, plain_empty);
        const int64_t cap = std::min<int64_t>(boxes / threshold, int64_t{1} << (dimension * (level - 1)));
        int32_t tried = 0;
        for (int32_t k = 2; k < f; ++k) {
            const int32_t n = (prev + k - 1) / k;  // decreasing in k
            if (n <= plain) break;
            if (n > cap || n == tried) continue;
            tried = n;
            double imbalance = 0.0;
            bool empty = false;
            std::vector<int32_t> c = layout(n, imbalance, empty);
            if (imbalance < 0.0 || (empty && !plain_empty)) continue;
            const double limit = plain_imbalance < 0.0 ? region_imbalance_limit()
                                                       : std::max(region_imbalance_limit(), plain_imbalance);
            if (imbalance <= limit) {
                chosen = c;
                break;
            }
        }
        count = chosen;
    }
    return count;
}

// The layout of every level: counts, starts, and the threshold that made them
// (the caller's, or a multiple of 2^d of it - see the top of this file).
struct RegionPlan {
    std::vector<int32_t> count;
    std::vector<std::vector<int64_t>> start;
    int64_t threshold = 0;
    double leaf_imbalance = 1.0;   // by the weights
    bool empty_region = false;     // a leaf region with no load
    bool bricks = false;           // the equal cubes of (2^d)^k ranks
};
inline RegionPlan plan_regions(int dimension, int num_levels, int num_procs, int64_t threshold,
                               const RegionWeights& weights = RegionWeights{}, bool keep_ranks = false) {
    auto counts_at = [&](int64_t t) {
        return keep_ranks ? region_counts_keeping_ranks(dimension, num_levels, num_procs, t, weights)
                          : region_counts(dimension, num_levels, num_procs, t);
    };
    const int64_t f = int64_t{1} << dimension;
    // the better of two layouts: a leaf region with no load is worse than
    // any imbalance; then the smaller imbalance; the first on a tie
    auto better = [](const RegionPlan& a, const RegionPlan& b) {
        if (a.empty_region != b.empty_region) return !a.empty_region;
        return a.leaf_imbalance < b.leaf_imbalance;
    };
    auto good = [](const RegionPlan& p) { return !p.empty_region && p.leaf_imbalance <= region_imbalance_limit(); };
    // the layout of these counts (no starts where whole parent boxes cannot be)
    auto build = [&](const std::vector<int32_t>& count, int64_t t) {
        RegionPlan p;
        p.count = count;
        p.threshold = t;
        p.start = region_starts(dimension, p.count, weights);
        if (p.start.empty()) return p;
        p.leaf_imbalance = region_leaf_imbalance(dimension, p.start, weights, &p.empty_region);
        p.bricks = region_count_is_cubic(dimension, num_procs);
        for (int level = 0; p.bricks && level < num_levels; ++level) {
            p.bricks = region_starts_are_bricks(dimension, level, p.start[static_cast<size_t>(level)]);
        }
        return p;
    };
    RegionPlan plan;
    bool have = false;
    if (region_count_is_cubic(dimension, num_procs)) {
        // (2^d)^k: the caller's threshold and the equal ranges, the bricks as
        // before (also where a brick is a single box) - kept unless the
        // weights make them uneven or leave one empty
        plan.count = region_counts(dimension, num_levels, num_procs, threshold);
        plan.threshold = threshold;
        plan.bricks = true;
        plan.start.resize(static_cast<size_t>(num_levels));
        for (int level = 0; level < num_levels; ++level) {
            const int64_t boxes = int64_t{1} << (dimension * level);
            const int64_t n = plan.count[static_cast<size_t>(level)];
            auto& t = plan.start[static_cast<size_t>(level)];
            t.resize(static_cast<size_t>(n) + 1);
            for (int64_t i = 0; i <= n; ++i) t[static_cast<size_t>(i)] = i * (boxes / n);
        }
        if (weights.empty() && !keep_ranks) return plan;
        plan.leaf_imbalance = region_leaf_imbalance(dimension, plan.start, weights, &plan.empty_region);
        if (good(plan)) {
            // keeping ranks: a layout with more on the coarse levels at this
            // threshold where it is as even (other thresholds keep fewer)
            if (keep_ranks) {
                const std::vector<int32_t> count = counts_at(threshold);
                if (count != plan.count) {
                    RegionPlan p = build(count, threshold);
                    if (!p.start.empty() && good(p)) return p;
                }
            }
            return plan;
        }
        have = true;
    }
    // the caller's threshold while the leaf regions are within
    // kRegionImbalance of the mean, else the first reduction a level finer
    // (the threshold times 2^d), until they are or the threshold passes the
    // leaf boxes; then the best layout seen
    const int64_t leaf_boxes = int64_t{1} << (dimension * (num_levels - 1));
    std::vector<int32_t> seen;
    for (int64_t t = threshold; t <= std::max(threshold, leaf_boxes); t *= f) {
        const std::vector<int32_t> count = counts_at(t);
        if (count == seen) continue;
        seen = count;
        RegionPlan p = build(count, t);
        if (p.start.empty()) continue;
        if (!have || better(p, plan)) plan = p;
        have = true;
        if (good(p)) break;
    }
    if (!have) {
        throw std::invalid_argument("H2 tree: no region layout for " + std::to_string(num_procs) + " ranks over " +
                                    std::to_string(num_levels) + " levels");
    }
    return plan;
}

// The region of box `morton` on a level with these starts (count + 1).
inline int32_t region_of_box(const std::vector<int64_t>& start, int64_t morton) {
    const auto it = std::upper_bound(start.begin(), start.end(), morton);
    if (it == start.begin() || it == start.end()) {
        throw std::out_of_range("region_of_box: box " + std::to_string(morton) + " outside the level");
    }
    return static_cast<int32_t>(it - start.begin()) - 1;
}

}  // namespace fmm
