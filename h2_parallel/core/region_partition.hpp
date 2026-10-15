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
// Storage: per level, count + 1 range starts - the size of the process count,
// like morton_to_rank; nothing grows with the number of boxes.
#include <algorithm>
#include <cstdint>
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

// The region starts of every level: start[level] has count[level] + 1
// entries, region i holding boxes [start[i], start[i + 1]) of the level.
// Empty when a reduction's parent region has fewer boxes than child regions
// (no split in whole parent boxes).
inline std::vector<std::vector<int64_t>> region_starts(int dimension, const std::vector<int32_t>& count) {
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
        if (n != (below + f - 1) / f) {
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
            // its children: whole boxes of this level, the boundaries at the
            // nearest one to their share of the leaf regions, one at least each
            const int32_t c0 = region_first_child(i, n, below);
            const int32_t c1 = region_first_child(i + 1, n, below);
            const int64_t m = c1 - c0;
            if (m < 1 || boxes < m) return {};
            int64_t total = 0;
            for (int32_t c = c0; c < c1; ++c) total += weight[static_cast<size_t>(level + 1)][static_cast<size_t>(c)];
            int64_t cum = 0, prev = -1;
            for (int32_t c = c0; c < c1; ++c) {
                // (2 boxes cum + total) / (2 total): boxes cum / total rounded
                int64_t at = static_cast<int64_t>((static_cast<__int128>(2) * boxes * cum + total) / (2 * total));
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

// The largest leaf region over the mean leaf region.
inline double region_leaf_imbalance(const std::vector<std::vector<int64_t>>& start) {
    if (start.empty()) return 0.0;
    const auto& t = start.back();
    const int64_t n = static_cast<int64_t>(t.size()) - 1;
    int64_t largest = 0;
    for (int64_t i = 0; i < n; ++i) largest = std::max(largest, t[static_cast<size_t>(i) + 1] - t[static_cast<size_t>(i)]);
    return static_cast<double>(largest) * static_cast<double>(n) / static_cast<double>(t.back());
}

// The largest leaf region a layout may have against the mean before
// plan_regions takes the first reduction a level finer.
inline constexpr double kRegionImbalance = 1.15;

// The layout of every level: counts, starts, and the threshold that made them
// (the caller's, or a multiple of 2^d of it - see the top of this file).
struct RegionPlan {
    std::vector<int32_t> count;
    std::vector<std::vector<int64_t>> start;
    int64_t threshold = 0;
    double leaf_imbalance = 1.0;
};
inline RegionPlan plan_regions(int dimension, int num_levels, int num_procs, int64_t threshold) {
    const int64_t f = int64_t{1} << dimension;
    RegionPlan plan;
    if (region_count_is_cubic(dimension, num_procs)) {
        // (2^d)^k: the caller's threshold and the equal ranges, the bricks as
        // before (also where a brick is a single box)
        plan.count = region_counts(dimension, num_levels, num_procs, threshold);
        plan.threshold = threshold;
        plan.start.resize(static_cast<size_t>(num_levels));
        for (int level = 0; level < num_levels; ++level) {
            const int64_t boxes = int64_t{1} << (dimension * level);
            const int64_t n = plan.count[static_cast<size_t>(level)];
            auto& t = plan.start[static_cast<size_t>(level)];
            t.resize(static_cast<size_t>(n) + 1);
            for (int64_t i = 0; i <= n; ++i) t[static_cast<size_t>(i)] = i * (boxes / n);
        }
        return plan;
    }
    // any other count: the caller's threshold while the leaf regions are within
    // kRegionImbalance of the mean, else the first reduction a level finer
    // (the threshold times 2^d), until it is or the threshold passes the leaf
    // boxes; then the most even layout seen
    const int64_t leaf_boxes = int64_t{1} << (dimension * (num_levels - 1));
    plan.leaf_imbalance = -1.0;
    std::vector<int32_t> seen;
    for (int64_t t = threshold; t <= std::max(threshold, leaf_boxes); t *= f) {
        RegionPlan p;
        p.count = region_counts(dimension, num_levels, num_procs, t);
        if (p.count == seen) continue;
        seen = p.count;
        p.threshold = t;
        p.start = region_starts(dimension, p.count);
        if (p.start.empty()) continue;
        p.leaf_imbalance = region_leaf_imbalance(p.start);
        if (plan.leaf_imbalance < 0.0 || p.leaf_imbalance < plan.leaf_imbalance) plan = p;
        if (p.leaf_imbalance <= kRegionImbalance) break;
    }
    if (plan.leaf_imbalance < 0.0) {
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
