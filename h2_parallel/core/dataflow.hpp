#pragma once
// ---------------------------------------------------------------------------
// dataflow.hpp — component-owner dataflow engine (MPI-free core).
//
// Design: concise-algorithm-ref/component-owner-dataflow-design.md.
//
// The boundary work of one level is a small static DAG:
//   blue clusters (corner, incl. orange)  ->  purple edge runs  ->  green face
//   patches  ->  each rank's interior.
// Every component is eliminated exactly once, by one *owner* chosen among the
// ranks it straddles; the owner ships generators (+ the factor section for
// boxes whose home is another rank).  Execution is counter driven: a
// component becomes ready when the generators of all its predecessors are
// available on the owner; the ready queue is ordered by priority.
//
// This header contains everything that is a pure function of the process
// grid and therefore identical on every rank without communication:
//   ProcessGrid      rank <-> brick geometry (same mapping as morton::
//                    assign_to_processes_*: pid = morton of the brick coords)
//   ComponentGraph   components (26-connected flood fill per color class over
//                    the classification the CA code already uses), sharer and
//                    adjacent-rank sets, predecessors / successors, per-rank
//                    interior join node, canonical ids
//   assign_owners    deterministic owner per component (pluggable policy)
//   FetchPlan        per owned component: which remote boxes to fetch from
//                    whom (staged by color), and the mirror "serve" plan
//   Runtime          per-rank counters, priority ready queue, expected
//                    arrivals, per-component time stamps, serial oracle order
//
// Nothing numerical lives here.  MPI transport (arrival loop, push) sits in
// the driver/serialization layer and calls into Runtime.
//
// Classification is injected (ColorOf callback) so that this header does not
// depend on tree_impl.hpp and can be unit-tested; production wraps
// is_blue_box / is_orange_box / is_purple_box (see make_ca_classifier).
// ---------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "morton.hpp"

namespace fmm {
namespace dataflow {

// ===========================================================================
// Colors
// ===========================================================================

/// Elimination classes in CA order.  BLUE includes the orange layer (they are
/// eliminated together as one group today: `blue_with_orange`).
enum class Color : uint8_t { BLUE = 0, PURPLE = 1, GREEN = 2, INTERIOR = 3, NONE = 255 };

inline const char* color_name(Color c) {
    switch (c) {
        case Color::BLUE: return "blue";
        case Color::PURPLE: return "purple";
        case Color::GREEN: return "green";
        case Color::INTERIOR: return "interior";
        default: return "none";
    }
}

/// Number of boundary colors for a dimension (blue, [purple], green).
inline int num_boundary_colors(int dimension) { return dimension == 3 ? 3 : 2; }

/// Boundary colors in elimination order for a dimension.
inline std::vector<Color> boundary_colors(int dimension) {
    if (dimension == 3) return {Color::BLUE, Color::PURPLE, Color::GREEN};
    return {Color::BLUE, Color::GREEN};
}

// ===========================================================================
// Process grid
// ===========================================================================

/// Geometry of one level's process decomposition.  Ranks are identified by
/// their *morton process id* (the value `TreeLevel::my_morton_id`), which is
/// the Morton code of the brick coordinates — exactly what
/// morton::assign_to_processes_{2d,3d} returns.
struct ProcessGrid {
    int dimension = 3;
    uint32_t grid_size = 0;      ///< boxes per dimension at this level (1 << level)
    uint32_t num_procs = 1;      ///< active processes (4^k in 2D, 8^k in 3D)
    uint32_t procs_per_dim = 1;
    uint32_t brick = 0;          ///< boxes per dimension per process

    ProcessGrid() = default;
    ProcessGrid(int dim, uint32_t grid, uint32_t nprocs)
        : dimension(dim), grid_size(grid), num_procs(nprocs) {
        if (dim != 2 && dim != 3) throw std::invalid_argument("ProcessGrid: dimension must be 2 or 3");
        if (nprocs == 0 || (nprocs & (nprocs - 1)) != 0)
            throw std::invalid_argument("ProcessGrid: num_procs must be a power of two");
        const uint32_t k = static_cast<uint32_t>(__builtin_ctz(nprocs)) / static_cast<uint32_t>(dim);
        if ((1u << (k * dim)) != nprocs)
            throw std::invalid_argument("ProcessGrid: num_procs must be 4^k (2D) or 8^k (3D)");
        procs_per_dim = 1u << k;
        if (grid_size % procs_per_dim != 0)
            throw std::invalid_argument("ProcessGrid: grid_size not divisible by procs_per_dim");
        brick = grid_size / procs_per_dim;
    }

    int64_t boxes_per_proc() const {
        int64_t b = brick;
        return dimension == 3 ? b * b * b : b * b;
    }

    std::array<uint32_t, 3> decode(int64_t morton) const {
        uint32_t x = 0, y = 0, z = 0;
        if (dimension == 2) morton::decode_2d(static_cast<uint64_t>(morton), x, y);
        else morton::decode_3d(static_cast<uint64_t>(morton), x, y, z);
        return {x, y, z};
    }

    int64_t encode(uint32_t x, uint32_t y, uint32_t z) const {
        return dimension == 2 ? static_cast<int64_t>(morton::encode_2d(x, y))
                              : static_cast<int64_t>(morton::encode_3d(x, y, z));
    }

    /// Morton process id of the brick containing a box.
    uint32_t proc_of_box(int64_t morton) const {
        const auto c = decode(morton);
        const uint32_t px = c[0] / brick, py = c[1] / brick, pz = c[2] / brick;
        return static_cast<uint32_t>(dimension == 2 ? morton::encode_2d(px, py)
                                                    : morton::encode_3d(px, py, pz));
    }

    /// Brick coordinates (in units of bricks) of a process id.
    std::array<uint32_t, 3> proc_coords(uint32_t pid) const {
        uint32_t x = 0, y = 0, z = 0;
        if (dimension == 2) morton::decode_2d(pid, x, y);
        else morton::decode_3d(pid, x, y, z);
        return {x, y, z};
    }

    /// Global box offset of a process's brick (== `local_offset` in tree_impl).
    std::array<uint32_t, 3> proc_offset(uint32_t pid) const {
        auto c = proc_coords(pid);
        return {c[0] * brick, c[1] * brick, c[2] * brick};
    }

    /// Same as tree_impl's `get_box_owner_offset`.
    std::array<uint32_t, 3> box_owner_offset(int64_t morton) const {
        const auto c = decode(morton);
        return {(c[0] / brick) * brick, (c[1] / brick) * brick, (c[2] / brick) * brick};
    }

    /// 26- (or 8-) neighbors inside the global grid, excluding self.
    std::vector<uint64_t> neighbors(int64_t morton) const {
        return dimension == 2 ? morton::neighbors_2d(static_cast<uint64_t>(morton), grid_size)
                              : morton::neighbors_3d(static_cast<uint64_t>(morton), grid_size);
    }

    /// Chebyshev distance (in boxes) from a box to a process's brick; 0 inside.
    uint32_t brick_distance(uint32_t pid, int64_t morton) const {
        const auto c = decode(morton);
        const auto o = proc_offset(pid);
        uint32_t d = 0;
        for (int k = 0; k < dimension; ++k) {
            const int64_t lo = o[static_cast<size_t>(k)], hi = lo + brick - 1, x = c[static_cast<size_t>(k)];
            const int64_t dk = x < lo ? lo - x : (x > hi ? x - hi : 0);
            d = std::max(d, static_cast<uint32_t>(dk));
        }
        return d;
    }

    /// How far the CA halo reaches beyond a brick: boundary boxes within this
    /// Chebyshev distance of the brick are ghosts on that rank (tree_impl
    /// compute_ghost_boxes_and_colors: 1-hop always, 2-hop boundary boxes in 3D).
    uint32_t halo_reach() const { return dimension == 3 ? 2u : 1u; }

    /// Whether a boundary box is held (local or ghost) by a process.
    bool holds_boundary_box(uint32_t pid, int64_t morton) const {
        return brick_distance(pid, morton) <= halo_reach();
    }

    /// Boxes within two hops (5^d region minus self), inside the global grid.
    std::vector<uint64_t> neighbors_2hop(int64_t morton) const {
        return dimension == 2 ? morton::neighbors_2hop_2d(static_cast<uint64_t>(morton), grid_size)
                              : morton::neighbors_2hop_3d(static_cast<uint64_t>(morton), grid_size);
    }

    /// Whether a box lies on the boundary of its brick.
    bool on_brick_boundary(int64_t morton) const {
        const auto c = decode(morton);
        for (int d = 0; d < dimension; ++d) {
            const uint32_t l = c[static_cast<size_t>(d)] % brick;
            if (l == 0 || l == brick - 1) return true;
        }
        return false;
    }

    /// All boundary boxes of one brick, in increasing Morton order.
    std::vector<int64_t> brick_boundary_boxes(uint32_t pid) const {
        std::vector<int64_t> out;
        const auto o = proc_offset(pid);
        const uint32_t zmax = dimension == 3 ? brick : 1;
        for (uint32_t z = 0; z < zmax; ++z)
            for (uint32_t y = 0; y < brick; ++y)
                for (uint32_t x = 0; x < brick; ++x) {
                    const bool bx = (x == 0 || x == brick - 1);
                    const bool by = (y == 0 || y == brick - 1);
                    const bool bz = dimension == 3 && (z == 0 || z == brick - 1);
                    if (!(bx || by || bz)) continue;
                    out.push_back(encode(o[0] + x, o[1] + y, dimension == 3 ? o[2] + z : 0));
                }
        std::sort(out.begin(), out.end());
        return out;
    }
};

// ===========================================================================
// Classification callback
// ===========================================================================

/// Boundary-box color class as seen from the box's *home* brick.  Must return
/// BLUE for blue and orange boxes, PURPLE for purple, GREEN for the rest of
/// the boundary, NONE for interior boxes.
using ColorOf = std::function<Color(int64_t morton)>;

/// Build a ColorOf from the three CA predicates (signature of
/// is_blue_box / is_orange_box in tree_impl.hpp:
///   bool(int64_t morton, int32_t level, int32_t dim, uint32_t grid,
///        uint32_t local_grid_size, const uint32_t local_offset[3])
/// and is_purple_box without the dimension argument).  Kept as a template so
/// this header never includes tree_impl.hpp.
template <typename IsBlue, typename IsOrange, typename IsPurple>
ColorOf make_ca_classifier(const ProcessGrid& grid, int32_t level,
                           IsBlue is_blue, IsOrange is_orange, IsPurple is_purple) {
    return [=](int64_t morton) -> Color {
        if (!grid.on_brick_boundary(morton)) return Color::NONE;
        const auto off = grid.box_owner_offset(morton);
        const uint32_t o[3] = {off[0], off[1], off[2]};
        if (is_blue(morton, level, grid.dimension, grid.grid_size, grid.brick, o)) return Color::BLUE;
        if (is_orange(morton, level, grid.dimension, grid.grid_size, grid.brick, o)) return Color::BLUE;
        if (grid.dimension == 3 && is_purple(morton, level, grid.grid_size, grid.brick, o)) return Color::PURPLE;
        return Color::GREEN;
    };
}

// ===========================================================================
// Components and the graph
// ===========================================================================

struct Component {
    int32_t id = -1;                    ///< canonical id: sorted by (color, min morton)
    Color color = Color::NONE;
    std::vector<int64_t> boxes;         ///< sorted Morton indices
    std::vector<uint32_t> sharers;      ///< home process ids of the boxes (sorted)
    std::vector<int64_t> home_counts;   ///< boxes homed at each sharer (aligned with sharers; both constructions)
    std::vector<uint32_t> adjacent;     ///< pids holding a box within 1 hop of the component (sorted, ⊇ sharers)
    std::vector<uint32_t> visible;      ///< pids whose halo holds at least one box of the component (sorted, ⊇ adjacent)
    std::vector<int32_t> preds;         ///< lower-color components adjacent (1 hop) — sorted by id
    std::vector<int32_t> succs;         ///< higher-color components that list us as pred
    /// Lower-color components with a box within TWO hops but not adjacent
    /// (sorted). Their boxes appear as far endpoints in this component's
    /// sketches, so their skeletons must be known before it starts.
    std::vector<int32_t> preds2;
    std::vector<int32_t> succs2;
    double cost = 0.0;                  ///< scheduling cost (box count by default; measured later)
    double upward_rank = 0.0;           ///< HEFT priority: cost + max over succs (upward_rank)
    int32_t owner = -1;                 ///< pid; set by assign_owners
    /// The local construction (build_component_graph_local) attaches the box
    /// lists of the rank's neighbourhood only: `first_morton` is the
    /// component's smallest Morton index over ALL its boxes, `boxes_complete`
    /// whether `boxes` holds them all (always, in the global build).
    int64_t first_morton = std::numeric_limits<int64_t>::max();
    bool boxes_complete = true;

    int64_t min_morton() const {
        if (first_morton != std::numeric_limits<int64_t>::max()) return first_morton;
        return boxes.empty() ? first_morton : boxes.front();
    }
    bool has_sharer(uint32_t pid) const { return std::binary_search(sharers.begin(), sharers.end(), pid); }
    /// Boxes of the component homed at `pid` (0 for a non-sharer).
    int64_t home_boxes(uint32_t pid) const {
        const auto it = std::lower_bound(sharers.begin(), sharers.end(), pid);
        if (it == sharers.end() || *it != pid) return 0;
        return home_counts[static_cast<size_t>(it - sharers.begin())];
    }
    bool is_adjacent(uint32_t pid) const { return std::binary_search(adjacent.begin(), adjacent.end(), pid); }
    bool is_visible(uint32_t pid) const { return std::binary_search(visible.begin(), visible.end(), pid); }
};

/// Per-rank interior join node: the components whose generators the rank's
/// interior elimination (and its level-end state) depends on.
struct InteriorNode {
    uint32_t pid = 0;
    std::vector<int32_t> preds;         ///< sorted component ids
};

struct ComponentGraph {
    ProcessGrid grid;
    std::vector<Component> comps;                       ///< indexed by id
    std::unordered_map<int64_t, int32_t> box_to_comp;   ///< boundary box -> component id
    std::vector<InteriorNode> interior;                 ///< indexed by pid
    /// The local construction: the bricks whose boxes are in box_to_comp
    /// (sorted); empty for the global build (every brick).
    std::vector<uint32_t> local_bricks;

    const Component& comp_of_box(int64_t morton) const {
        auto it = box_to_comp.find(morton);
        if (it == box_to_comp.end())
            throw std::runtime_error("ComponentGraph: box " + std::to_string(morton) + " is not in any component");
        return comps[static_cast<size_t>(it->second)];
    }
    int32_t comp_id_of_box(int64_t morton) const {
        auto it = box_to_comp.find(morton);
        return it == box_to_comp.end() ? -1 : it->second;
    }

    std::vector<int32_t> components_of_color(Color c) const {
        std::vector<int32_t> out;
        for (const auto& k : comps) if (k.color == c) out.push_back(k.id);
        return out;
    }

    /// Components owned by / shared with / adjacent to a rank.
    std::vector<int32_t> owned_by(uint32_t pid) const {
        std::vector<int32_t> out;
        for (const auto& k : comps) if (k.owner == static_cast<int32_t>(pid)) out.push_back(k.id);
        return out;
    }
    std::vector<int32_t> shared_with(uint32_t pid) const {
        std::vector<int32_t> out;
        for (const auto& k : comps) if (k.has_sharer(pid)) out.push_back(k.id);
        return out;
    }
    std::vector<int32_t> adjacent_to(uint32_t pid) const {
        std::vector<int32_t> out;
        for (const auto& k : comps) if (k.is_adjacent(pid)) out.push_back(k.id);
        return out;
    }

    size_t num_boundary_boxes() const { return box_to_comp.size(); }
};

namespace detail {

inline int color_rank(Color c) { return static_cast<int>(c); }

/// Flood fill over the boundary boxes of *all* bricks, one class at a time,
/// 26-connectivity (8 in 2D).  Two boxes of the same class that are neighbors
/// belong to the same component; same-class components are therefore never
/// adjacent, which is the independence property the elimination relies on.
inline std::vector<Component> flood_components(const ProcessGrid& grid,
                                               const std::unordered_map<int64_t, Color>& cls) {
    std::vector<Component> out;
    std::unordered_map<int64_t, int32_t> seen;   // morton -> provisional component index
    seen.reserve(cls.size());

    // Deterministic seed order: increasing morton.
    std::vector<int64_t> keys;
    keys.reserve(cls.size());
    for (const auto& kv : cls) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());

    std::vector<int64_t> stack;
    for (int64_t seed : keys) {
        if (seen.count(seed)) continue;
        const Color c = cls.at(seed);
        Component comp;
        comp.color = c;
        stack.clear();
        stack.push_back(seed);
        seen[seed] = static_cast<int32_t>(out.size());
        while (!stack.empty()) {
            const int64_t m = stack.back();
            stack.pop_back();
            comp.boxes.push_back(m);
            for (uint64_t nb_u : grid.neighbors(m)) {
                const int64_t nb = static_cast<int64_t>(nb_u);
                auto it = cls.find(nb);
                if (it == cls.end() || it->second != c) continue;
                if (seen.count(nb)) continue;
                seen[nb] = static_cast<int32_t>(out.size());
                stack.push_back(nb);
            }
        }
        std::sort(comp.boxes.begin(), comp.boxes.end());
        out.push_back(std::move(comp));
    }
    return out;
}

}  // namespace detail

/// Build the component graph for one level.  Pure function of (grid, color_of);
/// every rank obtains the identical graph.
inline ComponentGraph build_component_graph(const ProcessGrid& grid, const ColorOf& color_of) {
    ComponentGraph g;
    g.grid = grid;

    // 1. Classify every boundary box of every brick from its home brick's view.
    std::unordered_map<int64_t, Color> cls;
    for (uint32_t pid = 0; pid < grid.num_procs; ++pid) {
        for (int64_t m : grid.brick_boundary_boxes(pid)) {
            const Color c = color_of(m);
            if (c == Color::NONE || c == Color::INTERIOR)
                throw std::runtime_error("build_component_graph: boundary box " + std::to_string(m) +
                                         " classified as non-boundary");
            if (c == Color::PURPLE && grid.dimension != 3)
                throw std::runtime_error("build_component_graph: purple box in 2D");
            cls.emplace(m, c);
        }
    }

    // 2. Connected components per class; canonical ids by (color, min morton).
    g.comps = detail::flood_components(grid, cls);
    std::sort(g.comps.begin(), g.comps.end(), [](const Component& a, const Component& b) {
        if (a.color != b.color) return detail::color_rank(a.color) < detail::color_rank(b.color);
        return a.min_morton() < b.min_morton();
    });
    for (size_t i = 0; i < g.comps.size(); ++i) {
        g.comps[i].id = static_cast<int32_t>(i);
        g.comps[i].first_morton = g.comps[i].boxes.front();
        for (int64_t m : g.comps[i].boxes) g.box_to_comp[m] = g.comps[i].id;
    }

    // 3. Sharers, adjacent ranks, predecessors.
    for (auto& k : g.comps) {
        std::set<uint32_t> sharers, adjacent;
        std::set<int32_t> preds;
        for (int64_t m : k.boxes) {
            sharers.insert(grid.proc_of_box(m));
            adjacent.insert(grid.proc_of_box(m));
            for (uint64_t nb_u : grid.neighbors(m)) {
                const int64_t nb = static_cast<int64_t>(nb_u);
                adjacent.insert(grid.proc_of_box(nb));
                auto it = g.box_to_comp.find(nb);
                if (it == g.box_to_comp.end()) continue;          // interior neighbor
                const Component& other = g.comps[static_cast<size_t>(it->second)];
                if (other.id == k.id) continue;
                if (detail::color_rank(other.color) < detail::color_rank(k.color)) preds.insert(other.id);
                else if (other.color == k.color)
                    throw std::runtime_error("build_component_graph: adjacent same-color components " +
                                             std::to_string(k.id) + " and " + std::to_string(other.id));
                // higher color: it will list us as its pred
            }
        }
        k.sharers.assign(sharers.begin(), sharers.end());
        k.home_counts.assign(k.sharers.size(), 0);
        for (int64_t m : k.boxes) {
            const auto it = std::lower_bound(k.sharers.begin(), k.sharers.end(), grid.proc_of_box(m));
            ++k.home_counts[static_cast<size_t>(it - k.sharers.begin())];
        }
        k.adjacent.assign(adjacent.begin(), adjacent.end());
        {
            // visible: every rank whose brick is within halo reach of a box.
            // Candidate ranks: bricks within one brick of the sharers.
            std::set<uint32_t> visible(adjacent.begin(), adjacent.end());
            for (uint32_t s : sharers) {
                const auto pc = grid.proc_coords(s);
                const int zr = grid.dimension == 3 ? 1 : 0;
                for (int dz = -zr; dz <= zr; ++dz)
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx) {
                            const int64_t px = static_cast<int64_t>(pc[0]) + dx, py = static_cast<int64_t>(pc[1]) + dy,
                                          pz = static_cast<int64_t>(pc[2]) + dz;
                            if (px < 0 || py < 0 || pz < 0 || px >= grid.procs_per_dim || py >= grid.procs_per_dim ||
                                pz >= grid.procs_per_dim) continue;
                            const uint32_t q = static_cast<uint32_t>(grid.dimension == 2
                                ? morton::encode_2d(static_cast<uint32_t>(px), static_cast<uint32_t>(py))
                                : morton::encode_3d(static_cast<uint32_t>(px), static_cast<uint32_t>(py), static_cast<uint32_t>(pz)));
                            if (visible.count(q)) continue;
                            for (int64_t m : k.boxes)
                                if (grid.holds_boundary_box(q, m)) { visible.insert(q); break; }
                        }
            }
            k.visible.assign(visible.begin(), visible.end());
        }
        k.preds.assign(preds.begin(), preds.end());
        k.cost = static_cast<double>(k.boxes.size());
    }
    for (const auto& k : g.comps)
        for (int32_t p : k.preds) g.comps[static_cast<size_t>(p)].succs.push_back(k.id);
    for (auto& k : g.comps) std::sort(k.succs.begin(), k.succs.end());
    // 2-hop lower-color predecessors (far endpoints of the sketch)
    for (auto& k : g.comps) {
        std::set<int32_t> p2;
        for (int64_t m : k.boxes)
            for (uint64_t nb_u : grid.neighbors_2hop(m)) {
                auto it = g.box_to_comp.find(static_cast<int64_t>(nb_u));
                if (it == g.box_to_comp.end()) continue;
                const Component& other = g.comps[static_cast<size_t>(it->second)];
                if (other.id == k.id) continue;
                if (detail::color_rank(other.color) >= detail::color_rank(k.color)) continue;
                if (std::binary_search(k.preds.begin(), k.preds.end(), other.id)) continue;
                p2.insert(other.id);
            }
        k.preds2.assign(p2.begin(), p2.end());
    }
    for (const auto& k : g.comps)
        for (int32_t p : k.preds2) g.comps[static_cast<size_t>(p)].succs2.push_back(k.id);
    for (auto& k : g.comps) std::sort(k.succs2.begin(), k.succs2.end());

    // 4. Interior join nodes: every component holding a box of the brick or a
    //    box within one hop of it (the second set is contained in the first
    //    for bricks >= 8, but computed directly so we never rely on that).
    g.interior.resize(grid.num_procs);
    for (uint32_t pid = 0; pid < grid.num_procs; ++pid) {
        std::set<int32_t> preds;
        for (int64_t m : grid.brick_boundary_boxes(pid)) {
            preds.insert(g.box_to_comp.at(m));
            for (uint64_t nb_u : grid.neighbors(m)) {
                auto it = g.box_to_comp.find(static_cast<int64_t>(nb_u));
                if (it != g.box_to_comp.end()) preds.insert(it->second);
            }
        }
        g.interior[pid].pid = pid;
        g.interior[pid].preds.assign(preds.begin(), preds.end());
    }

    // 5. HEFT upward ranks (transfer cost folded into a unit latency per hop).
    //    Components are already in (color, morton) order, so successors have
    //    larger ids: one reverse pass suffices.
    for (size_t i = g.comps.size(); i-- > 0;) {
        Component& k = g.comps[i];
        double best = 0.0;
        for (int32_t s : k.succs) best = std::max(best, g.comps[static_cast<size_t>(s)].upward_rank);
        k.upward_rank = k.cost + best;
    }
    return g;
}

// ===========================================================================
// Local construction: the graph from the rank's neighbourhood
// ===========================================================================
//
// build_component_graph() flood-fills the boundary boxes of every brick of
// the grid, on every rank: 0.36 s at 64 bricks and 3.6 s at 512 per level
// (h2_gpu/GPU_CA_PLAN.md, M5e).  Every component a rank uses lies in the 3^d
// bricks around it (a corner cluster spans the bricks of its corner, an edge
// run those of its edge, a face its two; the two-hop readers' boxes are
// within two boxes of the brick), so the rank builds the components' PARTS
// (a component restricted to one brick, connected within it) of its own
// brick and its neighbours, describes its own brick's parts (box count, the
// ranks adjacent and visible, the links to same-colour parts across the
// brick's faces, the lower-colour parts within one and two hops) in a
// stream of int64 that the caller allgathers over the level's ranks, and
// unites all ranks' parts into the global components: the same ids
// (colour, smallest Morton index), costs, sharers, adjacency, predecessors
// and HEFT ranks as the global build, hence the same owners from
// assign_owners().  Box lists are attached for the parts of the
// neighbourhood: complete for every component touching the rank's brick
// (the wave partition, the fetch plan, the payload recipients and the
// schedule's filters read those), partial for the others (whose boxes near
// the rank are what the halo and the sketches see; `boxes_complete`).
//
//   LocalParts lp = local_component_parts(grid, color_of, my_pid);
//   std::vector<int64_t> all = allgather(lp.stream);   // over the level's ranks, any order
//   ComponentGraph g = assemble_component_graph(grid, my_pid, lp, all);
//
// or build_component_graph_local(grid, color_of, my_pid, allgather) for the
// three steps.  No MPI here: the test concatenates every brick's stream.

struct LocalParts {
    struct Part {
        uint32_t brick = 0;
        Color color = Color::NONE;
        std::vector<int64_t> boxes;   ///< sorted; the part's key is (brick, boxes.front())
    };
    std::vector<uint32_t> bricks;                       ///< the neighbourhood, sorted
    std::vector<Part> parts;                            ///< of the neighbourhood's bricks
    std::vector<int64_t> stream;                        ///< the descriptors of my brick's parts
};

namespace detail {

/// The bricks within one brick of `pid` in the process grid, `pid` included, sorted.
inline std::vector<uint32_t> brick_neighbourhood(const ProcessGrid& grid, uint32_t pid) {
    std::vector<uint32_t> out;
    const auto c = grid.proc_coords(pid);
    const int zr = grid.dimension == 3 ? 1 : 0;
    const int64_t n = grid.procs_per_dim;
    for (int dz = -zr; dz <= zr; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
                const int64_t x = static_cast<int64_t>(c[0]) + dx, y = static_cast<int64_t>(c[1]) + dy,
                              z = static_cast<int64_t>(c[2]) + dz;
                if (x < 0 || y < 0 || z < 0 || x >= n || y >= n || z >= n) continue;
                out.push_back(static_cast<uint32_t>(
                    grid.dimension == 2 ? morton::encode_2d(static_cast<uint32_t>(x), static_cast<uint32_t>(y))
                                        : morton::encode_3d(static_cast<uint32_t>(x), static_cast<uint32_t>(y),
                                                            static_cast<uint32_t>(z))));
            }
    std::sort(out.begin(), out.end());
    return out;
}

/// A part's key in the streams: the brick and the smallest Morton index.
inline uint64_t part_key(uint32_t brick, int64_t min_morton) {
    return (static_cast<uint64_t>(brick) << 40) ^ static_cast<uint64_t>(min_morton);
}

}  // namespace detail

/// Step 1: the parts of the bricks around `my_pid` and the descriptors of my
/// brick's parts.  Dense arrays over the neighbourhood's cells (3^d bricks:
/// 110K cells for 16^3 bricks) instead of hash maps and neighbour vectors:
/// a few milliseconds.
inline LocalParts local_component_parts(const ProcessGrid& grid, const ColorOf& color_of, uint32_t my_pid) {
    LocalParts lp;
    lp.bricks = detail::brick_neighbourhood(grid, my_pid);
    const int dim = grid.dimension;
    const uint32_t B = grid.brick;
    const uint32_t G = grid.grid_size;
    // the neighbourhood's cells: [lo, lo + n) per dimension, in boxes
    const auto mc = grid.proc_coords(my_pid);
    uint32_t lo[3] = {0, 0, 0}, n[3] = {1, 1, 1};
    for (int d = 0; d < dim; ++d) {
        const uint32_t b0 = mc[static_cast<size_t>(d)] > 0 ? mc[static_cast<size_t>(d)] - 1 : 0;
        const uint32_t b1 = std::min(mc[static_cast<size_t>(d)] + 2, grid.procs_per_dim);
        lo[d] = b0 * B;
        n[d] = (b1 - b0) * B;
    }
    const size_t cells = static_cast<size_t>(n[0]) * n[1] * n[2];
    std::vector<uint8_t> cls(cells, 255);    // Color, 255: not a boundary box
    std::vector<int32_t> part(cells, -1);    // index in lp.parts
    auto cell = [&](uint32_t x, uint32_t y, uint32_t z) -> size_t {
        return (static_cast<size_t>(z - lo[2]) * n[1] + (y - lo[1])) * n[0] + (x - lo[0]);
    };
    auto inside = [&](int64_t x, int64_t y, int64_t z) {   // in the grid and in the neighbourhood's cells
        return x >= lo[0] && y >= lo[1] && z >= lo[2] && x < static_cast<int64_t>(lo[0] + n[0]) &&
               y < static_cast<int64_t>(lo[1] + n[1]) && z < static_cast<int64_t>(lo[2] + n[2]) && x < G && y < G &&
               (dim == 2 || z < G);
    };
    auto brick_of = [&](uint32_t x, uint32_t y, uint32_t z) -> uint32_t {
        return static_cast<uint32_t>(dim == 2 ? morton::encode_2d(x / B, y / B) : morton::encode_3d(x / B, y / B, z / B));
    };
    auto on_boundary = [&](uint32_t x, uint32_t y, uint32_t z) {
        const uint32_t lx = x % B, ly = y % B, lz = z % B;
        return lx == 0 || lx == B - 1 || ly == 0 || ly == B - 1 || (dim == 3 && (lz == 0 || lz == B - 1));
    };
    // classification of every brick's boundary boxes
    for (uint32_t q : lp.bricks) {
        const auto o = grid.proc_offset(q);
        const uint32_t zmax = dim == 3 ? B : 1;
        for (uint32_t z = 0; z < zmax; ++z)
            for (uint32_t y = 0; y < B; ++y)
                for (uint32_t x = 0; x < B; ++x) {
                    const uint32_t gx = o[0] + x, gy = o[1] + y, gz = dim == 3 ? o[2] + z : 0;
                    if (!on_boundary(gx, gy, gz)) continue;
                    const int64_t m = grid.encode(gx, gy, gz);
                    const Color c = color_of(m);
                    if (c == Color::NONE || c == Color::INTERIOR)
                        throw std::runtime_error("local_component_parts: boundary box " + std::to_string(m) +
                                                 " classified as non-boundary");
                    if (c == Color::PURPLE && dim != 3)
                        throw std::runtime_error("local_component_parts: purple box in 2D");
                    cls[cell(gx, gy, gz)] = static_cast<uint8_t>(c);
                }
    }
    // the parts: flood fill within each brick (seeds in Morton order, so the
    // part's first box is its smallest Morton index), 26-connectivity
    const int zr = dim == 3 ? 1 : 0;
    std::vector<size_t> stack;
    for (uint32_t q : lp.bricks) {
        const auto o = grid.proc_offset(q);
        const int64_t qlo[3] = {o[0], o[1], dim == 3 ? static_cast<int64_t>(o[2]) : 0};
        const int64_t qhi[3] = {qlo[0] + B, qlo[1] + B, dim == 3 ? qlo[2] + B : 1};
        for (int64_t seed : grid.brick_boundary_boxes(q)) {
            const auto sc = grid.decode(seed);
            const size_t sidx = cell(sc[0], sc[1], dim == 3 ? sc[2] : 0);
            if (part[sidx] >= 0) continue;
            const uint8_t c = cls[sidx];
            const int32_t index = static_cast<int32_t>(lp.parts.size());
            LocalParts::Part p;
            p.brick = q;
            p.color = static_cast<Color>(c);
            part[sidx] = index;
            stack.clear();
            stack.push_back(sidx);
            while (!stack.empty()) {
                const size_t i = stack.back();
                stack.pop_back();
                const uint32_t x = lo[0] + static_cast<uint32_t>(i % n[0]);
                const uint32_t y = lo[1] + static_cast<uint32_t>((i / n[0]) % n[1]);
                const uint32_t z = lo[2] + static_cast<uint32_t>(i / (static_cast<size_t>(n[0]) * n[1]));
                p.boxes.push_back(grid.encode(x, y, dim == 3 ? z : 0));
                for (int dz = -zr; dz <= zr; ++dz)
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx) {
                            if (dx == 0 && dy == 0 && dz == 0) continue;
                            const int64_t nx = static_cast<int64_t>(x) + dx, ny = static_cast<int64_t>(y) + dy,
                                          nz = static_cast<int64_t>(z) + dz;
                            if (nx < qlo[0] || ny < qlo[1] || nz < qlo[2] || nx >= qhi[0] || ny >= qhi[1] || nz >= qhi[2]) continue;
                            const size_t j = cell(static_cast<uint32_t>(nx), static_cast<uint32_t>(ny), static_cast<uint32_t>(nz));
                            if (cls[j] != c || part[j] >= 0) continue;
                            part[j] = index;
                            stack.push_back(j);
                        }
            }
            std::sort(p.boxes.begin(), p.boxes.end());
            lp.parts.push_back(std::move(p));
        }
    }
    // my brick's parts' descriptors
    auto key_of_cell = [&](size_t j) {
        const LocalParts::Part& p = lp.parts[static_cast<size_t>(part[j])];
        return std::make_pair(static_cast<int64_t>(p.brick), p.boxes.front());
    };
    auto put_keys = [&](const std::set<std::pair<int64_t, int64_t>>& keys) {
        lp.stream.push_back(static_cast<int64_t>(keys.size()));
        for (const auto& k : keys) {
            lp.stream.push_back(k.first);
            lp.stream.push_back(k.second);
        }
    };
    auto put_ranks = [&](const std::set<uint32_t>& ranks) {
        lp.stream.push_back(static_cast<int64_t>(ranks.size()));
        for (uint32_t r : ranks) lp.stream.push_back(static_cast<int64_t>(r));
    };
    const uint32_t reach = grid.halo_reach();
    for (const LocalParts::Part& p : lp.parts) {
        if (p.brick != my_pid) continue;
        std::set<uint32_t> adjacent{my_pid}, visible;
        std::set<std::pair<int64_t, int64_t>> links, preds, preds2;
        const int my_rank = detail::color_rank(p.color);
        std::vector<char> seen_visible(lp.bricks.size(), 0);
        for (int64_t m : p.boxes) {
            const auto c = grid.decode(m);
            const int64_t x = c[0], y = c[1], z = dim == 3 ? c[2] : 0;
            // the two-hop region: 1-hop neighbours (adjacency, links, preds)
            // and the 2-hop shell (preds2)
            const int zr2 = dim == 3 ? 2 : 0;
            for (int dz = -zr2; dz <= zr2; ++dz)
                for (int dy = -2; dy <= 2; ++dy)
                    for (int dx = -2; dx <= 2; ++dx) {
                        if (dx == 0 && dy == 0 && dz == 0) continue;
                        const int64_t nx = x + dx, ny = y + dy, nz = z + dz;
                        if (nx < 0 || ny < 0 || nz < 0 || !inside(nx, ny, nz)) continue;
                        const bool one_hop = std::abs(dx) <= 1 && std::abs(dy) <= 1 && std::abs(dz) <= 1;
                        const uint32_t ux = static_cast<uint32_t>(nx), uy = static_cast<uint32_t>(ny), uz = static_cast<uint32_t>(nz);
                        if (one_hop) adjacent.insert(brick_of(ux, uy, uz));
                        const size_t j = cell(ux, uy, uz);
                        if (cls[j] == 255) continue;   // an interior box
                        const int other = detail::color_rank(static_cast<Color>(cls[j]));
                        if (one_hop) {
                            if (other == my_rank) {
                                if (brick_of(ux, uy, uz) != my_pid) links.insert(key_of_cell(j));
                            } else if (other < my_rank) {
                                preds.insert(key_of_cell(j));
                            }
                        } else if (other < my_rank) {
                            preds2.insert(key_of_cell(j));
                        }
                    }
            // the ranks whose halo holds the box
            for (size_t b = 0; b < lp.bricks.size(); ++b) {
                if (seen_visible[b]) continue;
                const auto o = grid.proc_offset(lp.bricks[b]);
                uint32_t dist = 0;
                for (int d = 0; d < dim; ++d) {
                    const int64_t l = o[static_cast<size_t>(d)], h = l + B - 1, v = c[static_cast<size_t>(d)];
                    const int64_t dk = v < l ? l - v : (v > h ? v - h : 0);
                    dist = std::max(dist, static_cast<uint32_t>(dk));
                }
                if (dist <= reach) seen_visible[b] = 1;
            }
        }
        visible = adjacent;
        for (size_t b = 0; b < lp.bricks.size(); ++b)
            if (seen_visible[b]) visible.insert(lp.bricks[b]);
        lp.stream.push_back(static_cast<int64_t>(p.brick));
        lp.stream.push_back(p.boxes.front());
        lp.stream.push_back(static_cast<int64_t>(p.color));
        lp.stream.push_back(static_cast<int64_t>(p.boxes.size()));
        put_ranks(adjacent);
        put_ranks(visible);
        put_keys(links);
        put_keys(preds);
        put_keys(preds2);
    }
    return lp;
}

/// Step 3: the global components from every rank's stream (`all`: the
/// concatenation of the streams, in any order), with the neighbourhood's
/// box lists attached.
inline ComponentGraph assemble_component_graph(const ProcessGrid& grid, uint32_t my_pid, const LocalParts& lp,
                                               const std::vector<int64_t>& all) {
    struct PartDesc {
        uint32_t brick = 0;
        int64_t min = 0;
        Color color = Color::NONE;
        int64_t count = 0;
        std::vector<uint32_t> adjacent, visible;
        std::vector<uint64_t> links, preds, preds2;   // keys
    };
    std::vector<PartDesc> parts;
    std::unordered_map<uint64_t, size_t> index_of_key;
    {
        size_t at = 0;
        auto next = [&]() -> int64_t {
            if (at >= all.size()) throw std::runtime_error("assemble_component_graph: truncated part stream");
            return all[at++];
        };
        auto ranks = [&](std::vector<uint32_t>& v) {
            v.resize(static_cast<size_t>(next()));
            for (auto& r : v) r = static_cast<uint32_t>(next());
        };
        auto keys = [&](std::vector<uint64_t>& v) {
            v.resize(static_cast<size_t>(next()));
            for (auto& k : v) {
                const int64_t brick = next();
                const int64_t min = next();
                k = detail::part_key(static_cast<uint32_t>(brick), min);
            }
        };
        while (at < all.size()) {
            PartDesc d;
            d.brick = static_cast<uint32_t>(next());
            d.min = next();
            d.color = static_cast<Color>(next());
            d.count = next();
            ranks(d.adjacent);
            ranks(d.visible);
            keys(d.links);
            keys(d.preds);
            keys(d.preds2);
            const uint64_t key = detail::part_key(d.brick, d.min);
            if (!index_of_key.emplace(key, parts.size()).second)
                throw std::runtime_error("assemble_component_graph: part " + std::to_string(d.min) + " of brick " +
                                         std::to_string(d.brick) + " described twice");
            parts.push_back(std::move(d));
        }
    }
    auto index_of = [&](uint64_t key, const char* what) -> size_t {
        auto it = index_of_key.find(key);
        if (it == index_of_key.end())
            throw std::runtime_error(std::string("assemble_component_graph: unknown part in ") + what +
                                     " (a brick's stream is missing)");
        return it->second;
    };
    // union-find over the parts by the links
    std::vector<size_t> parent(parts.size());
    for (size_t i = 0; i < parts.size(); ++i) parent[i] = i;
    std::function<size_t(size_t)> find = [&](size_t i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    for (size_t i = 0; i < parts.size(); ++i)
        for (uint64_t l : parts[i].links) {
            const size_t j = index_of(l, "links");
            if (parts[j].color != parts[i].color)
                throw std::runtime_error("assemble_component_graph: a link between parts of different colours");
            const size_t a = find(i), b = find(j);
            if (a != b) parent[std::max(a, b)] = std::min(a, b);
        }
    // the components: parts by root, ordered by (colour, smallest Morton index)
    std::unordered_map<size_t, std::vector<size_t>> members;
    for (size_t i = 0; i < parts.size(); ++i) members[find(i)].push_back(i);
    struct Proto { size_t root; int color; int64_t min; };
    std::vector<Proto> protos;
    protos.reserve(members.size());
    for (const auto& kv : members) {
        int64_t mn = std::numeric_limits<int64_t>::max();
        for (size_t i : kv.second) mn = std::min(mn, parts[i].min);
        protos.push_back({kv.first, detail::color_rank(parts[kv.first].color), mn});
    }
    std::sort(protos.begin(), protos.end(), [](const Proto& a, const Proto& b) {
        if (a.color != b.color) return a.color < b.color;
        return a.min < b.min;
    });
    ComponentGraph g;
    g.grid = grid;
    g.local_bricks = lp.bricks;
    g.comps.resize(protos.size());
    std::vector<int32_t> comp_of_part(parts.size(), -1);
    for (size_t c = 0; c < protos.size(); ++c) {
        Component& k = g.comps[c];
        k.id = static_cast<int32_t>(c);
        k.color = parts[protos[c].root].color;
        k.first_morton = protos[c].min;
        std::set<uint32_t> sharers, adjacent, visible;
        std::map<uint32_t, int64_t> home;
        for (size_t i : members.at(protos[c].root)) {
            comp_of_part[i] = k.id;
            sharers.insert(parts[i].brick);
            home[parts[i].brick] += parts[i].count;
            adjacent.insert(parts[i].adjacent.begin(), parts[i].adjacent.end());
            visible.insert(parts[i].visible.begin(), parts[i].visible.end());
            k.cost += static_cast<double>(parts[i].count);
        }
        k.sharers.assign(sharers.begin(), sharers.end());
        for (uint32_t s : k.sharers) k.home_counts.push_back(home.at(s));
        k.adjacent.assign(adjacent.begin(), adjacent.end());
        k.visible.assign(visible.begin(), visible.end());
    }
    for (size_t c = 0; c < protos.size(); ++c) {
        Component& k = g.comps[c];
        std::set<int32_t> preds, preds2;
        for (size_t i : members.at(protos[c].root)) {
            for (uint64_t key : parts[i].preds) preds.insert(comp_of_part[index_of(key, "predecessors")]);
            for (uint64_t key : parts[i].preds2) preds2.insert(comp_of_part[index_of(key, "2-hop predecessors")]);
        }
        preds.erase(k.id);
        for (int32_t p : preds) {
            if (detail::color_rank(g.comps[static_cast<size_t>(p)].color) >= detail::color_rank(k.color))
                throw std::runtime_error("assemble_component_graph: a predecessor of no lower colour");
        }
        k.preds.assign(preds.begin(), preds.end());
        for (int32_t p : preds) preds2.erase(p);
        preds2.erase(k.id);
        k.preds2.assign(preds2.begin(), preds2.end());
    }
    for (const auto& k : g.comps) {
        for (int32_t p : k.preds) g.comps[static_cast<size_t>(p)].succs.push_back(k.id);
        for (int32_t p : k.preds2) g.comps[static_cast<size_t>(p)].succs2.push_back(k.id);
    }
    for (auto& k : g.comps) {
        std::sort(k.succs.begin(), k.succs.end());
        std::sort(k.succs2.begin(), k.succs2.end());
    }
    // the neighbourhood's boxes
    {
        size_t total = 0;
        for (const LocalParts::Part& p : lp.parts) total += p.boxes.size();
        g.box_to_comp.reserve(total);
    }
    for (const LocalParts::Part& p : lp.parts) {
        const int32_t c = comp_of_part[index_of(detail::part_key(p.brick, p.boxes.front()), "the neighbourhood")];
        Component& k = g.comps[static_cast<size_t>(c)];
        k.boxes.insert(k.boxes.end(), p.boxes.begin(), p.boxes.end());
        for (int64_t m : p.boxes) g.box_to_comp[m] = c;
    }
    for (auto& k : g.comps) {
        std::sort(k.boxes.begin(), k.boxes.end());
        k.boxes_complete = true;
        for (uint32_t s : k.sharers)
            if (!std::binary_search(lp.bricks.begin(), lp.bricks.end(), s)) { k.boxes_complete = false; break; }
        if (k.boxes_complete && static_cast<double>(k.boxes.size()) != k.cost)
            throw std::runtime_error("assemble_component_graph: a complete component's boxes do not add up");
    }
    // the rank's interior join node (the other ranks' stay empty)
    g.interior.resize(grid.num_procs);
    for (uint32_t pid = 0; pid < grid.num_procs; ++pid) g.interior[pid].pid = pid;
    {
        std::set<int32_t> preds;
        for (int64_t m : grid.brick_boundary_boxes(my_pid)) {
            preds.insert(g.box_to_comp.at(m));
            for (uint64_t nb_u : grid.neighbors(m)) {
                auto it = g.box_to_comp.find(static_cast<int64_t>(nb_u));
                if (it != g.box_to_comp.end()) preds.insert(it->second);
            }
        }
        g.interior[my_pid].preds.assign(preds.begin(), preds.end());
    }
    // HEFT upward ranks, as the global build
    for (size_t i = g.comps.size(); i-- > 0;) {
        Component& k = g.comps[i];
        double best = 0.0;
        for (int32_t s : k.succs) best = std::max(best, g.comps[static_cast<size_t>(s)].upward_rank);
        k.upward_rank = k.cost + best;
    }
    return g;
}

/// The three steps: `allgather` returns the concatenation of every rank's
/// stream (over the level's ranks, in any order).
inline ComponentGraph build_component_graph_local(
    const ProcessGrid& grid, const ColorOf& color_of, uint32_t my_pid,
    const std::function<std::vector<int64_t>(const std::vector<int64_t>&)>& allgather) {
    const LocalParts lp = local_component_parts(grid, color_of, my_pid);
    const std::vector<int64_t> all = allgather(lp.stream);
    return assemble_component_graph(grid, my_pid, lp, all);
}

// ===========================================================================
// Owner assignment
// ===========================================================================

enum class AssignPolicy {
    LOWEST_SHARER,   ///< debug: lowest pid among sharers (reproduces "home of min morton" style ownership)
    ROTATION,        ///< sharer[(id) % |sharers|] — spreads work, ignores cost
    LPT_PER_COLOR,   ///< default: per color, largest component first to the least-loaded sharer
    /// Blue as LPT_PER_COLOR (every rank gets its blue), then purple and green
    /// largest-first to the sharer with the smallest TOTAL load over the
    /// colours so far; relief and polish on the totals.  For the arrival-
    /// driven schedule: a rank heavy in green is made light in blue/purple and
    /// starts its green early, so the level ends near the mean total instead
    /// of the sum of the per-colour maxima (ledger 32/33).
    LPT_CUMULATIVE,
};

namespace detail {

/// Deterministic local search on one color's assignment: move a component to
/// another eligible sharer, or swap two components between their owners, when
/// that lowers the sum of squared loads (a strictly decreasing objective, so
/// the loop terminates; it minimises the makespan proxy while preferring
/// even loads).  Same input -> same output on every rank.
inline void polish_assignment_reference(ComponentGraph& g, const std::vector<int32_t>& ids, std::vector<double>& load) {
    auto sq = [](double x) { return x * x; };
    const uint32_t P = static_cast<uint32_t>(load.size());

    // ---- 1. makespan relief by augmenting paths ---------------------------
    // From a rank r at the maximum load, BFS over ranks: an edge u -> s moves
    // one component of u to another of its sharers s.  Intermediate ranks
    // receive one component and pass another on; every rank on the chain,
    // including the last, must end strictly below load[r].  Applying the chain
    // lowers (max load, number of ranks at max) lexicographically, so the loop
    // terminates.  Deterministic iteration order (ranks, components, sharers
    // ascending) -> identical result on every rank.
    std::vector<std::vector<int32_t>> on_rank(P);
    auto rebuild_on_rank = [&]() {
        for (auto& v : on_rank) v.clear();
        for (int32_t id : ids) on_rank[static_cast<size_t>(g.comps[static_cast<size_t>(id)].owner)].push_back(id);
    };
    auto relieve = [&](uint32_t r) -> bool {
        const double cap = load[r];
        std::vector<int32_t> parent_rank(P, -1), in_comp(P, -1);   // in_comp[s] = component moved INTO s
        std::vector<double> in_cost(P, 0.0);
        std::vector<char> visited(P, 0);
        std::deque<uint32_t> q;
        visited[r] = 1;
        q.push_back(r);
        while (!q.empty()) {
            const uint32_t u = q.front();
            q.pop_front();
            const double gained = (u == r) ? 0.0 : in_cost[u];
            for (int32_t id : on_rank[u]) {
                const Component& k = g.comps[static_cast<size_t>(id)];
                if (u != r && load[u] + gained - k.cost >= cap - 1e-9) continue;   // u would stay at the max
                for (uint32_t s : k.sharers) {
                    if (visited[s]) continue;
                    if (load[s] + k.cost < cap - 1e-9) {
                        // terminal: apply the chain s <- u <- ... <- r
                        int32_t moved = id;
                        uint32_t dst = s, src = u;
                        while (true) {
                            g.comps[static_cast<size_t>(moved)].owner = static_cast<int32_t>(dst);
                            load[dst] += g.comps[static_cast<size_t>(moved)].cost;
                            load[src] -= g.comps[static_cast<size_t>(moved)].cost;
                            if (src == r) break;
                            moved = in_comp[src];
                            dst = src;
                            src = static_cast<uint32_t>(parent_rank[src]);
                        }
                        return true;
                    }
                    visited[s] = 1;
                    parent_rank[s] = static_cast<int32_t>(u);
                    in_comp[s] = id;
                    in_cost[s] = k.cost;
                    q.push_back(s);
                }
            }
        }
        return false;
    };
    for (int iter = 0; iter < 100000; ++iter) {
        rebuild_on_rank();
        const double mx = *std::max_element(load.begin(), load.end());
        bool relieved = false;
        for (uint32_t r = 0; r < P && !relieved; ++r)
            if (load[r] >= mx - 1e-9 && relieve(r)) relieved = true;
        if (!relieved) break;
    }

    // ---- 2. evenness: sum-of-squares descent (never raises the max) --------
    for (int pass = 0; pass < 200; ++pass) {
        bool changed = false;
        // moves
        for (int32_t id : ids) {
            Component& k = g.comps[static_cast<size_t>(id)];
            const uint32_t r = static_cast<uint32_t>(k.owner);
            uint32_t best = r;
            double best_gain = 0.0;
            for (uint32_t s : k.sharers) {
                if (s == r) continue;
                // gain = old - new of sq(load[r]) + sq(load[s])
                const double gain = sq(load[r]) + sq(load[s]) - sq(load[r] - k.cost) - sq(load[s] + k.cost);
                if (gain > best_gain + 1e-9) { best_gain = gain; best = s; }
            }
            if (best != r) {
                load[r] -= k.cost;
                load[best] += k.cost;
                k.owner = static_cast<int32_t>(best);
                changed = true;
            }
        }
        // ejection chains of length 2: move k r->s and push some k' s->t
        for (int32_t id : ids) {
            Component& k = g.comps[static_cast<size_t>(id)];
            const uint32_t r = static_cast<uint32_t>(k.owner);
            bool moved = false;
            for (uint32_t s : k.sharers) {
                if (moved || s == r) continue;
                for (int32_t id2 : ids) {
                    Component& k2 = g.comps[static_cast<size_t>(id2)];
                    if (moved || k2.owner != static_cast<int32_t>(s) || id2 == id) continue;
                    for (uint32_t t : k2.sharers) {
                        if (t == s || t == r) continue;
                        const double gain = sq(load[r]) + sq(load[s]) + sq(load[t])
                                          - sq(load[r] - k.cost) - sq(load[s] + k.cost - k2.cost) - sq(load[t] + k2.cost);
                        if (gain > 1e-9) {
                            load[r] -= k.cost;
                            load[s] += k.cost - k2.cost;
                            load[t] += k2.cost;
                            k.owner = static_cast<int32_t>(s);
                            k2.owner = static_cast<int32_t>(t);
                            changed = moved = true;
                            break;
                        }
                    }
                }
            }
        }
        // swaps (a on r, b on s, a eligible on s and b eligible on r)
        for (size_t i = 0; i < ids.size(); ++i) {
            Component& a = g.comps[static_cast<size_t>(ids[i])];
            for (size_t j = i + 1; j < ids.size(); ++j) {
                Component& b = g.comps[static_cast<size_t>(ids[j])];
                if (a.owner == b.owner || a.cost == b.cost) continue;
                const uint32_t r = static_cast<uint32_t>(a.owner), s = static_cast<uint32_t>(b.owner);
                if (!a.has_sharer(s) || !b.has_sharer(r)) continue;
                const double d = b.cost - a.cost;   // r gains d, s loses d
                const double gain = sq(load[r]) + sq(load[s]) - sq(load[r] + d) - sq(load[s] - d);
                if (gain > 1e-9) {
                    load[r] += d;
                    load[s] -= d;
                    a.owner = static_cast<int32_t>(s);
                    b.owner = static_cast<int32_t>(r);
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
}


/// polish_assignment_reference with the components of every rank indexed
/// (in `ids` order, the order the reference scans) so that the search for
/// a partner on a rank runs over that rank's components instead of all
/// `ids`: the same moves and swaps in the same order, hence the same
/// assignment (owner_graph_test.cpp checks it); 4401 components on 512
/// ranks take milliseconds instead of ~90 ms.
inline void polish_assignment(ComponentGraph& g, const std::vector<int32_t>& ids, std::vector<double>& load) {
    auto sq = [](double x) { return x * x; };
    const uint32_t P = static_cast<uint32_t>(load.size());
    std::vector<int32_t> pos(g.comps.size(), -1);   // position in ids
    for (size_t i = 0; i < ids.size(); ++i) pos[static_cast<size_t>(ids[i])] = static_cast<int32_t>(i);
    std::vector<std::vector<int32_t>> on_rank(P);   // each rank's components, by position
    auto rebuild_on_rank = [&]() {
        for (auto& v : on_rank) v.clear();
        for (int32_t id : ids) on_rank[static_cast<size_t>(g.comps[static_cast<size_t>(id)].owner)].push_back(id);
    };
    auto set_owner = [&](int32_t id, uint32_t to) {
        Component& k = g.comps[static_cast<size_t>(id)];
        std::vector<int32_t>& from = on_rank[static_cast<size_t>(k.owner)];
        from.erase(std::find(from.begin(), from.end(), id));
        std::vector<int32_t>& dst = on_rank[to];
        dst.insert(std::upper_bound(dst.begin(), dst.end(), id,
                                    [&](int32_t x, int32_t y) { return pos[static_cast<size_t>(x)] < pos[static_cast<size_t>(y)]; }),
                   id);
        k.owner = static_cast<int32_t>(to);
    };

    // ---- 1. makespan relief by augmenting paths (as the reference) --------
    auto relieve = [&](uint32_t r) -> bool {
        const double cap = load[r];
        std::vector<int32_t> parent_rank(P, -1), in_comp(P, -1);
        std::vector<double> in_cost(P, 0.0);
        std::vector<char> visited(P, 0);
        std::deque<uint32_t> q;
        visited[r] = 1;
        q.push_back(r);
        while (!q.empty()) {
            const uint32_t u = q.front();
            q.pop_front();
            const double gained = (u == r) ? 0.0 : in_cost[u];
            for (int32_t id : on_rank[u]) {
                const Component& k = g.comps[static_cast<size_t>(id)];
                if (u != r && load[u] + gained - k.cost >= cap - 1e-9) continue;
                for (uint32_t s : k.sharers) {
                    if (visited[s]) continue;
                    if (load[s] + k.cost < cap - 1e-9) {
                        int32_t moved = id;
                        uint32_t dst = s, src = u;
                        while (true) {
                            g.comps[static_cast<size_t>(moved)].owner = static_cast<int32_t>(dst);
                            load[dst] += g.comps[static_cast<size_t>(moved)].cost;
                            load[src] -= g.comps[static_cast<size_t>(moved)].cost;
                            if (src == r) break;
                            moved = in_comp[src];
                            dst = src;
                            src = static_cast<uint32_t>(parent_rank[src]);
                        }
                        return true;
                    }
                    visited[s] = 1;
                    parent_rank[s] = static_cast<int32_t>(u);
                    in_comp[s] = id;
                    in_cost[s] = k.cost;
                    q.push_back(s);
                }
            }
        }
        return false;
    };
    for (int iter = 0; iter < 100000; ++iter) {
        rebuild_on_rank();
        const double mx = *std::max_element(load.begin(), load.end());
        bool relieved = false;
        for (uint32_t r = 0; r < P && !relieved; ++r)
            if (load[r] >= mx - 1e-9 && relieve(r)) relieved = true;
        if (!relieved) break;
    }
    rebuild_on_rank();

    // ---- 2. local search on the sum of squared loads (as the reference) ----
    for (int pass = 0; pass < 200; ++pass) {
        bool changed = false;
        // moves
        for (int32_t id : ids) {
            Component& k = g.comps[static_cast<size_t>(id)];
            const uint32_t r = static_cast<uint32_t>(k.owner);
            uint32_t best = r;
            double best_gain = 0.0;
            for (uint32_t s : k.sharers) {
                if (s == r) continue;
                const double gain = sq(load[r]) + sq(load[s]) - sq(load[r] - k.cost) - sq(load[s] + k.cost);
                if (gain > best_gain + 1e-9) { best_gain = gain; best = s; }
            }
            if (best != r) {
                load[r] -= k.cost;
                load[best] += k.cost;
                set_owner(id, best);
                changed = true;
            }
        }
        // swaps: k (on r) to a sharer s, a component k2 of s to a third rank t
        for (int32_t id : ids) {
            Component& k = g.comps[static_cast<size_t>(id)];
            const uint32_t r = static_cast<uint32_t>(k.owner);
            bool moved = false;
            for (uint32_t s : k.sharers) {
                if (moved || s == r) continue;
                const std::vector<int32_t> candidates = on_rank[s];   // (a copy: the swap edits the lists)
                for (int32_t id2 : candidates) {
                    if (moved || id2 == id) continue;
                    Component& k2 = g.comps[static_cast<size_t>(id2)];
                    for (uint32_t t : k2.sharers) {
                        if (t == s || t == r) continue;
                        const double gain = sq(load[r]) + sq(load[s]) + sq(load[t])
                                          - sq(load[r] - k.cost) - sq(load[s] + k.cost - k2.cost) - sq(load[t] + k2.cost);
                        if (gain > 1e-9) {
                            load[r] -= k.cost;
                            load[s] += k.cost - k2.cost;
                            load[t] += k2.cost;
                            set_owner(id, s);
                            set_owner(id2, t);
                            changed = moved = true;
                            break;
                        }
                    }
                }
            }
        }
        // exchanges of two components between their owners, later ones in
        // `ids` order (the candidates: the components of a's sharers, by
        // position; rebuilt after an exchange, which changes a's owner)
        for (size_t i = 0; i < ids.size(); ++i) {
            Component& a = g.comps[static_cast<size_t>(ids[i])];
            int32_t after = static_cast<int32_t>(i);
            while (true) {
                const uint32_t r = static_cast<uint32_t>(a.owner);
                std::vector<int32_t> candidates;
                for (uint32_t s : a.sharers) {
                    if (s == r) continue;
                    for (int32_t id : on_rank[s])
                        if (pos[static_cast<size_t>(id)] > after) candidates.push_back(id);
                }
                std::sort(candidates.begin(), candidates.end(),
                          [&](int32_t x, int32_t y) { return pos[static_cast<size_t>(x)] < pos[static_cast<size_t>(y)]; });
                bool exchanged = false;
                for (int32_t idb : candidates) {
                    Component& b = g.comps[static_cast<size_t>(idb)];
                    if (a.cost == b.cost) continue;
                    const uint32_t s = static_cast<uint32_t>(b.owner);
                    if (!b.has_sharer(r)) continue;
                    const double d = b.cost - a.cost;
                    const double gain = sq(load[r]) + sq(load[s]) - sq(load[r] + d) - sq(load[s] - d);
                    if (gain > 1e-9) {
                        load[r] += d;
                        load[s] -= d;
                        set_owner(ids[i], s);
                        set_owner(idb, r);
                        changed = true;
                        after = pos[static_cast<size_t>(idb)];
                        exchanged = true;
                        break;
                    }
                }
                if (!exchanged) break;
            }
        }
        if (!changed) break;
    }
}

}  // namespace detail

/// Assign one owner per component among its sharers.  Deterministic; must be
/// called identically on every rank.  `fixed_load[pid]` (optional) seeds the
/// per-rank load (e.g. the interior work or measured costs).
/// `reference_polish`: the unindexed local search (the test compares them).
inline void assign_owners(ComponentGraph& g, AssignPolicy policy = AssignPolicy::LPT_PER_COLOR,
                          const std::vector<double>* fixed_load = nullptr, bool reference_polish = false) {
    auto polish = [&](const std::vector<int32_t>& ids, std::vector<double>& load) {
        if (reference_polish) detail::polish_assignment_reference(g, ids, load);
        else detail::polish_assignment(g, ids, load);
    };
    const uint32_t P = g.grid.num_procs;
    if (policy == AssignPolicy::LOWEST_SHARER) {
        for (auto& k : g.comps) k.owner = static_cast<int32_t>(k.sharers.front());
        return;
    }
    if (policy == AssignPolicy::ROTATION) {
        for (auto& k : g.comps)
            k.owner = static_cast<int32_t>(k.sharers[static_cast<size_t>(k.id) % k.sharers.size()]);
        return;
    }
    // LPT per color: colors are sequential phases, so balance each phase's
    // makespan independently.  Single-sharer components (domain boundary) are
    // forced and go first as pre-load; the rest largest-first to the
    // least-loaded eligible sharer (ties -> smallest cumulative load over the
    // colors so far, then lowest pid); then a deterministic local search.
    // LPT_CUMULATIVE: the first colour as above, the later colours by the
    // total (ties -> most boxes at home on that sharer, then lowest pid).
    const bool cumulative = policy == AssignPolicy::LPT_CUMULATIVE;
    // (from the per-sharer counts, not the box list: on a local graph the
    // list of a component away from the rank is partial)
    auto home_boxes = [&](const Component& k, uint32_t s) { return static_cast<int>(k.home_boxes(s)); };
    std::vector<double> total(P, 0.0);
    if (fixed_load) for (uint32_t p = 0; p < P && p < fixed_load->size(); ++p) total[p] = (*fixed_load)[p];
    std::vector<int32_t> later_ids;   // cumulative: components balanced on the totals
    bool first_color = true;
    for (Color c : boundary_colors(g.grid.dimension)) {
        std::vector<double> load(P, 0.0);
        std::vector<int32_t> ids = g.components_of_color(c);
        const bool by_total = cumulative && !first_color;
        first_color = false;
        std::stable_sort(ids.begin(), ids.end(), [&](int32_t a, int32_t b) {
            const Component& A = g.comps[static_cast<size_t>(a)];
            const Component& B = g.comps[static_cast<size_t>(b)];
            const bool fa = A.sharers.size() == 1, fb = B.sharers.size() == 1;
            if (fa != fb) return fa;                       // forced components first (pre-load)
            if (A.cost != B.cost) return A.cost > B.cost;  // then LPT: largest first
            if (A.sharers.size() != B.sharers.size()) return A.sharers.size() < B.sharers.size();
            return A.id < B.id;
        });
        for (int32_t id : ids) {
            Component& k = g.comps[static_cast<size_t>(id)];
            uint32_t best = k.sharers.front();
            if (by_total) {
                // least total so far; ties: keep the component nearer home
                // (fewer boxes shipped), then lowest pid
                int best_home = home_boxes(k, best);
                for (uint32_t s : k.sharers) {
                    const int h = home_boxes(k, s);
                    if (total[s] < total[best] ||
                        (total[s] == total[best] && (h > best_home || (h == best_home && s < best)))) {
                        best = s;
                        best_home = h;
                    }
                }
                k.owner = static_cast<int32_t>(best);
                total[best] += k.cost;
                later_ids.push_back(id);
                continue;
            }
            for (uint32_t s : k.sharers) {
                if (load[s] < load[best] ||
                    (load[s] == load[best] && (total[s] < total[best] ||
                                               (total[s] == total[best] && s < best))))
                    best = s;
            }
            k.owner = static_cast<int32_t>(best);
            load[best] += k.cost;
        }
        if (by_total) continue;
        polish(ids, load);
        for (int32_t id : ids) total[static_cast<size_t>(g.comps[static_cast<size_t>(id)].owner)] += g.comps[static_cast<size_t>(id)].cost;
    }
    // cumulative: relief and polish of the later colours on the totals (blue
    // stays as assigned: it is the pipeline's head and every rank keeps its share)
    if (cumulative && !later_ids.empty()) polish(later_ids, total);
}


/// Cheap consistency hash of the assignment (fold into an Allreduce at level
/// start so a divergent view fails loudly instead of deadlocking).
inline uint64_t assignment_hash(const ComponentGraph& g) {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&](uint64_t v) { h ^= v; h *= 1099511628211ull; };
    mix(g.comps.size());
    for (const auto& k : g.comps) {
        mix(static_cast<uint64_t>(k.id));
        mix(static_cast<uint64_t>(k.owner) + 7);
        mix(static_cast<uint64_t>(k.cost));
        mix(static_cast<uint64_t>(k.min_morton()));
    }
    return h;
}

/// Per-rank load summary for one color (box counts), for the phase report.
struct LoadSummary { double min = 0, max = 0, mean = 0; uint32_t argmax = 0; };
inline LoadSummary owned_load(const ComponentGraph& g, Color c) {
    std::vector<double> load(g.grid.num_procs, 0.0);
    for (const auto& k : g.comps) if (k.color == c && k.owner >= 0) load[static_cast<size_t>(k.owner)] += k.cost;
    LoadSummary s;
    s.min = *std::min_element(load.begin(), load.end());
    s.argmax = static_cast<uint32_t>(std::max_element(load.begin(), load.end()) - load.begin());
    s.max = load[s.argmax];
    s.mean = std::accumulate(load.begin(), load.end(), 0.0) / static_cast<double>(load.size());
    return s;
}
/// Total owned load per rank over all colours.
inline LoadSummary owned_load_total(const ComponentGraph& g) {
    std::vector<double> load(g.grid.num_procs, 0.0);
    for (const auto& k : g.comps) if (k.owner >= 0) load[static_cast<size_t>(k.owner)] += k.cost;
    LoadSummary s;
    s.min = *std::min_element(load.begin(), load.end());
    s.argmax = static_cast<uint32_t>(std::max_element(load.begin(), load.end()) - load.begin());
    s.max = load[s.argmax];
    s.mean = std::accumulate(load.begin(), load.end(), 0.0) / static_cast<double>(load.size());
    return s;
}

// ===========================================================================
// Fetch / serve plans (entry-state exchange, per owned component)
// ===========================================================================

/// What one rank must fetch: for each owned component, its boxes whose home
/// is another rank, grouped by home pid.  Mirror: what it must serve.
struct FetchPlan {
    /// fetch[color][home_pid] = sorted mortons this rank must receive from home_pid
    std::map<Color, std::map<uint32_t, std::vector<int64_t>>> fetch;
    /// serve[color][owner_pid] = sorted mortons (mine) this rank must send to owner_pid
    std::map<Color, std::map<uint32_t, std::vector<int64_t>>> serve;
    /// my boxes eliminated by someone else (factor section expected back): morton -> owner
    std::map<int64_t, uint32_t> foreign_eliminated;
};

inline FetchPlan build_fetch_plan(const ComponentGraph& g, uint32_t my_pid) {
    FetchPlan plan;
    for (const auto& k : g.comps) {
        if (k.owner < 0) throw std::runtime_error("build_fetch_plan: unassigned component");
        const uint32_t owner = static_cast<uint32_t>(k.owner);
        if (owner == my_pid) {
            for (int64_t m : k.boxes) {
                const uint32_t home = g.grid.proc_of_box(m);
                if (home != my_pid) plan.fetch[k.color][home].push_back(m);
            }
        } else if (k.has_sharer(my_pid)) {
            for (int64_t m : k.boxes) {
                if (g.grid.proc_of_box(m) == my_pid) {
                    plan.serve[k.color][owner].push_back(m);
                    plan.foreign_eliminated[m] = owner;
                }
            }
        }
    }
    return plan;
}

// ===========================================================================
// Runtime: counters, ready queue, expected arrivals, stamps
// ===========================================================================

/// Kind of payload a rank expects from another rank for one component.
enum class Arrival : uint8_t {
    ENTRY_STATE,   ///< remote boxes of an owned component (from their home)
    GENERATORS,    ///< generator section of a component eliminated elsewhere (from its owner)
    FACTORS,       ///< factor section for my boxes eliminated elsewhere (from the owner)
    SKELETON,      ///< auxiliary skeleton sections of a non-shared 2-hop predecessor (from its owner)
};

inline const char* arrival_name(Arrival a) {
    switch (a) {
        case Arrival::ENTRY_STATE: return "entry-state";
        case Arrival::GENERATORS: return "generators";
        case Arrival::FACTORS: return "factors";
        case Arrival::SKELETON: return "skeleton";
    }
    return "?";
}

struct ComponentStamps {
    double ready = -1, start = -1, finish = -1, sent = -1;
    double wait_after_ready = 0;   ///< ready -> start (queueing)
};

/// One rank's execution state for one level.  Touched only by the master
/// thread between team-wide units (see design §9).
/// Transport variants the runtime can be configured for.
struct RuntimeOptions {
    /// Entry state of owned components is already present (today's halo
    /// holds every box of every shared component): no ENTRY_STATE arrivals.
    bool entry_state_prefetched = false;
    /// One payload per component carries generators and factors together:
    /// no separate FACTORS arrivals.
    bool factors_with_generators = false;
    /// Expect generator payloads only for components this rank *shares*
    /// (has a local box in), not for every adjacent component.
    bool generators_for_shared_only = false;
    /// Owned components also wait for their 2-hop lower-color predecessors
    /// (`preds2`) — needed when sketches read far endpoints' skeletons.
    bool two_hop_dependencies = false;
    /// Non-shared components whose skeleton sections arrive as auxiliary
    /// payloads (see auxiliary_components()).
    std::vector<int32_t> aux_comps;
};

class Runtime {
public:
    Runtime() = default;
    Runtime(const ComponentGraph* graph, uint32_t my_pid, RuntimeOptions opts = {}) { reset(graph, my_pid, opts); }

    void reset(const ComponentGraph* graph, uint32_t my_pid, RuntimeOptions opts = {}) {
        g_ = graph;
        me_ = my_pid;
        opts_ = opts;
        const size_t n = g_->comps.size();
        available_.assign(n, 0);
        arrived_.assign(n, 0);
        aux_.assign(n, 0);
        entry_missing_.assign(n, 0);
        pred_missing_.assign(n, 0);
        pred2_missing_.assign(n, 0);
        done_.assign(n, 0);
        stamps_.assign(n, ComponentStamps{});
        ready_.clear();
        installable_.clear();
        owned_.clear();
        shared_not_owned_.clear();
        expected_.clear();
        finished_owned_ = 0;

        for (const auto& k : g_->comps) {
            const size_t i = static_cast<size_t>(k.id);
            if (k.owner == static_cast<int32_t>(me_)) {
                owned_.push_back(k.id);
                pred_missing_[i] = static_cast<int32_t>(k.preds.size());
                if (opts_.two_hop_dependencies) pred2_missing_[i] = static_cast<int32_t>(k.preds2.size());
                if (!opts_.entry_state_prefetched) {
                    // entry state: one arrival per remote home rank
                    std::set<uint32_t> homes;
                    for (int64_t m : k.boxes) {
                        const uint32_t h = g_->grid.proc_of_box(m);
                        if (h != me_) homes.insert(h);
                    }
                    entry_missing_[i] = static_cast<int32_t>(homes.size());
                    for (uint32_t h : homes) expected_.push_back({k.id, Arrival::ENTRY_STATE, h});
                }
                if (pred_missing_[i] == 0 && pred2_missing_[i] == 0 && entry_missing_[i] == 0) push_ready(k.id);
            } else if (opts_.generators_for_shared_only ? k.has_sharer(me_) : k.is_adjacent(me_)) {
                expected_.push_back({k.id, Arrival::GENERATORS, static_cast<uint32_t>(k.owner)});
                pred_missing_[i] = static_cast<int32_t>(k.preds.size());
                if (k.has_sharer(me_)) shared_not_owned_.push_back(k.id);
                if (!opts_.factors_with_generators && k.has_sharer(me_)) {
                    bool mine = false;
                    for (int64_t m : k.boxes) if (g_->grid.proc_of_box(m) == me_) { mine = true; break; }
                    if (mine) expected_.push_back({k.id, Arrival::FACTORS, static_cast<uint32_t>(k.owner)});
                }
            }
        }
        for (int32_t id : opts_.aux_comps) {
            const Component& k = g_->comps[static_cast<size_t>(id)];
            if (k.owner == static_cast<int32_t>(me_) || k.has_sharer(me_))
                throw std::runtime_error("dataflow: auxiliary component is owned or shared");
            aux_[static_cast<size_t>(id)] = 1;
            expected_.push_back({id, Arrival::SKELETON, static_cast<uint32_t>(k.owner)});
        }
        const auto& in = g_->interior[me_];
        interior_missing_ = static_cast<int32_t>(in.preds.size());
        std::sort(expected_.begin(), expected_.end());
        expected_consumed_.assign(expected_.size(), 0);
    }

    // ---- queries ----------------------------------------------------------
    uint32_t my_pid() const { return me_; }
    const ComponentGraph& graph() const { return *g_; }
    const std::vector<int32_t>& owned() const { return owned_; }
    /// Non-owned components this rank shares (their payloads are installed here).
    const std::vector<int32_t>& shared_not_owned() const { return shared_not_owned_; }
    /// Whether component `id`'s payload has arrived (not necessarily installed).
    bool has_arrived(int32_t id) const { return arrived_[static_cast<size_t>(id)] != 0; }
    bool has_installable() const { return !installable_.empty(); }
    /// Lowest-id arrived non-owned component whose predecessors are all
    /// available here (so its consumption replay is deterministic), or -1.
    int32_t pop_installable() {
        if (installable_.empty()) return -1;
        const int32_t id = *installable_.begin();
        installable_.erase(installable_.begin());
        return id;
    }
    bool has_ready() const { return !ready_.empty(); }
    bool all_owned_done() const { return finished_owned_ == owned_.size(); }
    bool interior_ready() const { return interior_missing_ == 0; }
    int32_t interior_missing() const { return interior_missing_; }
    bool is_done(int32_t id) const { return done_[static_cast<size_t>(id)] != 0; }
    bool is_available(int32_t id) const { return available_[static_cast<size_t>(id)] != 0; }
    const ComponentStamps& stamps(int32_t id) const { return stamps_[static_cast<size_t>(id)]; }
    ComponentStamps& stamps(int32_t id) { return stamps_[static_cast<size_t>(id)]; }

    struct Expected {
        int32_t comp; Arrival kind; uint32_t from;
        bool operator<(const Expected& o) const {
            if (comp != o.comp) return comp < o.comp;
            if (kind != o.kind) return kind < o.kind;
            return from < o.from;
        }
        bool operator==(const Expected& o) const { return comp == o.comp && kind == o.kind && from == o.from; }
    };
    const std::vector<Expected>& expected() const { return expected_; }

    /// True once every expected arrival has been consumed and every owned
    /// component is done: the level's boundary phase is over *locally*.
    bool boundary_complete() const {
        if (!all_owned_done()) return false;
        for (char c : expected_consumed_) if (!c) return false;
        return true;
    }

    // ---- transitions (master thread only) ----------------------------------

    /// Highest-priority ready owned component, or -1.  Priority: upward rank
    /// descending, then id ascending (deterministic).
    int32_t pop_ready(double now = 0.0) {
        if (ready_.empty()) return -1;
        auto it = ready_.begin();
        const int32_t id = it->second;
        ready_.erase(it);
        auto& s = stamps_[static_cast<size_t>(id)];
        s.start = now;
        if (s.ready >= 0) s.wait_after_ready = now - s.ready;
        return id;
    }
    int32_t peek_ready() const { return ready_.empty() ? -1 : ready_.begin()->second; }
    /// Stamp every currently-ready component with the caller's clock (reset()
    /// stamps them at 0, which is meaningless against a wall clock).
    void restamp_ready(double now) {
        for (const auto& e : ready_) stamps_[static_cast<size_t>(e.second)].ready = now;
    }
    /// Put a popped-but-not-started component back (serial oracle only).
    void requeue(int32_t id) {
        const Component& k = g_->comps[static_cast<size_t>(id)];
        ready_.insert({-k.upward_rank, id});
    }
    void requeue_installable(int32_t id) { installable_.insert(id); }

    /// An owned component's elimination finished on this rank.  Its generators
    /// are now available locally (for successors owned here and the interior).
    void on_owned_finished(int32_t id, double now = 0.0) {
        check_owned(id, "on_owned_finished");
        if (done_[static_cast<size_t>(id)]) throw std::runtime_error("dataflow: component finished twice");
        done_[static_cast<size_t>(id)] = 1;
        ++finished_owned_;
        stamps_[static_cast<size_t>(id)].finish = now;
        mark_available(id, now);
    }

    /// An entry-state payload for owned component `id` arrived from `from`.
    void on_entry_state(int32_t id, uint32_t from, double now = 0.0) {
        check_owned(id, "on_entry_state");
        consume({id, Arrival::ENTRY_STATE, from});
        if (--entry_missing_[static_cast<size_t>(id)] == 0 && pred_missing_[static_cast<size_t>(id)] == 0 &&
            pred2_missing_[static_cast<size_t>(id)] == 0)
            push_ready(id, now);
    }

    /// The generator section of component `id` (eliminated by `from`) arrived.
    /// In the two-phase form (payload_arrived / on_installed) availability is
    /// deferred until the receiver has replayed the component; the one-shot
    /// form keeps the original behaviour.
    void on_generators(int32_t id, uint32_t from, double now = 0.0) {
        consume({id, Arrival::GENERATORS, from});
        mark_available(id, now);
    }
    /// Two-phase: payload of non-owned component `id` arrived from `from`.
    void on_payload_arrived(int32_t id, uint32_t from) {
        consume({id, Arrival::GENERATORS, from});
        arrived_[static_cast<size_t>(id)] = 1;
        if (pred_missing_[static_cast<size_t>(id)] == 0) installable_.insert(id);
    }
    /// Two-phase: the receiver finished replaying component `id` (its
    /// consumption event); its generators now count as available here.
    void on_installed(int32_t id, double now = 0.0) {
        if (!arrived_[static_cast<size_t>(id)]) throw std::runtime_error("dataflow: installed before arrival");
        if (done_[static_cast<size_t>(id)]) throw std::runtime_error("dataflow: component installed twice");
        done_[static_cast<size_t>(id)] = 1;
        mark_available(id, now);
    }

    /// The factor section for my boxes of component `id` arrived from `from`.
    void on_factors(int32_t id, uint32_t from) { consume({id, Arrival::FACTORS, from}); }
    /// Auxiliary skeleton sections of non-shared component `id` arrived: its
    /// skeletons are known here, so owned 2-hop successors may proceed.
    void on_aux_arrived(int32_t id, uint32_t from, double now = 0.0) {
        consume({id, Arrival::SKELETON, from});
        arrived_[static_cast<size_t>(id)] = 1;
        mark_available(id, now);
    }
    bool is_aux(int32_t id) const { return aux_[static_cast<size_t>(id)] != 0; }

    /// Owned components in the debug (serial oracle) order: component id.
    std::vector<int32_t> serial_order() const { return owned_; }

    /// Passive diagnostic: what this rank is still waiting for.
    void dump_pending(FILE* f) const {
        std::fprintf(f, "[dataflow] pid %u: owned %zu finished %zu, interior missing %d, ready %zu\n",
                     me_, owned_.size(), finished_owned_, interior_missing_, ready_.size());
        for (size_t i = 0; i < expected_.size(); ++i) {
            if (expected_consumed_[i]) continue;
            const auto& e = expected_[i];
            std::fprintf(f, "  waiting: comp %d (%s) %s from pid %u\n", e.comp,
                         color_name(g_->comps[static_cast<size_t>(e.comp)].color), arrival_name(e.kind), e.from);
        }
        for (int32_t id : owned_) {
            const size_t i = static_cast<size_t>(id);
            if (done_[i]) continue;
            std::fprintf(f, "  owned comp %d (%s): entry missing %d, preds missing %d, 2-hop preds missing %d\n", id,
                         color_name(g_->comps[i].color), entry_missing_[i], pred_missing_[i], pred2_missing_[i]);
        }
    }

    /// Realized critical-path summary for the phase report.
    struct Summary { double max_wait_after_ready = 0; double sum_wait_after_ready = 0; int32_t n = 0; };
    Summary summary() const {
        Summary s;
        for (int32_t id : owned_) {
            const auto& st = stamps_[static_cast<size_t>(id)];
            if (st.start < 0) continue;
            s.max_wait_after_ready = std::max(s.max_wait_after_ready, st.wait_after_ready);
            s.sum_wait_after_ready += st.wait_after_ready;
            ++s.n;
        }
        return s;
    }

private:
    void check_owned(int32_t id, const char* where) const {
        if (id < 0 || static_cast<size_t>(id) >= g_->comps.size() ||
            g_->comps[static_cast<size_t>(id)].owner != static_cast<int32_t>(me_))
            throw std::runtime_error(std::string("dataflow: ") + where + ": component " +
                                     std::to_string(id) + " is not owned by pid " + std::to_string(me_));
    }
    void push_ready(int32_t id, double now = 0.0) {
        const Component& k = g_->comps[static_cast<size_t>(id)];
        ready_.insert({-k.upward_rank, id});   // set orders ascending: most negative rank first
        stamps_[static_cast<size_t>(id)].ready = now;
    }
    void consume(const Expected& e) {
        auto it = std::lower_bound(expected_.begin(), expected_.end(), e);
        if (it == expected_.end() || !(*it == e))
            throw std::runtime_error("dataflow: unexpected arrival: comp " + std::to_string(e.comp) + " " +
                                     arrival_name(e.kind) + " from pid " + std::to_string(e.from) +
                                     " on pid " + std::to_string(me_));
        const size_t i = static_cast<size_t>(it - expected_.begin());
        if (expected_consumed_[i])
            throw std::runtime_error("dataflow: duplicate arrival: comp " + std::to_string(e.comp) + " " +
                                     arrival_name(e.kind) + " from pid " + std::to_string(e.from));
        expected_consumed_[i] = 1;
    }
    /// Generators of `id` are available on this rank: decrement successors
    /// owned here and the interior counter.
    void mark_available(int32_t id, double now) {
        const size_t i = static_cast<size_t>(id);
        if (available_[i]) throw std::runtime_error("dataflow: component available twice");
        available_[i] = 1;
        const Component& k = g_->comps[i];
        for (int32_t s : k.succs) {
            const Component& S = g_->comps[static_cast<size_t>(s)];
            const size_t si = static_cast<size_t>(s);
            if (S.owner == static_cast<int32_t>(me_)) {
                if (--pred_missing_[si] == 0 && entry_missing_[si] == 0 && pred2_missing_[si] == 0) push_ready(s, now);
            } else if (pred_missing_[si] > 0) {
                // non-owned component tracked here (expected payload): it
                // becomes installable once its payload is in and preds are done
                if (--pred_missing_[si] == 0 && arrived_[si] && !done_[si]) installable_.insert(s);
            }
        }
        if (opts_.two_hop_dependencies)
            for (int32_t s : k.succs2) {
                const Component& S = g_->comps[static_cast<size_t>(s)];
                const size_t si = static_cast<size_t>(s);
                if (S.owner != static_cast<int32_t>(me_)) continue;
                if (--pred2_missing_[si] == 0 && pred_missing_[si] == 0 && entry_missing_[si] == 0) push_ready(s, now);
            }
        const auto& in = g_->interior[me_];
        if (std::binary_search(in.preds.begin(), in.preds.end(), id)) --interior_missing_;
    }

    const ComponentGraph* g_ = nullptr;
    uint32_t me_ = 0;
    RuntimeOptions opts_;
    std::vector<int32_t> owned_, shared_not_owned_;
    std::vector<char> available_, done_, arrived_, aux_;
    std::set<int32_t> installable_;
    std::vector<int32_t> entry_missing_, pred_missing_, pred2_missing_;
    std::vector<ComponentStamps> stamps_;
    std::set<std::pair<double, int32_t>> ready_;
    std::vector<Expected> expected_;
    std::vector<char> expected_consumed_;
    size_t finished_owned_ = 0;
    int32_t interior_missing_ = 0;
};

// ===========================================================================
// Validation against the existing CA color lists
// ===========================================================================

/// Check that the graph's classification agrees with the lists the CA level
/// carries (`blue`, `orange`, `purple`, `green` — local and ghost boxes) and
/// that every listed box is in a component.  Returns an empty string on
/// success, otherwise a description of the first mismatch.
inline std::string validate_against_lists(const ComponentGraph& g,
                                          const std::vector<int64_t>& blue,
                                          const std::vector<int64_t>& orange,
                                          const std::vector<int64_t>& purple,
                                          const std::vector<int64_t>& green) {
    auto check = [&](const std::vector<int64_t>& list, Color want, const char* name) -> std::string {
        for (int64_t m : list) {
            const int32_t id = g.comp_id_of_box(m);
            if (id < 0) return std::string("box ") + std::to_string(m) + " in list '" + name + "' has no component";
            const Color got = g.comps[static_cast<size_t>(id)].color;
            if (got != want)
                return std::string("box ") + std::to_string(m) + " in list '" + name + "' is " + color_name(got) +
                       " in the graph";
        }
        return {};
    };
    std::string e;
    if (!(e = check(blue, Color::BLUE, "blue")).empty()) return e;
    if (!(e = check(orange, Color::BLUE, "orange")).empty()) return e;
    if (!(e = check(purple, Color::PURPLE, "purple")).empty()) return e;
    if (!(e = check(green, Color::GREEN, "green")).empty()) return e;
    return {};
}

/// Structural invariants of a graph (used by the unit test and, under an env
/// flag, by the driver).  Returns an empty string on success.
inline std::string validate_graph(const ComponentGraph& g) {
    const ProcessGrid& grid = g.grid;
    // every boundary box of every brick (a local graph: of the neighbourhood's
    // bricks) is in exactly one component
    const bool local = !g.local_bricks.empty();
    size_t count = 0;
    std::vector<uint32_t> bricks = g.local_bricks;
    if (!local) for (uint32_t pid = 0; pid < grid.num_procs; ++pid) bricks.push_back(pid);
    for (uint32_t pid : bricks)
        for (int64_t m : grid.brick_boundary_boxes(pid)) {
            ++count;
            if (g.comp_id_of_box(m) < 0) return "boundary box " + std::to_string(m) + " has no component";
        }
    if (count != g.box_to_comp.size()) return "box_to_comp has extra entries";
    for (const auto& k : g.comps) {
        if (k.boxes.empty() && (!local || k.boxes_complete)) return "empty component " + std::to_string(k.id);
        if (!std::is_sorted(k.boxes.begin(), k.boxes.end())) return "unsorted boxes in " + std::to_string(k.id);
        if (k.boxes_complete && static_cast<double>(k.boxes.size()) != k.cost) return "cost not the box count in " + std::to_string(k.id);
        if (!k.boxes.empty() && k.boxes.front() < k.first_morton) return "first Morton index not the smallest in " + std::to_string(k.id);
        if (k.boxes_complete && !k.boxes.empty() && k.boxes.front() != k.first_morton) return "first Morton index wrong in " + std::to_string(k.id);
        if (k.sharers.empty()) return "component without sharers " + std::to_string(k.id);
        for (uint32_t s : k.sharers) if (!k.is_adjacent(s)) return "sharer not adjacent in " + std::to_string(k.id);
        for (uint32_t s : k.adjacent) if (!k.is_visible(s)) return "adjacent rank not visible in " + std::to_string(k.id);
        // preds: lower color and adjacent; no adjacent same-color component
        // (a local graph: among the boxes it holds)
        for (int64_t m : k.boxes)
            for (uint64_t nb_u : grid.neighbors(m)) {
                const int32_t o = g.comp_id_of_box(static_cast<int64_t>(nb_u));
                if (o < 0 || o == k.id) continue;
                const Component& O = g.comps[static_cast<size_t>(o)];
                if (O.color == k.color) return "adjacent same-color components " + std::to_string(k.id) + "," + std::to_string(o);
                if (detail::color_rank(O.color) < detail::color_rank(k.color) &&
                    !std::binary_search(k.preds.begin(), k.preds.end(), o))
                    return "missing pred " + std::to_string(o) + " of " + std::to_string(k.id);
            }
        for (int32_t p : k.preds) {
            const Component& P = g.comps[static_cast<size_t>(p)];
            if (detail::color_rank(P.color) >= detail::color_rank(k.color)) return "pred not lower color";
            if (!std::binary_search(P.succs.begin(), P.succs.end(), k.id)) return "succ list inconsistent";
            if (p >= k.id) return "pred id not smaller than successor id (canonical order broken)";
        }
        if (k.owner >= 0 && !k.has_sharer(static_cast<uint32_t>(k.owner)))
            return "owner not a sharer in " + std::to_string(k.id);
        for (int32_t p : k.preds2) {
            const Component& P = g.comps[static_cast<size_t>(p)];
            if (detail::color_rank(P.color) >= detail::color_rank(k.color)) return "2-hop pred not lower color";
            if (std::binary_search(k.preds.begin(), k.preds.end(), p)) return "2-hop pred is also a 1-hop pred";
        }
    }
    return {};
}

/// Boxes of interest to a rank: its local boundary boxes and every box of
/// the components it owns (interior boxes are local by definition).
inline std::unordered_set<int64_t> interest_boundary_boxes(const ComponentGraph& g, uint32_t pid) {
    std::unordered_set<int64_t> out;
    for (int64_t m : g.grid.brick_boundary_boxes(pid)) out.insert(m);
    for (const auto& k : g.comps)
        if (k.owner == static_cast<int32_t>(pid)) for (int64_t m : k.boxes) out.insert(m);
    return out;
}

/// Components a rank neither owns nor shares but whose boxes are far
/// endpoints (within two hops) of one of its boxes — their skeleton sections
/// must reach the rank as auxiliary payloads.  Two kinds of reader: a LOCAL
/// box needs every 2-hop box's skeleton by the level's end (transition), any
/// color; a box of an OWNED component needs, for its sketch, only 2-hop boxes
/// of LOWER color (same/higher color is treated as un-eliminated there).
inline std::vector<int32_t> auxiliary_components(const ComponentGraph& g, uint32_t pid) {
    std::set<int32_t> out;
    auto consider = [&](int64_t n, int reader_color /* -1: any */) {
        for (uint64_t mb : g.grid.neighbors_2hop(n)) {
            const int32_t c = g.comp_id_of_box(static_cast<int64_t>(mb));
            if (c < 0) {
                if (g.grid.on_brick_boundary(static_cast<int64_t>(mb)))
                    throw std::runtime_error("auxiliary_components: boundary box " + std::to_string(mb) +
                                             " outside the graph's neighbourhood");
                continue;
            }
            const Component& k = g.comps[static_cast<size_t>(c)];
            if (reader_color >= 0 && detail::color_rank(k.color) >= reader_color) continue;
            if (k.owner == static_cast<int32_t>(pid) || k.has_sharer(pid)) continue;
            out.insert(c);
        }
    };
    // (a box deeper than one layer inside the brick sees, within two hops,
    // only this brick's boundary, whose components the rank shares)
    const auto o = g.grid.proc_offset(pid);
    const uint32_t b = g.grid.brick, zmax = g.grid.dimension == 3 ? b : 1;
    auto outer = [b](uint32_t v) { return v < 2 || v + 2 >= b; };
    for (uint32_t z = 0; z < zmax; ++z)
        for (uint32_t y = 0; y < b; ++y)
            for (uint32_t x = 0; x < b; ++x) {
                if (!(outer(x) || outer(y) || (g.grid.dimension == 3 && outer(z)))) continue;
                consider(g.grid.encode(o[0] + x, o[1] + y, g.grid.dimension == 3 ? o[2] + z : 0), -1);
            }
    for (const auto& k : g.comps)
        if (k.owner == static_cast<int32_t>(pid))
            for (int64_t m : k.boxes) consider(m, detail::color_rank(k.color));
    return std::vector<int32_t>(out.begin(), out.end());
}

/// What a recipient gets of a box's post-elimination state.
enum class PayloadKind : int8_t { FULL = 0, COMPACT = 1, SKELETON = 2 };

/// Recipients of box m's state (pid, kind), excluding `sender`:
///  - FULL     to m's home rank (transition and solve read everything);
///  - COMPACT  to home ranks of m's 1-hop neighbours (their local copies of
///             the pairs) and to owners of HIGHER-colour components adjacent
///             to m (their elimination reads the pairs); lower-colour owners
///             finished before m and never read them again;
///  - SKELETON to home ranks of boxes within two hops (the transition
///             regenerates far blocks with both endpoints eliminated) and to
///             owners of HIGHER-colour components within two hops (their
///             sketches read lower-colour far endpoints' skeletons).
/// Consistency (checked by the unit test): every COMPACT/FULL recipient is a
/// sharer of m's component; every SKELETON recipient is a sharer or lists the
/// component in auxiliary_components().
inline std::vector<std::pair<uint32_t, PayloadKind>> payload_recipients(const ComponentGraph& g, uint32_t sender,
                                                                        int64_t m) {
    std::map<uint32_t, PayloadKind> out;
    auto add = [&](uint32_t p, PayloadKind kind) {
        auto it = out.find(p);
        if (it == out.end()) out.emplace(p, kind);
        else if (static_cast<int8_t>(kind) < static_cast<int8_t>(it->second)) it->second = kind;
    };
    const int32_t cm = g.comp_id_of_box(m);
    const int my_color = cm >= 0 ? detail::color_rank(g.comps[static_cast<size_t>(cm)].color) : -1;
    auto owner_if_higher = [&](int64_t n, PayloadKind kind) {
        const int32_t c = g.comp_id_of_box(n);
        if (c < 0) {
            if (g.grid.on_brick_boundary(n))
                throw std::runtime_error("payload_recipients: boundary box " + std::to_string(n) +
                                         " outside the graph's neighbourhood");
            return;
        }
        const Component& k = g.comps[static_cast<size_t>(c)];
        if (detail::color_rank(k.color) > my_color) add(static_cast<uint32_t>(k.owner), kind);
    };
    for (uint64_t nb : g.grid.neighbors_2hop(m)) {
        const int64_t n = static_cast<int64_t>(nb);
        add(g.grid.proc_of_box(n), PayloadKind::SKELETON);
        owner_if_higher(n, PayloadKind::SKELETON);
    }
    for (uint64_t nb : g.grid.neighbors(m)) {
        const int64_t n = static_cast<int64_t>(nb);
        add(g.grid.proc_of_box(n), PayloadKind::COMPACT);
        owner_if_higher(n, PayloadKind::COMPACT);
    }
    out[g.grid.proc_of_box(m)] = PayloadKind::FULL;
    out.erase(sender);
    return std::vector<std::pair<uint32_t, PayloadKind>>(out.begin(), out.end());
}

/// One-line description for the phase report.
inline std::string describe(const ComponentGraph& g) {
    std::map<Color, std::pair<size_t, size_t>> n;   // color -> (components, boxes)
    for (const auto& k : g.comps) { n[k.color].first++; n[k.color].second += static_cast<size_t>(k.cost); }
    std::string s = "components:";
    for (const auto& kv : n)
        s += std::string(" ") + color_name(kv.first) + " " + std::to_string(kv.second.first) + " (" +
             std::to_string(kv.second.second) + " boxes)";
    if (!g.comps.empty() && g.comps.front().owner >= 0) {
        for (Color c : boundary_colors(g.grid.dimension)) {
            const LoadSummary l = owned_load(g, c);
            s += std::string(" | ") + color_name(c) + " owned load min/mean/max " + std::to_string(static_cast<int>(l.min)) +
                 "/" + std::to_string(static_cast<int>(l.mean)) + "/" + std::to_string(static_cast<int>(l.max));
        }
        const LoadSummary t = owned_load_total(g);
        s += " | total owned load min/mean/max " + std::to_string(static_cast<int>(t.min)) + "/" +
             std::to_string(static_cast<int>(t.mean)) + "/" + std::to_string(static_cast<int>(t.max));
    }
    return s;
}

/// A local graph against the global one, for brick `pid`: the first
/// difference as text, empty when they agree (the owner_graph_test checks,
/// also run in situ under BPACK_CHECK=graphcheck).  Both graphs assigned.
inline std::string compare_component_graphs(const ComponentGraph& G, const ComponentGraph& L, uint32_t pid) {
    auto s = [](auto v) { return std::to_string(v); };
    if (G.comps.size() != L.comps.size()) return "component count " + s(L.comps.size()) + " against " + s(G.comps.size());
    for (size_t i = 0; i < G.comps.size(); ++i) {
        const Component& g = G.comps[i];
        const Component& l = L.comps[i];
        const std::string who = "component " + s(i) + ": ";
        if (l.id != g.id) return who + "id";
        if (l.color != g.color) return who + "colour";
        if (l.first_morton != g.first_morton) return who + "first Morton index";
        if (l.cost != g.cost) return who + "cost";
        if (l.owner != g.owner) return who + "owner " + s(l.owner) + " against " + s(g.owner);
        if (l.upward_rank != g.upward_rank) return who + "upward rank";
        if (l.sharers != g.sharers) return who + "sharers";
        if (l.home_counts != g.home_counts) return who + "home counts";
        if (l.adjacent != g.adjacent) return who + "adjacent";
        if (l.visible != g.visible) return who + "visible";
        if (l.preds != g.preds) return who + "preds";
        if (l.succs != g.succs) return who + "succs";
        if (l.preds2 != g.preds2) return who + "preds2";
        if (l.succs2 != g.succs2) return who + "succs2";
        const bool touches = g.has_sharer(pid) || g.is_adjacent(pid);
        if (touches && !l.boxes_complete) return who + "touches the brick but its boxes are not complete";
        if (l.boxes_complete) {
            if (l.boxes != g.boxes) return who + "boxes";
        } else {
            for (int64_t m : l.boxes)
                if (!std::binary_search(g.boxes.begin(), g.boxes.end(), m)) return who + "a box not in the component";
        }
    }
    size_t n = 0;
    for (uint32_t q : L.local_bricks) {
        for (int64_t m : G.grid.brick_boundary_boxes(q)) {
            ++n;
            if (L.comp_id_of_box(m) != G.comp_id_of_box(m)) return "box " + s(m) + " of brick " + s(q) + " in another component";
        }
    }
    if (L.box_to_comp.size() != n) return "box_to_comp size";
    if (L.interior[pid].preds != G.interior[pid].preds) return "interior node";
    if (auxiliary_components(L, pid) != auxiliary_components(G, pid)) return "auxiliary components";
    const FetchPlan fl = build_fetch_plan(L, pid), fg = build_fetch_plan(G, pid);
    if (fl.fetch != fg.fetch || fl.serve != fg.serve || fl.foreign_eliminated != fg.foreign_eliminated) return "fetch plan";
    for (const auto& k : G.comps) {
        if (k.owner != static_cast<int32_t>(pid)) continue;
        for (int64_t m : k.boxes)
            if (payload_recipients(L, pid, m) != payload_recipients(G, pid, m)) return "payload recipients of box " + s(m);
    }
    if (assignment_hash(L) != assignment_hash(G)) return "assignment hash";
    return {};
}

}  // namespace dataflow
}  // namespace fmm
