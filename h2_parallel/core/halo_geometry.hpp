#ifndef HALO_GEOMETRY_HPP
#define HALO_GEOMETRY_HPP
// Geometry of the replicated CA levels' halo: the box colors of the CA
// schedule (blue, orange, purple, green, from a box's position in its
// region) and the ghost and assisting sets of a region, which depend on the
// Morton layout of the level alone.  The tree computes them for the rank's
// own region (tree_impl.hpp compute_ghost_boxes_and_colors,
// compute_assisting_boxes) and the level-entry halo gather for the
// neighboring regions (serialization.hpp gather_CA_factorization_data), so
// no rank has to ask another what it needs.

#include "morton.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <set>
#include <stdexcept>
#include <vector>

namespace fmm {

/**
 * @brief Check if a box is a blue box (corner neighborhood) using LOCAL grid coordinates
 * 
 * 2D: Box is one of the 4 corners of the LOCAL process grid
 * 3D: Box is in a 2×2×2 neighborhood around any LOCAL cube corner, excluding diagonal interior
 * 
 * @param morton_index Morton index of the box (global)
 * @param level Tree level
 * @param dimension 2 or 3
 * @param grid_size Global grid size
 * @param local_grid_size Size of local grid per process
 * @param local_offset Starting global coordinates of this process's region
 * @return true if box is blue
 */
inline bool is_blue_box(int64_t morton_index, int32_t level, int32_t dimension, uint32_t grid_size,
                        uint32_t local_grid_size, const uint32_t local_offset[3]) {
    uint32_t x_global, y_global, z_global = 0;
    
    if (dimension == 2) {
        morton::decode_2d(morton_index, x_global, y_global);
        
        // Convert to local coordinates
        int32_t x = x_global - local_offset[0];
        int32_t y = y_global - local_offset[1];
        
        // Check if within this process's local region
        if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
            y < 0 || y >= static_cast<int32_t>(local_grid_size)) {
            return false;
        }
        
        // In 2D, blue boxes are the 4 corners of the LOCAL grid
        return (x == 0 || x == static_cast<int32_t>(local_grid_size) - 1) && 
               (y == 0 || y == static_cast<int32_t>(local_grid_size) - 1);
        
    } else {
        morton::decode_3d(morton_index, x_global, y_global, z_global);
        
        // Convert to local coordinates
        int32_t x = x_global - local_offset[0];
        int32_t y = y_global - local_offset[1];
        int32_t z = z_global - local_offset[2];
        
        // Check if within this process's local region
        if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
            y < 0 || y >= static_cast<int32_t>(local_grid_size) ||
            z < 0 || z >= static_cast<int32_t>(local_grid_size)) {
            return false;
        }
        
        // Local corners: {0, local_grid_size-1} for each dimension
        uint32_t corners[8][3] = {
            {0, 0, 0}, {local_grid_size-1, 0, 0},
            {0, local_grid_size-1, 0}, {local_grid_size-1, local_grid_size-1, 0},
            {0, 0, local_grid_size-1}, {local_grid_size-1, 0, local_grid_size-1},
            {0, local_grid_size-1, local_grid_size-1}, {local_grid_size-1, local_grid_size-1, local_grid_size-1}
        };
        
        for (int c = 0; c < 8; ++c) {
            uint32_t cx = corners[c][0];
            uint32_t cy = corners[c][1];
            uint32_t cz = corners[c][2];
            
            // Check if box is in 2×2×2 neighborhood (using local coordinates)
            int32_t dx = std::abs(x - static_cast<int32_t>(cx));
            int32_t dy = std::abs(y - static_cast<int32_t>(cy));
            int32_t dz = std::abs(z - static_cast<int32_t>(cz));
            
            if (dx <= 1 && dy <= 1 && dz <= 1) {
                // Exclude diagonal interior box (all distances = 1)
                if (dx == 1 && dy == 1 && dz == 1) {
                    continue;
                }
                return true;
            }
        }
        
        return false;
    }
}


/**
 * @brief Get Manhattan distance-1 neighbors (face-adjacent only)
 * Returns up to 6 neighbors in 3D, up to 4 neighbors in 2D
 */
inline std::vector<uint64_t> get_manhattan_neighbors_3d(uint64_t morton_index, uint32_t grid_size) {
    uint32_t x, y, z;
    morton::decode_3d(morton_index, x, y, z);
    
    std::vector<uint64_t> neighbors;
    neighbors.reserve(6);
    
    // Check all 6 face directions
    if (x > 0) neighbors.push_back(morton::encode_3d(x-1, y, z));
    if (x < grid_size-1) neighbors.push_back(morton::encode_3d(x+1, y, z));
    if (y > 0) neighbors.push_back(morton::encode_3d(x, y-1, z));
    if (y < grid_size-1) neighbors.push_back(morton::encode_3d(x, y+1, z));
    if (z > 0) neighbors.push_back(morton::encode_3d(x, y, z-1));
    if (z < grid_size-1) neighbors.push_back(morton::encode_3d(x, y, z+1));
    
    return neighbors;
}

inline std::vector<uint64_t> get_manhattan_neighbors_2d(uint64_t morton_index, uint32_t grid_size) {
    uint32_t x, y;
    morton::decode_2d(morton_index, x, y);
    
    std::vector<uint64_t> neighbors;
    neighbors.reserve(4);
    
    if (x > 0) neighbors.push_back(morton::encode_2d(x-1, y));
    if (x < grid_size-1) neighbors.push_back(morton::encode_2d(x+1, y));
    if (y > 0) neighbors.push_back(morton::encode_2d(x, y-1));
    if (y < grid_size-1) neighbors.push_back(morton::encode_2d(x, y+1));
    
    return neighbors;
}

/**
 * @brief Check if a box is an orange box (boundary box adjacent to blue) using LOCAL grid
 * 
 * 2D: LOCAL boundary box at Manhattan distance 1 from any LOCAL corner
 * 3D: ANY boundary box (edge or face) adjacent to at least one blue box
 * 
 * @param morton_index Morton index of the box (global)
 * @param level Tree level
 * @param dimension 2 or 3
 * @param grid_size Global grid size
 * @param local_grid_size Size of local grid per process
 * @param local_offset Starting global coordinates of this process's region
 * @return true if box is orange
 */
inline bool is_orange_box(int64_t morton_index, int32_t level, int32_t dimension, uint32_t grid_size,
                          uint32_t local_grid_size, const uint32_t local_offset[3]) {
    uint32_t x_global, y_global, z_global = 0;
    
    if (dimension == 2) {
        morton::decode_2d(morton_index, x_global, y_global);
        
        // Convert to local coordinates
        int32_t x = x_global - local_offset[0];
        int32_t y = y_global - local_offset[1];
        
        // Check if within this process's local region
        if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
            y < 0 || y >= static_cast<int32_t>(local_grid_size)) {
            return false;
        }
        
        // Count LOCAL boundary coordinates
        int boundary_count = 0;
        if (x == 0 || x == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        if (y == 0 || y == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        
        // Must be on boundary
        if (boundary_count == 0) return false;
        
        // Must not be a corner (blue)
        if (boundary_count == 2) return false;
        
        // Check Manhattan distance to each LOCAL corner
        uint32_t corners[4][2] = {
            {0, 0}, {local_grid_size-1, 0}, 
            {0, local_grid_size-1}, {local_grid_size-1, local_grid_size-1}
        };
        
        for (int c = 0; c < 4; ++c) {
            int32_t dist = std::abs(x - static_cast<int32_t>(corners[c][0])) +
                          std::abs(y - static_cast<int32_t>(corners[c][1]));
            if (dist == 1) return true;
        }
        
        return false;
        
    } else {
        morton::decode_3d(morton_index, x_global, y_global, z_global);
        
        // Convert to local coordinates
        int32_t x = x_global - local_offset[0];
        int32_t y = y_global - local_offset[1];
        int32_t z = z_global - local_offset[2];
        
        // Check if within this process's local region
        if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
            y < 0 || y >= static_cast<int32_t>(local_grid_size) ||
            z < 0 || z >= static_cast<int32_t>(local_grid_size)) {
            return false;
        }
        
        // Count LOCAL boundary coordinates
        int boundary_count = 0;
        if (x == 0 || x == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        if (y == 0 || y == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        if (z == 0 || z == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        
        // Must be on boundary (edge or face)
        if (boundary_count == 0) return false;
        
        // Must be adjacent (Manhattan distance 1) to at least one blue box
        // Check only face-adjacent neighbors (Manhattan distance 1)
        std::vector<uint64_t> neighbors = get_manhattan_neighbors_3d(morton_index, grid_size);
        
        for (uint64_t n_u : neighbors) {
            int64_t n = static_cast<int64_t>(n_u);
            if (is_blue_box(n, level, dimension, grid_size, local_grid_size, local_offset)) {
                return true;
            }
        }
        
        return false;
    }
}

/**
 * @brief Check if a box is a purple box (3D only: edges + face 2nd layer, not orange) using LOCAL grid
 * 
 * 3D: Edge boxes (b=2) OR face 2nd layer boxes (b=1, min_face_dist=1), 
 *     but NOT adjacent to any blue box
 * 
 * @param morton_index Morton index of the box (global)
 * @param level Tree level
 * @param grid_size Global grid size
 * @param local_grid_size Size of local grid per process
 * @param local_offset Starting global coordinates of this process's region
 * @return true if box is purple
 */
inline bool is_purple_box(int64_t morton_index, int32_t level, uint32_t grid_size,
                          uint32_t local_grid_size, const uint32_t local_offset[3]) {
    uint32_t x_global, y_global, z_global;
    morton::decode_3d(morton_index, x_global, y_global, z_global);
    
    // Convert to local coordinates
    int32_t x = x_global - local_offset[0];
    int32_t y = y_global - local_offset[1];
    int32_t z = z_global - local_offset[2];
    
    // Check if within this process's local region
    if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
        y < 0 || y >= static_cast<int32_t>(local_grid_size) ||
        z < 0 || z >= static_cast<int32_t>(local_grid_size)) {
        return false;
    }
    
    // Count LOCAL boundary coordinates
    int boundary_count = 0;
    if (x == 0 || x == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
    if (y == 0 || y == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
    if (z == 0 || z == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
    
    bool is_candidate = false;
    
    // Case 1: On an edge (boundary_count == 2)
    if (boundary_count == 2) {
        is_candidate = true;
    }
    // Case 2: On a face in the second-to-last layer (boundary_count == 1)
    else if (boundary_count == 1) {
        // Calculate minimum distance to edge within the face (LOCAL coordinates)
        auto dist_to_boundary = [local_grid_size](int32_t coord) -> uint32_t {
            return std::min(coord, static_cast<int32_t>(local_grid_size) - 1 - coord);
        };
        
        uint32_t min_face_dist = local_grid_size;
        if (x != 0 && x != static_cast<int32_t>(local_grid_size) - 1) {
            min_face_dist = std::min(min_face_dist, dist_to_boundary(x));
        }
        if (y != 0 && y != static_cast<int32_t>(local_grid_size) - 1) {
            min_face_dist = std::min(min_face_dist, dist_to_boundary(y));
        }
        if (z != 0 && z != static_cast<int32_t>(local_grid_size) - 1) {
            min_face_dist = std::min(min_face_dist, dist_to_boundary(z));
        }
        
        // Must be in second-to-last layer
        if (min_face_dist == 1) {
            is_candidate = true;
        }
    }
    
    if (!is_candidate) return false;
    
    // Must NOT be adjacent to any blue box (otherwise it would be orange)
    std::vector<uint64_t> neighbors = get_manhattan_neighbors_3d(morton_index, grid_size);
    
    for (uint64_t n_u : neighbors) {
        int64_t n = static_cast<int64_t>(n_u);
        if (is_blue_box(n, level, 3, grid_size, local_grid_size, local_offset)) {
            return false;  // Adjacent to blue means it's orange, not purple
        }
    }
    
    return true;
}

// ---------------------------------------------------------------------------
// The regions of a level and the halo sets of a region.
//
// With num_active_processes = (2^dimension)^k, the 2^level grid of boxes is
// tiled by cubes of `width` boxes, one per region (process); region q holds
// the Morton range [q * boxes_per_region, (q + 1) * boxes_per_region): the
// top bits of a box's Morton index are its region's (assign_to_processes_nd).
// ---------------------------------------------------------------------------
struct RegionLayout {
    int dimension = 3;
    int32_t level = 0;
    uint32_t grid = 1;             // boxes per dimension
    uint32_t procs_per_dim = 1;    // regions per dimension
    uint32_t width = 1;            // boxes per dimension of a region
    int64_t boxes_per_region = 1;
    uint32_t regions = 1;
    bool uniform = false;          // the layout above holds

    static RegionLayout of(int dimension, int32_t level, int num_active_processes) {
        RegionLayout L;
        L.dimension = dimension;
        L.level = level;
        L.grid = uint32_t{1} << level;
        if (num_active_processes <= 0 || dimension < 2 || dimension > 3 ||
            !morton::is_valid_process_count(dimension, static_cast<uint32_t>(num_active_processes))) {
            return L;
        }
        L.regions = static_cast<uint32_t>(num_active_processes);
        L.procs_per_dim = morton::processes_per_dimension(dimension, L.regions);
        if (L.procs_per_dim > L.grid) return L;
        L.width = L.grid / L.procs_per_dim;
        L.boxes_per_region = 1;
        for (int d = 0; d < dimension; ++d) L.boxes_per_region *= L.width;
        L.uniform = static_cast<int64_t>(L.regions) * L.boxes_per_region == (int64_t{1} << (dimension * level));
        return L;
    }
    uint32_t region_of(int64_t morton) const { return static_cast<uint32_t>(morton / boxes_per_region); }
    int64_t first(uint32_t q) const { return static_cast<int64_t>(q) * boxes_per_region; }
    bool in_region(int64_t morton, uint32_t q) const { return region_of(morton) == q; }
    // grid coordinates of a box, or of a region's first box
    std::array<int32_t, 3> coords(int64_t morton) const {
        uint32_t x, y, z;
        morton::decode_nd(dimension, static_cast<uint64_t>(morton), x, y, z);
        return {static_cast<int32_t>(x), static_cast<int32_t>(y), static_cast<int32_t>(z)};
    }
    std::array<uint32_t, 3> offset(uint32_t q) const {
        uint32_t x, y, z;
        morton::decode_nd(dimension, q, x, y, z);
        return {x * width, y * width, z * width};
    }
    std::array<uint32_t, 3> offset_of(int64_t morton) const { return offset(region_of(morton)); }
    int64_t encode(int32_t x, int32_t y, int32_t z) const {
        return static_cast<int64_t>(morton::encode_nd(dimension, static_cast<uint32_t>(x), static_cast<uint32_t>(y),
                                                      static_cast<uint32_t>(z)));
    }
    bool in_grid(int32_t x, int32_t y, int32_t z) const {
        const int32_t g = static_cast<int32_t>(grid);
        return x >= 0 && x < g && y >= 0 && y < g && (dimension < 3 || (z >= 0 && z < g));
    }
    // Chebyshev distance of a box to region q's cube (0 inside)
    int32_t distance_to(const std::array<int32_t, 3>& c, uint32_t q) const {
        const auto off = offset(q);
        int32_t dist = 0;
        for (int d = 0; d < dimension; ++d) {
            const int32_t lo = static_cast<int32_t>(off[static_cast<size_t>(d)]);
            const int32_t hi = lo + static_cast<int32_t>(width) - 1;
            dist = std::max({dist, lo - c[static_cast<size_t>(d)], c[static_cast<size_t>(d)] - hi});
        }
        return dist;
    }
    // Chebyshev distance between regions q and r on the grid of regions
    uint32_t region_distance(uint32_t q, uint32_t r) const {
        uint32_t qx, qy, qz, rx, ry, rz;
        morton::decode_nd(dimension, q, qx, qy, qz);
        morton::decode_nd(dimension, r, rx, ry, rz);
        auto ad = [](uint32_t a, uint32_t b) { return a > b ? a - b : b - a; };
        return std::max({ad(qx, rx), ad(qy, ry), dimension == 3 ? ad(qz, rz) : 0u});
    }
    // whether the box lies on a face of its own region
    bool on_region_face(const std::array<int32_t, 3>& c) const {
        for (int d = 0; d < dimension; ++d) {
            const uint32_t l = static_cast<uint32_t>(c[static_cast<size_t>(d)]) % width;
            if (l == 0 || l == width - 1) return true;
        }
        return false;
    }
    // The boxes of region r within Chebyshev distance `radius` of region q's
    // cube (r may be q), in Morton order.
    std::vector<int64_t> boxes_near(uint32_t r, uint32_t q, int32_t radius) const {
        const auto oq = offset(q), orr = offset(r);
        int32_t lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
        for (int d = 0; d < dimension; ++d) {
            const int32_t a = static_cast<int32_t>(oq[static_cast<size_t>(d)]) - radius;
            const int32_t b = static_cast<int32_t>(oq[static_cast<size_t>(d)]) + static_cast<int32_t>(width) - 1 + radius;
            lo[d] = std::max(a, static_cast<int32_t>(orr[static_cast<size_t>(d)]));
            hi[d] = std::min(b, static_cast<int32_t>(orr[static_cast<size_t>(d)]) + static_cast<int32_t>(width) - 1);
            if (lo[d] > hi[d]) return {};
        }
        std::vector<int64_t> out;
        for (int32_t z = lo[2]; z <= hi[2]; ++z) {
            for (int32_t y = lo[1]; y <= hi[1]; ++y) {
                for (int32_t x = lo[0]; x <= hi[0]; ++x) out.push_back(encode(x, y, z));
            }
        }
        std::sort(out.begin(), out.end());
        return out;
    }
    // The boxes (of any region) within Chebyshev distance `radius` of box c,
    // the box itself excluded, within the grid.
    template<typename Fn>
    void for_each_within(const std::array<int32_t, 3>& c, int32_t radius, Fn&& fn) const {
        const int32_t zr = dimension == 3 ? radius : 0;
        for (int32_t dz = -zr; dz <= zr; ++dz) {
            for (int32_t dy = -radius; dy <= radius; ++dy) {
                for (int32_t dx = -radius; dx <= radius; ++dx) {
                    if (dx == 0 && dy == 0 && dz == 0) continue;
                    const int32_t x = c[0] + dx, y = c[1] + dy, z = c[2] + dz;
                    if (in_grid(x, y, z)) fn(encode(x, y, z), std::array<int32_t, 3>{x, y, z});
                }
            }
        }
    }
    // The boundary boxes of region q: those on its faces (equivalently,
    // check_boundary_condition: with a neighbor in another region or on the
    // domain boundary, which coincides with a region face), sorted.
    std::vector<int64_t> boundary(uint32_t q) const {
        const auto off = offset(q);
        const int32_t w = static_cast<int32_t>(width);
        std::vector<int64_t> out;
        const int32_t zlo = dimension == 3 ? static_cast<int32_t>(off[2]) : 0;
        const int32_t zhi = dimension == 3 ? zlo + w - 1 : 0;
        for (int32_t z = zlo; z <= zhi; ++z) {
            for (int32_t y = static_cast<int32_t>(off[1]); y < static_cast<int32_t>(off[1]) + w; ++y) {
                for (int32_t x = static_cast<int32_t>(off[0]); x < static_cast<int32_t>(off[0]) + w; ++x) {
                    const bool face = x == static_cast<int32_t>(off[0]) || x == static_cast<int32_t>(off[0]) + w - 1 ||
                                      y == static_cast<int32_t>(off[1]) || y == static_cast<int32_t>(off[1]) + w - 1 ||
                                      (dimension == 3 && (z == zlo || z == zhi));
                    if (face) out.push_back(encode(x, y, z));
                }
            }
        }
        std::sort(out.begin(), out.end());
        return out;
    }
};

// ---------------------------------------------------------------------------
// The ghost boxes of a region (compute_ghost_boxes_and_colors), as the sets
//   S1: the one-hop shell of the region's cube (the neighbors of its
//       boundary boxes in other regions);
//   S2 (3D): the two-hop neighbors of its boundary boxes that lie on a face
//       of their own region, i.e. the boxes at distance 2 of the cube on a
//       face of their region (a box at distance 2 has a boundary box of the
//       region at distance exactly 2);
//   S3: the orange completion: the neighbors of the blue ghosts (blue in
//       their own region) that are orange in theirs and not ghosts already.
// Every rank computes these for itself (the tree) and, restricted to its own
// boxes, for the neighboring regions (the senders of the level-entry halo
// gather), from the layout alone: both sides agree by construction.
// ---------------------------------------------------------------------------
struct RegionGhosts {
    std::vector<int64_t> ghosts;   // sorted
    std::vector<int64_t> blue;     // the blue boxes among S1 and S2 (sorted)
    std::vector<int64_t> orange;   // S3 (sorted)
};

namespace halo_detail {

// The ghosts of q beyond its one-hop shell: S2 \ S1 and S3, sorted, with
// the blue ghosts (for the colors).
inline void extra_ghosts_of_region(const RegionLayout& L, uint32_t q, std::vector<int64_t>& extra,
                                   std::vector<int64_t>& blue, std::vector<int64_t>& orange) {
    extra.clear();
    blue.clear();
    orange.clear();
    const int dim = L.dimension;
    const int32_t w = static_cast<int32_t>(L.width);
    auto in_s1s2 = [&](const std::array<int32_t, 3>& c) {
        const int32_t dist = L.distance_to(c, q);
        return dist == 1 || (dim == 3 && dist == 2 && L.on_region_face(c));
    };
    // S2 \ S1: the second shell, on faces of their regions
    if (dim == 3) {
        const auto oq = L.offset(q);
        for (int32_t z = static_cast<int32_t>(oq[2]) - 2; z <= static_cast<int32_t>(oq[2]) + w + 1; ++z) {
            for (int32_t y = static_cast<int32_t>(oq[1]) - 2; y <= static_cast<int32_t>(oq[1]) + w + 1; ++y) {
                for (int32_t x = static_cast<int32_t>(oq[0]) - 2; x <= static_cast<int32_t>(oq[0]) + w + 1; ++x) {
                    const std::array<int32_t, 3> c{x, y, z};
                    if (!L.in_grid(x, y, z) || L.distance_to(c, q) != 2 || !L.on_region_face(c)) continue;
                    extra.push_back(L.encode(x, y, z));
                }
            }
        }
    }
    // the blue ghosts: blue boxes of the regions around q within S1 or S2
    // (a region's blue boxes lie next to its corners)
    std::set<int64_t> blue_set;
    uint32_t qx, qy, qz;
    morton::decode_nd(dim, q, qx, qy, qz);
    const int32_t P = static_cast<int32_t>(L.procs_per_dim);
    for (int32_t rz = (dim == 3 ? static_cast<int32_t>(qz) - 1 : 0); rz <= (dim == 3 ? static_cast<int32_t>(qz) + 1 : 0); ++rz) {
        for (int32_t ry = static_cast<int32_t>(qy) - 1; ry <= static_cast<int32_t>(qy) + 1; ++ry) {
            for (int32_t rx = static_cast<int32_t>(qx) - 1; rx <= static_cast<int32_t>(qx) + 1; ++rx) {
                if (rx < 0 || ry < 0 || rz < 0 || rx >= P || ry >= P || rz >= P) continue;
                const uint32_t r = static_cast<uint32_t>(morton::encode_nd(dim, static_cast<uint32_t>(rx),
                                                                           static_cast<uint32_t>(ry),
                                                                           static_cast<uint32_t>(rz)));
                if (r == q) continue;
                const auto orr = L.offset(r);
                const uint32_t off3[3] = {orr[0], orr[1], orr[2]};
                // the boxes within one of a corner of r
                for (int cz = 0; cz < (dim == 3 ? 2 : 1); ++cz) {
                    for (int cy = 0; cy < 2; ++cy) {
                        for (int cx = 0; cx < 2; ++cx) {
                            const int32_t kx = static_cast<int32_t>(orr[0]) + (cx ? w - 1 : 0);
                            const int32_t ky = static_cast<int32_t>(orr[1]) + (cy ? w - 1 : 0);
                            const int32_t kz = dim == 3 ? static_cast<int32_t>(orr[2]) + (cz ? w - 1 : 0) : 0;
                            const int32_t zr = dim == 3 ? 1 : 0;
                            for (int32_t dz = -zr; dz <= zr; ++dz) {
                                for (int32_t dy = -1; dy <= 1; ++dy) {
                                    for (int32_t dx = -1; dx <= 1; ++dx) {
                                        const std::array<int32_t, 3> c{kx + dx, ky + dy, kz + dz};
                                        if (!L.in_grid(c[0], c[1], c[2]) || L.distance_to(c, r) != 0) continue;
                                        const int64_t g = L.encode(c[0], c[1], c[2]);
                                        if (!in_s1s2(c)) continue;
                                        if (is_blue_box(g, L.level, dim, L.grid, L.width, off3)) blue_set.insert(g);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    blue.assign(blue_set.begin(), blue_set.end());
    // S3: their neighbors that are orange in their own region and not ghosts
    std::set<int64_t> orange_set;
    for (int64_t g : blue) {
        L.for_each_within(L.coords(g), 1, [&](int64_t n, const std::array<int32_t, 3>& c) {
            if (L.region_of(n) == q || in_s1s2(c)) return;
            const auto on = L.offset_of(n);
            const uint32_t off3[3] = {on[0], on[1], on[2]};
            if (is_orange_box(n, L.level, dim, L.grid, L.width, off3)) orange_set.insert(n);
        });
    }
    orange.assign(orange_set.begin(), orange_set.end());
    extra.insert(extra.end(), orange.begin(), orange.end());
    std::sort(extra.begin(), extra.end());
    extra.erase(std::unique(extra.begin(), extra.end()), extra.end());
}

}  // namespace halo_detail

// The ghosts of a region beyond its one-hop shell (S2 \ S1 and S3, sorted)
// with its blue ghosts and orange completions: computed once per region and
// shared by the functions below (the orange completion is their main cost).
struct RegionExtra {
    std::vector<int64_t> extra, blue, orange;
};
inline RegionExtra region_extra(const RegionLayout& L, uint32_t q) {
    if (!L.uniform) throw std::runtime_error("region_extra: the level's regions are not uniform cubes");
    RegionExtra e;
    halo_detail::extra_ghosts_of_region(L, q, e.extra, e.blue, e.orange);
    return e;
}

// The ghosts of region q (all regions), with the blue and orange lists the
// color classification needs.
inline RegionGhosts ghosts_of_region(const RegionLayout& L, uint32_t q, const RegionExtra& e) {
    if (!L.uniform) throw std::runtime_error("ghosts_of_region: the level's regions are not uniform cubes");
    RegionGhosts out;
    out.blue = e.blue;
    out.orange = e.orange;
    // S1: the one-hop shell
    const auto oq = L.offset(q);
    const int32_t w = static_cast<int32_t>(L.width);
    const int32_t zlo = L.dimension == 3 ? static_cast<int32_t>(oq[2]) - 1 : 0;
    const int32_t zhi = L.dimension == 3 ? static_cast<int32_t>(oq[2]) + w : 0;
    for (int32_t z = zlo; z <= zhi; ++z) {
        for (int32_t y = static_cast<int32_t>(oq[1]) - 1; y <= static_cast<int32_t>(oq[1]) + w; ++y) {
            for (int32_t x = static_cast<int32_t>(oq[0]) - 1; x <= static_cast<int32_t>(oq[0]) + w; ++x) {
                const std::array<int32_t, 3> c{x, y, z};
                if (!L.in_grid(x, y, z) || L.distance_to(c, q) != 1) continue;
                out.ghosts.push_back(L.encode(x, y, z));
            }
        }
    }
    out.ghosts.insert(out.ghosts.end(), e.extra.begin(), e.extra.end());
    std::sort(out.ghosts.begin(), out.ghosts.end());
    out.ghosts.erase(std::unique(out.ghosts.begin(), out.ghosts.end()), out.ghosts.end());
    return out;
}
inline RegionGhosts ghosts_of_region(const RegionLayout& L, uint32_t q) {
    return ghosts_of_region(L, q, region_extra(L, q));
}

// The ghosts of region q that lie in region r (r != q), sorted.
inline std::vector<int64_t> ghosts_of_region_in(const RegionLayout& L, uint32_t q, uint32_t r, const RegionExtra& e) {
    if (!L.uniform) throw std::runtime_error("ghosts_of_region_in: the level's regions are not uniform cubes");
    std::vector<int64_t> out;
    if (r == q) return out;
    for (int64_t m : L.boxes_near(r, q, 1)) out.push_back(m);  // (all at distance 1: r != q)
    for (int64_t m : e.extra) {
        if (L.region_of(m) == r) out.push_back(m);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}
inline std::vector<int64_t> ghosts_of_region_in(const RegionLayout& L, uint32_t q, uint32_t r) {
    if (r == q) return {};
    return ghosts_of_region_in(L, q, r, region_extra(L, q));
}

// The assisting boxes of a CA region q that lie in region r (r != q): the
// boxes within two hops of a ghost of q that are neither q's nor ghosts
// (compute_assisting_boxes with the ghosts as seeds), sorted.  A box is
// within two hops of the one-hop shell exactly when it lies within three of
// the cube; the ghosts beyond the shell (S2, S3) reach two further.
inline std::vector<int64_t> assisting_of_region_in(const RegionLayout& L, uint32_t q, uint32_t r,
                                                   const RegionExtra& e) {
    if (!L.uniform) throw std::runtime_error("assisting_of_region_in: the level's regions are not uniform cubes");
    std::vector<int64_t> out;
    if (r == q) return out;
    auto is_ghost = [&](int64_t m, const std::array<int32_t, 3>& c) {
        return L.distance_to(c, q) == 1 || std::binary_search(e.extra.begin(), e.extra.end(), m);
    };
    for (int64_t m : L.boxes_near(r, q, 3)) {
        if (!is_ghost(m, L.coords(m))) out.push_back(m);
    }
    for (int64_t g : e.extra) {
        // (only the ghosts within two hops of r's cube reach it)
        const std::array<int32_t, 3> gc = L.coords(g);
        if (L.distance_to(gc, r) > 2) continue;
        L.for_each_within(gc, 2, [&](int64_t n, const std::array<int32_t, 3>& c) {
            if (L.region_of(n) != r || is_ghost(n, c)) return;
            out.push_back(n);
        });
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}
inline std::vector<int64_t> assisting_of_region_in(const RegionLayout& L, uint32_t q, uint32_t r) {
    if (r == q) return {};
    return assisting_of_region_in(L, q, r, region_extra(L, q));
}

// Regions whose halo (ghosts, assisting boxes) may reach into region q, or
// that q's halo may reach: ghosts lie within 3 box hops of their region's
// cube, assisting boxes within 5.
inline uint32_t halo_region_radius(const RegionLayout& L) { return (5 + L.width - 1) / L.width; }

template<typename Fn>
void for_each_region_within(const RegionLayout& L, uint32_t q, uint32_t radius, Fn&& fn) {
    uint32_t qx, qy, qz;
    morton::decode_nd(L.dimension, q, qx, qy, qz);
    const int32_t P = static_cast<int32_t>(L.procs_per_dim), R = static_cast<int32_t>(radius);
    const int32_t zr = L.dimension == 3 ? R : 0;
    for (int32_t dz = -zr; dz <= zr; ++dz) {
        for (int32_t dy = -R; dy <= R; ++dy) {
            for (int32_t dx = -R; dx <= R; ++dx) {
                const int32_t rx = static_cast<int32_t>(qx) + dx, ry = static_cast<int32_t>(qy) + dy,
                              rz = static_cast<int32_t>(qz) + dz;
                if (rx < 0 || ry < 0 || rz < 0 || rx >= P || ry >= P || (L.dimension == 3 && rz >= P)) continue;
                const uint32_t r = static_cast<uint32_t>(morton::encode_nd(L.dimension, static_cast<uint32_t>(rx),
                                                                           static_cast<uint32_t>(ry),
                                                                           static_cast<uint32_t>(rz)));
                if (r != q) fn(r);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The reference forms (the original enumeration from the region's boundary
// boxes), for the self-test (tests/halo_geometry_test.cpp) and the
// assisting sets of Color levels (seeds: the boundary boxes).
// ---------------------------------------------------------------------------
inline RegionGhosts ghosts_of_region_reference(const RegionLayout& L, uint32_t q) {
    const int dim = L.dimension;
    const uint32_t grid = L.grid;
    const std::vector<int64_t> boundary = L.boundary(q);
    std::set<int64_t> ghosts;
    for (int64_t b : boundary) {
        for (uint64_t n : morton::neighbors_nd(dim, static_cast<uint64_t>(b), grid)) {
            if (L.region_of(static_cast<int64_t>(n)) != q) ghosts.insert(static_cast<int64_t>(n));
        }
    }
    if (dim == 3) {
        for (int64_t b : boundary) {
            for (uint64_t n : morton::neighbors_2hop_3d(static_cast<uint64_t>(b), grid)) {
                const int64_t m = static_cast<int64_t>(n);
                if (L.region_of(m) == q) continue;
                if (L.on_region_face(L.coords(m))) ghosts.insert(m);
            }
        }
    }
    RegionGhosts out;
    for (int64_t g : ghosts) {
        const auto og = L.offset_of(g);
        const uint32_t off3[3] = {og[0], og[1], og[2]};
        if (is_blue_box(g, L.level, dim, grid, L.width, off3)) out.blue.push_back(g);
    }
    std::set<int64_t> orange;
    for (int64_t g : out.blue) {
        for (uint64_t nu : morton::neighbors_nd(dim, static_cast<uint64_t>(g), grid)) {
            const int64_t n = static_cast<int64_t>(nu);
            if (ghosts.count(n) != 0 || L.region_of(n) == q) continue;
            const auto on = L.offset_of(n);
            const uint32_t off3[3] = {on[0], on[1], on[2]};
            if (is_orange_box(n, L.level, dim, grid, L.width, off3)) {
                ghosts.insert(n);
                orange.insert(n);
            }
        }
    }
    out.ghosts.assign(ghosts.begin(), ghosts.end());
    out.orange.assign(orange.begin(), orange.end());
    return out;
}

// The boxes within one hop (`one_hop`) and within two hops (`all`) of the
// seeds, in other regions and not ghosts, sorted (compute_assisting_boxes).
struct RegionAssisting {
    std::vector<int64_t> all, one_hop;
};
inline RegionAssisting assisting_of_region_reference(const RegionLayout& L, uint32_t q,
                                                     const std::vector<int64_t>& seeds,
                                                     const std::vector<int64_t>& ghosts) {
    const int dim = L.dimension;
    const uint32_t grid = L.grid;
    std::set<int64_t> all, one;
    auto wanted = [&](int64_t n) { return L.region_of(n) != q && !std::binary_search(ghosts.begin(), ghosts.end(), n); };
    for (int64_t s : seeds) {
        for (uint64_t nu : morton::neighbors_nd(dim, static_cast<uint64_t>(s), grid)) {
            const int64_t n = static_cast<int64_t>(nu);
            if (wanted(n)) {
                all.insert(n);
                one.insert(n);
            }
        }
        for (uint64_t nu : morton::neighbors_2hop_nd(dim, static_cast<uint64_t>(s), grid)) {
            const int64_t n = static_cast<int64_t>(nu);
            if (wanted(n)) all.insert(n);
        }
    }
    RegionAssisting out;
    out.all.assign(all.begin(), all.end());
    out.one_hop.assign(one.begin(), one.end());
    return out;
}

}  // namespace fmm

#endif  // HALO_GEOMETRY_HPP
