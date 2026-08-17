#ifndef MORTON_HPP
#define MORTON_HPP

#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace morton {

/**
 * @brief Spread bits of a 32-bit value for 2D Morton encoding
 * Inserts a zero between each bit: abc -> a0b0c0
 */
inline uint64_t spread_bits_2d(uint32_t v) {
    uint64_t x = v;
    x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
    x = (x | (x << 8))  & 0x00FF00FF00FF00FF;
    x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0F;
    x = (x | (x << 2))  & 0x3333333333333333; 
    x = (x | (x << 1))  & 0x5555555555555555;
    return x;
}

/**
 * @brief Spread bits of a 32-bit value for 3D Morton encoding
 * Inserts two zeros between each bit: abc -> a00b00c00
 */
inline uint64_t spread_bits_3d(uint32_t v) {
    uint64_t x = v;
    x = (x | (x << 32)) & 0x001F00000000FFFF;
    x = (x | (x << 16)) & 0x001F0000FF0000FF;
    x = (x | (x << 8))  & 0x100F00F00F00F00F;
    x = (x | (x << 4))  & 0x10C30C30C30C30C3;
    x = (x | (x << 2))  & 0x1249249249249249;
    return x;
}

/**
 * @brief Compact bits for 2D Morton decoding
 * Extracts every other bit: a0b0c0 -> abc
 */
inline uint32_t compact_bits_2d(uint64_t x) {
    x = x & 0x5555555555555555;
    x = (x | (x >> 1))  & 0x3333333333333333;
    x = (x | (x >> 2))  & 0x0F0F0F0F0F0F0F0F;
    x = (x | (x >> 4))  & 0x00FF00FF00FF00FF;
    x = (x | (x >> 8))  & 0x0000FFFF0000FFFF;
    x = (x | (x >> 16)) & 0x00000000FFFFFFFF;
    return static_cast<uint32_t>(x);
}

/**
 * @brief Compact bits for 3D Morton decoding
 * Extracts every third bit: a00b00c00 -> abc
 */
inline uint32_t compact_bits_3d(uint64_t x) {
    x = x & 0x1249249249249249;
    x = (x | (x >> 2))  & 0x10C30C30C30C30C3;
    x = (x | (x >> 4))  & 0x100F00F00F00F00F;
    x = (x | (x >> 8))  & 0x001F0000FF0000FF;
    x = (x | (x >> 16)) & 0x001F00000000FFFF;
    x = (x | (x >> 32)) & 0x00000000001FFFFF;
    return static_cast<uint32_t>(x);
}

/**
 * @brief Convert 2D coordinates to Morton Z-order index
 * @param x X coordinate
 * @param y Y coordinate
 * @return Morton index
 */
inline uint64_t encode_2d(uint32_t x, uint32_t y) {
    return spread_bits_2d(x) | (spread_bits_2d(y) << 1);
}

/**
 * @brief Convert 3D coordinates to Morton Z-order index
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Morton index
 */
inline uint64_t encode_3d(uint32_t x, uint32_t y, uint32_t z) {
    return spread_bits_3d(x) | (spread_bits_3d(y) << 1) | (spread_bits_3d(z) << 2);
}

/**
 * @brief Decode Morton index to 2D coordinates
 * @param index Morton index
 * @param x Output X coordinate
 * @param y Output Y coordinate
 */
inline void decode_2d(uint64_t index, uint32_t& x, uint32_t& y) {
    x = compact_bits_2d(index);
    y = compact_bits_2d(index >> 1);
}

/**
 * @brief Decode Morton index to 3D coordinates
 * @param index Morton index
 * @param x Output X coordinate
 * @param y Output Y coordinate
 * @param z Output Z coordinate
 */
inline void decode_3d(uint64_t index, uint32_t& x, uint32_t& y, uint32_t& z) {
    x = compact_bits_3d(index);
    y = compact_bits_3d(index >> 1);
    z = compact_bits_3d(index >> 2);
}

/**
 * @brief Get all 8-connected neighbors in 2D (includes diagonals)
 * @param index Morton index of the cell
 * @param grid_size Size of the grid (assumes square grid)
 * @return Vector of Morton indices of neighbors
 */
inline std::vector<uint64_t> neighbors_2d(uint64_t index, uint32_t grid_size) {
    uint32_t x, y;
    decode_2d(index, x, y);
    
    std::vector<uint64_t> neighbors;
    neighbors.reserve(8);
    
    // Check all 8 neighbors (including diagonals)
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) continue; // Skip self
            
            int32_t nx = static_cast<int32_t>(x) + dx;
            int32_t ny = static_cast<int32_t>(y) + dy;
            
            // Check bounds
            if (nx >= 0 && nx < static_cast<int32_t>(grid_size) &&
                ny >= 0 && ny < static_cast<int32_t>(grid_size)) {
                neighbors.push_back(encode_2d(nx, ny));
            }
        }
    }
    
    return neighbors;
}

/**
 * @brief Get all 26-connected neighbors in 3D (includes diagonals)
 * @param index Morton index of the cell
 * @param grid_size Size of the grid (assumes cubic grid)
 * @return Vector of Morton indices of neighbors
 */
inline std::vector<uint64_t> neighbors_3d(uint64_t index, uint32_t grid_size) {
    uint32_t x, y, z;
    decode_3d(index, x, y, z);
    
    std::vector<uint64_t> neighbors;
    neighbors.reserve(26);
    
    // Check all 26 neighbors (3x3x3 - self)
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) continue; // Skip self
                
                int32_t nx = static_cast<int32_t>(x) + dx;
                int32_t ny = static_cast<int32_t>(y) + dy;
                int32_t nz = static_cast<int32_t>(z) + dz;
                
                // Check bounds
                if (nx >= 0 && nx < static_cast<int32_t>(grid_size) &&
                    ny >= 0 && ny < static_cast<int32_t>(grid_size) &&
                    nz >= 0 && nz < static_cast<int32_t>(grid_size)) {
                    neighbors.push_back(encode_3d(nx, ny, nz));
                }
            }
        }
    }
    
    return neighbors;
}

/**
 * @brief Get all neighbors within 2-hop distance in 2D (5x5 region)
 * @param index Morton index of the cell
 * @param grid_size Size of the grid (assumes square grid)
 * @return Vector of Morton indices of 2-hop neighbors
 */
inline std::vector<uint64_t> neighbors_2hop_2d(uint64_t index, uint32_t grid_size) {
    uint32_t x, y;
    decode_2d(index, x, y);
    
    std::vector<uint64_t> neighbors;
    neighbors.reserve(16); // 5x5 - 1
    
    // Check all cells in 5x5 region
    for (int dy = -2; dy <= 2; ++dy) {
        for (int dx = -2; dx <= 2; ++dx) {
            if (std::abs(dx) <= 1 && std::abs(dy) <= 1) continue; // Skip 1-hop and self
            
            int32_t nx = static_cast<int32_t>(x) + dx;
            int32_t ny = static_cast<int32_t>(y) + dy;
            
            // Check bounds
            if (nx >= 0 && nx < static_cast<int32_t>(grid_size) &&
                ny >= 0 && ny < static_cast<int32_t>(grid_size)) {
                neighbors.push_back(encode_2d(nx, ny));
            }
        }
    }
    
    return neighbors;
}

/**
 * @brief Get all neighbors within 2-hop distance in 3D (5x5x5 region)
 * @param index Morton index of the cell
 * @param grid_size Size of the grid (assumes cubic grid)
 * @return Vector of Morton indices of 2-hop neighbors
 */
inline std::vector<uint64_t> neighbors_2hop_3d(uint64_t index, uint32_t grid_size) {
    uint32_t x, y, z;
    decode_3d(index, x, y, z);
    
    std::vector<uint64_t> neighbors;
    neighbors.reserve(98); // 5x5x5 - 1
    
    // Check all cells in 5x5x5 region
    for (int dz = -2; dz <= 2; ++dz) {
        for (int dy = -2; dy <= 2; ++dy) {
            for (int dx = -2; dx <= 2; ++dx) {
                if (std::abs(dx) <= 1 && std::abs(dy) <= 1 && std::abs(dz) <= 1) continue; // Skip 1-hop and self
                
                int32_t nx = static_cast<int32_t>(x) + dx;
                int32_t ny = static_cast<int32_t>(y) + dy;
                int32_t nz = static_cast<int32_t>(z) + dz;
                
                // Check bounds
                if (nx >= 0 && nx < static_cast<int32_t>(grid_size) &&
                    ny >= 0 && ny < static_cast<int32_t>(grid_size) &&
                    nz >= 0 && nz < static_cast<int32_t>(grid_size)) {
                    neighbors.push_back(encode_3d(nx, ny, nz));
                }
            }
        }
    }
    
    return neighbors;
}

/**
 * @brief Check if a number is a power of 4
 */
inline bool is_power_of_4(uint32_t n) {
    if (n == 0) return false;
    if (n == 1) return true;
    
    // n must be a power of 2 first
    if ((n & (n - 1)) != 0) return false;
    
    // Count trailing zeros - must be even for power of 4
    int zeros = __builtin_ctz(n);
    return (zeros % 2) == 0;
}

/**
 * @brief Check if a number is a power of 8
 */
inline bool is_power_of_8(uint32_t n) {
    if (n == 0) return false;
    if (n == 1) return true;
    
    // n must be a power of 2 first
    if ((n & (n - 1)) != 0) return false;
    
    // Count trailing zeros - must be divisible by 3 for power of 8
    int zeros = __builtin_ctz(n);
    return (zeros % 3) == 0;
}
/**
 * @brief Assign 2D Morton indices to processes
 * @param indices Vector of Morton indices
 * @param num_processes Number of processes (must be 4^k for some k >= 0)
 * @param grid_size Size of the grid per dimension
 * @return Vector of process IDs corresponding to each index
 * @throws std::invalid_argument if num_processes is not a power of 4
 */
inline std::vector<uint32_t> assign_to_processes_2d(
    const std::vector<uint64_t>& indices, 
    uint32_t num_processes,
    uint32_t grid_size) {
    
    if (!is_power_of_4(num_processes)) {
        throw std::invalid_argument("num_processes must be 4^k for some k >= 0");
    }
    
    // Calculate number of processes per dimension: 2^k where num_processes = 4^k
    uint32_t k = __builtin_ctz(num_processes) / 2;
    uint32_t procs_per_dim = 1 << k;  // 2^k
    
    // Calculate cells per process per dimension
    uint32_t cells_per_proc = grid_size / procs_per_dim;
    
    // Calculate shift amount
    uint32_t shift = __builtin_ctz(cells_per_proc);
    
    std::vector<uint32_t> process_ids;
    process_ids.reserve(indices.size());
    
    for (uint64_t index : indices) {
        uint32_t x, y;
        decode_2d(index, x, y);
        
        // Get process coordinates
        uint32_t proc_x = x >> shift;
        uint32_t proc_y = y >> shift;
        
        // Encode process coordinates to get process ID
        uint32_t proc_id = static_cast<uint32_t>(encode_2d(proc_x, proc_y));
        process_ids.push_back(proc_id);
    }
    
    return process_ids;
}

/**
 * @brief Assign 3D Morton indices to processes
 * @param indices Vector of Morton indices
 * @param num_processes Number of processes (must be 8^k for some k >= 0)
 * @param grid_size Size of the grid per dimension
 * @return Vector of process IDs corresponding to each index
 * @throws std::invalid_argument if num_processes is not a power of 8
 */
inline std::vector<uint32_t> assign_to_processes_3d(
    const std::vector<uint64_t>& indices, 
    uint32_t num_processes,
    uint32_t grid_size) {
    
    if (!is_power_of_8(num_processes)) {
        throw std::invalid_argument("num_processes must be 8^k for some k >= 0");
    }
    
    // Calculate number of processes per dimension: 2^k where num_processes = 8^k
    uint32_t k = __builtin_ctz(num_processes) / 3;
    uint32_t procs_per_dim = 1 << k;  // 2^k
    
    // Calculate cells per process per dimension
    uint32_t cells_per_proc = grid_size / procs_per_dim;
    
    // Calculate shift amount
    uint32_t shift = __builtin_ctz(cells_per_proc);
    
    std::vector<uint32_t> process_ids;
    process_ids.reserve(indices.size());
    
    for (uint64_t index : indices) {
        uint32_t x, y, z;
        decode_3d(index, x, y, z);
        
        // Get process coordinates
        uint32_t proc_x = x >> shift;
        uint32_t proc_y = y >> shift;
        uint32_t proc_z = z >> shift;
        
        // Encode process coordinates to get process ID
        uint32_t proc_id = static_cast<uint32_t>(encode_3d(proc_x, proc_y, proc_z));
        process_ids.push_back(proc_id);
    }
    
    return process_ids;
}

} // namespace morton

#endif // MORTON_HPP