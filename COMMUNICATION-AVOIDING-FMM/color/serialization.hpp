#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include "tree.hpp"
#include <cstring>
#include <vector>
#include <chrono>
#include <climits>
#include <atomic>
#include <exception>
#include <mutex>


namespace ser {

template <class T>
inline char* write_pod(char* ptr, const T& v) {
    static_assert(std::is_trivially_copyable_v<T>, "write_pod requires trivially copyable");
    std::memcpy(ptr, &v, sizeof(T));
    return ptr + sizeof(T);
}

template <class T>
inline const char* read_pod(const char* ptr, T& v) {
    static_assert(std::is_trivially_copyable_v<T>, "read_pod requires trivially copyable");
    std::memcpy(&v, ptr, sizeof(T));
    return ptr + sizeof(T);
}

template <class T>
inline size_t bytes_vector(const std::vector<T>& v) {
    static_assert(std::is_trivially_copyable_v<T>, "bytes_vector requires trivially copyable elements");
    return sizeof(uint64_t) + sizeof(T) * v.size();
}

template <class T>
inline char* write_vector(char* ptr, const std::vector<T>& v) {
    static_assert(std::is_trivially_copyable_v<T>, "write_vector requires trivially copyable elements");
    uint64_t n = static_cast<uint64_t>(v.size());
    ptr = write_pod(ptr, n);
    if (n) {
        std::memcpy(ptr, v.data(), sizeof(T) * static_cast<size_t>(n));
        ptr += sizeof(T) * static_cast<size_t>(n);
    }
    return ptr;
}

template <class T>
inline const char* read_vector(const char* ptr, std::vector<T>& v) {
    static_assert(std::is_trivially_copyable_v<T>, "read_vector requires trivially copyable elements");
    uint64_t n = 0;
    ptr = read_pod(ptr, n);
    v.resize(static_cast<size_t>(n));
    if (n) {
        std::memcpy(v.data(), ptr, sizeof(T) * static_cast<size_t>(n));
        ptr += sizeof(T) * static_cast<size_t>(n);
    }
    return ptr;
}

}

namespace fmm {

/**
 * @brief Get serialized size of MatrixStorage
 */
template<typename DataType>
size_t get_serialized_size(const MatrixStorage<DataType>& mat);

/**
 * @brief Serialize MatrixStorage into buffer
 * @return Pointer to position after serialized data
 */
template<typename DataType>
char* serialize(const MatrixStorage<DataType>& mat, char* buffer);

/**
 * @brief Deserialize MatrixStorage from buffer
 * @return Pointer to position after deserialized data
 */
template<typename DataType>
const char* deserialize(MatrixStorage<DataType>& mat, const char* buffer);

/**
 * @brief Get serialized size of ModifiedBlock
 */
template<typename DataType>
size_t get_serialized_size(const ModifiedBlock<DataType>& block);

/**
 * @brief Serialize ModifiedBlock into buffer
 * @return Pointer to position after serialized data
 */
template<typename DataType>
char* serialize(const ModifiedBlock<DataType>& block, char* buffer);

/**
 * @brief Deserialize ModifiedBlock from buffer
 * @return Pointer to position after deserialized data
 */
template<typename DataType>
const char* deserialize(ModifiedBlock<DataType>& block, const char* buffer);

/**
 * @brief Get serialized size of BoxData
 */
template<typename CoordType, typename DataType>
size_t get_serialized_size(const BoxData<CoordType, DataType>& box);

/**
 * @brief Serialize BoxData into buffer
 * @param box BoxData to serialize
 * @param buffer Buffer to write to (must be large enough)
 * @return Pointer to position after serialized data
 */
template<typename CoordType, typename DataType>
char* serialize(const BoxData<CoordType, DataType>& box, char* buffer);

/**
 * @brief Deserialize BoxData from buffer
 * @param box BoxData to deserialize into
 * @param buffer Buffer to read from
 * @return Pointer to position after deserialized data
 */
template<typename CoordType, typename DataType>
const char* deserialize(BoxData<CoordType, DataType>& box, const char* buffer);


/**
 * @brief Get serialized size of SolveDataRequest
 */
template<typename CoordType, typename DataType>
size_t get_serialized_size(const SolveDataRequest<CoordType, DataType>& request);

/**
 * @brief Serialize SolveDataRequest into buffer
 * @return Pointer to position after serialized data
 */
template<typename CoordType, typename DataType>
char* serialize(const SolveDataRequest<CoordType, DataType>& request, char* buffer);

/**
 * @brief Deserialize SolveDataRequest from buffer
 * @return Pointer to position after deserialized data
 */
template<typename CoordType, typename DataType>
const char* deserialize(SolveDataRequest<CoordType, DataType>& request, const char* buffer);

// ============================================================================
// Implementation
// ============================================================================



// -------- size in bytes (for MPI packing) --------
template <typename DataType>
static inline uint64_t bytes_pending(const PendingSolveUpdates<DataType>& p)
{
    static_assert(std::is_trivially_copyable_v<DataType>,
                  "PendingSolveUpdates serialization requires trivially copyable DataType");

    uint64_t bytes = 0;
    bytes += sizeof(uint64_t); // num full entries
    for (const auto& kv : p.full_updates) {
        bytes += sizeof(int64_t);   // morton
        bytes += sizeof(uint64_t);  // len
        bytes += static_cast<uint64_t>(kv.second.size()) * sizeof(DataType);
    }

    bytes += sizeof(uint64_t); // num skel entries
    for (const auto& kv : p.skel_updates) {
        bytes += sizeof(int64_t);   // morton
        bytes += sizeof(uint64_t);  // len
        bytes += static_cast<uint64_t>(kv.second.size()) * sizeof(DataType);
    }
    return bytes;
}

// -------- serialization helpers --------
namespace detail {
template <typename T>
static inline char* write_pod(char* out, const T& v) {
    std::memcpy(out, &v, sizeof(T));
    return out + sizeof(T);
}

template <typename T>
static inline const char* read_pod(const char* in, T& v) {
    std::memcpy(&v, in, sizeof(T));
    return in + sizeof(T);
}
} // namespace detail

/**
 * @brief Serialize PendingSolveUpdates into a caller-provided buffer of size bytes_pending(p).
 * @return pointer one-past-the-last written byte.
 */
template <typename DataType>
static inline char* serialize(const PendingSolveUpdates<DataType>& p, char* out)
{
    static_assert(std::is_trivially_copyable_v<DataType>,
                  "PendingSolveUpdates serialization requires trivially copyable DataType");

    // full_updates
    uint64_t nfull = static_cast<uint64_t>(p.full_updates.size());
    out = detail::write_pod(out, nfull);
    for (const auto& kv : p.full_updates) {
        const int64_t morton = kv.first;
        const uint64_t len   = static_cast<uint64_t>(kv.second.size());
        out = detail::write_pod(out, morton);
        out = detail::write_pod(out, len);
        if (len) {
            std::memcpy(out, kv.second.data(), static_cast<size_t>(len) * sizeof(DataType));
            out += static_cast<size_t>(len) * sizeof(DataType);
        }
    }

    // skel_updates
    uint64_t nskel = static_cast<uint64_t>(p.skel_updates.size());
    out = detail::write_pod(out, nskel);
    for (const auto& kv : p.skel_updates) {
        const int64_t morton = kv.first;
        const uint64_t len   = static_cast<uint64_t>(kv.second.size());
        out = detail::write_pod(out, morton);
        out = detail::write_pod(out, len);
        if (len) {
            std::memcpy(out, kv.second.data(), static_cast<size_t>(len) * sizeof(DataType));
            out += static_cast<size_t>(len) * sizeof(DataType);
        }
    }

    return out;
}

/**
 * @brief Deserialize PendingSolveUpdates from a buffer produced by serialize().
 *
 * This overwrites the contents of @p p (clears maps first).
 * @return pointer one-past-the-last read byte.
 */
template <typename DataType>
static inline const char* deserialize(PendingSolveUpdates<DataType>& p, const char* in)
{
    static_assert(std::is_trivially_copyable_v<DataType>,
                  "PendingSolveUpdates serialization requires trivially copyable DataType");

    p.full_updates.clear();
    p.skel_updates.clear();

    // full_updates
    uint64_t nfull = 0;
    in = detail::read_pod(in, nfull);
    for (uint64_t t = 0; t < nfull; ++t) {
        int64_t morton = -1;
        uint64_t len   = 0;
        in = detail::read_pod(in, morton);
        in = detail::read_pod(in, len);

        auto& v = p.full_updates[morton];
        v.resize(static_cast<size_t>(len));
        if (len) {
            std::memcpy(v.data(), in, static_cast<size_t>(len) * sizeof(DataType));
            in += static_cast<size_t>(len) * sizeof(DataType);
        }
    }

    // skel_updates
    uint64_t nskel = 0;
    in = detail::read_pod(in, nskel);
    for (uint64_t t = 0; t < nskel; ++t) {
        int64_t morton = -1;
        uint64_t len   = 0;
        in = detail::read_pod(in, morton);
        in = detail::read_pod(in, len);

        auto& v = p.skel_updates[morton];
        v.resize(static_cast<size_t>(len));
        if (len) {
            std::memcpy(v.data(), in, static_cast<size_t>(len) * sizeof(DataType));
            in += static_cast<size_t>(len) * sizeof(DataType);
        }
    }

    return in;
}

template <typename DataType>
inline size_t bytes_denseblock(const DenseBlock<DataType>& b) {
    static_assert(std::is_trivially_copyable_v<DataType>, "DenseBlock<DataType>: DataType must be trivially copyable");
    return sizeof(int64_t) + sizeof(int64_t) + ser::bytes_vector(b.data);
}

template <typename DataType>
inline char* serialize(const DenseBlock<DataType>& b, char* ptr) {
    static_assert(std::is_trivially_copyable_v<DataType>, "serialize(DenseBlock): DataType must be trivially copyable");
    ptr = ser::write_pod(ptr, b.rows);
    ptr = ser::write_pod(ptr, b.cols);
    ptr = ser::write_vector(ptr, b.data);
    return ptr;
}

template <typename DataType>
inline const char* deserialize(DenseBlock<DataType>& b, const char* ptr) {
    static_assert(std::is_trivially_copyable_v<DataType>, "deserialize(DenseBlock): DataType must be trivially copyable");
    ptr = ser::read_pod(ptr, b.rows);
    ptr = ser::read_pod(ptr, b.cols);
    ptr = ser::read_vector(ptr, b.data);

    // sanity (optional but recommended)
    const auto expect = static_cast<size_t>(b.rows * b.cols);
    if (b.rows < 0 || b.cols < 0) throw std::runtime_error("DenseBlock: negative dims");
    if (b.data.size() != expect) throw std::runtime_error("DenseBlock: data.size != rows*cols");
    return ptr;
}


inline size_t bytes_key(const ReplaceKey&) { return sizeof(int64_t) + sizeof(int64_t); }

inline char* serialize(const ReplaceKey& k, char* ptr) {
    ptr = ser::write_pod(ptr, k.target_box);
    ptr = ser::write_pod(ptr, k.neighbor_box);
    return ptr;
}

inline const char* deserialize(ReplaceKey& k, const char* ptr) {
    ptr = ser::read_pod(ptr, k.target_box);
    ptr = ser::read_pod(ptr, k.neighbor_box);
    return ptr;
}

inline size_t bytes_key(const EdgeKey&) { return sizeof(int64_t) + sizeof(int64_t) + sizeof(uint8_t); }

inline char* serialize(const EdgeKey& k, char* ptr) {
    ptr = ser::write_pod(ptr, k.lo);
    ptr = ser::write_pod(ptr, k.hi);
    uint8_t kind = static_cast<uint8_t>(k.kind);
    ptr = ser::write_pod(ptr, kind);
    return ptr;
}

inline const char* deserialize(EdgeKey& k, const char* ptr) {
    ptr = ser::read_pod(ptr, k.lo);
    ptr = ser::read_pod(ptr, k.hi);
    uint8_t kind = 0;
    ptr = ser::read_pod(ptr, kind);
    k.kind = static_cast<EdgeKind>(kind);
    return ptr;
}


template <typename DataType>
size_t bytes_pending(const PendingFactorUpdates<DataType>& p) {
    size_t n = 0;

    // counts
    n += sizeof(uint64_t); // replace_blocks count
    for (const auto& [k, b] : p.replace_blocks) {
        n += bytes_key(k);
        n += bytes_denseblock(b);
    }

    n += sizeof(uint64_t); // accumulated_deltas count
    for (const auto& [k, b] : p.accumulated_deltas) {
        n += bytes_key(k);
        n += bytes_denseblock(b);
    }

    return n;
}

template <typename DataType>
char* serialize(const PendingFactorUpdates<DataType>& p, char* ptr) {
    static_assert(std::is_trivially_copyable_v<DataType>, "PendingFactorUpdates: DataType must be trivially copyable");

    // replace_blocks
    {
        uint64_t count = static_cast<uint64_t>(p.replace_blocks.size());
        ptr = ser::write_pod(ptr, count);
        for (const auto& [k, b] : p.replace_blocks) {
            ptr = serialize(k, ptr);
            ptr = serialize(b, ptr);
        }
    }

    // accumulated_deltas
    {
        uint64_t count = static_cast<uint64_t>(p.accumulated_deltas.size());
        ptr = ser::write_pod(ptr, count);
        for (const auto& [k, b] : p.accumulated_deltas) {
            ptr = serialize(k, ptr);
            ptr = serialize(b, ptr);
        }
    }

    return ptr;
}

template <typename DataType>
const char* deserialize(PendingFactorUpdates<DataType>& p, const char* ptr) {
    static_assert(std::is_trivially_copyable_v<DataType>, "PendingFactorUpdates: DataType must be trivially copyable");

    p.replace_blocks.clear();
    p.accumulated_deltas.clear();

    // replace_blocks
    {
        uint64_t count = 0;
        ptr = ser::read_pod(ptr, count);
        p.replace_blocks.reserve(static_cast<size_t>(count));
        for (uint64_t i = 0; i < count; ++i) {
            ReplaceKey k{};
            DenseBlock<DataType> b{};
            ptr = deserialize(k, ptr);
            ptr = deserialize(b, ptr);
            p.replace_blocks.emplace(std::move(k), std::move(b));
        }
    }

    // accumulated_deltas
    {
        uint64_t count = 0;
        ptr = ser::read_pod(ptr, count);
        p.accumulated_deltas.reserve(static_cast<size_t>(count));
        for (uint64_t i = 0; i < count; ++i) {
            EdgeKey k{};
            DenseBlock<DataType> b{};
            ptr = deserialize(k, ptr);
            ptr = deserialize(b, ptr);
            p.accumulated_deltas.emplace(std::move(k), std::move(b));
        }
    }

    return ptr;
}

template<typename DataType>
size_t get_serialized_size(const MatrixStorage<DataType>& mat) {
    size_t size = 0;
    size += sizeof(int64_t);  // rows
    size += sizeof(int64_t);  // cols
    size += sizeof(int64_t);  // lda
    size += sizeof(typename MatrixStorage<DataType>::Format);  // format
    size += sizeof(size_t);   // data.size()
    size += mat.data.size() * sizeof(DataType);  // data elements
    return size;
}

template<typename DataType>
char* serialize(const MatrixStorage<DataType>& mat, char* buffer) {
    // Pack metadata
    std::memcpy(buffer, &mat.rows, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    std::memcpy(buffer, &mat.cols, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    std::memcpy(buffer, &mat.lda, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    std::memcpy(buffer, &mat.format, sizeof(typename MatrixStorage<DataType>::Format));
    buffer += sizeof(typename MatrixStorage<DataType>::Format);
    
    // Pack data vector
    size_t data_size = mat.data.size();
    std::memcpy(buffer, &data_size, sizeof(size_t));
    buffer += sizeof(size_t);
    
    if (data_size > 0) {
        std::memcpy(buffer, mat.data.data(), data_size * sizeof(DataType));
        buffer += data_size * sizeof(DataType);
    }
    
    return buffer;
}

template<typename DataType>
const char* deserialize(MatrixStorage<DataType>& mat, const char* buffer) {
    // Unpack metadata
    std::memcpy(&mat.rows, buffer, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    std::memcpy(&mat.cols, buffer, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    std::memcpy(&mat.lda, buffer, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    std::memcpy(&mat.format, buffer, sizeof(typename MatrixStorage<DataType>::Format));
    buffer += sizeof(typename MatrixStorage<DataType>::Format);
    
    // Unpack data vector
    size_t data_size;
    std::memcpy(&data_size, buffer, sizeof(size_t));
    buffer += sizeof(size_t);
    
    if (data_size > 0) {
        mat.data.resize(data_size);
        std::memcpy(mat.data.data(), buffer, data_size * sizeof(DataType));
        buffer += data_size * sizeof(DataType);
    } else {
        mat.data.clear();
    }
    
    return buffer;
}

template<typename DataType>
size_t get_serialized_size(const ModifiedBlock<DataType>& block) {
    size_t size = 0;
    size += sizeof(int64_t);  // neighbor_morton
    size += get_serialized_size(block.A_SN);
    size += get_serialized_size(block.A_NS);
    return size;
}

template<typename DataType>
char* serialize(const ModifiedBlock<DataType>& block, char* buffer) {
    std::memcpy(buffer, &block.neighbor_morton, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    buffer = serialize(block.A_SN, buffer);
    buffer = serialize(block.A_NS, buffer);
    
    return buffer;
}

template<typename DataType>
const char* deserialize(ModifiedBlock<DataType>& block, const char* buffer) {
    std::memcpy(&block.neighbor_morton, buffer, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    buffer = deserialize(block.A_SN, buffer);
    buffer = deserialize(block.A_NS, buffer);
    
    return buffer;
}

template<typename CoordType, typename DataType>
size_t get_serialized_size(const BoxData<CoordType, DataType>& box) {
    size_t size = 0;
    
    // ===== POD fields =====
    size += sizeof(int64_t);    // morton_index
    size += sizeof(int32_t);    // level
    size += sizeof(int32_t) * 3;  // grid_coords
    size += sizeof(CoordType) * 6;  // bounds
    size += sizeof(CoordType) * 3;  // center
    size += sizeof(CoordType);    // size
    size += sizeof(int64_t);    // parent_morton
    size += sizeof(int64_t) * 8;  // children_morton
    size += sizeof(int32_t);    // num_children
    size += sizeof(int64_t);    // num_points
    // REMOVED: num_redundant and num_skeleton
    size += sizeof(int64_t);    // num_near_field_interactions
    size += sizeof(int64_t);    // num_far_field_interactions
    size += sizeof(bool);    // on_boundary
    
    // ===== Vectors =====
    size += sizeof(size_t) + box.point_indices.size() * sizeof(int64_t);
    size += sizeof(size_t) + box.point_coords.size() * sizeof(CoordType);
    size += sizeof(size_t) + box.redundant_indices.size() * sizeof(int64_t);
    size += sizeof(size_t) + box.skeleton_indices.size() * sizeof(int64_t);
    size += sizeof(size_t) + box.one_hop.size() * sizeof(int64_t);
    size += sizeof(size_t) + box.use_full_set.size() * sizeof(int64_t);
    size += sizeof(size_t) + box.two_hop.size() * sizeof(int64_t);
    size += sizeof(size_t) + box.X_RR_pivots.size() * sizeof(int);
    
    // ===== MatrixStorage objects =====
    size += get_serialized_size(box.interpolation_matrix);
    size += get_serialized_size(box.X_RR);
    size += get_serialized_size(box.X_RS);
    size += get_serialized_size(box.X_SR);
    size += get_serialized_size(box.schur_complement);
    size += get_serialized_size(box.X_RN);
    size += get_serialized_size(box.X_NR);
    
    // ===== Interaction maps =====
    size += sizeof(size_t);  // near_field_interaction_map.size()
    size += box.near_field_interaction_map.size() * (sizeof(int64_t) + sizeof(int64_t));
    
    size += sizeof(size_t);  // far_field_interaction_map.size()
    size += box.far_field_interaction_map.size() * (sizeof(int64_t) + sizeof(int64_t));

    size += sizeof(size_t);  // near_field_interaction_map_nonsymmetry.size()
    size += box.near_field_interaction_map_nonsymmetry.size() * (sizeof(int64_t) + sizeof(int64_t));

    size += sizeof(size_t);  // near_field_interaction_map_nonsymmetry.size()
    size += box.far_field_interaction_map_nonsymmetry.size() * (sizeof(int64_t) + sizeof(int64_t));
    
    // ===== ModifiedBlock vectors =====
    size += sizeof(size_t);  // near_field_modified_interactions.size()
    for (const auto& block : box.near_field_modified_interactions) {
        size += get_serialized_size(block);
    }
    
    size += sizeof(size_t);  // far_field_modified_interactions.size()
    for (const auto& block : box.far_field_modified_interactions) {
        size += get_serialized_size(block);
    }
    
    return size;
}


template<typename CoordType, typename DataType>
char* serialize(const BoxData<CoordType, DataType>& box, char* buffer) {
    char* ptr = buffer;
    
    // ===== POD fields =====
    std::memcpy(ptr, &box.morton_index, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    std::memcpy(ptr, &box.level, sizeof(int32_t));
    ptr += sizeof(int32_t);
    
    std::memcpy(ptr, box.grid_coords, sizeof(int32_t) * 3);
    ptr += sizeof(int32_t) * 3;
    
    std::memcpy(ptr, box.bounds, sizeof(CoordType) * 6);
    ptr += sizeof(CoordType) * 6;
    
    std::memcpy(ptr, box.center, sizeof(CoordType) * 3);
    ptr += sizeof(CoordType) * 3;
    
    std::memcpy(ptr, &box.size, sizeof(CoordType));
    ptr += sizeof(CoordType);
    
    std::memcpy(ptr, &box.parent_morton, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    std::memcpy(ptr, box.children_morton, sizeof(int64_t) * 8);
    ptr += sizeof(int64_t) * 8;
    
    std::memcpy(ptr, &box.num_children, sizeof(int32_t));
    ptr += sizeof(int32_t);
    
    std::memcpy(ptr, &box.num_points, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    // REMOVED: num_redundant and num_skeleton serialization
    
    std::memcpy(ptr, &box.num_near_field_interactions, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    std::memcpy(ptr, &box.num_far_field_interactions, sizeof(int64_t));
    ptr += sizeof(int64_t);

    std::memcpy(ptr, &box.on_boundary, sizeof(bool));
    ptr += sizeof(bool);
    
    // ===== Helper lambda for vector serialization =====
    auto serialize_vector = [&ptr](const auto& vec) {
        size_t vec_size = vec.size();
        std::memcpy(ptr, &vec_size, sizeof(size_t));
        ptr += sizeof(size_t);
        
        if (vec_size > 0) {
            std::memcpy(ptr, vec.data(), vec_size * sizeof(typename std::decay_t<decltype(vec)>::value_type));
            ptr += vec_size * sizeof(typename std::decay_t<decltype(vec)>::value_type);
        }
    };
    
    // ===== Vectors =====
    serialize_vector(box.point_indices);
    serialize_vector(box.point_coords);
    serialize_vector(box.redundant_indices);
    serialize_vector(box.skeleton_indices);
    serialize_vector(box.one_hop);
    serialize_vector(box.use_full_set);
    serialize_vector(box.two_hop);
    serialize_vector(box.X_RR_pivots);
    
    // ===== MatrixStorage objects =====
    ptr = serialize(box.interpolation_matrix, ptr);
    ptr = serialize(box.X_RR, ptr);
    ptr = serialize(box.X_RS, ptr);
    ptr = serialize(box.X_SR, ptr);
    ptr = serialize(box.schur_complement, ptr);
    ptr = serialize(box.X_RN, ptr);
    ptr = serialize(box.X_NR, ptr);
    
    // ===== Interaction maps =====
    size_t near_map_size = box.near_field_interaction_map.size();
    std::memcpy(ptr, &near_map_size, sizeof(size_t));
    ptr += sizeof(size_t);
    
    for (const auto& [key, value] : box.near_field_interaction_map) {
        std::memcpy(ptr, &key, sizeof(int64_t));
        ptr += sizeof(int64_t);
        std::memcpy(ptr, &value, sizeof(int64_t));
        ptr += sizeof(int64_t);
    }

    size_t near_map_size_nonsymmetry = box.near_field_interaction_map_nonsymmetry.size();
    std::memcpy(ptr, &near_map_size_nonsymmetry, sizeof(size_t));
    ptr += sizeof(size_t);
    
    for (const auto& [key, value] : box.near_field_interaction_map_nonsymmetry) {
        std::memcpy(ptr, &key, sizeof(int64_t));
        ptr += sizeof(int64_t);
        std::memcpy(ptr, &value, sizeof(int64_t));
        ptr += sizeof(int64_t);
    }
    
    size_t far_map_size = box.far_field_interaction_map.size();
    std::memcpy(ptr, &far_map_size, sizeof(size_t));
    ptr += sizeof(size_t);
    
    for (const auto& [key, value] : box.far_field_interaction_map) {
        std::memcpy(ptr, &key, sizeof(int64_t));
        ptr += sizeof(int64_t);
        std::memcpy(ptr, &value, sizeof(int64_t));
        ptr += sizeof(int64_t);
    }

    size_t far_map_size_nonsymmetry = box.far_field_interaction_map_nonsymmetry.size();
    std::memcpy(ptr, &far_map_size_nonsymmetry, sizeof(size_t));
    ptr += sizeof(size_t);
    
    for (const auto& [key, value] : box.far_field_interaction_map_nonsymmetry) {
        std::memcpy(ptr, &key, sizeof(int64_t));
        ptr += sizeof(int64_t);
        std::memcpy(ptr, &value, sizeof(int64_t));
        ptr += sizeof(int64_t);
    }
    
    // ===== ModifiedBlock vectors =====
    size_t near_blocks_size = box.near_field_modified_interactions.size();
    std::memcpy(ptr, &near_blocks_size, sizeof(size_t));
    ptr += sizeof(size_t);
    
    for (const auto& block : box.near_field_modified_interactions) {
        ptr = serialize(block, ptr);
    }
    
    size_t far_blocks_size = box.far_field_modified_interactions.size();
    std::memcpy(ptr, &far_blocks_size, sizeof(size_t));
    ptr += sizeof(size_t);
    
    for (const auto& block : box.far_field_modified_interactions) {
        ptr = serialize(block, ptr);
    }
    
    return ptr;
}

template<typename CoordType, typename DataType>
const char* deserialize(BoxData<CoordType, DataType>& box, const char* buffer) {
    const char* ptr = buffer;
    
    // ===== POD fields =====
    std::memcpy(&box.morton_index, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    std::memcpy(&box.level, ptr, sizeof(int32_t));
    ptr += sizeof(int32_t);
    
    std::memcpy(box.grid_coords, ptr, sizeof(int32_t) * 3);
    ptr += sizeof(int32_t) * 3;
    
    std::memcpy(box.bounds, ptr, sizeof(CoordType) * 6);
    ptr += sizeof(CoordType) * 6;
    
    std::memcpy(box.center, ptr, sizeof(CoordType) * 3);
    ptr += sizeof(CoordType) * 3;
    
    std::memcpy(&box.size, ptr, sizeof(CoordType));
    ptr += sizeof(CoordType);
    
    std::memcpy(&box.parent_morton, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    std::memcpy(box.children_morton, ptr, sizeof(int64_t) * 8);
    ptr += sizeof(int64_t) * 8;
    
    std::memcpy(&box.num_children, ptr, sizeof(int32_t));
    ptr += sizeof(int32_t);
    
    std::memcpy(&box.num_points, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    // REMOVED: num_redundant and num_skeleton deserialization
    
    std::memcpy(&box.num_near_field_interactions, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    std::memcpy(&box.num_far_field_interactions, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);

    std::memcpy(&box.on_boundary, ptr, sizeof(bool));
    ptr += sizeof(bool);
    
    // ===== Helper lambda for vector deserialization =====
    auto deserialize_vector = [&ptr](auto& vec) {
        size_t vec_size;
        std::memcpy(&vec_size, ptr, sizeof(size_t));
        ptr += sizeof(size_t);
        
        if (vec_size > 0) {
            vec.resize(vec_size);
            std::memcpy(vec.data(), ptr, vec_size * sizeof(typename std::decay_t<decltype(vec)>::value_type));
            ptr += vec_size * sizeof(typename std::decay_t<decltype(vec)>::value_type);
        } else {
            vec.clear();
        }
    };
    
    // ===== Vectors =====
    deserialize_vector(box.point_indices);
    deserialize_vector(box.point_coords);
    deserialize_vector(box.redundant_indices);
    deserialize_vector(box.skeleton_indices);
    deserialize_vector(box.one_hop);
    deserialize_vector(box.use_full_set);
    deserialize_vector(box.two_hop);
    deserialize_vector(box.X_RR_pivots);
    
    // ===== MatrixStorage objects =====
    ptr = deserialize(box.interpolation_matrix, ptr);
    ptr = deserialize(box.X_RR, ptr);
    ptr = deserialize(box.X_RS, ptr);
    ptr = deserialize(box.X_SR, ptr);
    ptr = deserialize(box.schur_complement, ptr);
    ptr = deserialize(box.X_RN, ptr);
    ptr = deserialize(box.X_NR, ptr);
    
    // ===== Interaction maps =====
    size_t near_map_size;
    std::memcpy(&near_map_size, ptr, sizeof(size_t));
    ptr += sizeof(size_t);
    
    box.near_field_interaction_map.clear();
    for (size_t i = 0; i < near_map_size; ++i) {
        int64_t key, value;
        std::memcpy(&key, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        std::memcpy(&value, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        box.near_field_interaction_map[key] = value;
    }

    size_t near_map_size_nonsymmetry;
    std::memcpy(&near_map_size_nonsymmetry, ptr, sizeof(size_t));
    ptr += sizeof(size_t);
    
    box.near_field_interaction_map_nonsymmetry.clear();
    for (size_t i = 0; i < near_map_size_nonsymmetry; ++i) {
        int64_t key, value;
        std::memcpy(&key, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        std::memcpy(&value, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        box.near_field_interaction_map_nonsymmetry[key] = value;
    }
    
    size_t far_map_size;
    std::memcpy(&far_map_size, ptr, sizeof(size_t));
    ptr += sizeof(size_t);
    
    box.far_field_interaction_map.clear();
    for (size_t i = 0; i < far_map_size; ++i) {
        int64_t key, value;
        std::memcpy(&key, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        std::memcpy(&value, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        box.far_field_interaction_map[key] = value;
    }

    size_t far_map_size_nonsymmetry;
    std::memcpy(&far_map_size_nonsymmetry, ptr, sizeof(size_t));
    ptr += sizeof(size_t);
    
    box.far_field_interaction_map_nonsymmetry.clear();
    for (size_t i = 0; i < far_map_size_nonsymmetry; ++i) {
        int64_t key, value;
        std::memcpy(&key, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        std::memcpy(&value, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        box.far_field_interaction_map_nonsymmetry[key] = value;
    }
    
    // ===== ModifiedBlock vectors =====
    size_t near_blocks_size;
    std::memcpy(&near_blocks_size, ptr, sizeof(size_t));
    ptr += sizeof(size_t);
    
    box.near_field_modified_interactions.resize(near_blocks_size);
    for (auto& block : box.near_field_modified_interactions) {
        ptr = deserialize(block, ptr);
    }
    
    size_t far_blocks_size;
    std::memcpy(&far_blocks_size, ptr, sizeof(size_t));
    ptr += sizeof(size_t);
    
    box.far_field_modified_interactions.resize(far_blocks_size);
    for (auto& block : box.far_field_modified_interactions) {
        ptr = deserialize(block, ptr);
    }
    
    return ptr;
}

/**
 * @brief Get serialized size of PointDataRequest
 */
template<typename CoordType>
size_t get_serialized_size(const PointDataRequest<CoordType>& request) {
    size_t size = 0;
    size += sizeof(int64_t);  // morton_index
    size += sizeof(int);      // source_rank
    size += sizeof(size_t);   // coords.size()
    size += request.coords.size() * sizeof(CoordType);
    size += sizeof(size_t);   // skel_indices.size()
    size += request.skel_indices.size() * sizeof(int64_t);
    size += sizeof(bool);     // on_boundary
    return size;
}

/**
 * @brief Serialize PointDataRequest into buffer
 * @return Pointer to position after serialized data
 */
template<typename CoordType>
char* serialize(const PointDataRequest<CoordType>& request, char* buffer) {
    std::memcpy(buffer, &request.morton_index, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    std::memcpy(buffer, &request.source_rank, sizeof(int));
    buffer += sizeof(int);
    
    size_t coords_size = request.coords.size();
    std::memcpy(buffer, &coords_size, sizeof(size_t));
    buffer += sizeof(size_t);
    
    if (coords_size > 0) {
        std::memcpy(buffer, request.coords.data(), coords_size * sizeof(CoordType));
        buffer += coords_size * sizeof(CoordType);
    }

    size_t skel_size = request.skel_indices.size();
    std::memcpy(buffer, &skel_size, sizeof(size_t));
    buffer += sizeof(size_t);

    if (skel_size > 0) {
        std::memcpy(buffer, request.skel_indices.data(), skel_size * sizeof(int64_t));
        buffer += skel_size * sizeof(int64_t);
    }

    std::memcpy(buffer, &request.on_boundary, sizeof(bool));
    buffer += sizeof(bool);
    
    return buffer;
}

/**
 * @brief Deserialize PointDataRequest from buffer
 * @return Pointer to position after deserialized data
 */
template<typename CoordType>
const char* deserialize(PointDataRequest<CoordType>& request, const char* buffer) {
    std::memcpy(&request.morton_index, buffer, sizeof(int64_t));
    buffer += sizeof(int64_t);
    
    std::memcpy(&request.source_rank, buffer, sizeof(int));
    buffer += sizeof(int);
    
    size_t coords_size;
    std::memcpy(&coords_size, buffer, sizeof(size_t));
    buffer += sizeof(size_t);
    
    if (coords_size > 0) {
        request.coords.resize(coords_size);
        std::memcpy(request.coords.data(), buffer, coords_size * sizeof(CoordType));
        buffer += coords_size * sizeof(CoordType);
    } else {
        request.coords.clear();
    }

    size_t skel_size;
    std::memcpy(&skel_size, buffer, sizeof(size_t));
    buffer += sizeof(size_t);
    if (skel_size > 0) {
        request.skel_indices.resize(skel_size);
        std::memcpy(request.skel_indices.data(), buffer, skel_size * sizeof(int64_t));
        buffer += skel_size * sizeof(int64_t);
    } else {
        request.skel_indices.clear();
    }

    std::memcpy(&request.on_boundary, buffer, sizeof(bool));
    buffer += sizeof(bool);
    
    return buffer;
}


/**
 * @brief Get serialized size of SolveDataRequest
 */
template<typename CoordType, typename DataType>
size_t get_serialized_size(const SolveDataRequest<CoordType, DataType>& request) {
    size_t size = 0;
    
    // POD fields
    size += sizeof(int64_t);  // morton_index
    size += sizeof(int);      // source_rank
    
    // Vector fields with DataType elements
    size += sizeof(size_t) + request.right_side.size() * sizeof(DataType);
    size += sizeof(size_t) + request.left_side.size() * sizeof(DataType);
    
    // Vector fields with int64_t elements
    size += sizeof(size_t) + request.redundant_indices.size() * sizeof(int64_t);
    size += sizeof(size_t) + request.skeleton_indices.size() * sizeof(int64_t);
    size += sizeof(size_t) + request.one_hop.size() * sizeof(int64_t);
    size += sizeof(size_t) + request.use_full_set.size() * sizeof(int64_t);
    size += sizeof(size_t) + request.X_RR_pivots.size() * sizeof(int);
    
    // MatrixStorage fields
    size += get_serialized_size(request.interpolation_matrix);
    size += get_serialized_size(request.X_RR);
    size += get_serialized_size(request.X_RS);
    size += get_serialized_size(request.X_SR);
    size += get_serialized_size(request.schur_complement);
    size += get_serialized_size(request.X_RN);
    size += get_serialized_size(request.X_NR);
    
    return size;
}

/**
 * @brief Serialize SolveDataRequest into buffer
 * @return Pointer to position after serialized data
 */
template<typename CoordType, typename DataType>
char* serialize(const SolveDataRequest<CoordType, DataType>& request, char* buffer) {
    char* ptr = buffer;
    
    // POD fields
    std::memcpy(ptr, &request.morton_index, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    std::memcpy(ptr, &request.source_rank, sizeof(int));
    ptr += sizeof(int);
    
    // Helper lambda for vector serialization
    auto serialize_vector = [&ptr](const auto& vec) {
        size_t vec_size = vec.size();
        std::memcpy(ptr, &vec_size, sizeof(size_t));
        ptr += sizeof(size_t);
        
        if (vec_size > 0) {
            using ValueType = typename std::decay_t<decltype(vec)>::value_type;
            std::memcpy(ptr, vec.data(), vec_size * sizeof(ValueType));
            ptr += vec_size * sizeof(ValueType);
        }
    };
    
    // Serialize vectors
    serialize_vector(request.right_side);
    serialize_vector(request.left_side);
    serialize_vector(request.redundant_indices);
    serialize_vector(request.skeleton_indices);
    serialize_vector(request.one_hop);
    serialize_vector(request.use_full_set);
    serialize_vector(request.X_RR_pivots);
    
    // Serialize MatrixStorage objects
    ptr = serialize(request.interpolation_matrix, ptr);
    ptr = serialize(request.X_RR, ptr);
    ptr = serialize(request.X_RS, ptr);
    ptr = serialize(request.X_SR, ptr);
    ptr = serialize(request.schur_complement, ptr);
    ptr = serialize(request.X_RN, ptr);
    ptr = serialize(request.X_NR, ptr);
    
    return ptr;
}

/**
 * @brief Deserialize SolveDataRequest from buffer
 * @return Pointer to position after deserialized data
 */
template<typename CoordType, typename DataType>
const char* deserialize(SolveDataRequest<CoordType, DataType>& request, const char* buffer) {
    const char* ptr = buffer;
    
    // POD fields
    std::memcpy(&request.morton_index, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    std::memcpy(&request.source_rank, ptr, sizeof(int));
    ptr += sizeof(int);
    
    // Helper lambda for vector deserialization
    auto deserialize_vector = [&ptr](auto& vec) {
        size_t vec_size;
        std::memcpy(&vec_size, ptr, sizeof(size_t));
        ptr += sizeof(size_t);
        
        if (vec_size > 0) {
            vec.resize(vec_size);
            using ValueType = typename std::decay_t<decltype(vec)>::value_type;
            std::memcpy(vec.data(), ptr, vec_size * sizeof(ValueType));
            ptr += vec_size * sizeof(ValueType);
        } else {
            vec.clear();
        }
    };
    
    // Deserialize vectors
    deserialize_vector(request.right_side);
    deserialize_vector(request.left_side);
    deserialize_vector(request.redundant_indices);
    deserialize_vector(request.skeleton_indices);
    deserialize_vector(request.one_hop);
    deserialize_vector(request.use_full_set);
    deserialize_vector(request.X_RR_pivots);
    
    // Deserialize MatrixStorage objects
    ptr = deserialize(request.interpolation_matrix, ptr);
    ptr = deserialize(request.X_RR, ptr);
    ptr = deserialize(request.X_RS, ptr);
    ptr = deserialize(request.X_SR, ptr);
    ptr = deserialize(request.schur_complement, ptr);
    ptr = deserialize(request.X_RN, ptr);
    ptr = deserialize(request.X_NR, ptr);
    
    return ptr;
}


/**
 * @brief MPI_Send wrapper that handles message sizes exceeding INT_MAX
 * 
 * Splits large messages into chunks of INT_MAX bytes and sends them sequentially.
 * For messages <= INT_MAX, behaves identically to MPI_Send.
 */
inline int MPI_Send_large(const void* buf, size_t count, MPI_Datatype datatype,
                          int dest, int tag, MPI_Comm comm) {
    const char* ptr = static_cast<const char*>(buf);
    while (count > 0) {
        int chunk = (count > static_cast<size_t>(INT_MAX)) ? INT_MAX : static_cast<int>(count);
        int err = MPI_Send(ptr, chunk, datatype, dest, tag, comm);
        if (err != MPI_SUCCESS) return err;
        ptr += chunk;
        count -= chunk;
    }
    return MPI_SUCCESS;
}

/**
 * @brief MPI_Isend wrapper that handles message sizes exceeding INT_MAX
 * 
 * For messages <= INT_MAX, posts a single Isend.
 * For larger messages, sends all but the last chunk with blocking MPI_Send,
 * then posts an Isend for the final chunk so the caller gets a valid request.
 */
inline int MPI_Isend_large(const void* buf, size_t count, MPI_Datatype datatype,
                           int dest, int tag, MPI_Comm comm,
                           std::vector<MPI_Request>& requests) {
    const char* ptr = static_cast<const char*>(buf);
    while (count > 0) {
        int chunk = std::min(count, (size_t)INT_MAX);
        MPI_Request req;
        int err = MPI_Isend(ptr, chunk, datatype, dest, tag, comm, &req);
        if (err != MPI_SUCCESS) return err;
        requests.push_back(req);
        ptr += chunk;
        count -= chunk;
    }
    return MPI_SUCCESS;
}

/**
 * @brief MPI_Irecv wrapper that handles message sizes exceeding INT_MAX
 * 
 * For messages <= INT_MAX, posts a single Irecv.
 * For larger messages, receives all but the last chunk with blocking MPI_Recv,
 * then posts an Irecv for the final chunk so the caller gets a valid request.
 */
inline int MPI_Irecv_large(void* buf, size_t count, MPI_Datatype datatype,
                           int source, int tag, MPI_Comm comm,
                           std::vector<MPI_Request>& requests) {
    char* ptr = static_cast<char*>(buf);
    while (count > 0) {
        int chunk = std::min(count, (size_t)INT_MAX);
        MPI_Request req;
        int err = MPI_Irecv(ptr, chunk, datatype, source, tag, comm, &req);
        if (err != MPI_SUCCESS) return err;
        requests.push_back(req);
        ptr += chunk;
        count -= chunk;
    }
    return MPI_SUCCESS;
}

/**
 * @brief MPI_Sendrecv wrapper that handles message sizes exceeding INT_MAX.
 *
 * Sends and receives the message in lockstep chunks so callers can safely use
 * deterministic pairwise exchanges without many-peer progress assumptions.
 */
inline int MPI_Sendrecv_large(const void* sendbuf, size_t sendcount, MPI_Datatype sendtype,
                              int dest, int sendtag,
                              void* recvbuf, size_t recvcount, MPI_Datatype recvtype,
                              int source, int recvtag,
                              MPI_Comm comm, MPI_Status* status) {
    const char* send_ptr = static_cast<const char*>(sendbuf);
    char* recv_ptr = static_cast<char*>(recvbuf);
    MPI_Status local_status;
    MPI_Status* active_status = (status == nullptr || status == MPI_STATUS_IGNORE)
        ? &local_status : status;

    uint64_t local_sizes[2] = {
        static_cast<uint64_t>(sendcount),
        static_cast<uint64_t>(recvcount)
    };
    uint64_t remote_sizes[2] = {0, 0};
    int size_err = MPI_Sendrecv(local_sizes,
                                2,
                                MPI_UINT64_T,
                                dest,
                                sendtag,
                                remote_sizes,
                                2,
                                MPI_UINT64_T,
                                source,
                                recvtag,
                                comm,
                                MPI_STATUS_IGNORE);
    if (size_err != MPI_SUCCESS) return size_err;
    if (remote_sizes[0] != local_sizes[1] || remote_sizes[1] != local_sizes[0]) {
        int rank = -1;
        MPI_Comm_rank(comm, &rank);
        std::fprintf(stderr,
                     "Rank %d <-> %d: local_send=%zu local_recv=%zu remote_send=%zu remote_recv=%zu\n",
                     rank,
                     dest,
                     sendcount,
                     recvcount,
                     static_cast<size_t>(remote_sizes[0]),
                     static_cast<size_t>(remote_sizes[1]));
        std::fflush(stderr);
        return MPI_ERR_TRUNCATE;
    }

    while (sendcount > 0 || recvcount > 0) {
        const int send_chunk = static_cast<int>(std::min(sendcount, static_cast<size_t>(INT_MAX)));
        const int recv_chunk = static_cast<int>(std::min(recvcount, static_cast<size_t>(INT_MAX)));
        const void* active_send_ptr = send_chunk > 0 ? send_ptr : nullptr;
        void* active_recv_ptr = recv_chunk > 0 ? recv_ptr : nullptr;
        int err = MPI_Sendrecv(active_send_ptr, send_chunk, sendtype, dest, sendtag,
                               active_recv_ptr, recv_chunk, recvtype, source, recvtag,
                               comm, active_status);
        if (err != MPI_SUCCESS) return err;
        if (send_chunk > 0) {
            send_ptr += send_chunk;
            sendcount -= static_cast<size_t>(send_chunk);
        }
        if (recv_chunk > 0) {
            recv_ptr += recv_chunk;
            recvcount -= static_cast<size_t>(recv_chunk);
        }
    }
    return MPI_SUCCESS;
}

/**
 * @brief MPI_Recv wrapper that handles message sizes exceeding INT_MAX
 * 
 * Splits large messages into chunks of INT_MAX bytes and receives them sequentially.
 * For messages <= INT_MAX, behaves identically to MPI_Recv.
 */
inline int MPI_Recv_large(void* buf, size_t count, MPI_Datatype datatype,
                          int source, int tag, MPI_Comm comm, MPI_Status* status) {
    char* ptr = static_cast<char*>(buf);
    while (count > 0) {
        int chunk = (count > static_cast<size_t>(INT_MAX)) ? INT_MAX : static_cast<int>(count);
        int err = MPI_Recv(ptr, chunk, datatype, source, tag, comm, status);
        if (err != MPI_SUCCESS) return err;
        ptr += chunk;
        count -= chunk;
    }
    return MPI_SUCCESS;
}

/**
 * @brief Check whether two active process regions are spatial neighbors.
 */
inline bool process_morton_regions_are_neighbors(int lhs_morton_id,
                                                 int rhs_morton_id,
                                                 int dimension) {
    if (lhs_morton_id == rhs_morton_id) {
        return false;
    }

    if (dimension == 2) {
        uint32_t lx = 0, ly = 0, rx = 0, ry = 0;
        morton::decode_2d(static_cast<uint64_t>(lhs_morton_id), lx, ly);
        morton::decode_2d(static_cast<uint64_t>(rhs_morton_id), rx, ry);
        return std::abs(static_cast<int>(lx) - static_cast<int>(rx)) <= 1 &&
               std::abs(static_cast<int>(ly) - static_cast<int>(ry)) <= 1;
    }

    if (dimension == 3) {
        uint32_t lx = 0, ly = 0, lz = 0, rx = 0, ry = 0, rz = 0;
        morton::decode_3d(static_cast<uint64_t>(lhs_morton_id), lx, ly, lz);
        morton::decode_3d(static_cast<uint64_t>(rhs_morton_id), rx, ry, rz);
        return std::abs(static_cast<int>(lx) - static_cast<int>(rx)) <= 1 &&
               std::abs(static_cast<int>(ly) - static_cast<int>(ry)) <= 1 &&
               std::abs(static_cast<int>(lz) - static_cast<int>(rz)) <= 1;
    }

    throw std::runtime_error(
        "process_morton_regions_are_neighbors: unsupported dimension " +
        std::to_string(dimension));
}

/**
 * @brief Build and cache a deterministic one-peer-per-phase schedule.
 */
template<typename CoordType, typename DataType>
ProcessPairCommunicationSchedule& get_or_build_process_pair_communication_schedule(
    TreeLevel<CoordType, DataType>& lvl) {

    if (lvl.process_pair_schedule.initialized) {
        return lvl.process_pair_schedule;
    }

    auto& schedule = lvl.process_pair_schedule;
    schedule.reset();

    std::vector<int> active_ranks;
    active_ranks.reserve(lvl.rank_to_morton.size());
    for (const auto& entry : lvl.rank_to_morton) {
        active_ranks.push_back(entry.first);
    }
    std::sort(active_ranks.begin(), active_ranks.end());

    std::vector<std::pair<int, int>> remaining_edges;
    for (size_t i = 0; i < active_ranks.size(); ++i) {
        for (size_t j = i + 1; j < active_ranks.size(); ++j) {
            const int lhs_rank = active_ranks[i];
            const int rhs_rank = active_ranks[j];
            const int lhs_morton = lvl.rank_to_morton.at(lhs_rank);
            const int rhs_morton = lvl.rank_to_morton.at(rhs_rank);
            if (process_morton_regions_are_neighbors(lhs_morton, rhs_morton, lvl.dimension)) {
                remaining_edges.emplace_back(lhs_rank, rhs_rank);
            }
        }
    }

    while (!remaining_edges.empty()) {
        std::unordered_set<int> used_ranks;
        std::vector<std::pair<int, int>> phase_edges;
        std::vector<std::pair<int, int>> deferred_edges;
        phase_edges.reserve(remaining_edges.size());
        deferred_edges.reserve(remaining_edges.size());

        for (const auto& edge : remaining_edges) {
            if (used_ranks.find(edge.first) != used_ranks.end() ||
                used_ranks.find(edge.second) != used_ranks.end()) {
                deferred_edges.push_back(edge);
                continue;
            }
            phase_edges.push_back(edge);
            used_ranks.insert(edge.first);
            used_ranks.insert(edge.second);
        }

        if (phase_edges.empty()) {
            throw std::runtime_error(
                "get_or_build_process_pair_communication_schedule: failed to build a non-empty phase");
        }

        schedule.phase_edges.push_back(std::move(phase_edges));
        remaining_edges.swap(deferred_edges);
    }

    int my_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    schedule.local_phase_peers.assign(schedule.phase_edges.size(), -1);
    for (size_t phase_idx = 0; phase_idx < schedule.phase_edges.size(); ++phase_idx) {
        for (const auto& edge : schedule.phase_edges[phase_idx]) {
            int peer = -1;
            if (edge.first == my_rank) peer = edge.second;
            else if (edge.second == my_rank) peer = edge.first;
            else continue;
            if (schedule.local_phase_peers[phase_idx] != -1) {
                throw std::runtime_error(
                    "get_or_build_process_pair_communication_schedule: rank " +
                    std::to_string(my_rank) + " appears multiple times in phase " +
                    std::to_string(phase_idx));
            }
            schedule.local_phase_peers[phase_idx] = peer;
        }
    }

    schedule.initialized = true;
    return schedule;
}

// ---------------- small math helpers ----------------
template <class T>
static inline T maybe_conj(const T& x, bool) { return x; }

template <class T>
static inline std::complex<T> maybe_conj(const std::complex<T>& x, bool do_conj) {
    return do_conj ? std::conj(x) : x;
}

// column-major transpose (rows x cols) -> (cols x rows)
template <class DataType>
static std::vector<DataType> transpose_colmajor(
    const std::vector<DataType>& A, int64_t rows, int64_t cols, bool do_conj)
{
    if (A.size() != static_cast<size_t>(rows * cols))
        throw std::runtime_error("transpose_colmajor: size mismatch");

    std::vector<DataType> AT(static_cast<size_t>(cols * rows));
    for (int64_t j = 0; j < cols; ++j)
        for (int64_t i = 0; i < rows; ++i)
            AT[static_cast<size_t>(j + i * cols)] =
                maybe_conj(A[static_cast<size_t>(i + j * rows)], do_conj);
    return AT;
}
#include <cassert>
// accumulate DenseBlock (ADD)
template <class DataType>
static void accumulate_denseblock(DenseBlock<DataType>& dst,
                                  const DenseBlock<DataType>& add)
{
    if (dst.data.empty()) {
        dst.rows = add.rows;
        dst.cols = add.cols;
        dst.data = add.data;
        return;
    }
    if (dst.rows != add.rows || dst.cols != add.cols || dst.data.size() != add.data.size())
        throw std::runtime_error("accumulate_denseblock: dimension mismatch");
    
    for (size_t i = 0; i < dst.data.size(); ++i) dst.data[i] += add.data[i];
}

// merge received packets into one (so we can gather assisting first)
template <class DataType>
static void merge_pending(PendingFactorUpdates<DataType>& dst,
                          const PendingFactorUpdates<DataType>& src)
{
    // Step 7 REPLACE: last-wins overwrite is fine (should be unique anyway)
    for (const auto& [k, b] : src.replace_blocks) dst.replace_blocks[k] = b;

    // Step 9 ADD: sum if same EdgeKey appears from multiple senders
    for (const auto& [k, b] : src.accumulated_deltas) {
        auto& out = dst.accumulated_deltas[k];
        DenseBlock<DataType> tmp = b;
        accumulate_denseblock(out, tmp);
    }
}



template <typename CoordType, typename DataType>
static ModifiedBlock<DataType>& upsert_block_symmetric(BoxData<CoordType, DataType>&, int64_t, EdgeKind); // not used

template <typename CoordType, typename DataType>
static ModifiedBlock<DataType>& upsert_block_symmetric(
    BoxData<CoordType, DataType>& box, int64_t neighbor_morton, EdgeKind kind)
{
    if (kind == EdgeKind::Near) {
        auto it = box.near_field_interaction_map.find(neighbor_morton);
        if (it != box.near_field_interaction_map.end())
            return box.near_field_modified_interactions[static_cast<size_t>(it->second)];

        ModifiedBlock<DataType> nb; nb.neighbor_morton = neighbor_morton;
        const int64_t new_idx = static_cast<int64_t>(box.near_field_modified_interactions.size());
        box.near_field_modified_interactions.push_back(std::move(nb));
        box.near_field_interaction_map[neighbor_morton] = new_idx;
        box.num_near_field_interactions = static_cast<int64_t>(box.near_field_modified_interactions.size());
        return box.near_field_modified_interactions.back();
    }

    if (kind == EdgeKind::Far) {
        auto it = box.far_field_interaction_map.find(neighbor_morton);
        if (it != box.far_field_interaction_map.end())
            return box.far_field_modified_interactions[static_cast<size_t>(it->second)];

        ModifiedBlock<DataType> nb; nb.neighbor_morton = neighbor_morton;
        const int64_t new_idx = static_cast<int64_t>(box.far_field_modified_interactions.size());
        box.far_field_modified_interactions.push_back(std::move(nb));
        box.far_field_interaction_map[neighbor_morton] = new_idx;
        box.num_far_field_interactions = static_cast<int64_t>(box.far_field_modified_interactions.size());
        return box.far_field_modified_interactions.back();
    }

    throw std::runtime_error("upsert_block_symmetric: kind must be Near/Far for off-diagonal");
}

template <typename CoordType, typename DataType>
static ModifiedBlock<DataType>* find_block_symmetric(
    BoxData<CoordType, DataType>& box, int64_t neighbor_morton, EdgeKind kind)
{
    if (kind == EdgeKind::Near) {
        auto it = box.near_field_interaction_map.find(neighbor_morton);
        if (it == box.near_field_interaction_map.end()) return nullptr;
        return &box.near_field_modified_interactions[static_cast<size_t>(it->second)];
    }
    if (kind == EdgeKind::Far) {
        auto it = box.far_field_interaction_map.find(neighbor_morton);
        if (it == box.far_field_interaction_map.end()) return nullptr;
        return &box.far_field_modified_interactions[static_cast<size_t>(it->second)];
    }
    return nullptr;
}

// Extract exactly n points worth of coords (point-major: [x0 y0 z0 x1 y1 z1 ...])
template <typename CoordType, typename DataType>
static std::vector<CoordType> coords_for_count_from_box(
    const BoxData<CoordType, DataType>& b, int64_t n_needed, int dim)
{
    if (n_needed < 0) throw std::runtime_error("coords_for_count_from_box: n_needed < 0");

    // Use skeleton indices if they match
    if (!b.skeleton_indices.empty() && static_cast<int64_t>(b.skeleton_indices.size()) == n_needed) {
        std::vector<CoordType> out(static_cast<size_t>(n_needed * dim));
        for (int64_t i = 0; i < n_needed; ++i) {
            const int64_t src = b.skeleton_indices[static_cast<size_t>(i)];
            for (int d = 0; d < dim; ++d)
                out[static_cast<size_t>(i * dim + d)] =
                    b.point_coords[static_cast<size_t>(src * dim + d)];
        }
        return out;
    }

    // Otherwise, require full coords to match
    if (b.num_points == n_needed) {
        return b.point_coords;
    }

    throw std::runtime_error("coords_for_count_from_box: cannot produce requested count=" +
                             std::to_string(n_needed) + " from box morton=" +
                             std::to_string(b.morton_index) +
                             " (num_points=" + std::to_string(b.num_points) +
                             ", skel_size=" + std::to_string(b.skeleton_indices.size()) + ")");
}

template <typename CoordType>
static std::vector<CoordType> coords_for_count_from_assist(
    const PointDataRequest<CoordType>& a, int64_t n_needed, int dim)
{
    // prefer skeleton indices if present and matching
    if (!a.skel_indices.empty() && static_cast<int64_t>(a.skel_indices.size()) == n_needed) {
        // a.coords is full coords in point-major; extract skeleton coords
        std::vector<CoordType> out(static_cast<size_t>(n_needed * dim));
        for (int64_t i = 0; i < n_needed; ++i) {
            const int64_t src = a.skel_indices[static_cast<size_t>(i)];
            for (int d = 0; d < dim; ++d)
                out[static_cast<size_t>(i * dim + d)] =
                    a.coords[static_cast<size_t>(src * dim + d)];
        }
        return out;
    }

    // else, require full coords length to match
    const int64_t n_full = static_cast<int64_t>(a.coords.size()) / dim;
    if (n_full == n_needed) return a.coords;

    throw std::runtime_error("coords_for_count_from_assist: cannot produce requested count=" +
                             std::to_string(n_needed) + " from assist morton=" +
                             std::to_string(a.morton_index) +
                             " (full_n=" + std::to_string(n_full) +
                             ", skel_size=" + std::to_string(a.skel_indices.size()) + ")");
}

// upsert assisting box entry in lvl.assisting_boxes + map
template <typename CoordType, typename DataType>
static void upsert_assisting(TreeLevel<CoordType, DataType>& lvl,
                             PointDataRequest<CoordType>&& req)
{
    auto it = lvl.assisting_box_points_for_kernel_evaluation.find(req.morton_index);
    if (it == lvl.assisting_box_points_for_kernel_evaluation.end()) {
        const int64_t idx = static_cast<int64_t>(lvl.assisting_boxes.size());
        lvl.assisting_boxes.push_back(std::move(req));
        lvl.assisting_box_points_for_kernel_evaluation[req.morton_index] = idx;
    } else {
        lvl.assisting_boxes[static_cast<size_t>(it->second)] = std::move(req);
    }
}


template <typename CoordType, typename DataType>
static int owner_rank_of_morton(
    ParallelTree<CoordType, DataType>* tree,
    TreeLevel<CoordType, DataType>& lvl,
    int level,
    int64_t morton)
{
    const uint32_t grid_size = 1u << level;
    std::vector<uint64_t> single = { static_cast<uint64_t>(morton) };
    std::vector<uint32_t> region =
        (tree->dimension == 2)
        ? morton::assign_to_processes_2d(single, lvl.num_active_processes, grid_size)
        : morton::assign_to_processes_3d(single, lvl.num_active_processes, grid_size);
    const int region_id = static_cast<int>(region[0]);
    return lvl.morton_to_rank.at(region_id);
}

template <typename CoordType, typename DataType>
std::vector<int> compute_one_hop_neighbor_ranks(
    ParallelTree<CoordType, DataType>* tree,
    const TreeLevel<CoordType, DataType>& lvl,
    int level)
{
    if (!lvl.is_process_active) return {};
    if (lvl.my_morton_id < 0) throw std::runtime_error("my_morton_id < 0 on active process");

    const uint32_t P = static_cast<uint32_t>(lvl.num_active_processes);

    uint32_t k = 0;
    uint32_t procs_per_dim = 0;
    if (tree->dimension == 2) {
        // P = 4^k => procs_per_dim = 2^k
        k = __builtin_ctz(P) / 2;
        procs_per_dim = 1u << k;
    } else {
        // P = 8^k => procs_per_dim = 2^k
        k = __builtin_ctz(P) / 3;
        procs_per_dim = 1u << k;
    }

    const uint64_t my_region = static_cast<uint64_t>(lvl.my_morton_id);

    std::vector<uint64_t> neigh_regions =
        (tree->dimension == 2)
        ? morton::neighbors_2d(my_region, procs_per_dim)      // 8 neighbors
        : morton::neighbors_3d(my_region, procs_per_dim);     // 26 neighbors

    std::vector<int> neigh_ranks;
    neigh_ranks.reserve(neigh_regions.size());

    for (uint64_t rid64 : neigh_regions) {
        const int rid = static_cast<int>(rid64);
        auto it = lvl.morton_to_rank.find(rid);
        if (it != lvl.morton_to_rank.end()) {
            neigh_ranks.push_back(it->second);
        }
    }

    std::sort(neigh_ranks.begin(), neigh_ranks.end());
    neigh_ranks.erase(std::unique(neigh_ranks.begin(), neigh_ranks.end()), neigh_ranks.end());
    return neigh_ranks;
}


/**
 * @brief Request and receive assisting box point data from 1-hop neighboring ranks.
 *
 * This routine is used to support on-demand kernel evaluation of missing interaction blocks
 * during update application. If a local endpoint needs to evaluate a base interaction block
 * involving a remote endpoint box, it must have access to the remote box's point coordinates
 * and (optionally) skeleton indices.
 *
 * Protocol (neighbor-only point-to-point; no collectives):
 *  1) For each needed remote morton index, compute its owning rank and group requests per owner.
 *     (Enforces that each owner is within the provided 1-hop neighbor set; otherwise throws.)
 *  2) Exchange request counts with each neighbor (tag 600).
 *  3) Exchange requested morton-index lists (tag 601).
 *  4) For each neighbor's requests, pack a response buffer containing one or more serialized
 *     PointDataRequest<CoordType> records, prefixed by per-record byte size.
 *     Exchange total response sizes (tag 602).
 *  5) Exchange response payload buffers (tag 603).
 *  6) Deserialize responses and populate:
 *      - lvl.assisting_boxes
 *      - lvl.assisting_box_points_for_kernel_evaluation (morton -> index in assisting_boxes)
 *
 * Contents returned per assisting box:
 *  - full point coordinates (BoxData::point_coords)
 *  - skeleton indices (BoxData::skeleton_indices), enabling skeleton-sized kernel evaluation
 *  - on_boundary flag
 *
 * Error handling:
 *  - Throws if a requested morton index is not locally owned when serving requests.
 *  - Throws if a needed owner is not in the 1-hop neighbor set (assumption violation).
 *  - Throws on buffer parse mismatches.
 *
 * @return Communication-only duration spent in MPI send/recv/wait operations.
 */
template <typename CoordType, typename DataType>
std::chrono::high_resolution_clock::duration exchange_assisting_for_mortons_onehop(
    ParallelTree<CoordType, DataType>* tree,
    TreeLevel<CoordType, DataType>& lvl,
    int level,
    const std::vector<int>& neighbor_ranks,
    const std::vector<int64_t>& needed_remote_mortons)
{
    using clock = std::chrono::high_resolution_clock;
    clock::duration communication_time{};

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Step 1: Group requested assisting mortons by owner rank.
    std::unordered_map<int, std::vector<int64_t>> req_to_send;
    req_to_send.reserve(neighbor_ranks.size());

    const uint32_t grid_size = 1u << level;

    auto owner_of_morton = [&](int64_t morton_idx) -> int {
        std::vector<uint64_t> single = { static_cast<uint64_t>(morton_idx) };
        std::vector<uint32_t> region =
            (tree->dimension == 2)
            ? morton::assign_to_processes_2d(single, lvl.num_active_processes, grid_size)
            : morton::assign_to_processes_3d(single, lvl.num_active_processes, grid_size);
        const int rid = static_cast<int>(region[0]);
        return lvl.morton_to_rank.at(rid);
    };

    std::unordered_set<int> neigh_set(neighbor_ranks.begin(), neighbor_ranks.end());

    for (int64_t morton_idx : needed_remote_mortons) {
        const int owner = owner_of_morton(morton_idx);
        if (owner == rank) continue;

        // Enforce the “1-hop only” assumption:
        if (!neigh_set.count(owner)) {
            throw std::runtime_error(
                "Assisting request targets non-1hop rank " + std::to_string(owner) +
                " from rank " + std::to_string(rank) +
                ". If this is valid, expand neighbor list to 2-hop process regions.");
        }

        req_to_send[owner].push_back(morton_idx);
    }

    // Dedup per-dest
    for (auto& [dst, v] : req_to_send) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    auto throw_mpi_error = [&](int err, const char* op, int peer) {
        if (err == MPI_SUCCESS) return;
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            std::string("exchange_assisting_for_mortons_onehop: ") + op +
            " failed for rank " + std::to_string(peer) +
            " with error " + std::string(errbuf, errlen));
    };

    // Step 2: Exchange per-neighbor assisting request counts.
    std::vector<int> recv_counts(neighbor_ranks.size(), 0);
    std::vector<int> send_counts(neighbor_ranks.size(), 0);
    std::vector<MPI_Request> requests;

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        auto it = req_to_send.find(neighbor_ranks[i]);
        send_counts[i] = (it == req_to_send.end()) ? 0 : static_cast<int>(it->second.size());
    }

    auto comm_start = clock::now();
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Request req;
        int err = MPI_Irecv(&recv_counts[i], 1, MPI_INT, neighbor_ranks[i], 600, MPI_COMM_WORLD, &req);
        throw_mpi_error(err, "count Irecv", neighbor_ranks[i]);
        requests.push_back(req);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        int err = MPI_Send(&send_counts[i], 1, MPI_INT, neighbor_ranks[i], 600, MPI_COMM_WORLD);
        throw_mpi_error(err, "count Send", neighbor_ranks[i]);
    }
    if (!requests.empty()) {
        int err = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        throw_mpi_error(err, "count Waitall", -1);
        requests.clear();
    }
    communication_time += (clock::now() - comm_start);

    // Step 3: Exchange the requested Morton lists.
    std::vector<std::vector<int64_t>> req_received(neighbor_ranks.size());
    comm_start = clock::now();
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_counts[i] > 0) {
            req_received[i].resize(static_cast<size_t>(recv_counts[i]));
            MPI_Request req;
            int err = MPI_Irecv(req_received[i].data(), recv_counts[i], MPI_INT64_T,
                                neighbor_ranks[i], 601, MPI_COMM_WORLD, &req);
            throw_mpi_error(err, "request-list Irecv", neighbor_ranks[i]);
            requests.push_back(req);
        }
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_counts[i] > 0) {
            int err = MPI_Send(req_to_send[neighbor_ranks[i]].data(),
                               send_counts[i], MPI_INT64_T, neighbor_ranks[i], 601, MPI_COMM_WORLD);
            throw_mpi_error(err, "request-list Send", neighbor_ranks[i]);
        }
    }
    if (!requests.empty()) {
        int err = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        throw_mpi_error(err, "request-list Waitall", -1);
        requests.clear();
    }
    communication_time += (clock::now() - comm_start);

    // Step 4: Build packed responses and exchange payload sizes.
    std::vector<uint64_t> recv_sizes(neighbor_ranks.size(), 0);
    std::vector<uint64_t> send_sizes(neighbor_ranks.size(), 0);
    std::vector<std::vector<char>> send_buffers(neighbor_ranks.size());

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (req_received[i].empty()) {
            continue;
        }

        size_t total = 0;
        for (int64_t morton_idx : req_received[i]) {
            auto* b = lvl.find_local_box(morton_idx);
            if (!b) throw std::runtime_error("assist: requested morton not local " + std::to_string(morton_idx));

            PointDataRequest<CoordType> p;
            p.morton_index = morton_idx;
            p.source_rank  = rank;
            p.coords       = b->point_coords;
            p.skel_indices = b->skeleton_indices;
            p.on_boundary  = b->on_boundary;

            total += sizeof(uint64_t) + get_serialized_size(p);
        }

        send_buffers[i].resize(total);
        char* ptr = send_buffers[i].data();
        for (int64_t morton_idx : req_received[i]) {
            auto* b = lvl.find_local_box(morton_idx);

            PointDataRequest<CoordType> p;
            p.morton_index = morton_idx;
            p.source_rank  = rank;
            p.coords       = b->point_coords;
            p.skel_indices = b->skeleton_indices;
            p.on_boundary  = b->on_boundary;

            const uint64_t nbytes = (uint64_t)get_serialized_size(p);
            std::memcpy(ptr, &nbytes, sizeof(uint64_t));
            ptr += sizeof(uint64_t);
            ptr = serialize(p, ptr);
        }
        send_sizes[i] = static_cast<uint64_t>(total);
    }

    comm_start = clock::now();
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Request req;
        int err = MPI_Irecv(&recv_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i], 602, MPI_COMM_WORLD, &req);
        throw_mpi_error(err, "size Irecv", neighbor_ranks[i]);
        requests.push_back(req);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        int err = MPI_Send(&send_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i], 602, MPI_COMM_WORLD);
        throw_mpi_error(err, "size Send", neighbor_ranks[i]);
    }
    if (!requests.empty()) {
        int err = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        throw_mpi_error(err, "size Waitall", -1);
        requests.clear();
    }
    communication_time += (clock::now() - comm_start);

    // Step 5: Post all payload receives, then send the serialized payload buffers.
    std::vector<std::vector<char>> recv_buffers(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] > 0) {
            recv_buffers[i].resize((size_t)recv_sizes[i]);
        }
    }

    comm_start = clock::now();
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] > 0) {
            int err = MPI_Irecv_large(recv_buffers[i].data(), (size_t)recv_sizes[i], MPI_CHAR,
                                      neighbor_ranks[i], 603, MPI_COMM_WORLD, requests);
            throw_mpi_error(err, "payload Irecv_large", neighbor_ranks[i]);
        }
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_sizes[i] > 0) {
            int err = MPI_Send_large(send_buffers[i].data(), send_buffers[i].size(), MPI_CHAR,
                                     neighbor_ranks[i], 603, MPI_COMM_WORLD);
            throw_mpi_error(err, "payload Send_large", neighbor_ranks[i]);
        }
    }
    if (!requests.empty()) {
        int err = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        throw_mpi_error(err, "payload Waitall", -1);
        requests.clear();
    }
    communication_time += (clock::now() - comm_start);

    // Step 6: Deserialize received assisting boxes into the level caches.
    
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] == 0) continue;
        const char* ptr = recv_buffers[i].data();
        const char* end = ptr + recv_buffers[i].size();

        while (ptr < end) {
            uint64_t rec_bytes = 0;
            std::memcpy(&rec_bytes, ptr, sizeof(uint64_t));
            ptr += sizeof(uint64_t);

            PointDataRequest<CoordType> req;
            ptr = deserialize(req, ptr);

            // populate lvl.assisting_boxes and map
            auto it = lvl.assisting_box_points_for_kernel_evaluation.find(req.morton_index);
            if (it == lvl.assisting_box_points_for_kernel_evaluation.end()) {
                std::runtime_error("assist recv: unexpected new assisting morton " + std::to_string(req.morton_index));
                // const int64_t idx = (int64_t)lvl.assisting_boxes.size();
                // lvl.assisting_boxes.push_back(std::move(req));
                // lvl.assisting_box_points_for_kernel_evaluation[req.morton_index] = idx;
            } else {
                lvl.assisting_boxes[(size_t)it->second] = std::move(req);
            }
        }
        if (ptr != end) throw std::runtime_error("assist recv buffer parse mismatch");
    }

    return communication_time;
}


template <typename DataType>
static void copy_contiguous_into(MatrixStorage<DataType>& dst,
                                 const std::vector<DataType>& src,
                                 int64_t rows, int64_t cols)
{
    if (rows <= 0 || cols <= 0) throw std::runtime_error("copy_contiguous_into: bad dims");
    if ((int64_t)src.size() != rows * cols) throw std::runtime_error("copy_contiguous_into: src size mismatch");

    if (!dst.is_allocated() || dst.rows != rows || dst.cols != cols)
        throw std::runtime_error("copy_contiguous_into: dst not allocated or dim mismatch");

    // src is column-major with lda=rows
    for (int64_t j = 0; j < cols; ++j)
        for (int64_t i = 0; i < rows; ++i)
            dst(i, j) = src[(size_t)(i + j * rows)];
}

template <typename DataType>
static void add_contiguous_into(MatrixStorage<DataType>& dst,
                                const std::vector<DataType>& add,
                                int64_t rows, int64_t cols)
{
    if (rows <= 0 || cols <= 0) throw std::runtime_error("add_contiguous_into: bad dims");
    if ((int64_t)add.size() != rows * cols) throw std::runtime_error("add_contiguous_into: add size mismatch");

    if (!dst.is_allocated() || dst.rows != rows || dst.cols != cols)
        throw std::runtime_error("add_contiguous_into: dst not allocated or dim mismatch");
    
    // add is column-major with lda=rows
    for (int64_t j = 0; j < cols; ++j)
        for (int64_t i = 0; i < rows; ++i)
            dst(i, j) += (add[(size_t)(i + j * rows)]);
}

template <typename DataType>
static void add_transposed_contiguous_into(MatrixStorage<DataType>& dst,
                                           const std::vector<DataType>& add,
                                           int64_t src_rows,
                                           int64_t src_cols,
                                           bool do_conj)
{
    if (src_rows <= 0 || src_cols <= 0)
        throw std::runtime_error("add_transposed_contiguous_into: bad dims");
    if ((int64_t)add.size() != src_rows * src_cols)
        throw std::runtime_error("add_transposed_contiguous_into: add size mismatch");

    const int64_t rows = src_cols;
    const int64_t cols = src_rows;
    if (!dst.is_allocated() || dst.rows != rows || dst.cols != cols)
        throw std::runtime_error("add_transposed_contiguous_into: dst not allocated or dim mismatch");

    for (int64_t j = 0; j < cols; ++j) {
        for (int64_t i = 0; i < rows; ++i) {
            dst(i, j) += maybe_conj(add[(size_t)(j + i * src_rows)], do_conj);
        }
    }
}


enum class IncomingApplyTaskKind : uint8_t {
    Replace = 0,
    Diag = 1,
    Offdiag = 2
};

template <typename DataType>
struct IncomingApplyTask {
    IncomingApplyTaskKind kind = IncomingApplyTaskKind::Replace;
    int64_t target_morton = -1;
    int64_t neighbor_morton = -1;
    EdgeKind edge_kind = EdgeKind::Near;
    const DenseBlock<DataType>* block = nullptr;
    bool transpose_delta = false;
};

/**
 * @brief Apply received factorization updates to locally-owned boxes (symmetric/Hermitian storage),
 *        allocating missing base matrix blocks via on-demand kernel evaluation.
 *
 * This routine consumes the merged incoming update packets and applies them to local storage:
 *
 *  - Step 7 (REPLACE): directed near-field replacement blocks.
 *    For symmetric storage, we overwrite the locally-stored A_NS for (target, neighbor).
 *    The payload is assumed to already be in the target's storage orientation. (For Hermitian,
 *    the sender must have packed the conjugate-transpose as required by the storage convention.)
 *
 *  - Step 9 (ADD): accumulated deltas keyed by canonical (lo, hi, kind).
 *    * Diagonal (EdgeKind::Diag): adds into the local Schur complement. If the Schur complement
 *      is not allocated yet (delta can arrive before allocation), we evaluate the base diagonal
 *      block using the kernel and then add the delta.
 *    * Off-diagonal (Near/Far): updates one or both endpoint views (if endpoints are local).
 *      If the required interaction block is missing or not allocated, we evaluate the base
 *      interaction block using the kernel, using either local point coordinates or previously
 *      received assisting box data (lvl.assisting_boxes). Then we add the delta.
 *
 * Error handling:
 *  - Throws if a required local target is missing.
 *  - Throws if received block payload sizes are inconsistent.
 *  - Throws if block dimensions are inconsistent with existing storage (coloring should prevent
 *    size changes; we check and fail fast if violated).
 *  - Throws if assisting point data is required but not present.
 *
 * Coordinate layout assumptions:
 *  - coords_for_count_from_box / coords_for_count_from_assist must produce coordinates in the
 *    format expected by kernel->evaluate_block(...).
 *
 * Threading:
 *  - Incoming updates are first rewritten into per-target task lists.
 *  - A static OpenMP loop then assigns disjoint local boxes to threads, so interaction-map
 *    inserts and matrix writes stay lock-free.
 */
template <typename KernelType, typename CoordType, typename DataType>
static void apply_updates_with_kernel_symmetric(
    ParallelTree<CoordType, DataType>* tree,
    TreeLevel<CoordType, DataType>& lvl,
    KernelType* kernel,
    const PendingFactorUpdates<DataType>& incoming,
    bool is_hermitian)
{
    const int dim = tree->dimension;

    auto get_assist = [&](int64_t morton) -> const PointDataRequest<CoordType>* {
        auto it = lvl.assisting_box_points_for_kernel_evaluation.find(morton);
        if (it == lvl.assisting_box_points_for_kernel_evaluation.end()) return nullptr;
        return &lvl.assisting_boxes[static_cast<size_t>(it->second)];
    };

    auto local_box_index = [&](int64_t morton) -> int64_t {
        if (morton < lvl.local_morton_start || morton > lvl.local_morton_end) {
            return -1;
        }
        return morton - lvl.local_morton_start;
    };

    auto local_or_ghost_box = [&](int64_t morton) -> const BoxData<CoordType, DataType>* {
        if (auto* b = lvl.find_local_box(morton)) {
            return b;
        }
        auto git = lvl.ghost_id_to_index.find(morton);
        if (git != lvl.ghost_id_to_index.end()) {
            return &lvl.ghost_boxes[static_cast<size_t>(git->second)];
        }
        return nullptr;
    };

    auto endpoint_skeleton = [&](int64_t morton) -> const std::vector<int64_t>* {
        if (const auto* b = local_or_ghost_box(morton)) {
            return b->skeleton_indices.empty() ? nullptr : &b->skeleton_indices;
        }
        if (const auto* assist = get_assist(morton)) {
            return assist->skel_indices.empty() ? nullptr : &assist->skel_indices;
        }
        return nullptr;
    };

    auto current_point_count = [&](int64_t morton) -> int64_t {
        const bool eliminated =
            lvl.eliminated_boxes.find(morton) != lvl.eliminated_boxes.end();
        if (const auto* b = local_or_ghost_box(morton)) {
            if (eliminated && !b->skeleton_indices.empty()) {
                return static_cast<int64_t>(b->skeleton_indices.size());
            }
            return b->num_points;
        }
        if (const auto* assist = get_assist(morton)) {
            if (eliminated && !assist->skel_indices.empty()) {
                return static_cast<int64_t>(assist->skel_indices.size());
            }
            return static_cast<int64_t>(assist->coords.size()) / dim;
        }
        throw std::runtime_error(
            "apply task: cannot determine current point count for morton " +
            std::to_string(morton));
    };

    std::vector<size_t> task_counts(lvl.local_boxes.size(), 0);

    for (const auto& [rk, block] : incoming.replace_blocks) {
        (void)block;
        const int64_t local_idx = local_box_index(rk.target_box);
        if (local_idx >= 0) {
            ++task_counts[static_cast<size_t>(local_idx)];
        }
    }

    for (const auto& [ek, delta_lo_hi] : incoming.accumulated_deltas) {
        (void)delta_lo_hi;
        if (ek.kind == EdgeKind::Diag) {
            const int64_t local_idx = local_box_index(ek.lo);
            if (local_idx >= 0) {
                ++task_counts[static_cast<size_t>(local_idx)];
            }
            continue;
        }

        const int64_t local_lo = local_box_index(ek.lo);
        const int64_t local_hi = local_box_index(ek.hi);
        if (local_lo >= 0) {
            ++task_counts[static_cast<size_t>(local_lo)];
        }
        if (local_hi >= 0) {
            ++task_counts[static_cast<size_t>(local_hi)];
        }
    }

    std::vector<std::vector<IncomingApplyTask<DataType>>> tasks_per_box(lvl.local_boxes.size());
    for (size_t idx = 0; idx < tasks_per_box.size(); ++idx) {
        tasks_per_box[idx].reserve(task_counts[idx]);
    }

    for (const auto& [rk, block] : incoming.replace_blocks) {
        const int64_t local_idx = local_box_index(rk.target_box);
        if (local_idx < 0) {
            continue;
        }

        IncomingApplyTask<DataType> task;
        task.kind = IncomingApplyTaskKind::Replace;
        task.target_morton = rk.target_box;
        task.neighbor_morton = rk.neighbor_box;
        task.edge_kind = EdgeKind::Near;
        task.block = &block;
        tasks_per_box[static_cast<size_t>(local_idx)].push_back(task);
    }

    for (const auto& [ek, delta_lo_hi] : incoming.accumulated_deltas) {
        if (ek.kind == EdgeKind::Diag) {
            const int64_t local_idx = local_box_index(ek.lo);
            if (local_idx < 0) {
                continue;
            }

            IncomingApplyTask<DataType> task;
            task.kind = IncomingApplyTaskKind::Diag;
            task.target_morton = ek.lo;
            task.neighbor_morton = ek.hi;
            task.edge_kind = ek.kind;
            task.block = &delta_lo_hi;
            tasks_per_box[static_cast<size_t>(local_idx)].push_back(task);
            continue;
        }

        const int64_t local_lo = local_box_index(ek.lo);
        if (local_lo >= 0) {
            IncomingApplyTask<DataType> task;
            task.kind = IncomingApplyTaskKind::Offdiag;
            task.target_morton = ek.lo;
            task.neighbor_morton = ek.hi;
            task.edge_kind = ek.kind;
            task.block = &delta_lo_hi;
            task.transpose_delta = false;
            tasks_per_box[static_cast<size_t>(local_lo)].push_back(task);
        }

        const int64_t local_hi = local_box_index(ek.hi);
        if (local_hi >= 0) {
            IncomingApplyTask<DataType> task;
            task.kind = IncomingApplyTaskKind::Offdiag;
            task.target_morton = ek.hi;
            task.neighbor_morton = ek.lo;
            task.edge_kind = ek.kind;
            task.block = &delta_lo_hi;
            task.transpose_delta = true;
            tasks_per_box[static_cast<size_t>(local_hi)].push_back(task);
        }
    }

    std::exception_ptr apply_exception;
    std::mutex apply_exception_mutex;
    std::atomic<bool> apply_failed{false};

    // Each local box is owned by exactly one loop iteration, so the apply phase
    // can mutate interaction maps and matrices without locks.
    #pragma omp parallel for schedule(static) if (lvl.local_boxes.size() > 1)
    for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(lvl.local_boxes.size()); ++local_idx) {
        if (apply_failed.load(std::memory_order_relaxed)) {
            continue;
        }

        try {
            auto& target_box = lvl.local_boxes[static_cast<size_t>(local_idx)];
            const auto& tasks = tasks_per_box[static_cast<size_t>(local_idx)];

            for (const auto& task : tasks) {
                if (apply_failed.load(std::memory_order_relaxed)) {
                    break;
                }

                const auto* payload = task.block;
                if (payload == nullptr) {
                    throw std::runtime_error("apply task: missing payload");
                }
                if (payload->rows <= 0 || payload->cols <= 0 ||
                    static_cast<int64_t>(payload->data.size()) != payload->rows * payload->cols) {
                    throw std::runtime_error(
                        "apply task: invalid DenseBlock payload for target " +
                        std::to_string(task.target_morton));
                }

                if (task.kind == IncomingApplyTaskKind::Replace) {
                    auto& mb = upsert_block_symmetric(target_box, task.neighbor_morton, EdgeKind::Near);
                    mb.A_NS.allocate(payload->rows, payload->cols, MatrixStorage<DataType>::FULL);
                    copy_contiguous_into(mb.A_NS, payload->data, payload->rows, payload->cols);
                    continue;
                }

                if (task.kind == IncomingApplyTaskKind::Diag) {
                    if (task.target_morton != task.neighbor_morton) {
                        throw std::runtime_error("apply diag task: lo!=hi");
                    }
                    if (!target_box.schur_complement.is_allocated()) {
                        if (payload->rows != payload->cols) {
                            throw std::runtime_error(
                                "apply diag task: nonsquare delta at morton " +
                                std::to_string(task.target_morton));
                        }

                        auto coords = coords_for_count_from_box(target_box, payload->rows, dim);
                        std::vector<DataType> A(static_cast<size_t>(payload->rows * payload->rows));
                        kernel->evaluate_block(
                            coords.data(), payload->rows,
                            coords.data(), payload->rows,
                            A.data(), payload->rows);

                        target_box.schur_complement.allocate(
                            payload->rows, payload->cols, MatrixStorage<DataType>::FULL);
                        copy_contiguous_into(
                            target_box.schur_complement, A, payload->rows, payload->cols);
                    }

                    if (!target_box.schur_complement.is_allocated() ||
                        target_box.schur_complement.rows != payload->rows ||
                        target_box.schur_complement.cols != payload->cols) {
                        throw std::runtime_error(
                            "apply diag task: dim mismatch at morton " +
                            std::to_string(task.target_morton));
                    }

                    add_contiguous_into(
                        target_box.schur_complement,
                        payload->data,
                        payload->rows,
                        payload->cols);
                    continue;
                }

                const int64_t current_self =
                    (lvl.eliminated_boxes.find(task.target_morton) != lvl.eliminated_boxes.end() &&
                     !target_box.skeleton_indices.empty())
                        ? static_cast<int64_t>(target_box.skeleton_indices.size())
                        : target_box.num_points;
                const int64_t current_neighbor = current_point_count(task.neighbor_morton);

                const std::vector<DataType>* payload_data = &payload->data;
                std::vector<DataType> sliced_payload;
                int64_t payload_rows = payload->rows;
                int64_t payload_cols = payload->cols;

                if (task.transpose_delta) {
                    const bool need_row_slice = (payload_rows != current_self);
                    const bool need_col_slice = (payload_cols != current_neighbor);
                    if (need_row_slice || need_col_slice) {
                        if (payload_rows < current_self || payload_cols < current_neighbor) {
                            std::ostringstream oss;
                            oss << "apply offdiag task: payload smaller than current view"
                                << " target=" << task.target_morton
                                << " neighbor=" << task.neighbor_morton
                                << " payload_rows=" << payload_rows
                                << " payload_cols=" << payload_cols
                                << " current_self=" << current_self
                                << " current_neighbor=" << current_neighbor
                                << " transpose_delta=" << task.transpose_delta;
                            throw std::runtime_error(oss.str());
                        }
                        const auto* neighbor_skeleton = need_col_slice ? endpoint_skeleton(task.neighbor_morton) : nullptr;
                        if (need_row_slice && static_cast<int64_t>(target_box.skeleton_indices.size()) < current_self) {
                            throw std::runtime_error(
                                "apply offdiag task: target skeleton too small for transpose slice target=" +
                                std::to_string(task.target_morton));
                        }
                        if (need_col_slice &&
                            (neighbor_skeleton == nullptr ||
                             static_cast<int64_t>(neighbor_skeleton->size()) < current_neighbor)) {
                            throw std::runtime_error(
                                "apply offdiag task: neighbor skeleton too small for transpose slice neighbor=" +
                                std::to_string(task.neighbor_morton));
                        }

                        sliced_payload.resize(static_cast<size_t>(current_self * current_neighbor));
                        for (int64_t j = 0; j < current_neighbor; ++j) {
                            const int64_t src_col = need_col_slice
                                ? (*neighbor_skeleton)[static_cast<size_t>(j)]
                                : j;
                            if (src_col < 0 || src_col >= payload_cols) {
                                throw std::runtime_error(
                                    "apply offdiag task: transpose slice neighbor index out of bounds");
                            }
                            for (int64_t i = 0; i < current_self; ++i) {
                                const int64_t src_row = need_row_slice
                                    ? target_box.skeleton_indices[static_cast<size_t>(i)]
                                    : i;
                                if (src_row < 0 || src_row >= payload_rows) {
                                    throw std::runtime_error(
                                        "apply offdiag task: transpose slice target index out of bounds");
                                }
                                sliced_payload[static_cast<size_t>(i + j * current_self)] =
                                    payload->data[static_cast<size_t>(src_row + src_col * payload_rows)];
                            }
                        }
                        payload_data = &sliced_payload;
                        payload_rows = current_self;
                        payload_cols = current_neighbor;
                    }
                } else {
                    const bool need_row_slice = (payload_rows != current_neighbor);
                    const bool need_col_slice = (payload_cols != current_self);
                    if (need_row_slice || need_col_slice) {
                        if (payload_rows < current_neighbor || payload_cols < current_self) {
                            std::ostringstream oss;
                            oss << "apply offdiag task: payload smaller than current view"
                                << " target=" << task.target_morton
                                << " neighbor=" << task.neighbor_morton
                                << " payload_rows=" << payload_rows
                                << " payload_cols=" << payload_cols
                                << " current_neighbor=" << current_neighbor
                                << " current_self=" << current_self
                                << " transpose_delta=" << task.transpose_delta;
                            throw std::runtime_error(oss.str());
                        }
                        const auto* neighbor_skeleton = need_row_slice ? endpoint_skeleton(task.neighbor_morton) : nullptr;
                        if (need_col_slice && static_cast<int64_t>(target_box.skeleton_indices.size()) < current_self) {
                            throw std::runtime_error(
                                "apply offdiag task: target skeleton too small for direct slice target=" +
                                std::to_string(task.target_morton));
                        }
                        if (need_row_slice &&
                            (neighbor_skeleton == nullptr ||
                             static_cast<int64_t>(neighbor_skeleton->size()) < current_neighbor)) {
                            throw std::runtime_error(
                                "apply offdiag task: neighbor skeleton too small for direct slice neighbor=" +
                                std::to_string(task.neighbor_morton));
                        }

                        sliced_payload.resize(static_cast<size_t>(current_neighbor * current_self));
                        for (int64_t j = 0; j < current_self; ++j) {
                            const int64_t src_col = need_col_slice
                                ? target_box.skeleton_indices[static_cast<size_t>(j)]
                                : j;
                            if (src_col < 0 || src_col >= payload_cols) {
                                throw std::runtime_error(
                                    "apply offdiag task: direct slice target index out of bounds");
                            }
                            for (int64_t i = 0; i < current_neighbor; ++i) {
                                const int64_t src_row = need_row_slice
                                    ? (*neighbor_skeleton)[static_cast<size_t>(i)]
                                    : i;
                                if (src_row < 0 || src_row >= payload_rows) {
                                    throw std::runtime_error(
                                        "apply offdiag task: direct slice neighbor index out of bounds");
                                }
                                sliced_payload[static_cast<size_t>(i + j * current_neighbor)] =
                                    payload->data[static_cast<size_t>(src_row + src_col * payload_rows)];
                            }
                        }
                        payload_data = &sliced_payload;
                        payload_rows = current_neighbor;
                        payload_cols = current_self;
                    }
                }

                const int64_t n_self = task.transpose_delta ? payload_rows : payload_cols;
                const int64_t n_nei = task.transpose_delta ? payload_cols : payload_rows;

                ModifiedBlock<DataType>* mb =
                    find_block_symmetric(target_box, task.neighbor_morton, task.edge_kind);

                if (!mb || !mb->A_NS.is_allocated()) {
                    auto& newb =
                        upsert_block_symmetric(target_box, task.neighbor_morton, task.edge_kind);

                    auto coords_self = coords_for_count_from_box(target_box, n_self, dim);

                    std::vector<CoordType> coords_nei;
                    const int64_t neighbor_local_idx = local_box_index(task.neighbor_morton);
                    if (neighbor_local_idx >= 0) {
                        coords_nei = coords_for_count_from_box(
                            lvl.local_boxes[static_cast<size_t>(neighbor_local_idx)], n_nei, dim);
                    } else {
                        const auto* assist = get_assist(task.neighbor_morton);
                        if (!assist) {
                            throw std::runtime_error(
                                "apply task: missing assisting coords for morton " +
                                std::to_string(task.neighbor_morton));
                        }
                        coords_nei = coords_for_count_from_assist(*assist, n_nei, dim);
                    }

                    std::vector<DataType> A_self_nei(static_cast<size_t>(n_self * n_nei));
                    kernel->evaluate_block(
                        coords_self.data(), n_self,
                        coords_nei.data(), n_nei,
                        A_self_nei.data(), n_self);

                    const auto stored = transpose_colmajor(
                        A_self_nei,
                        n_self,
                        n_nei,
                        /*do_conj=*/is_hermitian);

                    newb.A_NS.allocate(n_nei, n_self, MatrixStorage<DataType>::FULL);
                    copy_contiguous_into(newb.A_NS, stored, n_nei, n_self);
                    mb = &newb;
                }

                if (mb->A_NS.is_allocated() &&
                    (mb->A_NS.rows != n_nei || mb->A_NS.cols != n_self)) {
                    if (mb->A_NS.rows < n_nei || mb->A_NS.cols < n_self) {
                        std::ostringstream oss;
                        oss << "apply offdiag task: existing block smaller than current view"
                            << " target=" << task.target_morton
                            << " neighbor=" << task.neighbor_morton
                            << " existing_rows=" << mb->A_NS.rows
                            << " existing_cols=" << mb->A_NS.cols
                            << " expected_rows=" << n_nei
                            << " expected_cols=" << n_self
                            << " transpose_delta=" << task.transpose_delta;
                        throw std::runtime_error(oss.str());
                    }

                    const bool need_row_slice_existing = (mb->A_NS.rows != n_nei);
                    const bool need_col_slice_existing = (mb->A_NS.cols != n_self);
                    const auto* neighbor_skeleton =
                        need_row_slice_existing ? endpoint_skeleton(task.neighbor_morton) : nullptr;

                    if (need_col_slice_existing &&
                        static_cast<int64_t>(target_box.skeleton_indices.size()) < n_self) {
                        throw std::runtime_error(
                            "apply offdiag task: target skeleton too small for existing-block slice target=" +
                            std::to_string(task.target_morton));
                    }
                    if (need_row_slice_existing &&
                        (neighbor_skeleton == nullptr ||
                         static_cast<int64_t>(neighbor_skeleton->size()) < n_nei)) {
                        throw std::runtime_error(
                            "apply offdiag task: neighbor skeleton too small for existing-block slice neighbor=" +
                            std::to_string(task.neighbor_morton));
                    }

                    std::vector<DataType> sliced_existing(
                        static_cast<size_t>(n_nei * n_self));
                    for (int64_t j = 0; j < n_self; ++j) {
                        const int64_t src_col = need_col_slice_existing
                            ? target_box.skeleton_indices[static_cast<size_t>(j)]
                            : j;
                        if (src_col < 0 || src_col >= mb->A_NS.cols) {
                            throw std::runtime_error(
                                "apply offdiag task: existing-block target index out of bounds");
                        }
                        for (int64_t i = 0; i < n_nei; ++i) {
                            const int64_t src_row = need_row_slice_existing
                                ? (*neighbor_skeleton)[static_cast<size_t>(i)]
                                : i;
                            if (src_row < 0 || src_row >= mb->A_NS.rows) {
                                throw std::runtime_error(
                                    "apply offdiag task: existing-block neighbor index out of bounds");
                            }
                            sliced_existing[static_cast<size_t>(i + j * n_nei)] =
                                mb->A_NS(src_row, src_col);
                        }
                    }

                    mb->A_NS.allocate(n_nei, n_self, MatrixStorage<DataType>::FULL);
                    mb->A_NS.data = std::move(sliced_existing);
                }

                if (!mb->A_NS.is_allocated() ||
                    mb->A_NS.rows != n_nei || mb->A_NS.cols != n_self) {
                    std::ostringstream oss;
                    oss << "apply offdiag task: dim mismatch"
                        << " target=" << task.target_morton
                        << " neighbor=" << task.neighbor_morton
                        << " existing_rows=" << mb->A_NS.rows
                        << " existing_cols=" << mb->A_NS.cols
                        << " expected_rows=" << n_nei
                        << " expected_cols=" << n_self
                        << " transpose_delta=" << task.transpose_delta;
                    throw std::runtime_error(oss.str());
                }

                if (task.transpose_delta) {
                    add_transposed_contiguous_into(
                        mb->A_NS,
                        *payload_data,
                        payload_rows,
                        payload_cols,
                        is_hermitian);
                } else {
                    add_contiguous_into(mb->A_NS, *payload_data, n_nei, n_self);
                }
            }
        } catch (...) {
            if (!apply_failed.exchange(true, std::memory_order_relaxed)) {
                std::lock_guard<std::mutex> lock(apply_exception_mutex);
                apply_exception = std::current_exception();
            }
        }
    }

    if (apply_exception) {
        std::rethrow_exception(apply_exception);
    }
}

static inline const char* edgekind_name(EdgeKind k) {
    switch (k) {
        case EdgeKind::Near: return "Near";
        case EdgeKind::Far:  return "Far";
        case EdgeKind::Diag: return "Diag";
        default:             return "Unknown";
    }
}

/**
 * @brief Print all pending factor update entries (keys + matrix dims) for debugging.
 *
 * Prints:
 *  - replace_blocks: (target_box, neighbor_box) and DenseBlock rows×cols (+ data.size)
 *  - accumulated_deltas: (lo, hi, kind) and DenseBlock rows×cols (+ data.size)
 *
 * If sort_keys=true, output is deterministic (sorted by key fields).
 */
template <typename DataType>
void print_pending_factor_updates(
    const PendingFactorUpdates<DataType>& pending,
    std::ostream& os,
    const char* prefix = "",
    bool sort_keys = true)
{
    os << prefix << "PendingFactorUpdates:\n";
    os << prefix << "  replace_blocks: " << pending.replace_blocks.size() << " entries\n";
    os << prefix << "  accumulated_deltas: " << pending.accumulated_deltas.size() << " entries\n";

    // -------- replace_blocks --------
    if (!pending.replace_blocks.empty()) {
        if (!sort_keys) {
            for (const auto& [k, b] : pending.replace_blocks) {
                os << prefix << "  REPLACE (target=" << k.target_box
                   << ", neighbor=" << k.neighbor_box << ") : "
                   << b.rows << "x" << b.cols
                   << " data.size=" << b.data.size()
                   << " (expected " << (static_cast<int64_t>(b.rows) * b.cols) << ")\n";
            }
        } else {
            struct Item { ReplaceKey k; const DenseBlock<DataType>* b; };
            std::vector<Item> items;
            items.reserve(pending.replace_blocks.size());
            for (const auto& kv : pending.replace_blocks) items.push_back({kv.first, &kv.second});

            std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
                if (a.k.target_box != b.k.target_box) return a.k.target_box < b.k.target_box;
                return a.k.neighbor_box < b.k.neighbor_box;
            });

            for (const auto& it : items) {
                const auto& k = it.k;
                const auto& b = *it.b;
                os << prefix << "  REPLACE (target=" << k.target_box
                   << ", neighbor=" << k.neighbor_box << ") : "
                   << b.rows << "x" << b.cols
                   << " data.size=" << b.data.size()
                   << " (expected " << (static_cast<int64_t>(b.rows) * b.cols) << ")\n";
            }
        }
    }

    // -------- accumulated_deltas --------
    if (!pending.accumulated_deltas.empty()) {
        if (!sort_keys) {
            for (const auto& [k, b] : pending.accumulated_deltas) {
                os << prefix << "  DELTA (" << edgekind_name(k.kind)
                   << " lo=" << k.lo << ", hi=" << k.hi << ") : "
                   << b.rows << "x" << b.cols
                   << " data.size=" << b.data.size()
                   << " (expected " << (static_cast<int64_t>(b.rows) * b.cols) << ")\n";
            }
        } else {
            struct Item { EdgeKey k; const DenseBlock<DataType>* b; };
            std::vector<Item> items;
            items.reserve(pending.accumulated_deltas.size());
            for (const auto& kv : pending.accumulated_deltas) items.push_back({kv.first, &kv.second});

            std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
                if (a.k.lo != b.k.lo) return a.k.lo < b.k.lo;
                if (a.k.hi != b.k.hi) return a.k.hi < b.k.hi;
                return static_cast<uint8_t>(a.k.kind) < static_cast<uint8_t>(b.k.kind);
            });

            for (const auto& it : items) {
                const auto& k = it.k;
                const auto& b = *it.b;
                os << prefix << "  DELTA (" << edgekind_name(k.kind)
                   << " lo=" << k.lo << ", hi=" << k.hi << ") : "
                   << b.rows << "x" << b.cols
                   << " data.size=" << b.data.size()
                   << " (expected " << (static_cast<int64_t>(b.rows) * b.cols) << ")\n";
            }
        }
    }
}


/**
 * @brief Exchange pending factorization updates with 1-hop neighboring process-regions and apply them.
 *
 * This function performs sparse, neighbor-only communication (no MPI_Alltoall/Alltoallv):
 *  1) Determine the set of 1-hop neighboring MPI ranks using the process-region Morton id
 *     (lvl.my_morton_id) and morton::{neighbors_2d, neighbors_3d} on the process grid.
 *  2) Partition locally-produced updates into per-destination packets:
 *      - Step 7 REPLACE: sent only to owner(target_box)
 *      - Step 9 ADD: sent to owner(lo) and/or owner(hi); sent once if both endpoints share owner
 *     Ownership is derived from assign_to_processes_{2d,3d} followed by lvl.morton_to_rank.
 *     (This routine enforces the assumption that all destinations are within 1-hop neighbor ranks;
 *      it throws if violated.)
 *  3) Pairwise exchange with each neighbor:
 *      - exchange uint64 payload size
 *      - exchange payload bytes (serialized PendingFactorUpdates)
 *     The implementation uses a simple size handshake and large-message send/recv helpers.
 *  4) Merge all received packets into one incoming_total.
 *  5) Before applying, request any missing assisting box point data needed for on-demand kernel
 *     evaluations of missing interaction blocks (via exchange_assisting_for_mortons_onehop()).
 *  6) Apply all updates locally (apply_updates_with_kernel_symmetric()).
 *
 * Assumptions:
 *  - Cross-process dependencies for this phase are limited to 1-hop process-region neighbors.
 *    If this is not true for some configurations, replace the neighbor set by a 2-hop stencil
 *    (neighbors_2hop_{2d,3d}) or generalize the neighbor discovery.
 *
 * Side effects:
 *  - Clears pending.replace_blocks and pending.accumulated_deltas on completion.
 *  - Populates lvl.assisting_boxes and lvl.assisting_box_points_for_kernel_evaluation as needed.
 *
 * Error handling:
 *  - Throws if non-1-hop communication would be required (assumption violation).
 *  - Throws on serialization size mismatches or invalid payloads.
 *
 * @return Communication-only duration spent in MPI send/recv/wait operations
 *         during this exchange (including assisting-data exchange).
 */
template <typename KernelType, typename CoordType, typename DataType>
std::chrono::high_resolution_clock::duration transport_and_apply_factor_updates_symmetric_onehop(
    ParallelTree<CoordType, DataType>* tree,
    int level,
    KernelType* kernel,
    PendingFactorUpdates<DataType>& pending,
    bool is_hermitian=false)
{
    using clock = std::chrono::high_resolution_clock;
    clock::duration communication_time{};

    auto& lvl = tree->levels[level];
    if (!lvl.is_process_active) return communication_time;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Step 1: Compute the 1-hop process-neighbor list used by this exchange.
    const std::vector<int> neighbor_ranks = compute_one_hop_neighbor_ranks(tree, lvl, level);
    std::unordered_set<int> neigh_set(neighbor_ranks.begin(), neighbor_ranks.end());

    const uint32_t grid_size = 1u << level;

    auto owner_of_morton = [&](int64_t morton_idx) -> int {
        std::vector<uint64_t> single = { static_cast<uint64_t>(morton_idx) };
        std::vector<uint32_t> region =
            (tree->dimension == 2)
            ? morton::assign_to_processes_2d(single, lvl.num_active_processes, grid_size)
            : morton::assign_to_processes_3d(single, lvl.num_active_processes, grid_size);
        const int rid = static_cast<int>(region[0]);
        return lvl.morton_to_rank.at(rid);
    };

    // Build outgoing packets per destination rank, while retaining locally-owned
    // updates so they still flow through the unified apply path below.
    std::unordered_map<int, PendingFactorUpdates<DataType>> out;
    PendingFactorUpdates<DataType> local_apply;

    // Route replacement blocks to the owner of target_box.
    for (const auto& [k, b] : pending.replace_blocks) {
        const int dst = owner_of_morton(k.target_box);
        if (dst == rank) {
            local_apply.replace_blocks[k] = b;
            continue;
        }

        if (!neigh_set.count(dst)) {
            throw std::runtime_error(
                "Step7 REPLACE targets non-1hop rank " + std::to_string(dst) +
                ". Expand neighbor list or fix assumption.");
        }
        out[dst].replace_blocks.emplace(k, b);
    }

    // Route additive deltas to one or both owners, depending on the endpoints.
    for (const auto& [ek, b] : pending.accumulated_deltas) {
        const int owner_lo = owner_of_morton(ek.lo);
        const int owner_hi = owner_of_morton(ek.hi);

        auto add_to = [&](int dst) {
            if (dst == rank) {
                return;
            }
            if (!neigh_set.count(dst)) {
                throw std::runtime_error(
                    "Step9 ADD targets non-1hop rank " + std::to_string(dst) +
                    ". Expand neighbor list or fix assumption.");
            }
            out[dst].accumulated_deltas.emplace(ek, b);
        };

        if (owner_lo == owner_hi) add_to(owner_lo);
        else { add_to(owner_lo); add_to(owner_hi); }
    }

    // Step 2: Exchange per-neighbor payload sizes by posting receives, then sending local sizes.
    std::vector<uint64_t> recv_sizes(neighbor_ranks.size(), 0);
    std::vector<uint64_t> send_sizes(neighbor_ranks.size(), 0);
    std::vector<MPI_Request> requests;
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        const int peer = neighbor_ranks[i];
        auto it = out.find(peer);
        send_sizes[i] = (it == out.end()) ? 0ull : (uint64_t)bytes_pending(it->second);
    }

    auto comm_start = clock::now();
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Request req;
        int err = MPI_Irecv(&recv_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i], 700, MPI_COMM_WORLD, &req);
        if (err != MPI_SUCCESS) {
            char errbuf[MPI_MAX_ERROR_STRING];
            int errlen = 0;
            MPI_Error_string(err, errbuf, &errlen);
            throw std::runtime_error(
                "transport_and_apply_factor_updates_symmetric_onehop: size Irecv failed for rank " +
                std::to_string(neighbor_ranks[i]) + " with error " + std::string(errbuf, errlen));
        }
        requests.push_back(req);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        int err = MPI_Send(&send_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i], 700, MPI_COMM_WORLD);
        if (err != MPI_SUCCESS) {
            char errbuf[MPI_MAX_ERROR_STRING];
            int errlen = 0;
            MPI_Error_string(err, errbuf, &errlen);
            throw std::runtime_error(
                "transport_and_apply_factor_updates_symmetric_onehop: size Send failed for rank " +
                std::to_string(neighbor_ranks[i]) + " with error " + std::string(errbuf, errlen));
        }
    }
    if (!requests.empty()) {
        int err = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        if (err != MPI_SUCCESS) {
            char errbuf[MPI_MAX_ERROR_STRING];
            int errlen = 0;
            MPI_Error_string(err, errbuf, &errlen);
            throw std::runtime_error(
                "transport_and_apply_factor_updates_symmetric_onehop: size Waitall failed with error " +
                std::string(errbuf, errlen));
        }
        requests.clear();
    }
    communication_time += (clock::now() - comm_start);

    // Step 3: Post all payload receives, then send the serialized update payloads.
    std::vector<std::vector<char>> recv_bufs(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] > 0) {
            recv_bufs[i].resize((size_t)recv_sizes[i]);
        }
    }

    std::vector<char> send_buffer; // reusable
    comm_start = clock::now();
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] > 0) {
            int err = MPI_Irecv_large(recv_bufs[i].data(), (size_t)recv_sizes[i], MPI_CHAR,
                                      neighbor_ranks[i], 701, MPI_COMM_WORLD, requests);
            if (err != MPI_SUCCESS) {
                char errbuf[MPI_MAX_ERROR_STRING];
                int errlen = 0;
                MPI_Error_string(err, errbuf, &errlen);
                throw std::runtime_error(
                    "transport_and_apply_factor_updates_symmetric_onehop: payload Irecv_large failed for rank " +
                    std::to_string(neighbor_ranks[i]) + " with error " + std::string(errbuf, errlen));
            }
        }
    }

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_sizes[i] == 0) {
            continue;
        }
        const int peer = neighbor_ranks[i];
        auto it = out.find(peer);
        if (it == out.end()) {
            throw std::runtime_error(
                "transport_and_apply_factor_updates_symmetric_onehop: missing payload for rank " +
                std::to_string(peer));
        }
        send_buffer.resize((size_t)send_sizes[i]);
        char* end_ptr = serialize(it->second, send_buffer.data());
        if ((size_t)(end_ptr - send_buffer.data()) != (size_t)send_sizes[i]) {
            throw std::runtime_error("serialize pending: byte mismatch");
        }
        int err = MPI_Send_large(send_buffer.data(), send_buffer.size(), MPI_CHAR,
                                 peer, 701, MPI_COMM_WORLD);
        if (err != MPI_SUCCESS) {
            char errbuf[MPI_MAX_ERROR_STRING];
            int errlen = 0;
            MPI_Error_string(err, errbuf, &errlen);
            throw std::runtime_error(
                "transport_and_apply_factor_updates_symmetric_onehop: payload Send_large failed for rank " +
                std::to_string(peer) + " with error " + std::string(errbuf, errlen));
        }
    }
    if (!requests.empty()) {
        int err = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        if (err != MPI_SUCCESS) {
            char errbuf[MPI_MAX_ERROR_STRING];
            int errlen = 0;
            MPI_Error_string(err, errbuf, &errlen);
            throw std::runtime_error(
                "transport_and_apply_factor_updates_symmetric_onehop: payload Waitall failed with error " +
                std::string(errbuf, errlen));
        }
        requests.clear();
    }
    communication_time += (clock::now() - comm_start);

    // Step 4: Merge incoming updates from all neighbors.
    PendingFactorUpdates<DataType> incoming_total;
    merge_pending(incoming_total, local_apply);
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] == 0) continue;
        PendingFactorUpdates<DataType> in;
        deserialize(in, recv_bufs[i].data());
        merge_pending(incoming_total, in);
    }
    
    // Step 5: Refresh assisting-box data needed by the apply path.
    std::vector<int64_t> need_assist;
    for (const auto& kv : lvl.assisting_box_points_for_kernel_evaluation) {
        need_assist.push_back(kv.first);
    }
    if (!need_assist.empty()) {
        communication_time += exchange_assisting_for_mortons_onehop(
            tree, lvl, level, neighbor_ranks, need_assist);
    }
    

    // Step 6: Apply the merged updates to local boxes.
    apply_updates_with_kernel_symmetric(tree, lvl, kernel, incoming_total, is_hermitian);
    
    // if(rank == 0)
    //     print_pending_factor_updates(incoming_total, std::cout, ("rank " + std::to_string(rank) + ": ").c_str());

    // clear the updates after each color
    clear_pending_factor_updates_memory(pending);

    return communication_time;
}






template<typename CoordType, typename DataType>
void gather_boxes_solve(
    ParallelTree<CoordType, DataType>* tree,
    int level,
    const std::vector<SolveDataRequest<CoordType, DataType>>& local_solve_data, 
    bool minimal_only = false,
    bool DEBUG = false) {

    auto& lvl = tree->levels[level];

    if (!lvl.is_process_active) {
        return;
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint32_t grid_size = 1 << level;

    auto abort_on_mpi_error = [&](int err,
                                  const char* op,
                                  int peer_rank,
                                  size_t bytes) {
        if (err == MPI_SUCCESS) {
            return;
        }
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        std::fprintf(stderr,
                     "[gather-solve-mpi-error] rank=%d level=%d op=%s peer=%d bytes=%zu err=%.*s\n",
                     rank, level, op, peer_rank, bytes, errlen, errbuf);
        std::fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, err);
    };

    struct RequestInfo {
        int64_t morton;
        bool is_ghost;
    };

    const size_t empty_vector_serialized_size = sizeof(size_t);
    const size_t empty_matrix_serialized_size =
        get_serialized_size(MatrixStorage<DataType>{});

    auto uses_factorization_data = [&](const RequestInfo& req_info) {
        return req_info.is_ghost && !minimal_only;
    };

    auto get_local_idx_for_request = [&](int64_t morton) -> int64_t {
        const int64_t local_idx = morton - lvl.local_morton_start;

        if (local_idx < 0 || local_idx >= lvl.num_boxes_local) {
            throw std::runtime_error(
                "gather_boxes_solve: Requested box " + std::to_string(morton) +
                " not found locally on rank " + std::to_string(rank) +
                " (local range: [" + std::to_string(lvl.local_morton_start) +
                ", " + std::to_string(lvl.local_morton_start + lvl.num_boxes_local - 1) + "])");
        }

        if (local_idx >= static_cast<int64_t>(local_solve_data.size())) {
            throw std::runtime_error(
                "gather_boxes_solve: Solve data for box " + std::to_string(morton) +
                " (local idx " + std::to_string(local_idx) +
                ") not found in local_solve_data (size: " +
                std::to_string(local_solve_data.size()) + ")");
        }

        return local_idx;
    };

    auto get_serialized_request_size_for_box =
        [&](const RequestInfo& req_info, int64_t local_idx) -> size_t {
            if (uses_factorization_data(req_info)) {
                throw std::runtime_error("color algorithm doesn't have halo region");
            }

            const auto& box = lvl.local_boxes[local_idx];
            const auto& rhs_info = local_solve_data[local_idx];

            auto vector_serialized_size = [](const auto& vec) -> size_t {
                using ValueType = typename std::decay_t<decltype(vec)>::value_type;
                return sizeof(size_t) + vec.size() * sizeof(ValueType);
            };

            size_t size = 0;
            size += sizeof(int64_t);  // morton_index
            size += sizeof(int);      // source_rank
            size += vector_serialized_size(rhs_info.right_side);
            size += vector_serialized_size(rhs_info.left_side);
            size += vector_serialized_size(box.redundant_indices);
            size += vector_serialized_size(box.skeleton_indices);
            size += empty_vector_serialized_size;  // one_hop
            size += empty_vector_serialized_size;  // use_full_set
            size += empty_vector_serialized_size;  // X_RR_pivots
            size += 7 * empty_matrix_serialized_size;

            return size;
        };

    auto populate_solve_request =
        [&](SolveDataRequest<CoordType, DataType>& solve_req,
            const RequestInfo& req_info,
            int64_t local_idx) {
            solve_req.morton_index = req_info.morton;
            solve_req.source_rank = rank;

            if (uses_factorization_data(req_info)) {
                throw std::runtime_error("color algorithm doesn't have halo region");
            }

            solve_req.copy_minimal_from_box(
                lvl.local_boxes[local_idx],
                local_solve_data[local_idx]);
        };

    // ========================================================================
    // Step 1: Group requested solve data by owning process.
    // ========================================================================
    std::unordered_map<int, std::vector<RequestInfo>> solve_requests_to_send;

    for (const auto& [morton, idx] : lvl.ghost_and_assisting_box_points_for_solve_map) {
        std::vector<uint64_t> single_box = {static_cast<uint64_t>(morton)};
        std::vector<uint32_t> morton_ids;

        if (tree->dimension == 2) {
            morton_ids = morton::assign_to_processes_2d(single_box, lvl.num_active_processes, grid_size);
        } else {
            morton_ids = morton::assign_to_processes_3d(single_box, lvl.num_active_processes, grid_size);
        }

        int morton_region_id = static_cast<int>(morton_ids[0]);

        auto rank_it = lvl.morton_to_rank.find(morton_region_id);
        if (rank_it == lvl.morton_to_rank.end()) {
            throw std::runtime_error(
                "gather_boxes_solve: Morton region " + std::to_string(morton_region_id) +
                " not found in morton_to_rank map at level " + std::to_string(level));
        }

        int owner_rank = rank_it->second;

        if (owner_rank != rank) {
            RequestInfo req_info;
            req_info.morton = morton;
            req_info.is_ghost = lvl.is_ghost_solve[idx];
            solve_requests_to_send[owner_rank].push_back(req_info);
        }
    }

    std::vector<int> neighbor_ranks;
    neighbor_ranks.reserve(solve_requests_to_send.size());
    for (const auto& [dest_rank, _] : solve_requests_to_send) {
        neighbor_ranks.push_back(dest_rank);
    }
    std::sort(neighbor_ranks.begin(), neighbor_ranks.end());

    std::vector<MPI_Request> requests;
    std::vector<int> send_counts;
    send_counts.reserve(neighbor_ranks.size());
    std::vector<int> recv_counts(neighbor_ranks.size(), 0);

    // ========================================================================
    // Step 2: Exchange per-neighbor solve request counts.
    // ========================================================================
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Request req;
        int ierr = MPI_Irecv(&recv_counts[i], 1, MPI_INT,
                             neighbor_ranks[i], 300, MPI_COMM_WORLD, &req);
        abort_on_mpi_error(ierr, "MPI_Irecv(counts)", neighbor_ranks[i], sizeof(int));
        requests.push_back(req);
    }

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        send_counts.push_back(static_cast<int>(solve_requests_to_send[neighbor_ranks[i]].size()));
        int ierr = MPI_Send(&send_counts[i], 1, MPI_INT,
                            neighbor_ranks[i], 300, MPI_COMM_WORLD);
        abort_on_mpi_error(ierr, "MPI_Send(counts)", neighbor_ranks[i], sizeof(int));
    }

    if (!requests.empty()) {
        int ierr = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        abort_on_mpi_error(ierr, "MPI_Waitall(counts)", -1, requests.size() * sizeof(MPI_Request));
        requests.clear();
    }

    // ========================================================================
    // Step 3: Exchange the requested solve Morton lists plus ghost/full-data flags.
    // ========================================================================
    std::vector<std::vector<RequestInfo>> requests_received(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_counts[i] > 0) {
            requests_received[i].resize(static_cast<size_t>(recv_counts[i]));
            MPI_Request req;
            int ierr = MPI_Irecv(requests_received[i].data(),
                                 recv_counts[i] * static_cast<int>(sizeof(RequestInfo)),
                                 MPI_CHAR,
                                 neighbor_ranks[i],
                                 301,
                                 MPI_COMM_WORLD,
                                 &req);
            abort_on_mpi_error(ierr,
                               "MPI_Irecv(request-list)",
                               neighbor_ranks[i],
                               static_cast<size_t>(recv_counts[i]) * sizeof(RequestInfo));
            requests.push_back(req);
        }
    }

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_counts[i] > 0) {
            int ierr = MPI_Send(solve_requests_to_send[neighbor_ranks[i]].data(),
                                send_counts[i] * static_cast<int>(sizeof(RequestInfo)),
                                MPI_CHAR,
                                neighbor_ranks[i],
                                301,
                                MPI_COMM_WORLD);
            abort_on_mpi_error(ierr,
                               "MPI_Send(request-list)",
                               neighbor_ranks[i],
                               static_cast<size_t>(send_counts[i]) * sizeof(RequestInfo));
        }
    }

    if (!requests.empty()) {
        int ierr = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        abort_on_mpi_error(ierr, "MPI_Waitall(request-list)", -1, requests.size() * sizeof(MPI_Request));
        requests.clear();
    }

    std::vector<size_t> send_sizes(neighbor_ranks.size(), 0);
    std::vector<size_t> recv_sizes(neighbor_ranks.size(), 0);

    // ========================================================================
    // Step 4: Compute response sizes from the requested solve boxes.
    // ========================================================================
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        size_t total_size = 0;
        for (const auto& req_info : requests_received[i]) {
            const int64_t local_idx = get_local_idx_for_request(req_info.morton);
            total_size += sizeof(size_t) +
                          get_serialized_request_size_for_box(req_info, local_idx);
        }
        send_sizes[i] = total_size;
    }

    // ========================================================================
    // Step 5: Exchange payload sizes by posting receives, then sending local sizes.
    // ========================================================================
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Request req;
        int ierr = MPI_Irecv(&recv_sizes[i], 1, MPI_UINT64_T,
                             neighbor_ranks[i], 302, MPI_COMM_WORLD, &req);
        abort_on_mpi_error(ierr, "MPI_Irecv(sizes)", neighbor_ranks[i], sizeof(size_t));
        requests.push_back(req);
    }

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        int ierr = MPI_Send(&send_sizes[i], 1, MPI_UINT64_T,
                            neighbor_ranks[i], 302, MPI_COMM_WORLD);
        abort_on_mpi_error(ierr, "MPI_Send(sizes)", neighbor_ranks[i], sizeof(size_t));
    }

    if (!requests.empty()) {
        int ierr = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        abort_on_mpi_error(ierr, "MPI_Waitall(sizes)", -1, requests.size() * sizeof(MPI_Request));
        requests.clear();
    }

    // ========================================================================
    // Step 6: Post all solve payload receives, then send packed responses with one reusable buffer.
    // ========================================================================
    std::vector<std::vector<char>> recv_buffers(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] > 0) {
            recv_buffers[i].resize(recv_sizes[i]);
            int ierr = MPI_Irecv_large(recv_buffers[i].data(),
                                       recv_sizes[i],
                                       MPI_CHAR,
                                       neighbor_ranks[i],
                                       303,
                                       MPI_COMM_WORLD,
                                       requests);
            abort_on_mpi_error(ierr, "MPI_Irecv_large(payload)", neighbor_ranks[i], recv_sizes[i]);
        }
    }

    std::vector<char> reusable_send_buffer;
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_sizes[i] == 0) {
            continue;
        }

        const size_t total_size = send_sizes[i];
        reusable_send_buffer.resize(total_size);
        char* send_ptr = reusable_send_buffer.data();

        for (const auto& req_info : requests_received[i]) {
            const int64_t local_idx = get_local_idx_for_request(req_info.morton);
            const size_t req_size =
                get_serialized_request_size_for_box(req_info, local_idx);

            SolveDataRequest<CoordType, DataType> solve_req;
            populate_solve_request(solve_req, req_info, local_idx);

            std::memcpy(send_ptr, &req_size, sizeof(size_t));
            send_ptr += sizeof(size_t);
            send_ptr = serialize(solve_req, send_ptr);
        }

        if (send_ptr != reusable_send_buffer.data() + total_size) {
            throw std::runtime_error(
                "gather_boxes_solve: send buffer size mismatch for rank " +
                std::to_string(neighbor_ranks[i]));
        }

        int ierr = MPI_Send_large(reusable_send_buffer.data(),
                                  total_size,
                                  MPI_CHAR,
                                  neighbor_ranks[i],
                                  303,
                                  MPI_COMM_WORLD);
        abort_on_mpi_error(ierr, "MPI_Send_large(payload)", neighbor_ranks[i], total_size);
    }

    if (!requests.empty()) {
        int ierr = MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        abort_on_mpi_error(ierr, "MPI_Waitall(payload)", -1, requests.size() * sizeof(MPI_Request));
        requests.clear();
    }

    // ========================================================================
    // Step 7: Deserialize received solve data into the local solve cache.
    // ========================================================================
    size_t num_solve_boxes = lvl.ghost_and_assisting_box_points_for_solve_map.size();
    lvl.ghost_and_assisting_boxes_for_solve.resize(num_solve_boxes);

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] == 0) continue;

        const char* ptr = recv_buffers[i].data();
        const char* buffer_end = ptr + recv_sizes[i];
        int boxes_received = 0;

        for (const auto& req_info : solve_requests_to_send[neighbor_ranks[i]]) {
            int64_t morton = req_info.morton;

            size_t req_size;
            std::memcpy(&req_size, ptr, sizeof(size_t));
            ptr += sizeof(size_t);

            auto idx_it = lvl.ghost_and_assisting_box_points_for_solve_map.find(morton);
            if (idx_it == lvl.ghost_and_assisting_box_points_for_solve_map.end()) {
                throw std::runtime_error(
                    "gather_boxes_solve: Morton " + std::to_string(morton) +
                    " not found in ghost_and_assisting_box_points_for_solve_map");
            }

            int64_t idx = idx_it->second;

            if (minimal_only) {
                SolveDataRequest<CoordType, DataType> temp;
                ptr = deserialize(temp, ptr);
                lvl.ghost_and_assisting_boxes_for_solve[idx].morton_index = temp.morton_index;
                lvl.ghost_and_assisting_boxes_for_solve[idx].source_rank = temp.source_rank;
                lvl.ghost_and_assisting_boxes_for_solve[idx].right_side = std::move(temp.right_side);
                lvl.ghost_and_assisting_boxes_for_solve[idx].left_side = std::move(temp.left_side);
                lvl.ghost_and_assisting_boxes_for_solve[idx].redundant_indices = std::move(temp.redundant_indices);
                lvl.ghost_and_assisting_boxes_for_solve[idx].skeleton_indices = std::move(temp.skeleton_indices);
            } else {
                ptr = deserialize(lvl.ghost_and_assisting_boxes_for_solve[idx], ptr);
            }

            boxes_received++;
        }

        if (boxes_received != send_counts[i]) {
            throw std::runtime_error(
                "gather_boxes_solve: Solve data count mismatch from rank " +
                std::to_string(neighbor_ranks[i]) + ": expected " +
                std::to_string(send_counts[i]) + ", got " +
                std::to_string(boxes_received));
        }

        if (ptr != buffer_end) {
            throw std::runtime_error(
                "gather_boxes_solve: Buffer size mismatch from rank " +
                std::to_string(neighbor_ranks[i]) + ": " +
                std::to_string(ptr - recv_buffers[i].data()) +
                " bytes consumed, " + std::to_string(recv_sizes[i]) +
                " bytes received");
        }
    }
}


template <typename CoordType, typename DataType>
static inline SolveDataRequest<CoordType, DataType>&
get_local_solve_entry_checked(TreeLevel<CoordType, DataType>& lvl,
                             std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                             int64_t morton)
{
    // Sanity check: must be locally owned at this level.
    // Uses your O(1) bracket-indexing API.
    BoxData<CoordType, DataType>* b = lvl.find_local_box(morton);
    if (!b) {
        throw std::runtime_error("solve update targets nonlocal morton " + std::to_string(morton));
    }

    // SolveData is stored in the same contiguous Morton bracket.
    const int64_t offset = morton - lvl.local_morton_start;
    if (offset < 0 || offset >= (int64_t)level_solve_data.size()) {
        throw std::runtime_error("local solve_data offset out of range for morton " + std::to_string(morton));
    }
    return level_solve_data[(size_t)offset];
}

template <typename CoordType, typename DataType>
static inline void apply_full_update_local(
    TreeLevel<CoordType, DataType>& lvl,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    int64_t morton,
    const std::vector<DataType>& upd)
{
    auto& dst = get_local_solve_entry_checked(lvl, level_solve_data, morton);

    if ((int64_t)dst.left_side.size() != (int64_t)upd.size()) {
        throw std::runtime_error("apply_full_update_local: size mismatch morton=" + std::to_string(morton) +
                                 " left_side=" + std::to_string(dst.left_side.size()) +
                                 " upd=" + std::to_string(upd.size()));
    }
    for (int64_t i = 0; i < (int64_t)upd.size(); ++i)
        dst.left_side[(size_t)i] += upd[(size_t)i];
}

template <typename CoordType, typename DataType>
static inline void apply_skel_update_local(
    TreeLevel<CoordType, DataType>& lvl,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    int64_t morton,
    const std::vector<DataType>& upd)
{
    auto& dst = get_local_solve_entry_checked(lvl, level_solve_data, morton);

    if ((int64_t)dst.skeleton_indices.size() != (int64_t)upd.size()) {
        throw std::runtime_error("apply_skel_update_local: size mismatch morton=" + std::to_string(morton) +
                                 " skel=" + std::to_string(dst.skeleton_indices.size()) +
                                 " upd=" + std::to_string(upd.size()));
    }
    for (int64_t i = 0; i < (int64_t)upd.size(); ++i) {
        const int64_t dof = dst.skeleton_indices[(size_t)i];
        dst.left_side[(size_t)dof] += upd[(size_t)i];
    }
}


/**
 * @brief Exchange PendingSolveUpdates with 1-hop neighbor ranks and apply to local solve_data.
 *
 * - Outgoing updates are partitioned by destination owner rank (computed from morton index).
 * - Nonlocal updates are sent as PendingSolveUpdates packets.
 * - Received packets are deserialized and accumulated (added) into local solve_data[level_index].
 *
 * Notes:
 * - A rank may receive contributions for the same target morton from multiple neighbors; we
 *   accumulate them all by addition when applying.
 * - pending is cleared at the end (including maps).
 *
 * @return Communication-only duration spent in MPI send/recv/wait operations.
 */
template <typename CoordType, typename DataType>
std::chrono::high_resolution_clock::duration transport_and_apply_solve_updates_onehop(
    ParallelTree<CoordType, DataType>* tree,
    TreeLevel<CoordType, DataType>& lvl,
    int level_index,
    PendingSolveUpdates<DataType>& pending,
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& solve_data)
{
    using clock = std::chrono::high_resolution_clock;
    clock::duration communication_time{};

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (level_index < 0 || level_index >= (int)solve_data.size())
        throw std::runtime_error("transport_and_apply_solve_updates_onehop: level_index out of range");

    auto& level_solve_data = solve_data[(size_t)level_index];

    const int64_t local_count = lvl.local_morton_end - lvl.local_morton_start + 1;
    if (local_count < 0 || (int64_t)level_solve_data.size() != local_count)
        throw std::runtime_error("transport_and_apply_solve_updates_onehop: solve_data[level] size mismatch");

    const uint32_t grid_size = 1u << level_index;

    auto owner_of_morton = [&](int64_t morton_idx) -> int {
        std::vector<uint64_t> single_box = { (uint64_t)morton_idx };
        std::vector<uint32_t> region =
            (tree->dimension == 2)
                ? morton::assign_to_processes_2d(single_box, lvl.num_active_processes, grid_size)
                : morton::assign_to_processes_3d(single_box, lvl.num_active_processes, grid_size);

        const int morton_region_id = (int)region[0];

        auto it = lvl.morton_to_rank.find(morton_region_id);
        if (it == lvl.morton_to_rank.end())
            throw std::runtime_error("transport_and_apply_solve_updates_onehop: morton_region_id not in morton_to_rank");
        return it->second;
    };

    // ========================================================================
    // Step 1: Build the 1-hop neighbor rank list from the solve ghost/assist map.
    // ========================================================================
    std::unordered_set<int> neigh_set;
    neigh_set.reserve(lvl.ghost_and_assisting_box_points_for_solve_map.size());

    for (const auto& [morton, idx] : lvl.ghost_and_assisting_box_points_for_solve_map) {
        (void)idx;
        const int owner_rank = owner_of_morton(morton);
        if (owner_rank != rank) neigh_set.insert(owner_rank);
    }

    std::vector<int> neighbor_ranks;
    neighbor_ranks.reserve(neigh_set.size());
    for (int r : neigh_set) neighbor_ranks.push_back(r);
    std::sort(neighbor_ranks.begin(), neighbor_ranks.end());

    // ========================================================================
    // Step 2: Build per-destination solve-update packets and apply local ones immediately.
    // ========================================================================
    std::unordered_map<int, PendingSolveUpdates<DataType>> send_to;

    for (auto& [morton, upd] : pending.full_updates) {
        if (lvl.find_local_box(morton)) {
            apply_full_update_local(lvl, level_solve_data, morton, upd);
        } else {
            const int owner = owner_of_morton(morton);

            // must be in derived 1-hop neighbor set, else protocol would deadlock
            if (!neigh_set.count(owner)) {
                throw std::runtime_error(
                    "transport_and_apply_solve_updates_onehop: pending FULL update targets non-1hop rank " +
                    std::to_string(owner) + " morton=" + std::to_string(morton));
            }
            send_to[owner].full_updates.emplace(morton, std::move(upd));
        }
    }

    for (auto& [morton, upd] : pending.skel_updates) {
        if (lvl.find_local_box(morton)) {
            apply_skel_update_local(lvl, level_solve_data, morton, upd);
        } else {
            const int owner = owner_of_morton(morton);
            if (!neigh_set.count(owner)) {
                throw std::runtime_error(
                    "transport_and_apply_solve_updates_onehop: pending SKEL update targets non-1hop rank " +
                    std::to_string(owner) + " morton=" + std::to_string(morton));
            }
            send_to[owner].skel_updates.emplace(morton, std::move(upd));
        }
    }

    clear_pending_solve_updates_memory(pending);

    // ========================================================================
    // Step 3: Exchange per-neighbor payload sizes (tag 700).
    // ========================================================================
    std::vector<uint64_t> recv_sizes(neighbor_ranks.size(), 0);
    std::vector<MPI_Request> reqs;
    reqs.reserve(neighbor_ranks.size());

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Request r;
        MPI_Irecv(&recv_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i], 700, MPI_COMM_WORLD, &r);
        reqs.push_back(r);
    }

    std::vector<std::vector<char>> send_buffers(neighbor_ranks.size());

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        const int peer = neighbor_ranks[i];
        auto it = send_to.find(peer);

        uint64_t nbytes = 0;
        if (it != send_to.end()) {
            nbytes = bytes_pending(it->second);
            send_buffers[i].resize((size_t)nbytes);
            if (nbytes > 0) {
                char* end = serialize(it->second, send_buffers[i].data());
                if ((uint64_t)(end - send_buffers[i].data()) != nbytes)
                    throw std::runtime_error("transport_and_apply_solve_updates_onehop: serialize size mismatch");
            }
        }

        // Always send a size to every neighbor in the derived list.
        auto comm_start = clock::now();
        MPI_Send(&nbytes, 1, MPI_UINT64_T, peer, 700, MPI_COMM_WORLD);
        communication_time += (clock::now() - comm_start);
    }

    if (!reqs.empty()) {
        auto comm_start = clock::now();
        MPI_Waitall((int)reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
        communication_time += (clock::now() - comm_start);
    }
    reqs.clear();

    // ========================================================================
    // Step 4: Post all payload receives, then send the serialized solve-update payloads (tag 701).
    // ========================================================================
    std::vector<std::vector<char>> recv_buffers(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] > 0) {
            recv_buffers[i].resize((size_t)recv_sizes[i]);
        }
    }

    auto comm_start = clock::now();
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_sizes[i] > 0) {
            int err = MPI_Irecv_large(recv_buffers[i].data(),
                                      (size_t)recv_sizes[i],
                                      MPI_CHAR,
                                      neighbor_ranks[i],
                                      701,
                                      MPI_COMM_WORLD,
                                      reqs);
            if (err != MPI_SUCCESS) {
                char errbuf[MPI_MAX_ERROR_STRING];
                int errlen = 0;
                MPI_Error_string(err, errbuf, &errlen);
                throw std::runtime_error(
                    "transport_and_apply_solve_updates_onehop: payload Irecv_large failed for rank " +
                    std::to_string(neighbor_ranks[i]) + " with error " + std::string(errbuf, errlen));
            }
        }
    }

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_buffers[i].empty()) {
            continue;
        }

        int err = MPI_Send_large(send_buffers[i].data(),
                                 send_buffers[i].size(),
                                 MPI_CHAR,
                                 neighbor_ranks[i],
                                 701,
                                 MPI_COMM_WORLD);
        if (err != MPI_SUCCESS) {
            char errbuf[MPI_MAX_ERROR_STRING];
            int errlen = 0;
            MPI_Error_string(err, errbuf, &errlen);
            throw std::runtime_error(
                "transport_and_apply_solve_updates_onehop: payload Send_large failed for rank " +
                std::to_string(neighbor_ranks[i]) + " with error " + std::string(errbuf, errlen));
        }
    }

    if (!reqs.empty()) {
        int err = MPI_Waitall((int)reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
        if (err != MPI_SUCCESS) {
            char errbuf[MPI_MAX_ERROR_STRING];
            int errlen = 0;
            MPI_Error_string(err, errbuf, &errlen);
            throw std::runtime_error(
                "transport_and_apply_solve_updates_onehop: payload Waitall failed with error " +
                std::string(errbuf, errlen));
        }
        reqs.clear();
    }
    communication_time += (clock::now() - comm_start);

    // ========================================================================
    // Step 5: Deserialize incoming updates and apply them to local solve_data.
    // ========================================================================
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_buffers[i].empty()) continue;

        PendingSolveUpdates<DataType> incoming;
        const char* ptr = recv_buffers[i].data();
        const char* end = ptr + recv_buffers[i].size();

        ptr = deserialize(incoming, ptr);
        if (ptr != end)
            throw std::runtime_error("transport_and_apply_solve_updates_onehop: recv buffer parse mismatch");

        for (const auto& [morton, upd] : incoming.full_updates) {
            apply_full_update_local(lvl, level_solve_data, morton, upd);
        }
        for (const auto& [morton, upd] : incoming.skel_updates) {
            apply_skel_update_local(lvl, level_solve_data, morton, upd);
        }
    }

    return communication_time;
}

// void send_large_buffer(const std::vector<char>& buffer, int dest, int tag, MPI_Comm comm) {
//     const int64_t chunk_size = 1'000'000'000;  // 1GB chunks
//     int64_t total_size = buffer.size();
    
//     // Send metadata
//     MPI_Send(&total_size, 1, MPI_INT64_T, dest, tag, comm);
    
//     // Send in chunks
//     for (int64_t offset = 0; offset < total_size; offset += chunk_size) {
//         int current_size = std::min(chunk_size, total_size - offset);
//         MPI_Send(buffer.data() + offset, current_size, MPI_CHAR, dest, tag + 1, comm);
//     }
// }

// void recv_large_buffer(std::vector<char>& buffer, int src, int tag, MPI_Comm comm) {
//     const int64_t chunk_size = 1'000'000'000;
//     int64_t total_size;
//     MPI_Status status;
    
//     // Receive metadata
//     MPI_Recv(&total_size, 1, MPI_INT64_T, src, tag, comm, &status);
//     buffer.resize(total_size);
    
//     // Receive chunks
//     for (int64_t offset = 0; offset < total_size; offset += chunk_size) {
//         int current_size = std::min(chunk_size, total_size - offset);
//         MPI_Recv(buffer.data() + offset, current_size, MPI_CHAR, src, tag + 1, comm, &status);
//     }
// }

} // namespace fmm

#endif // SERIALIZATION_HPP
