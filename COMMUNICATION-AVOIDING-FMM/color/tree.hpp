#ifndef TREE_HPP
#define TREE_HPP

#include <vector>
#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <mpi.h>
#include "morton.hpp"
#include <unordered_map>
#include <unordered_set>
#include <omp.h>


namespace fmm {

// Forward declarations
template<typename CoordType, typename DataType> struct BoxData;
template<typename CoordType, typename DataType> struct SolveDataRequest;
template<typename CoordType, typename DataType> struct TreeLevel;
template<typename CoordType, typename DataType> struct ParallelTree;
template<typename DataType> struct PendingFactorUpdates;


enum class EdgeKind : uint8_t { Near, Far, Diag };

struct ReplaceKey {
    int64_t target_box;   // G (the box to be updated on its owner rank)
    int64_t neighbor_box; // B (the eliminated box whose skeleton interaction is stored)
    bool operator==(const ReplaceKey& o) const noexcept {
        return target_box == o.target_box && neighbor_box == o.neighbor_box;
    }
};

struct EdgeKey {
    int64_t lo;
    int64_t hi;
    EdgeKind kind;
    bool operator==(const EdgeKey& o) const noexcept {
        return lo == o.lo && hi == o.hi && kind == o.kind;
    }
};

struct ReplaceKeyHash {
    size_t operator()(const ReplaceKey& k) const noexcept {
        // simple combine
        size_t h1 = std::hash<int64_t>{}(k.target_box);
        size_t h2 = std::hash<int64_t>{}(k.neighbor_box);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1<<6) + (h1>>2));
    }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& k) const noexcept {
        size_t h1 = std::hash<int64_t>{}(k.lo);
        size_t h2 = std::hash<int64_t>{}(k.hi);
        size_t h3 = std::hash<uint8_t>{}(static_cast<uint8_t>(k.kind));
        size_t h = h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1<<6) + (h1>>2));
        return h ^ (h3 + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
    }
};

template<typename DataType>
struct DenseBlock {
    int64_t rows = 0;
    int64_t cols = 0;
    std::vector<DataType> data; // column-major, size rows*cols
};


template<typename DataType>
struct PendingFactorUpdates {
    // Step 7: REPLACE blocks in near-field of target_box for neighbor_box
    std::unordered_map<ReplaceKey, DenseBlock<DataType>, ReplaceKeyHash> replace_blocks;

    // Step 9: ADD deltas for canonical edge (lo,hi). Stored as ΔM_{lo→hi} with dims (n_hi x n_lo)
    std::unordered_map<EdgeKey, DenseBlock<DataType>, EdgeKeyHash> accumulated_deltas;
};


/**
 * @brief Pending solve-phase updates destined for (possibly remote) neighbor boxes.
 *
 * We separate two cases to match your forward-elimination Step 3 semantics:
 *  - full_updates[morton]: add into neighbor.left_side[0..n-1]
 *  - skel_updates[morton]: add into neighbor.left_side at the neighbor's skeleton DOFs
 *                          in the neighbor's skeleton_indices order.
 *
 * NOTE: For skel_updates to be unambiguous, sender and receiver must agree on the
 * ordering of skeleton_indices for that box. (current code uses that ordering.)
 */
template <typename DataType>
struct PendingSolveUpdates {
    std::unordered_map<int64_t, std::vector<DataType>> full_updates;
    std::unordered_map<int64_t, std::vector<DataType>> skel_updates;
};

template <typename CoordType, typename DataType>
struct FactorizationThreadScratch {
    std::vector<DataType> workspace;
    std::vector<DataType> x_nn_full;
    std::vector<DataType> sketch_storage;
    int64_t workspace_rows = 0;
    int64_t workspace_cols = 0;

    std::vector<DataType> x_bb;
    std::vector<int64_t> neighbor_point_counts;
    std::vector<DataType> a_ns_all;
    std::vector<DataType> a_sn_all;

    std::vector<DataType> temp1;
    std::vector<DataType> temp2;
    std::vector<DataType> temp3;
    std::vector<DataType> temp4;
    std::vector<DataType> x_rs_original;
    std::vector<DataType> x_ns_update;
    std::vector<DataType> x_sn_update;
    std::vector<DataType> update_buffer;
    std::vector<DataType> eval_buffer;
    std::vector<CoordType> coord_buffer;
};


/**
 * @brief Cached phased pair schedule for process-neighbor communication at one level.
 *
 * Each phase is a deterministic matching on the spatial process-neighbor graph,
 * so a rank appears in at most one communicating pair per phase. This is used
 * to run robust pairwise MPI_Sendrecv-style exchanges without many-peer
 * progress assumptions.
 */
struct ProcessPairCommunicationSchedule {
    std::vector<std::vector<std::pair<int, int>>> phase_edges;
    std::vector<int> local_phase_peers;
    bool initialized;

    ProcessPairCommunicationSchedule() : initialized(false) {}

    void reset() {
        phase_edges.clear();
        local_phase_peers.clear();
        initialized = false;
    }
};

/**
 * @brief Matrix storage with format metadata (FMM3D.pdf Section 2.4)
 * 
 * Encapsulates matrix data with information about storage format.
 * Supports different matrix representations (full, Cholesky factored, etc.)
 * 
 * Memory Layout:
 * - Column-major (Fortran/BLAS style): A[i,j] = data[i + j*lda]
 * - lda = leading dimension (typically == rows, but can be larger for padding)
 */
template<typename DataType>
struct MatrixStorage {
    std::vector<DataType> data;   ///< Matrix data (column-major)
    int64_t rows;                 ///< Number of rows
    int64_t cols;                 ///< Number of columns
    int64_t lda;                  ///< Leading dimension (stride between columns)
    
    /**
     * @brief Matrix format enumeration
     * 
     * Tracks the current state/format of the matrix data:
     * - NONE: Not allocated
     * - FULL: General dense matrix (no special structure)
     * - CHOLESKY_L: Lower triangular Cholesky factor L (A = L*L^T or L*L^H)
     * - CHOLESKY_U: Upper triangular Cholesky factor U (A = U^T*U or U^H*U)
     * - LU_FACTORED: LU factorization with partial pivoting
     * - INVERSE: Matrix inverse (explicitly computed)
     */
    enum Format {
        NONE,                 ///< Not allocated
        FULL,                 ///< General dense matrix
        CHOLESKY_L,           ///< Lower triangular Cholesky factor
        CHOLESKY_U,           ///< Upper triangular Cholesky factor
        LU_FACTORED,          ///< LU factorization with pivoting
        BUNCH_KAUFMAN,        ///< Bunch-Kaufman symmetric factorization (A = P*L*D*L^T*P^T)
        INVERSE               ///< Explicit inverse
    };
    
    Format format;            ///< Current format of the matrix
    
    // Constructor
    MatrixStorage() : rows(0), cols(0), lda(0), format(NONE) {}
    
    // Destructor - no manual cleanup needed with vector
    ~MatrixStorage() = default;
    
    // Delete copy constructor and assignment
    MatrixStorage(const MatrixStorage&) = delete;
    MatrixStorage& operator=(const MatrixStorage&) = delete;

    // Add move constructor and move assignment
    MatrixStorage(MatrixStorage&&) = default;
    MatrixStorage& operator=(MatrixStorage&&) = default;
    
    /**
     * @brief Allocate matrix storage
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param fmt Initial format (default: FULL)
     * @param leading_dim Leading dimension (default: nrows)
     */
    void allocate(int64_t nrows, int64_t ncols, Format fmt = FULL, int64_t leading_dim = 0) {
        // rows = nrows;
        // cols = ncols;
        // lda = (leading_dim > 0) ? leading_dim : nrows;
        // format = fmt;
        
        // if (rows > 0 && cols > 0) {
        //     data.resize(lda * cols, DataType{});  // Zero-initialize
        //     data.shrink_to_fit();
        // } else {
        //     data.clear();
        // }
        rows = nrows;
        cols = ncols;
        lda = (leading_dim > 0) ? leading_dim : nrows;
        format = fmt;
        
        if (rows > 0 && cols > 0) {
            int64_t new_size = lda * cols;
            
            if (new_size < static_cast<int64_t>(data.size())) {
                // Shrinking: swap idiom guarantees memory is released
                std::vector<DataType>(new_size, DataType{}).swap(data);
            } else {
                // Growing or same size: resize is fine
                data.resize(new_size, DataType{});
            }
        } else {
            // Free everything
            std::vector<DataType>().swap(data);
        }
    }

    /**
     * @brief Take ownership of an already-built matrix buffer.
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param owned_data Matrix buffer to adopt
     * @param fmt Matrix format
     * @param leading_dim Leading dimension (default: nrows)
     */
    void set_owned(int64_t nrows,
                   int64_t ncols,
                   std::vector<DataType>&& owned_data,
                   Format fmt = FULL,
                   int64_t leading_dim = 0) {
        rows = nrows;
        cols = ncols;
        lda = (leading_dim > 0) ? leading_dim : nrows;
        format = fmt;

        const int64_t expected_size = lda * cols;
        if (static_cast<int64_t>(owned_data.size()) != expected_size) {
            throw std::runtime_error("MatrixStorage::set_owned size mismatch");
        }

        data = std::move(owned_data);
    }
    
    /**
     * @brief Check if matrix is allocated
     * @return true if data is allocated
     */
    bool is_allocated() const {
        return !data.empty() && rows > 0 && cols > 0;
    }
    
    /**
     * @brief Get matrix element (column-major indexing)
     * @param i Row index
     * @param j Column index
     * @return Reference to matrix element
     */
    DataType& operator()(int64_t i, int64_t j) {
        return data[i + j * lda];
    }
    
    const DataType& operator()(int64_t i, int64_t j) const {
        return data[i + j * lda];
    }
    
    /**
     * @brief Get total allocated memory size
     * @return Number of DataType elements allocated
     */
    int64_t size() const {
        return lda * cols;
    }
};


/**
 * @brief Helper struct for MPI point data exchange
 */
template<typename CoordType>
struct PointDataRequest {
    int64_t morton_index;     ///< Box Morton index
    int source_rank;          ///< MPI rank to request from
    std::vector<CoordType> coords;  ///< Point coordinates (dim × num_points, column-major)
    std::vector<int64_t> skel_indices; ///< Skeleton indices
    bool on_boundary; ///< if this box is on the boundary
    
    PointDataRequest() : morton_index(-1), source_rank(-1) {}
};

/**
 * @brief MPI solve data exchange for assisting boxes
 * 
 * Stores both vector data (right_side, left_side) and factorization matrices
 * for boxes on other processes needed during the solve phase.
 * 
 * For ghost boxes: All fields populated during factorization phase
 * For assisting boxes: Only needed fields populated on-demand during solve
 * 
 * This enables distributed solve where:
 * - V^{-1} sweep: Apply L^{-1}, U_T^{-1} using T, X_SR, X_NR, X_RR
 * - W^{-1} sweep: Apply U^{-1}, L_T^{-1} using T, X_RS, X_RN, X_RR
 * - Diagonal solve: Use schur_complement at coarsest level
 */
template<typename CoordType, typename DataType>
struct SolveDataRequest {
    // ===== Identification =====
    int64_t morton_index;              ///< Box Morton index
    int source_rank;                   ///< MPI rank to request from
    
    // ===== Vector Data (always needed for solve) =====
    std::vector<DataType> right_side;  ///< RHS vector (modified during V^{-1})
    std::vector<DataType> left_side;   ///< Solution vector (modified during W^{-1})

    std::vector<int64_t> redundant_indices;  ///< Redundant DOF indices R (local to box)
    std::vector<int64_t> skeleton_indices;   ///< Skeleton DOF indices S (local to box)
    std::vector<int64_t> one_hop; ///< One-hop neighbors
    std::vector<int64_t> use_full_set; ///< whether to use full set or skeleton
    
    // ===== Factorization Matrices (populated on-demand) =====
    
    MatrixStorage<DataType> interpolation_matrix;
    

    MatrixStorage<DataType> X_RR;
    std::vector<int> X_RR_pivots;
    
 
    MatrixStorage<DataType> X_RS;
    

    MatrixStorage<DataType> X_SR;
    

    MatrixStorage<DataType> schur_complement;
    

    MatrixStorage<DataType> X_RN;
    

    MatrixStorage<DataType> X_NR;
    
    // ===== Constructors =====
    SolveDataRequest() : morton_index(-1), source_rank(-1) {}
    
    /**
     * @brief Initialize with vector data only (minimal mode)
     * @param box_morton Box Morton index
     * @param rank Source MPI rank
     * @param num_points Number of DOFs
     */
    void initialize(int64_t box_morton, int rank, int64_t num_points) {
        morton_index = box_morton;
        source_rank = rank;
        right_side.resize(num_points, DataType{0.0});
        left_side.resize(num_points, DataType{0.0});
    }

    /**
    * @brief Used by assisting boxes, copy only the necessary info from BoxData
    * @param box Source box with factorization data
    * 
    * Used during solve phase as an assisting box when only partial data is needed.
    */
    void copy_minimal_from_box(const BoxData<CoordType, DataType>& box, const SolveDataRequest<CoordType, DataType>& rhs_info) {
        // Copy index vectors
        redundant_indices = box.redundant_indices;
        skeleton_indices = box.skeleton_indices;
        right_side = rhs_info.right_side;
        left_side = rhs_info.left_side;
        
        // No matrix data copied in minimal mode
    }
    
    /**
    * @brief Copy factorization matrices and index data from BoxData
    * @param box Source box with factorization data
    * 
    * Used when we need full factorization data for a remote box.
    * Typically called for ghost boxes or critical assisting boxes.
    */
    void copy_factorization_from_box(const BoxData<CoordType, DataType>& box, const SolveDataRequest<CoordType, DataType>& rhs_info) {
        // Copy index vectors
        redundant_indices = box.redundant_indices;
        skeleton_indices = box.skeleton_indices;
        one_hop = box.one_hop;
        use_full_set = box.use_full_set;
        right_side = rhs_info.right_side;
        left_side = rhs_info.left_side;
        
        // Copy all factorization matrices
        if (box.interpolation_matrix.is_allocated()) {
            interpolation_matrix.allocate(
                box.interpolation_matrix.rows,
                box.interpolation_matrix.cols,
                box.interpolation_matrix.format);
            interpolation_matrix.data = box.interpolation_matrix.data;
        }
        
        if (box.X_RR.is_allocated()) {
            X_RR.allocate(box.X_RR.rows, box.X_RR.cols, box.X_RR.format);
            X_RR.data = box.X_RR.data;
        }
        X_RR_pivots = box.X_RR_pivots;
        
        if (box.X_RS.is_allocated()) {
            X_RS.allocate(box.X_RS.rows, box.X_RS.cols, box.X_RS.format);
            X_RS.data = box.X_RS.data;
        }
        
        if (box.X_SR.is_allocated()) {
            X_SR.allocate(box.X_SR.rows, box.X_SR.cols, box.X_SR.format);
            X_SR.data = box.X_SR.data;
        }
        
        if (box.schur_complement.is_allocated()) {
            schur_complement.allocate(
                box.schur_complement.rows,
                box.schur_complement.cols,
                box.schur_complement.format);
            schur_complement.data = box.schur_complement.data;
        }
        
        if (box.X_RN.is_allocated()) {
            X_RN.allocate(box.X_RN.rows, box.X_RN.cols, box.X_RN.format);
            X_RN.data = box.X_RN.data;
        }
        
        if (box.X_NR.is_allocated()) {
            X_NR.allocate(box.X_NR.rows, box.X_NR.cols, box.X_NR.format);
            X_NR.data = box.X_NR.data;
        }
    }
    
    /**
     * @brief Get number of DOFs
     */
    int64_t size() const {
        return right_side.size();
    }
    
    /**
     * @brief Check if initialized (vector data allocated)
     */
    bool is_initialized() const {
        return !right_side.empty();
    }
    
    /**
     * @brief Check if full factorization data is available
     */
    bool has_factorization_data() const {
        return X_RR.is_allocated();  // X_RR is required for all boxes
    }
};



/**
 * @brief Modified interaction block (FMM3D.pdf Section 2.4.4)
 * 
 * Stores the updated near-field interaction matrices after redundant DOF elimination.
 * These are the "tilda" matrices from equations (7), (8):
 * 
 *   tilda A_SN = A_SN - X_SR * X_RR^{-1} * X_RN     (if stored)
 *   tilda A_NS = A_NS - X_NR * X_RR^{-1} * X_RS     (default for symmetric)
 * 
 * Storage Convention (as per user specification):
 * - For symmetric matrices: Store A_NS only (A_SN = nullptr, use A_NS^T)
 * - For nonsymmetric matrices: Store both A_NS and A_SN
 * 
 * These modified blocks accumulate during the upward pass and are used
 * in the solve phase to compute interactions with reduced DOFs.
 */
template<typename DataType>
struct ModifiedBlock {
    int64_t neighbor_morton;           ///< Morton index of neighbor box
    MatrixStorage<DataType> A_SN;      ///< Skeleton-to-neighbor (k × n_neighbor), null if symmetric
    MatrixStorage<DataType> A_NS;      ///< Neighbor-to-skeleton (n_neighbor × k), default storage
    
    ModifiedBlock() : neighbor_morton(-1) {}
    ~ModifiedBlock() = default;
    
    // Delete copy constructor and assignment
    ModifiedBlock(const ModifiedBlock&) = delete;
    ModifiedBlock& operator=(const ModifiedBlock&) = delete;

    // Add move constructor and move assignment
    ModifiedBlock(ModifiedBlock&&) = default;
    ModifiedBlock& operator=(ModifiedBlock&&) = default;
};

/**
 * @brief Box data structure (FMM3D.pdf Section 4.1)
 * 
 * Contains all geometric, hierarchical, and matrix data for a single box.
 * Supports both local boxes (owned by process) and ghost boxes (neighbors from other processes).
 * 
 * Data Organization (per FMM3D.pdf):
 * - Section 4.1: Basic box geometry and hierarchy
 * - Section 2.3: First-level elimination data (X blocks)
 * - Section 2.4: Second-level elimination data (Schur complement, modified blocks)
 * - Section 4.2: Ghost box data structure
 * 
 * @tparam CoordType Floating-point type for geometric data
 * @tparam DataType Type for matrix data (can be complex)
 */
template<typename CoordType, typename DataType = CoordType>
struct BoxData {
    // ===== Box Identification (Section 4.1) =====
    int64_t morton_index;        ///< Morton index (64-bit global identifier)
    int32_t level;               ///< Tree level (0 = root, L = leaf)
    int32_t grid_coords[3];      ///< Grid coordinates [x, y, z]
    
    // ===== Box Geometry =====
    CoordType bounds[6];         ///< Bounds [xmin, xmax, ymin, ymax, zmin, zmax]
    CoordType center[3];         ///< Box center coordinates
    CoordType size;              ///< Box size (same in all dimensions for cubic boxes)
    
    // ===== Hierarchical Information =====
    int64_t parent_morton;       ///< Parent box Morton index (-1 if root)
    int64_t children_morton[8];  ///< Children Morton indices (-1 if not present)
    int32_t num_children;        ///< Number of children (0 for leaf boxes)
    
    // ===== Point Information =====
    std::vector<int64_t> point_indices;  ///< Global indices of points (int64_t per Section 3.4.1)
    int64_t num_points;                  ///< Number of points in this box
    std::vector<CoordType> point_coords; ///< Point coordinates (column-major: dim × num_points)
    
    // ===== DOF Information (populated during compression) =====
    std::vector<int64_t> redundant_indices;  ///< Redundant DOF indices R (local to box)
    std::vector<int64_t> skeleton_indices;   ///< Skeleton DOF indices S (local to box)

    
    // ===== Neighbor Information =====
    std::vector<int64_t> one_hop; ///< Morton indices of 1-hop neighbors
    std::vector<int64_t> use_full_set; ///< whether to use full set or skeleton for each neighbor
    std::vector<int64_t> two_hop; ///< Morton indices of 2-hop neighbors
    bool on_boundary; ///< if this box is on the boundary
    
    // ===== Matrix Data (Section 2.3-2.4) =====
    
    /**
     * Interpolation matrix T (Section 2.2):
     * Maps redundant DOFs to skeleton DOFs via proxy points
     * Size: k × |R|
     * Format: FULL (never factored)
     */
    MatrixStorage<DataType> interpolation_matrix;
    
    /**
     * X_RR block (Section 2.3.1):
     * Self-interaction of redundant DOFs
     * Size: |R| × |R|
     * Format: Initially FULL, then CHOLESKY_L after factorization
     * Usage: Solve X_RR^{-1} * b using the chosen factorization
     */
    MatrixStorage<DataType> X_RR;
    std::vector<int> X_RR_pivots;   ///< Pivot vector for LU_FACTORED X_RR
    
    /**
     * X_RS block (Section 2.3.3):
     * Redundant-to-skeleton interaction
     * Size: |R| × k
     * Format: FULL
     * Note: nullptr if matrix is symmetric (use X_SR^T instead)
     */
    MatrixStorage<DataType> X_RS;
    
    /**
     * X_SR block (Section 2.3.3):
     * Skeleton-to-redundant interaction
     * Size: k × |R|
     * Format: FULL
     */
    MatrixStorage<DataType> X_SR;
    
    /**
     * Schur complement S (Equation 6):
     * S = A_SS - X_SR * X_RR^{-1} * X_RS
     * Size: k × k
     * Format: Initially FULL, then CHOLESKY_L at higher levels
     */
    MatrixStorage<DataType> schur_complement;
    
    /**
     * X_RN blocks (Section 2.3.3):
     * Redundant-to-neighbor interactions for each near neighbor
     * Array of size num_near_neighbors
     * Each block: |R| × |N_i| where |N_i| is neighbor's skeleton size
     * Format: FULL
     * Note: nullptr entries if symmetric (use X_NR^T instead)
     */
    MatrixStorage<DataType> X_RN;
    
    /**
     * X_NR blocks (Section 2.3.3):
     * Neighbor-to-redundant interactions for each near neighbor
     * Array of size num_near_neighbors
     * Each block: |N_i| × |R|
     * Format: FULL
     */
    MatrixStorage<DataType> X_NR;
    std::vector<int64_t> deferred_xnn_neighbor_point_counts;
    std::vector<DataType> deferred_xnn_temp2;

    
    // ===== Near-Field and Far-Field Interactions =====
    
    // ===== Interaction Index Maps (for fast lookup) =====
    std::unordered_map<int64_t, int64_t> near_field_interaction_map;  ///< neighbor_morton → index in near_field_modified_interactions
    std::unordered_map<int64_t, int64_t> far_field_interaction_map;   ///< neighbor_morton → index in far_field_modified_interactions
    std::unordered_map<int64_t, int64_t> near_field_interaction_map_nonsymmetry; 
    std::unordered_map<int64_t, int64_t> far_field_interaction_map_nonsymmetry;   

    /**
     * Near-field blocks (Section 2.4.3):
     * Stores kernel evaluations with 1-hop neighbors
     * Initially contains unmodified A matrices, will be updated to Ã during elimination
     */
    std::vector<ModifiedBlock<DataType>> near_field_modified_interactions;
    int64_t num_near_field_interactions;  // Typically == num_near_neighbors
    
    /**
     * Far-field blocks (for proxy-based compression):
     * Stores kernel evaluations with 2-hop neighbors (used for ID computation)
     * Will be compressed via interpolative decomposition
     */
    std::vector<ModifiedBlock<DataType>> far_field_modified_interactions;
    int64_t num_far_field_interactions;
    
    // ===== Constructors and Destructor =====
    
    BoxData();
    ~BoxData() = default;  // Vector handles cleanup automatically
    
    // Delete copy constructor and assignment
    BoxData(const BoxData&) = delete;
    BoxData& operator=(const BoxData&) = delete;

    // Add move constructor and move assignment
    BoxData(BoxData&&) = default;
    BoxData& operator=(BoxData&&) = default;
};

/**
 * @brief Tree level data structure (FMM3D.pdf Section 4.4.1)
 * 
 * Each level contains local and ghost boxes, Morton range, and active processes.
 * 
 * @tparam CoordType Floating-point type for geometric data
 * @tparam DataType Type for matrix data
 */
template<typename CoordType, typename DataType = CoordType>
struct TreeLevel {
    int32_t level;                   ///< Level number (0 = root)
    int64_t num_boxes_global;        ///< Total number of boxes at this level (global)
    int64_t num_boxes_local;         ///< Number of boxes owned by this process
    int dimension;                   ///< Dimension of problem

    // ========== Ordering ==========
    std::vector<int64_t> blue;
    std::vector<int64_t> orange;
    std::vector<int64_t> purple;
    std::vector<int64_t> green;
    
    
    
    // ========== Box Storage (Section 4.4.1) ==========
    std::vector<int64_t> ghost_id;
    std::vector<int64_t> interior_id;
    std::vector<int64_t> boundary_id;

    std::unordered_map<int64_t, int64_t> assisting_box_points_for_kernel_evaluation;
    std::unordered_map<int64_t, int64_t> ghost_id_to_index;
    
    std::vector<BoxData<CoordType, DataType>> local_boxes;  ///< Local boxes (contiguous array)
    std::vector<std::vector<int64_t>> solve_neighbor_size; ///< solve neighbor size
    std::vector<BoxData<CoordType, DataType>> ghost_boxes;  ///< Ghost boxes (1-hop neighbors)
    std::vector<PointDataRequest<CoordType>> assisting_boxes;

    // for shared memory sync
    std::unordered_map<int64_t, omp_lock_t*> box_locks;

    std::unordered_set<int64_t> eliminated_boxes;

    /**
     * @brief Map from Morton index to assisting_boxes_for_solve index
     * 
     * Used to quickly find assisting box data during solve.
     * Similar to assisting_box_points_for_kernel_evaluation but for solve vectors.
     */
    std::unordered_map<int64_t, int64_t> ghost_and_assisting_box_points_for_solve_map;
     /**
     * @brief Assisting boxes for solve phase
     * 
     * Stores partial vector data (right_side, left_side) for boxes on other processes
     * that are needed for coupling during V^{-1} and W^{-1} sweeps.
     * 
     * Example: Box B eliminates redundant DOFs, which couples to neighbor N.
     * If N is on another process (ghost or assisting box) it goes here.
     */
    std::vector<SolveDataRequest<CoordType, DataType>> ghost_and_assisting_boxes_for_solve;
    /**
     * @brief indicate where it is a ghost box or an assisting box
     * 
     * 
     */
    std::vector<int> is_ghost_solve; // true for ghost, false for assisting box
    
    // ========== Morton Range (Section 4.4.2) ==========
    int64_t local_morton_start;      ///< Start of Morton range for this process
    int64_t local_morton_end;        ///< End of Morton range (inclusive)
    
    // ========== Process Reduction Support (Section 5.3) ==========
    int32_t num_active_processes;              ///< Number of active processes at this level
    std::vector<int32_t> active_process_ranks; ///< Array of active process ranks
    bool is_process_active;                    ///< Whether this process is active at this level
    
    // ========== Process Morton ID Mapping (Spatial Abstraction Layer) ==========
    std::unordered_map<int, int> morton_to_rank; ///< Morton region ID → MPI rank
    std::unordered_map<int, int> rank_to_morton; ///< MPI rank → Morton region ID
    int my_morton_id;                  ///< This process's Morton region (-1 if inactive)
    ProcessPairCommunicationSchedule process_pair_schedule; ///< Cached one-peer-per-phase schedule for deterministic pairwise communication
    
    // ========== Inter-Level Communication for Reduction ==========
    int parent_level_owner;            ///< MPI rank that owns parent (if inactive at L-1)
    std::vector<int> children_senders; ///< MPI ranks sending child data (if active at L-1)
    
    // ========== Constructors and Destructor ==========
    TreeLevel();
    ~TreeLevel() = default;  // Vector handles cleanup automatically
    
    // Delete copy constructor and assignment
    TreeLevel(const TreeLevel&) = delete;
    TreeLevel& operator=(const TreeLevel&) = delete;

    // Add move constructor and move assignment
    TreeLevel(TreeLevel&&) = default;
    TreeLevel& operator=(TreeLevel&&) = default;
    
    // ========== Member Functions ==========
    
    /**
     * @brief Find local box by Morton index (O(1) lookup)
     * @param morton_index Global Morton index
     * @return Pointer to box, or nullptr if not found
     */
    BoxData<CoordType, DataType>* find_local_box(int64_t morton_index);
    
    /**
     * @brief Find ghost box by Morton index (O(n) search)
     * @param morton_index Global Morton index
     * @return Pointer to box, or nullptr if not found
     */
    BoxData<CoordType, DataType>* find_ghost_box(int64_t morton_index);

     /**
     * @brief check if the morton index is in process
     * @param morton_index Global Morton index
     * @return true if in process, else return false
     */
    bool is_box_on_process(int64_t morton_index)
    {
        return (morton_index >= local_morton_start && morton_index <= local_morton_end);
    }


};

/**
 * @brief Calculate the total memory usage of a BoxData instance in bytes.
 *
 * Accounts for all dynamically allocated storage (vectors, matrix data,
 * hash maps, modified blocks) as well as the fixed-size fields of the struct.
 *
 * @tparam CoordType Floating-point type for geometric data
 * @tparam DataType  Type for matrix data (can be complex)
 * @param box        The BoxData instance to measure
 * @return           Total estimated memory in bytes
 */
template<typename CoordType, typename DataType>
size_t calculate_box_data_size(const BoxData<CoordType, DataType>& box) {
    size_t total = 0;

    // ===== Fixed-size (inline) fields =====
    // morton_index, level, grid_coords, bounds, center, size,
    // parent_morton, children_morton, num_children, num_points,
    // num_near_field_interactions, num_far_field_interactions,
    // on_boundary, and the vector/map headers themselves.
    total += sizeof(BoxData<CoordType, DataType>);

    // ===== Helper lambdas =====

    // Heap memory used by a MatrixStorage (just the vector's heap allocation)
    auto matrix_heap = [](const MatrixStorage<DataType>& m) -> size_t {
        return m.data.capacity() * sizeof(DataType);
    };

    // Heap memory used by a ModifiedBlock
    auto modified_block_heap = [&matrix_heap](const ModifiedBlock<DataType>& mb) -> size_t {
        // sizeof(ModifiedBlock) is already counted when we account for the
        // vector<ModifiedBlock> capacity, so only count heap inside each block.
        return matrix_heap(mb.A_SN) + matrix_heap(mb.A_NS);
    };

    // Heap memory of a std::vector<T> (element storage only; header is inline)
    auto vec_heap = [](auto const& v) -> size_t {
        return v.capacity() * sizeof(typename std::remove_reference_t<decltype(v)>::value_type);
    };

    // Rough heap estimate for an unordered_map.
    // Each bucket is typically one pointer, and each element lives in a
    // separately allocated node containing the key-value pair plus a pointer.
    auto map_heap = [](auto const& m) -> size_t {
        using map_type = std::remove_reference_t<decltype(m)>;
        using value_type = typename map_type::value_type;
        // bucket array
        size_t bytes = m.bucket_count() * sizeof(void*);
        // per-element node overhead (pair + next pointer + possible padding)
        bytes += m.size() * (sizeof(value_type) + sizeof(void*) + sizeof(size_t));
        return bytes;
    };

    // ===== Vectors of indices / coordinates =====
    total += vec_heap(box.point_indices);
    total += vec_heap(box.point_coords);
    total += vec_heap(box.redundant_indices);
    total += vec_heap(box.skeleton_indices);
    total += vec_heap(box.one_hop);
    total += vec_heap(box.use_full_set);
    total += vec_heap(box.two_hop);
    total += vec_heap(box.X_RR_pivots);

    // ===== Matrix data =====
    total += matrix_heap(box.interpolation_matrix);
    total += matrix_heap(box.X_RR);
    total += matrix_heap(box.X_RS);
    total += matrix_heap(box.X_SR);
    total += matrix_heap(box.schur_complement);
    total += matrix_heap(box.X_RN);
    total += matrix_heap(box.X_NR);

    // ===== Interaction maps =====
    total += map_heap(box.near_field_interaction_map);
    total += map_heap(box.far_field_interaction_map);
    total += map_heap(box.near_field_interaction_map_nonsymmetry);
    total += map_heap(box.far_field_interaction_map_nonsymmetry);

    // ===== Modified interaction blocks =====
    // The vector stores ModifiedBlock objects inline; their sizeof is already
    // covered by vec_heap. We additionally need the heap inside each block.
    total += vec_heap(box.near_field_modified_interactions);
    for (const auto& mb : box.near_field_modified_interactions) {
        total += modified_block_heap(mb);
    }

    total += vec_heap(box.far_field_modified_interactions);
    for (const auto& mb : box.far_field_modified_interactions) {
        total += modified_block_heap(mb);
    }

    return total;
}

/**
 * @brief Process reduction strategy (FMM3D.pdf Section 5.4)
 * 
 * Determines how processes are reduced at coarser tree levels:
 * - UNIFORM: Reduce by fixed factor at each level (default)
 * - ADAPTIVE: Reduce based on work vs communication cost
 * - WORK_PRESERVING: Maintain constant work per process
 * - AGGRESSIVE: Reduce to single process as soon as possible
 */
enum class ReductionPattern {
    UNIFORM,           ///< Fixed reduction factor (4 for 2D, 8 for 3D)
    ADAPTIVE,          ///< Dynamic based on profiling (not yet implemented)
    WORK_PRESERVING,   ///< Maintain boxes/process ratio (not yet implemented)
    AGGRESSIVE         ///< Reduce to 1 process ASAP (not yet implemented)
};

/**
 * @brief Main parallel tree structure (FMM3D.pdf Section 4.4.1)
 * 
 * Represents the complete FMM tree distributed across MPI processes.
 * Boxes are distributed in contiguous Morton-ordered ranges for load balance.
 * 
 * Domain Handling:
 * - Supports shifted domains (e.g., [0,1]×[1,2]) via bounds offset
 * - Uses SQUARE (2D) or CUBIC (3D) boxes for uniform geometry
 * - Box size = max(dimension ranges) / grid_size
 * - Rectangular domains waste some boxes but simplify proxy geometry
 */
template<typename CoordType, typename DataType = CoordType>
struct ParallelTree {
    int dimension;            ///< 2 or 3
    int num_levels;           ///< Number of levels (L+1), leaf is at level L
    int64_t num_points;       ///< Total number of points globally (int64_t per Section 3.4.1)
    
    CoordType global_bounds[6];  ///< [xmin, xmax, ymin, ymax, zmin, zmax]
    
    std::vector<TreeLevel<CoordType, DataType>> levels;  ///< Array of levels [0..num_levels]
    
    // MPI information
    int mpi_rank;             ///< Current process rank
    int mpi_size;             ///< Total number of processes
    MPI_Comm comm;            ///< MPI communicator
    
    // Process reduction configuration (Section 5.3-5.4)
    ReductionPattern reduction_pattern;  ///< Reduction strategy
    int64_t reduction_threshold;         ///< Boxes/process threshold for reduction
    
    ParallelTree();
    ~ParallelTree() = default;  // Vector handles cleanup automatically
    
    ParallelTree(const ParallelTree&) = delete;
    ParallelTree& operator=(const ParallelTree&) = delete;

    // Add move constructor and move assignment
    ParallelTree(ParallelTree&&) = default;
    ParallelTree& operator=(ParallelTree&&) = default;


};



} // namespace fmm

#endif // TREE_HPP
