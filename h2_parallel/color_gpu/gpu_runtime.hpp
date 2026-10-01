#pragma once
// GPU runtime of the H2 Color backend: the CUDA device and MAGMA queue of this
// MPI rank, grow-only device and pinned-host buffers, and the batched kernels
// the backend calls.  Compiled only with H2_HAVE_GPU (CMake enable_h2_gpu).

#ifdef H2_HAVE_GPU

#include <cuda_runtime.h>
#include "magma_v2.h"
#include "device_kernels.hpp"
#include <mpi.h>

#include "../bpack_env.hpp"

#include <algorithm>
#include <complex>
#include <cstring>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace fmm {
namespace gpu {

inline void check_cuda(cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA error in ") + what + ": " +
                                 cudaGetErrorString(status));
    }
}

// One device per MPI rank: ranks on a node are spread round-robin over its
// GPUs.  CUDA keeps the current device per host thread, so every entry point
// of the backend calls activate() on the thread that issues GPU work.
class Context {
public:
    static Context& instance() {
        static Context context;
        return context;
    }

    void activate() const { check_cuda(cudaSetDevice(device_), "cudaSetDevice"); }
    int device() const { return device_; }
    magma_queue_t queue() const { return queue_; }
    cudaStream_t stream() const { return magma_queue_get_cuda_stream(queue_); }

    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

private:
    Context() {
        int count = 0;
        check_cuda(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
        if (count <= 0) {
            throw std::runtime_error("H2 GPU backend: no CUDA device visible to this rank");
        }
        MPI_Comm node_comm;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
        int local_rank = 0;
        MPI_Comm_rank(node_comm, &local_rank);
        MPI_Comm_free(&node_comm);
        device_ = local_rank % count;
        activate();
        magma_init();
        magma_queue_create(device_, &queue_);
    }
    // The queue and MAGMA state live until process exit: tearing them down
    // from a static destructor would race the CUDA runtime's own shutdown.
    ~Context() = default;

    int device_ = 0;
    magma_queue_t queue_ = nullptr;
};

// Grow-only device allocation; growing discards the contents.
class DeviceBuffer {
public:
    DeviceBuffer() = default;
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    ~DeviceBuffer() { release(); }

    void* reserve(size_t bytes) {
        if (bytes > capacity_) {
            release();
            const size_t grown = std::max(bytes, capacity_ + capacity_ / 4);
            check_cuda(cudaMalloc(&ptr_, grown), "cudaMalloc");
            capacity_ = grown;
        }
        return ptr_;
    }
    void release() {
        if (ptr_ != nullptr) {
            cudaFree(ptr_);  // may fail during process teardown; nothing to recover
            ptr_ = nullptr;
            capacity_ = 0;
        }
    }
    size_t capacity() const { return capacity_; }

private:
    void* ptr_ = nullptr;
    size_t capacity_ = 0;
};

// Grow-only page-locked host allocation for asynchronous transfers.  Pinning
// costs about a second per GB, so buffers grow geometrically and the
// backend's buffers live in pinned_pool() for the whole run.
class PinnedBuffer {
public:
    PinnedBuffer() = default;
    PinnedBuffer(const PinnedBuffer&) = delete;
    PinnedBuffer& operator=(const PinnedBuffer&) = delete;
    ~PinnedBuffer() { release(); }

    void* reserve(size_t bytes) {
        if (bytes > capacity_) {
            release();
            constexpr size_t granule = size_t{64} << 20;
            const size_t grown = (std::max(bytes, 2 * capacity_) + granule - 1) / granule * granule;
            check_cuda(cudaMallocHost(&ptr_, grown), "cudaMallocHost");
            capacity_ = grown;
        }
        return ptr_;
    }
    void release() {
        if (ptr_ != nullptr) {
            cudaFreeHost(ptr_);
            ptr_ = nullptr;
            capacity_ = 0;
        }
    }

private:
    void* ptr_ = nullptr;
    size_t capacity_ = 0;
};

inline size_t align_up(size_t bytes, size_t alignment = 256) {
    return (bytes + alignment - 1) / alignment * alignment;
}

// Scalar type MAGMA uses for DataType (std::complex<double> and
// magmaDoubleComplex share their layout).
template<typename DataType> struct MagmaScalar;
template<> struct MagmaScalar<double> { using type = double; };
template<> struct MagmaScalar<std::complex<double>> { using type = magmaDoubleComplex; };
template<> struct MagmaScalar<dcomplex> { using type = magmaDoubleComplex; };

template<typename DataType>
typename MagmaScalar<DataType>::type to_magma(DataType value) {
    if constexpr (std::is_same_v<DataType, double>) {
        return value;
    } else if constexpr (std::is_same_v<DataType, dcomplex>) {
        return MAGMA_Z_MAKE(value.re, value.im);
    } else {
        return MAGMA_Z_MAKE(value.real(), value.imag());
    }
}

// H2_use_gpu=2: the backend's double GEMM batches (VBatch::gemm) run on the
// FP64 tensor cores (launch_dgemm_vbatched_tc) instead of MAGMA.
inline bool& tensor_core_gemm() {
    static bool enabled = false;
    return enabled;
}

// One variable-size batch of C_i = alpha * op(A_i) * op(B_i) + beta * C_i
// (transposes are plain, also for complex data).  All arrays live on the
// device; the size arrays have batch + 1 entries (MAGMA uses the last one).
// The checked MAGMA entry points find the maximum sizes with a device
// reduction and a host synchronization per call; when the caller knows them
// (max_* > 0) the core routines run instead and nothing blocks.
template<typename DataType>
void gemm_vbatched(magma_trans_t trans_a, magma_trans_t trans_b,
                   magma_int_t* d_m, magma_int_t* d_n, magma_int_t* d_k,
                   DataType alpha,
                   const DataType* const* d_a, magma_int_t* d_lda,
                   const DataType* const* d_b, magma_int_t* d_ldb,
                   DataType beta,
                   DataType** d_c, magma_int_t* d_ldc,
                   magma_int_t batch, magma_queue_t queue,
                   magma_int_t max_m = 0, magma_int_t max_n = 0, magma_int_t max_k = 0) {
    if (batch == 0) {
        return;
    }
    using M = typename MagmaScalar<DataType>::type;
    const bool known = max_m > 0 && max_n > 0 && max_k > 0;
    if constexpr (std::is_same_v<DataType, double>) {
        if (known) {
            magmablas_dgemm_vbatched_core(trans_a, trans_b, max_m, max_n, max_k, d_m, d_n, d_k, alpha,
                                          d_a, 0, 0, d_lda, d_b, 0, 0, d_ldb, beta, d_c, 0, 0, d_ldc,
                                          batch, queue);
        } else {
            magmablas_dgemm_vbatched(trans_a, trans_b, d_m, d_n, d_k, alpha,
                                     d_a, d_lda, d_b, d_ldb, beta, d_c, d_ldc, batch, queue);
        }
    } else {
        if (known) {
            magmablas_zgemm_vbatched_core(trans_a, trans_b, max_m, max_n, max_k, d_m, d_n, d_k, to_magma(alpha),
                                          reinterpret_cast<const M* const*>(d_a), 0, 0, d_lda,
                                          reinterpret_cast<const M* const*>(d_b), 0, 0, d_ldb,
                                          to_magma(beta), reinterpret_cast<M**>(d_c), 0, 0, d_ldc,
                                          batch, queue);
        } else {
            magmablas_zgemm_vbatched(trans_a, trans_b, d_m, d_n, d_k, to_magma(alpha),
                                     reinterpret_cast<const M* const*>(d_a), d_lda,
                                     reinterpret_cast<const M* const*>(d_b), d_ldb,
                                     to_magma(beta), reinterpret_cast<M**>(d_c), d_ldc,
                                     batch, queue);
        }
    }
}

template<typename DataType>
void gemm_nt_vbatched(magma_int_t* d_m, magma_int_t* d_n, magma_int_t* d_k,
                      DataType alpha,
                      const DataType* const* d_a, magma_int_t* d_lda,
                      const DataType* const* d_b, magma_int_t* d_ldb,
                      DataType beta,
                      DataType** d_c, magma_int_t* d_ldc,
                      magma_int_t batch, magma_queue_t queue) {
    gemm_vbatched<DataType>(MagmaNoTrans, MagmaTrans, d_m, d_n, d_k, alpha, d_a, d_lda,
                            d_b, d_ldb, beta, d_c, d_ldc, batch, queue);
}

// B_i = alpha * B_i * op(A_i)^{-1} or alpha * op(A_i)^{-1} * B_i, batched,
// with the maximum m and n known on the host.
template<typename DataType>
void trsm_vbatched(magma_side_t side, magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
                   magma_int_t max_m, magma_int_t max_n,
                   magma_int_t* d_m, magma_int_t* d_n, DataType alpha,
                   DataType** d_a, magma_int_t* d_lda, DataType** d_b, magma_int_t* d_ldb,
                   magma_int_t batch, magma_queue_t queue) {
    if (batch == 0 || max_m <= 0 || max_n <= 0) {
        return;
    }
    using M = typename MagmaScalar<DataType>::type;
    if constexpr (std::is_same_v<DataType, double>) {
        magmablas_dtrsm_vbatched_max_nocheck(side, uplo, trans, diag, max_m, max_n, d_m, d_n, alpha,
                                             d_a, d_lda, d_b, d_ldb, batch, queue);
    } else {
        magmablas_ztrsm_vbatched_max_nocheck(side, uplo, trans, diag, max_m, max_n, d_m, d_n, to_magma(alpha),
                                             reinterpret_cast<M**>(d_a), d_lda,
                                             reinterpret_cast<M**>(d_b), d_ldb, batch, queue);
    }
}

// LU with partial pivoting of square matrices (1-based pivots), batched, with
// the maximum order known on the host; `work` grows as needed.
template<typename DataType>
void getrf_vbatched(magma_int_t max_n, magma_int_t* d_m, magma_int_t* d_n, DataType** d_a, magma_int_t* d_lda,
                    magma_int_t** d_ipiv, magma_int_t* d_info, magma_int_t batch,
                    DeviceBuffer& work, magma_queue_t queue) {
    if (batch == 0 || max_n <= 0) {
        return;
    }
    using M = typename MagmaScalar<DataType>::type;
    const magma_int_t max_minmn = max_n, max_mxn = max_n * max_n;
    magma_int_t lwork[1] = {-1};
    magma_int_t status = 0;
    if constexpr (std::is_same_v<DataType, double>) {
        magma_dgetrf_vbatched_max_nocheck_work(nullptr, nullptr, max_n, max_n, max_minmn, max_mxn, nullptr, nullptr,
                                               nullptr, nullptr, nullptr, lwork, batch, queue);
        void* w = work.reserve(static_cast<size_t>(lwork[0]));
        status = magma_dgetrf_vbatched_max_nocheck_work(d_m, d_n, max_n, max_n, max_minmn, max_mxn, d_a, d_lda,
                                                        d_ipiv, d_info, w, lwork, batch, queue);
    } else {
        magma_zgetrf_vbatched_max_nocheck_work(nullptr, nullptr, max_n, max_n, max_minmn, max_mxn, nullptr, nullptr,
                                               nullptr, nullptr, nullptr, lwork, batch, queue);
        void* w = work.reserve(static_cast<size_t>(lwork[0]));
        status = magma_zgetrf_vbatched_max_nocheck_work(d_m, d_n, max_n, max_n, max_minmn, max_mxn,
                                                        reinterpret_cast<M**>(d_a), d_lda, d_ipiv, d_info,
                                                        w, lwork, batch, queue);
    }
    if (status != 0) {
        throw std::runtime_error("magma getrf_vbatched failed with status " + std::to_string(status));
    }
}

// Size class of a batched LU (and of the solves with its factors): the next
// power of two, at least 32.  Launched per class with max_n = the class, a
// matrix's factors do not depend on its batch companions (MAGMA's blocking
// follows max_n; tests/batch_determinism.cpp).  (Classes under 32 would load
// MAGMA kernels of their own at first use, a one-time ~0.3 s.)
inline int lu_size_class(int n) {
    int c = 32;
    while (c < n) c *= 2;
    return c;
}

// Device time between consecutive marks on a stream, read after the stream
// has been synchronized (no extra synchronization).
class StreamMarks {
public:
    StreamMarks() = default;
    StreamMarks(const StreamMarks&) = delete;
    StreamMarks& operator=(const StreamMarks&) = delete;
    ~StreamMarks() {
        for (cudaEvent_t e : events_) cudaEventDestroy(e);
    }
    void mark(cudaStream_t stream) {
        cudaEvent_t e = nullptr;
        check_cuda(cudaEventCreate(&e), "cudaEventCreate");
        check_cuda(cudaEventRecord(e, stream), "cudaEventRecord");
        events_.push_back(e);
    }
    // wait for the last mark
    void synchronize() const {
        if (!events_.empty()) check_cuda(cudaEventSynchronize(events_.back()), "cudaEventSynchronize");
    }
    // seconds from mark i to mark i + 1
    double seconds(size_t i) const {
        float ms = 0.0f;
        check_cuda(cudaEventElapsedTime(&ms, events_[i], events_[i + 1]), "cudaEventElapsedTime");
        return 1e-3 * static_cast<double>(ms);
    }

private:
    std::vector<cudaEvent_t> events_;
};

// Host image of one wave's device metadata (index lists, launch items, MAGMA
// size and pointer arrays), built directly in pinned memory and uploaded with
// a single copy.  Entries are addressed by byte offsets, so the image can
// grow (geometrically, keeping its contents) while it is built.  Clearing or
// growing the image first waits for its last upload, so the host may build
// the next image while the device still runs the previous work.
class MetaBuilder {
public:
    MetaBuilder() = default;
    MetaBuilder(const MetaBuilder&) = delete;
    MetaBuilder& operator=(const MetaBuilder&) = delete;

    void clear() {
        wait_upload();
        size_ = 0;
    }
    size_t size() const { return size_; }

    // Reserve `bytes` at a 16-byte aligned offset, to be filled via host().
    size_t alloc(size_t bytes) {
        const size_t offset = align_up(size_, 16);
        ensure(offset + bytes);
        size_ = offset + bytes;
        return offset;
    }
    char* host(size_t offset) { return base_ + offset; }

    template<typename T>
    size_t append(const T* data, size_t count) {
        const size_t offset = alloc(count * sizeof(T));
        if (count > 0) std::memcpy(base_ + offset, data, count * sizeof(T));
        return offset;
    }
    template<typename T>
    size_t append(const std::vector<T>& values) {
        return append(values.data(), values.size());
    }

    // Copy the image to the device and return the device address of offset
    // 0.  The caller keeps the stream ordered.
    char* upload(DeviceBuffer& device, cudaStream_t stream) {
        char* d = static_cast<char*>(device.reserve(std::max<size_t>(size_, 16)));
        if (size_ > 0) {
            check_cuda(cudaMemcpyAsync(d, base_, size_, cudaMemcpyHostToDevice, stream), "metadata upload");
            if (uploaded_ == nullptr) {
                check_cuda(cudaEventCreateWithFlags(&uploaded_, cudaEventDisableTiming), "cudaEventCreate");
            }
            check_cuda(cudaEventRecord(uploaded_, stream), "cudaEventRecord");
            pending_ = true;
        }
        return d;
    }

private:
    void wait_upload() {
        if (!pending_) return;
        check_cuda(cudaEventSynchronize(uploaded_), "metadata upload");
        pending_ = false;
    }

    void ensure(size_t bytes) {
        if (bytes <= capacity_) return;
        wait_upload();
        constexpr size_t granule = size_t{64} << 20;
        const size_t grown = (std::max(bytes, 2 * capacity_) + granule - 1) / granule * granule;
        void* fresh = nullptr;
        check_cuda(cudaMallocHost(&fresh, grown), "cudaMallocHost");
        if (size_ > 0) std::memcpy(fresh, base_, size_);
        if (base_ != nullptr) cudaFreeHost(base_);
        base_ = static_cast<char*>(fresh);
        capacity_ = grown;
    }

    char* base_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;
    cudaEvent_t uploaded_ = nullptr;  // after the last upload (never destroyed: see PinnedPool)
    bool pending_ = false;
};

// MPI messages go between the ranks' device memories when MPI is GPU-aware
// (BPACK_GPU_AWARE_MPI: 1 yes, 0 no; unset: Cray MPICH's
// MPICH_GPU_SUPPORT_ENABLED=1); otherwise they are staged through the host.
// The H2 and HODLR GPU backends both ask here.
inline bool device_exchange_enabled() {
    static const bool enabled = [] {
        const char* v = std::getenv("BPACK_GPU_AWARE_MPI");
        if (v == nullptr) v = std::getenv("MPICH_GPU_SUPPORT_ENABLED");
        return v != nullptr && std::atoi(v) != 0;
    }();
    return enabled;
}

// BPACK_CHECK=replica: on a replicated CA level, every ghost copy of a box is
// compared (by a hash of its factors) with its owner's.
inline bool ca_replica_check_enabled() { return env::check("replica"); }

// Pinned staging buffers of the backend, kept for the whole run (never
// freed: releasing from a static destructor would race the CUDA runtime's
// own shutdown).
struct PinnedPool {
    MetaBuilder meta;          // per-wave metadata image
    MetaBuilder owner_meta;    // owner pass of a wave, built while its first part runs
    PinnedBuffer result;       // sketches, factors, level-end blocks
    PinnedBuffer ring[4];      // fixed-size slots of chunked downloads (two per copier)
    PinnedBuffer staging[2];   // alternating uploads (restored fill sources)
    static constexpr size_t kRingSlot = size_t{64} << 20;  // pinning costs ~1 s/GB
};

inline PinnedPool& pinned_pool() {
    static PinnedPool* pool = new PinnedPool;
    return *pool;
}

// A variable-size MAGMA batch described on the host; its six size arrays and
// three pointer arrays are placed in the metadata image by stage().
template<typename DataType>
struct VBatch {
    struct Entry {
        const DataType* a;
        const DataType* b;
        DataType* c;
        magma_int_t m, n, k, lda, ldb, ldc;
    };
    std::vector<Entry> entries;
    size_t sizes_offset = 0;     // m, n, k, lda, ldb, ldc: (count + 1) each
    size_t pointers_offset = 0;  // a, b, c: count each
    magma_int_t max_m = 0, max_n = 0, max_k = 0;
    // tensor-core schedule (H2_use_gpu=2, set by stage): the tile blocks of
    // every entry, the 48 x 48 ones first, as {entry, row << 16 | column}
    bool tc = false;
    size_t tc_blocks_offset = 0;
    int tc_blocks48 = 0, tc_blocks64 = 0;

    size_t count() const { return entries.size(); }


    void stage(MetaBuilder& meta) {
        const size_t n = entries.size();
        std::vector<magma_int_t> sizes(6 * (n + 1), 0);
        std::vector<const void*> pointers(3 * std::max<size_t>(n, 1), nullptr);
        max_m = max_n = max_k = 0;
        for (size_t i = 0; i < n; ++i) {
            const Entry& e = entries[i];
            max_m = std::max(max_m, e.m);
            max_n = std::max(max_n, e.n);
            max_k = std::max(max_k, e.k);
            sizes[0 * (n + 1) + i] = e.m;
            sizes[1 * (n + 1) + i] = e.n;
            sizes[2 * (n + 1) + i] = e.k;
            sizes[3 * (n + 1) + i] = e.lda;
            sizes[4 * (n + 1) + i] = e.ldb;
            sizes[5 * (n + 1) + i] = e.ldc;
            pointers[i] = e.a;
            pointers[n + i] = e.b;
            pointers[2 * n + i] = e.c;
        }
        sizes_offset = meta.append(sizes);
        pointers_offset = meta.append(pointers);
        tc = false;
        if constexpr (std::is_same_v<DataType, double>) {
            if (tensor_core_gemm() && n > 0) stage_tensor_core_blocks(meta);
        }
    }

    // Tile of an entry on the tensor cores: 64 x 64, unless 48 x 48 tiles
    // waste much less of it (measured on an A100: 40 x 40 x 176 and
    // 88 x 88 x 176 run faster with 48, 128 x 128 and above with 64).
    static int tc_tile_of(const Entry& e) {
        const double mn = static_cast<double>(e.m) * e.n;
        const double used48 = mn / (static_cast<double>((e.m + 47) / 48 * ((e.n + 47) / 48)) * 48.0 * 48.0);
        const double used64 = mn / (static_cast<double>((e.m + 63) / 64 * ((e.n + 63) / 64)) * 64.0 * 64.0);
        return used48 > 1.4 * used64 ? 48 : 64;
    }

    void stage_tensor_core_blocks(MetaBuilder& meta) {
        static_assert(sizeof(magma_int_t) == sizeof(int), "the tensor-core GEMM takes 32-bit sizes");
        std::vector<int2> blocks48, blocks64;
        for (size_t i = 0; i < entries.size(); ++i) {
            const Entry& e = entries[i];
            if (e.m <= 0 || e.n <= 0) continue;
            const int tile = tc_tile_of(e);
            std::vector<int2>& out = tile == 48 ? blocks48 : blocks64;
            const int rows = (e.m + tile - 1) / tile, cols = (e.n + tile - 1) / tile;
            if (rows > 0xffff || cols > 0xffff) {
                throw std::runtime_error("VBatch: matrix too large for the tensor-core GEMM");
            }
            for (int c = 0; c < cols; ++c)
                for (int r = 0; r < rows; ++r) out.push_back(make_int2(static_cast<int>(i), (r << 16) | c));
        }
        if (blocks48.size() + blocks64.size() > 0x7fffffffu) {
            throw std::runtime_error("VBatch: too many tensor-core blocks");
        }
        tc_blocks48 = static_cast<int>(blocks48.size());
        tc_blocks64 = static_cast<int>(blocks64.size());
        blocks48.insert(blocks48.end(), blocks64.begin(), blocks64.end());
        tc_blocks_offset = meta.append(blocks48);
        tc = true;
    }

    magma_int_t* size_array(char* meta_device, int which) const {
        return reinterpret_cast<magma_int_t*>(meta_device + sizes_offset) + which * (entries.size() + 1);
    }
    template<typename P>
    P* pointer_array(char* meta_device, int which) const {
        return reinterpret_cast<P*>(meta_device + pointers_offset) + which * entries.size();
    }

    // Launch as a GEMM: C = alpha * op(A) * op(B) + beta * C.  The entries
    // write distinct blocks, so the tensor-core path may run them in any order.
    void gemm(char* meta_device, magma_trans_t ta, magma_trans_t tb, DataType alpha, DataType beta,
              magma_queue_t queue) const {
        if constexpr (std::is_same_v<DataType, double>) {
            if (tc) {
                const int2* blocks = reinterpret_cast<const int2*>(meta_device + tc_blocks_offset);
                const cudaStream_t stream = magma_queue_get_cuda_stream(queue);
                for (int part = 0; part < 2; ++part) {
                    launch_dgemm_vbatched_tc(ta != MagmaNoTrans, tb != MagmaNoTrans, size_array(meta_device, 0),
                                             size_array(meta_device, 1), size_array(meta_device, 2), alpha,
                                             pointer_array<const double* const>(meta_device, 0),
                                             size_array(meta_device, 3),
                                             pointer_array<const double* const>(meta_device, 1),
                                             size_array(meta_device, 4), beta,
                                             pointer_array<double* const>(meta_device, 2), size_array(meta_device, 5),
                                             part == 0 ? blocks : blocks + tc_blocks48,
                                             part == 0 ? tc_blocks48 : tc_blocks64, part == 0 ? 48 : 64, stream);
                }
                return;
            }
        }
        gemm_vbatched<DataType>(ta, tb, size_array(meta_device, 0), size_array(meta_device, 1),
                                size_array(meta_device, 2),
                                alpha, pointer_array<const DataType* const>(meta_device, 0),
                                size_array(meta_device, 3),
                                pointer_array<const DataType* const>(meta_device, 1),
                                size_array(meta_device, 4), beta,
                                pointer_array<DataType*>(meta_device, 2), size_array(meta_device, 5),
                                static_cast<magma_int_t>(entries.size()), queue, max_m, max_n, max_k);
    }
};

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
