// The butterfly-specific glue around the H2 backend's batched interpolative
// decomposition.  See bf_id.hpp: the factorization is
// fmm::gpu::launch_qrcp; what is here is the load (with the 'C' half's
// transpose) and the assembly of X from (jpvt, rank, T).

#include "bf_id.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace bf {
namespace gpu {

using fmm::gpu::DeviceHeap;
using fmm::gpu::check_cuda;

namespace {

constexpr int kThreads = 256;

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

// One thread block per panel: copy the entries into the matrix the
// factorization destroys, transposing the panels the 'C' half works on.
__global__ void __launch_bounds__(kThreads)
id_load_kernel(const DataType* entries, const IdPanel* panels, DataType* mat) {
    const IdPanel p = panels[blockIdx.x];
    const int m = p.m(), n = p.n();
    const DataType* E = entries + p.entry_offset;
    DataType* A = mat + p.mat_offset;
    const int64_t lde = p.rows;
    for (int64_t t = threadIdx.x; t < static_cast<int64_t>(m) * n; t += kThreads) {
        const int i = static_cast<int>(t % m), j = static_cast<int>(t / m);
        A[t] = p.transposed ? E[j + static_cast<int64_t>(i) * lde] : E[i + static_cast<int64_t>(j) * lde];
    }
}

// X, as the butterfly stores it: rows 0..rank-1 of an m x n buffer, with
// X(:, jpvt[k]) = e_k on the skeleton columns and X(:, jpvt[rank+j]) the j-th
// column of T, which launch_qrcp left in rows 0..rank-1 of columns
// rank..n-1 of the factorized matrix.
__global__ void __launch_bounds__(kThreads)
id_scatter_kernel(const IdPanel* panels, const DataType* mat, const int* jpvt, const int* rank,
                  DataType* x) {
    const IdPanel p = panels[blockIdx.x];
    const int m = p.m(), n = p.n();
    const int k = rank[blockIdx.x];
    const DataType* A = mat + p.mat_offset;
    DataType* X = x + p.x_offset;
    const int* J = jpvt + p.jpvt_offset;

    // the skeleton columns: the identity
    for (int64_t t = threadIdx.x; t < static_cast<int64_t>(k) * k; t += kThreads) {
        const int row = static_cast<int>(t % k), col = static_cast<int>(t / k);
        X[row + static_cast<int64_t>(J[col]) * m] = (row == col) ? DataType(1) : DataType(0);
    }
    // the rest: T
    const int64_t rest = static_cast<int64_t>(k) * (n - k);
    for (int64_t t = threadIdx.x; t < rest; t += kThreads) {
        const int row = static_cast<int>(t % k), col = static_cast<int>(t / k);
        X[row + static_cast<int64_t>(J[k + col]) * m] = A[row + static_cast<int64_t>(k + col) * m];
    }
}

// launch_qrcp returns rank 0 for a panel it finds numerically empty; the CPU
// turns that into rank 1 with a zero X and the first pivot, which is what the
// level's bookkeeping expects.
__global__ void id_fixup_kernel(const IdPanel* panels, int count, int* jpvt, int* rank, DataType* x) {
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= count || rank[p] != 0) return;
    const IdPanel pan = panels[p];
    const int m = pan.m(), n = pan.n();
    DataType* X = x + pan.x_offset;
    for (int j = 0; j < n; ++j) X[static_cast<int64_t>(j) * m] = DataType(0);
    jpvt[pan.jpvt_offset] = 0;
    rank[p] = 1;
}

}  // namespace

bool plan_id_batch(IdBatch& batch) {
    int64_t mats = 0, jpvts = 0;
    int max_n = 0;
    for (IdPanel& p : batch.panels) {
        const int m = p.m(), n = p.n();
        p.mat_offset = mats;
        p.x_offset = mats;
        p.jpvt_offset = jpvts;
        mats += static_cast<int64_t>(m) * n;
        jpvts += n;
        max_n = std::max(max_n, n);
    }
    batch.mat_elements = mats;
    batch.jpvt_elements = jpvts;
    batch.max_n = max_n;
    if (batch.panels.empty() || mats == 0) return true;

    DeviceHeap& heap = DeviceHeap::instance();
    heap.ensure_initialized();
    const auto kind = fmm::gpu::ResidentKind::Adaptive;
    const int count = batch.count();
    const size_t work_bytes = fmm::gpu::qrcp_work_bytes<DataType>(count, max_n);
    std::vector<char*> taken;
    auto take = [&](size_t bytes) -> char* {
        char* q = bytes == 0 ? nullptr : heap.try_alloc(bytes, kind);
        if (bytes != 0) taken.push_back(q);
        return q;
    };
    char* panels = take(static_cast<size_t>(count) * sizeof(IdPanel));
    char* mat = take(static_cast<size_t>(mats) * sizeof(DataType));
    char* xbuf = take(static_cast<size_t>(mats) * sizeof(DataType));
    char* jpvt = take(static_cast<size_t>(jpvts) * sizeof(int));
    char* rank = take(static_cast<size_t>(count) * sizeof(int));
    char* norm = take(static_cast<size_t>(count) * sizeof(double));
    char* flag = take(static_cast<size_t>(count) * sizeof(int));
    char* items = take(static_cast<size_t>(count) * sizeof(fmm::gpu::QrcpItemT<DataType>));
    char* work = take(work_bytes);
    bool ok = true;
    for (size_t i = 0; i < taken.size(); ++i) {
        if (taken[i] == nullptr) ok = false;
    }
    if (!ok) {
        for (size_t i = taken.size(); i-- > 0;) {
            if (taken[i] != nullptr) heap.free(taken[i]);
        }
        return false;
    }
    batch.d_panels = reinterpret_cast<IdPanel*>(panels);
    batch.d_mat = reinterpret_cast<DataType*>(mat);
    batch.d_x = reinterpret_cast<DataType*>(xbuf);
    batch.d_jpvt = reinterpret_cast<int*>(jpvt);
    batch.d_rank = reinterpret_cast<int*>(rank);
    batch.d_norm = reinterpret_cast<double*>(norm);
    batch.d_flag = reinterpret_cast<int*>(flag);
    batch.d_items = reinterpret_cast<fmm::gpu::QrcpItemT<DataType>*>(items);
    batch.d_work = work;
    check_cuda(cudaMemcpy(batch.d_panels, batch.panels.data(), static_cast<size_t>(count) * sizeof(IdPanel),
                          cudaMemcpyHostToDevice), "bf id panels");

    // the factorization's descriptors, one per panel, over the same buffers
    std::vector<fmm::gpu::QrcpItemT<DataType>> host_items(static_cast<size_t>(count));
    for (int p = 0; p < count; ++p) {
        const IdPanel& pan = batch.panels[static_cast<size_t>(p)];
        fmm::gpu::QrcpItemT<DataType>& it = host_items[static_cast<size_t>(p)];
        it.a = batch.d_mat + pan.mat_offset;
        it.m = pan.m();
        it.n = pan.n();
        it.lda = pan.m();
        it.jpvt = batch.d_jpvt + pan.jpvt_offset;
        it.rank = batch.d_rank + p;
        it.norm = batch.d_norm + p;
        it.flag = batch.d_flag + p;
        it.tol = -1.0;   // the launch's tolerance
    }
    check_cuda(cudaMemcpy(batch.d_items, host_items.data(),
                          static_cast<size_t>(count) * sizeof(fmm::gpu::QrcpItemT<DataType>),
                          cudaMemcpyHostToDevice), "bf id items");
    return true;
}

void release_id_batch(IdBatch& batch) {
    DeviceHeap& heap = DeviceHeap::instance();
    for (void* q : {batch.d_work, static_cast<void*>(batch.d_items), static_cast<void*>(batch.d_flag),
                    static_cast<void*>(batch.d_norm), static_cast<void*>(batch.d_rank),
                    static_cast<void*>(batch.d_jpvt), static_cast<void*>(batch.d_x),
                    static_cast<void*>(batch.d_mat), static_cast<void*>(batch.d_panels)}) {
        if (q != nullptr) heap.free(reinterpret_cast<char*>(q));
    }
    batch.d_work = nullptr;
    batch.d_items = nullptr;
    batch.d_flag = nullptr;
    batch.d_norm = nullptr;
    batch.d_rank = nullptr;
    batch.d_jpvt = nullptr;
    batch.d_x = nullptr;
    batch.d_mat = nullptr;
    batch.d_panels = nullptr;
}

void launch_id_load(const DataType* entries, const IdPanel* panels, int count, DataType* mat,
                    cudaStream_t stream) {
    if (count <= 0) return;
    id_load_kernel<<<count, kThreads, 0, stream>>>(entries, panels, mat);
    check_launch("bf id load");
}

void launch_id_scatter(const IdPanel* panels, int count, const DataType* mat, const int* jpvt,
                       const int* rank, DataType* x, cudaStream_t stream) {
    if (count <= 0) return;
    id_scatter_kernel<<<count, kThreads, 0, stream>>>(panels, mat, jpvt, rank, x);
    check_launch("bf id scatter");
}

void run_id_batch(const DataType* entries, IdBatch& batch, double rtol, double atol, cudaStream_t stream) {
    const int count = batch.count();
    if (count == 0 || batch.mat_elements == 0) return;
    launch_id_load(entries, batch.d_panels, count, batch.d_mat, stream);
    // X is written only where the scatter puts something, so it starts zeroed:
    // a panel's X has n columns but only `rank` meaningful rows, and the
    // skeleton block is filled column by column
    check_cuda(cudaMemsetAsync(batch.d_x, 0, static_cast<size_t>(batch.mat_elements) * sizeof(DataType),
                               stream), "bf id zero X");
    fmm::gpu::launch_qrcp<DataType>(batch.d_items, count, batch.max_n, rtol, atol, batch.d_work, stream);
    launch_id_scatter(batch.d_panels, count, batch.d_mat, batch.d_jpvt, batch.d_rank, batch.d_x, stream);
    id_fixup_kernel<<<(count + kThreads - 1) / kThreads, kThreads, 0, stream>>>(
        batch.d_panels, count, batch.d_jpvt, batch.d_rank, batch.d_x);
    check_launch("bf id fixup");
}

}  // namespace gpu
}  // namespace bf
