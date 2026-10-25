// The two kernels of the butterfly's entry extraction (bf_extract.hpp): the
// gather of factor rows into the first partial products, and the scatter of
// the entries from leaf order to the caller's.  Everything between them is
// batched GEMMs.

#include "bf_extract.hpp"

#include <stdexcept>
#include <string>

namespace bf {
namespace gpu {

namespace {

constexpr int kThreads = 256;

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

// One thread block per gather: rows idx of the factor block, transposed
// (V: r x count) or as they are (U: count x r).
__global__ void __launch_bounds__(kThreads)
extract_gather_kernel(const ExtractGather* items, const int* idx, const DataType* values, DataType* dst) {
    const ExtractGather g = items[blockIdx.x];
    const DataType* S = values + g.src;
    const int* I = idx + g.idx;
    DataType* D = dst + g.dst;
    const int64_t total = static_cast<int64_t>(g.count) * g.r;
    for (int64_t t = threadIdx.x; t < total; t += kThreads) {
        int c, q;
        if (g.transposed) {
            q = static_cast<int>(t % g.r);
            c = static_cast<int>(t / g.r);
        } else {
            c = static_cast<int>(t % g.count);
            q = static_cast<int>(t / g.count);
        }
        D[t] = S[I[c] + static_cast<int64_t>(q) * g.ld_src];
    }
}

// One thread block per intersection: entry (p, q) of the caller's block is
// entry (pos_r[p], pos_c[q]) of the leaf-ordered one.
__global__ void __launch_bounds__(kThreads)
extract_scatter_kernel(const ExtractScatter* items, const int* idx, const DataType* src, DataType* dst) {
    const ExtractScatter s = items[blockIdx.x];
    const int* pr = idx + s.pos_r;
    const int* pc = idx + s.pos_c;
    const int64_t total = static_cast<int64_t>(s.m) * s.n;
    for (int64_t t = threadIdx.x; t < total; t += kThreads) {
        const int p = static_cast<int>(t % s.m), q = static_cast<int>(t / s.m);
        dst[s.dst + t] = src[s.src + pr[p] + static_cast<int64_t>(pc[q]) * s.m];
    }
}

}  // namespace

void launch_extract_gather(const ExtractGather* items, int count, const int* idx, const DataType* values,
                           DataType* dst, cudaStream_t stream) {
    if (count <= 0) return;
    extract_gather_kernel<<<count, kThreads, 0, stream>>>(items, idx, values, dst);
    check_launch("launch_extract_gather");
}

void launch_extract_scatter(const ExtractScatter* items, int count, const int* idx, const DataType* src,
                            DataType* dst, cudaStream_t stream) {
    if (count <= 0) return;
    extract_scatter_kernel<<<count, kThreads, 0, stream>>>(items, idx, src, dst);
    check_launch("launch_extract_scatter");
}

}  // namespace gpu
}  // namespace bf
