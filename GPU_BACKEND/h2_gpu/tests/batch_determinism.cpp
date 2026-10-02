// Batch-independence of the GPU backend's batched kernels (M1.1 of
// h2_gpu/GPU_CA_PLAN.md).  In replicated CA, a boundary box is eliminated
// by its owner and again, as a ghost, by its neighbor ranks, each time in a
// different batch; all copies must come out bitwise identical.  Here every
// test box runs through a kernel alone, then in batches of other sizes and
// companions (including counts across the kernels' launch thresholds), and
// each result is compared bitwise with the box's alone result.
//
//   make h2_gpu_batch_determinism
//   srun -n 1 --gpus 1 ./h2_gpu_batch_determinism [boxes]
//
// Kernels: the device ID (launch_qrcp), MAGMA getrf/trsm/gemm vbatched (the
// X_RR LU, the solves, the updates) and the FP64 tensor-core GEMM, in double
// and complex.

#include "../device_kernels.hpp"
#include "../gpu_runtime.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace fmm::gpu;

namespace {

template<typename T>
struct Scalar;
template<>
struct Scalar<double> {
    using host = double;
    static const char* name() { return "double"; }
};
template<>
struct Scalar<dcomplex> {
    using host = std::complex<double>;
    static const char* name() { return "complex"; }
};

template<typename T>
T random_value(std::mt19937_64& g) {
    std::normal_distribution<double> d(0.0, 1.0);
    if constexpr (std::is_same_v<T, double>) {
        return d(g);
    } else {
        const double re = d(g), im = d(g);
        return dcomplex(re, im);
    }
}

template<typename T>
std::vector<T> random_matrix(int m, int n, std::mt19937_64& g) {
    std::vector<T> a(static_cast<size_t>(m) * n);
    for (auto& x : a) x = random_value<T>(g);
    return a;
}

// m x n with a decaying spectrum (a sketch-like ID target): sum over j of
// 10^(-8 j / r) u_j v_j^T for r = n / 2, plus a few duplicated columns so
// that pivot norms tie
template<typename T>
std::vector<T> id_target(int m, int n, std::mt19937_64& g) {
    const int r = std::max(1, n / 2);
    std::vector<T> a(static_cast<size_t>(m) * n, T(0.0));
    for (int j = 0; j < r; ++j) {
        const double s = std::pow(10.0, -8.0 * j / r);
        std::vector<T> u = random_matrix<T>(m, 1, g), v = random_matrix<T>(n, 1, g);
        for (int c = 0; c < n; ++c)
            for (int i = 0; i < m; ++i) a[static_cast<size_t>(c) * m + i] += s * u[i] * v[c];
    }
    for (int c = 1; c < n; c += 7) {  // ties
        std::copy(a.begin() + static_cast<size_t>(c - 1) * m, a.begin() + static_cast<size_t>(c) * m,
                  a.begin() + static_cast<size_t>(c) * m);
    }
    return a;
}

// Groupings of the boxes: alone (the reference), then others.  Each grouping
// is a list of batches (box indices).
struct Grouping {
    std::string name;
    std::vector<std::vector<int>> batches;
};

std::vector<Grouping> groupings(int boxes, std::mt19937_64& g) {
    std::vector<Grouping> out;
    auto chunks = [&](const std::string& name, const std::vector<int>& order, int size) {
        Grouping gr{name, {}};
        for (size_t i = 0; i < order.size(); i += static_cast<size_t>(size)) {
            gr.batches.emplace_back(order.begin() + i,
                                    order.begin() + std::min(order.size(), i + static_cast<size_t>(size)));
        }
        out.push_back(std::move(gr));
    };
    std::vector<int> order(static_cast<size_t>(boxes));
    std::iota(order.begin(), order.end(), 0);
    chunks("alone", order, 1);
    chunks("all", order, boxes);
    chunks("by 150", order, 150);
    chunks("by 60", order, 60);
    chunks("by 7", order, 7);
    std::shuffle(order.begin(), order.end(), g);
    chunks("shuffled by 50", order, 50);
    chunks("shuffled by 250", order, 250);
    return out;
}

struct Report {
    int mismatched = 0;
    double max_rel = 0.0;
    std::string first;
};

template<typename T>
void compare(const std::vector<T>& ref, const std::vector<T>& got, int box, Report& r, const char* what) {
    if (ref.size() == got.size() && (ref.empty() || std::memcmp(ref.data(), got.data(), ref.size() * sizeof(T)) == 0)) {
        return;
    }
    double diff = 0.0, scale = 0.0;
    for (size_t i = 0; i < std::min(ref.size(), got.size()); ++i) {
        if constexpr (std::is_same_v<T, int>) {
            diff = std::max(diff, ref[i] == got[i] ? 0.0 : 1.0);
            scale = 1.0;
        } else {
            using H = typename Scalar<T>::host;
            const H a = *reinterpret_cast<const H*>(&ref[i]), b = *reinterpret_cast<const H*>(&got[i]);
            diff = std::max(diff, std::abs(a - b));
            scale = std::max(scale, std::abs(a));
        }
    }
    if (r.mismatched == 0) r.first = std::string(what) + " of box " + std::to_string(box);
    ++r.mismatched;
    r.max_rel = std::max(r.max_rel, scale > 0.0 ? diff / scale : diff);
}

void print(const char* kernel, const char* type, const std::string& grouping, const Report& r, int boxes) {
    if (r.mismatched == 0) {
        std::printf("  %-26s %-8s %-16s identical (%d boxes)\n", kernel, type, grouping.c_str(), boxes);
    } else {
        std::printf("  %-26s %-8s %-16s DIFFERENT: %d of %d boxes, max rel diff %.2e (first: %s)\n", kernel, type,
                    grouping.c_str(), r.mismatched, boxes, r.max_rel, r.first.c_str());
    }
    std::fflush(stdout);
}

template<typename T>
T* upload(const std::vector<T>& h) {
    T* d = nullptr;
    check_cuda(cudaMalloc(&d, std::max<size_t>(h.size(), 1) * sizeof(T)), "cudaMalloc");
    if (!h.empty()) check_cuda(cudaMemcpy(d, h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice), "upload");
    return d;
}
template<typename T>
std::vector<T> download(const T* d, size_t n) {
    std::vector<T> h(n);
    if (n > 0) check_cuda(cudaMemcpy(h.data(), d, n * sizeof(T), cudaMemcpyDeviceToHost), "download");
    return h;
}

// ---- device ID: column-pivoted QR, rank, T (launch_qrcp)
template<typename T>
void test_qrcp(int boxes, std::mt19937_64& g, double tol, bool batch_independent) {
    cudaStream_t stream = Context::instance().stream();
    std::uniform_int_distribution<int> sizes(16, 600);
    std::vector<int> n(static_cast<size_t>(boxes));
    std::vector<std::vector<T>> input(static_cast<size_t>(boxes));
    for (int b = 0; b < boxes; ++b) {
        n[b] = sizes(g);
        input[b] = id_target<T>(n[b], n[b], g);  // d = n rows, as the sketch
    }
    struct Out {
        std::vector<T> a;
        std::vector<int> jpvt, rank, flag;
        std::vector<double> norm;
    };
    auto run = [&](const Grouping& gr) {
        std::vector<Out> out(static_cast<size_t>(boxes));
        for (const auto& batch : gr.batches) {
            const int count = static_cast<int>(batch.size());
            int max_n = 0;
            std::vector<T*> a(batch.size());
            std::vector<int*> jpvt(batch.size());
            int* d_rank = nullptr;
            int* d_flag = nullptr;
            double* d_norm = nullptr;
            check_cuda(cudaMalloc(&d_rank, count * sizeof(int)), "cudaMalloc");
            check_cuda(cudaMalloc(&d_flag, count * sizeof(int)), "cudaMalloc");
            check_cuda(cudaMalloc(&d_norm, count * sizeof(double)), "cudaMalloc");
            std::vector<QrcpItemT<T>> items(batch.size());
            for (int i = 0; i < count; ++i) {
                const int b = batch[static_cast<size_t>(i)];
                max_n = std::max(max_n, n[b]);
                a[i] = upload(input[b]);
                check_cuda(cudaMalloc(&jpvt[i], n[b] * sizeof(int)), "cudaMalloc");
                items[i] = QrcpItemT<T>{a[i], n[b], n[b], n[b], jpvt[i], d_rank + i, d_norm + i, d_flag + i};
            }
            QrcpItemT<T>* d_items = upload(items);
            void* work = nullptr;
            if (const size_t w = qrcp_work_bytes<T>(count, max_n)) check_cuda(cudaMalloc(&work, w), "cudaMalloc");
            launch_qrcp<T>(d_items, count, max_n, tol, work, stream, batch_independent);
            check_cuda(cudaStreamSynchronize(stream), "qrcp");
            const std::vector<int> rank = download(d_rank, count), flag = download(d_flag, count);
            const std::vector<double> norm = download(d_norm, count);
            for (int i = 0; i < count; ++i) {
                const int b = batch[static_cast<size_t>(i)];
                out[b].a = download(a[i], static_cast<size_t>(n[b]) * n[b]);
                out[b].jpvt = download(jpvt[i], n[b]);
                out[b].rank = {rank[i]};
                out[b].flag = {flag[i]};
                out[b].norm = {norm[i]};
                cudaFree(a[i]);
                cudaFree(jpvt[i]);
            }
            cudaFree(d_items);
            cudaFree(work);
            cudaFree(d_rank);
            cudaFree(d_flag);
            cudaFree(d_norm);
        }
        return out;
    };
    const auto all = groupings(boxes, g);
    const std::vector<Out> ref = run(all[0]);
    for (size_t k = 1; k < all.size(); ++k) {
        const std::vector<Out> got = run(all[k]);
        Report factors, pivots, norms;
        for (int b = 0; b < boxes; ++b) {
            compare(ref[b].a, got[b].a, b, factors, "R/T");
            compare(ref[b].jpvt, got[b].jpvt, b, pivots, "pivots");
            compare(ref[b].rank, got[b].rank, b, pivots, "rank");
            compare(ref[b].norm, got[b].norm, b, norms, "traced norm");
        }
        print(batch_independent ? "ID: pivots and rank (fixed)" : "ID: pivots and rank", Scalar<T>::name(), all[k].name, pivots, boxes);
        print(batch_independent ? "ID: R and T (fixed)" : "ID: R and T", Scalar<T>::name(), all[k].name, factors, boxes);
        print(batch_independent ? "ID: traced sketch norm (fixed)" : "ID: traced sketch norm", Scalar<T>::name(), all[k].name, norms, boxes);
    }
}

// ---- a MAGMA batch through VBatch (sizes and pointers in a metadata image)
// One pinned metadata image for every batch: a MetaBuilder pins its memory
// for good (the backend keeps its builders for the run).
inline MetaBuilder& shared_meta() {
    static MetaBuilder* meta = new MetaBuilder;
    return *meta;
}
template<typename T>
struct Batch {
    VBatch<T> v;
    DeviceBuffer device;
    char* md = nullptr;
    void stage(cudaStream_t stream) {
        MetaBuilder& meta = shared_meta();
        meta.clear();
        v.stage(meta);
        md = meta.upload(device, stream);
    }
};

// ---- X_RR LU (getrf_vbatched)
template<typename T>
void test_getrf(int boxes, std::mt19937_64& g, bool by_class) {
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    std::uniform_int_distribution<int> sizes(4, 400);
    std::vector<int> n(static_cast<size_t>(boxes));
    std::vector<std::vector<T>> input(static_cast<size_t>(boxes));
    for (int b = 0; b < boxes; ++b) {
        n[b] = sizes(g);
        input[b] = random_matrix<T>(n[b], n[b], g);
    }
    struct Out {
        std::vector<T> lu;
        std::vector<magma_int_t> ipiv;
    };
    DeviceBuffer work;
    auto run = [&](const Grouping& gr) {
        std::vector<Out> out(static_cast<size_t>(boxes));
        // by_class: each batch splits into its size classes, launched with
        // max_n = the class size (a box's launch then depends on itself only)
        std::vector<std::vector<int>> launches;
        for (const auto& batch : gr.batches) {
            if (!by_class) {
                launches.push_back(batch);
                continue;
            }
            std::vector<int> classes;
            for (int b : batch) classes.push_back(lu_size_class(n[b]));
            std::sort(classes.begin(), classes.end());
            classes.erase(std::unique(classes.begin(), classes.end()), classes.end());
            for (int c : classes) {
                std::vector<int> part;
                for (int b : batch) {
                    if (lu_size_class(n[b]) == c) part.push_back(b);
                }
                launches.push_back(part);
            }
        }
        for (const auto& batch : launches) {
            Batch<T> bt;
            std::vector<T*> a;
            std::vector<magma_int_t*> ipiv;
            for (int b : batch) {
                a.push_back(upload(input[b]));
                magma_int_t* p = nullptr;
                check_cuda(cudaMalloc(&p, n[b] * sizeof(magma_int_t)), "cudaMalloc");
                ipiv.push_back(p);
                bt.v.entries.push_back({a.back(), nullptr, nullptr, n[b], n[b], 1, n[b], 1, 1});
            }
            bt.stage(stream);
            const magma_int_t max_n = by_class ? lu_size_class(bt.v.max_n) : bt.v.max_n;
            magma_int_t** d_ipiv = upload(ipiv);
            magma_int_t* d_info = nullptr;
            check_cuda(cudaMalloc(&d_info, batch.size() * sizeof(magma_int_t)), "cudaMalloc");
            getrf_vbatched<T>(max_n, bt.v.size_array(bt.md, 0), bt.v.size_array(bt.md, 1),
                              bt.v.template pointer_array<T*>(bt.md, 0), bt.v.size_array(bt.md, 3), d_ipiv, d_info,
                              static_cast<magma_int_t>(batch.size()), work, ctx.queue());
            check_cuda(cudaStreamSynchronize(stream), "getrf");
            for (size_t i = 0; i < batch.size(); ++i) {
                const int b = batch[i];
                out[b].lu = download(a[i], static_cast<size_t>(n[b]) * n[b]);
                out[b].ipiv = download(ipiv[i], n[b]);
                cudaFree(a[i]);
                cudaFree(ipiv[i]);
            }
            cudaFree(d_ipiv);
            cudaFree(d_info);
        }
        return out;
    };
    const auto all = groupings(boxes, g);
    const std::vector<Out> ref = run(all[0]);
    for (size_t k = 1; k < all.size(); ++k) {
        const std::vector<Out> got = run(all[k]);
        Report r;
        for (int b = 0; b < boxes; ++b) {
            compare(ref[b].lu, got[b].lu, b, r, "LU");
            std::vector<int> ri(ref[b].ipiv.begin(), ref[b].ipiv.end()), gi(got[b].ipiv.begin(), got[b].ipiv.end());
            compare(ri, gi, b, r, "pivots");
        }
        print(by_class ? "getrf by size class" : "getrf (X_RR LU)", Scalar<T>::name(), all[k].name, r, boxes);
    }
}

// ---- triangular solves B := B U^{-1} and B := B L^{-1} (unit L)
// (trsm_vbatched, as the temp solves).  by_class: square B (as the
// eliminator's W = U^{-1} L^{-1}), each batch split into its size classes,
// launched with max_m = max_n = the class size.
template<typename T>
void test_trsm(int boxes, std::mt19937_64& g, magma_uplo_t uplo, magma_diag_t diag, const char* label,
               bool by_class = false) {
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    std::uniform_int_distribution<int> rs(4, 400), ms(1, 600);
    std::vector<int> r(static_cast<size_t>(boxes)), m(static_cast<size_t>(boxes));
    std::vector<std::vector<T>> u(static_cast<size_t>(boxes)), bm(static_cast<size_t>(boxes));
    for (int b = 0; b < boxes; ++b) {
        r[b] = rs(g);
        m[b] = by_class ? r[b] : ms(g);
        u[b] = random_matrix<T>(r[b], r[b], g);
        for (int i = 0; i < r[b]; ++i) u[b][static_cast<size_t>(i) * r[b] + i] += T(2.0 * r[b]);  // well conditioned
        bm[b] = random_matrix<T>(m[b], r[b], g);
    }
    auto run = [&](const Grouping& gr) {
        std::vector<std::vector<T>> out(static_cast<size_t>(boxes));
        std::vector<std::vector<int>> launches;
        for (const auto& batch : gr.batches) {
            if (!by_class) {
                launches.push_back(batch);
                continue;
            }
            std::map<int, std::vector<int>> parts;
            for (int b : batch) parts[lu_size_class(r[b])].push_back(b);
            for (auto& kv : parts) launches.push_back(std::move(kv.second));
        }
        for (const auto& batch : launches) {
            Batch<T> bt;
            std::vector<T*> a, c;
            for (int b : batch) {
                a.push_back(upload(u[b]));
                c.push_back(upload(bm[b]));
                bt.v.entries.push_back({a.back(), nullptr, c.back(), m[b], r[b], 1, r[b], 1, m[b]});
            }
            bt.stage(stream);
            const magma_int_t max_m = by_class ? lu_size_class(bt.v.max_n) : bt.v.max_m;
            const magma_int_t max_n = by_class ? lu_size_class(bt.v.max_n) : bt.v.max_n;
            trsm_vbatched<T>(MagmaRight, uplo, MagmaNoTrans, diag, max_m, max_n,
                             bt.v.size_array(bt.md, 0), bt.v.size_array(bt.md, 1), T(1.0),
                             bt.v.template pointer_array<T*>(bt.md, 0), bt.v.size_array(bt.md, 3),
                             bt.v.template pointer_array<T*>(bt.md, 2), bt.v.size_array(bt.md, 5),
                             static_cast<magma_int_t>(batch.size()), ctx.queue());
            check_cuda(cudaStreamSynchronize(stream), "trsm");
            for (size_t i = 0; i < batch.size(); ++i) {
                const int b = batch[i];
                out[b] = download(c[i], static_cast<size_t>(m[b]) * r[b]);
                cudaFree(a[i]);
                cudaFree(c[i]);
            }
        }
        return out;
    };
    const auto all = groupings(boxes, g);
    const auto ref = run(all[0]);
    for (size_t k = 1; k < all.size(); ++k) {
        const auto got = run(all[k]);
        Report rep;
        for (int b = 0; b < boxes; ++b) compare(ref[b], got[b], b, rep, "B A^-1");
        print(label, Scalar<T>::name(), all[k].name, rep, boxes);
    }
}

// ---- C := C - op(A) op(B) (VBatch::gemm: MAGMA, or the tensor cores)
template<typename T>
void test_gemm(int boxes, std::mt19937_64& g, magma_trans_t ta, magma_trans_t tb, const char* label) {
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    std::uniform_int_distribution<int> dims(1, 600);
    struct Shape {
        int m, n, k;
    };
    std::vector<Shape> s(static_cast<size_t>(boxes));
    std::vector<std::vector<T>> a(static_cast<size_t>(boxes)), bm(static_cast<size_t>(boxes)),
        c(static_cast<size_t>(boxes));
    for (int b = 0; b < boxes; ++b) {
        s[b] = {dims(g), dims(g), dims(g)};
        a[b] = random_matrix<T>(s[b].m, s[b].k, g);  // op(A) is m x k either way (lda below)
        bm[b] = random_matrix<T>(s[b].k, s[b].n, g);
        c[b] = random_matrix<T>(s[b].m, s[b].n, g);
    }
    auto run = [&](const Grouping& gr) {
        std::vector<std::vector<T>> out(static_cast<size_t>(boxes));
        for (const auto& batch : gr.batches) {
            Batch<T> bt;
            std::vector<T*> da, db, dc;
            for (int b : batch) {
                da.push_back(upload(a[b]));
                db.push_back(upload(bm[b]));
                dc.push_back(upload(c[b]));
                const int lda = ta == MagmaNoTrans ? s[b].m : s[b].k;
                const int ldb = tb == MagmaNoTrans ? s[b].k : s[b].n;
                bt.v.entries.push_back({da.back(), db.back(), dc.back(), s[b].m, s[b].n, s[b].k, lda, ldb, s[b].m});
            }
            bt.stage(stream);
            bt.v.gemm(bt.md, ta, tb, T(-1.0), T(1.0), ctx.queue());
            check_cuda(cudaStreamSynchronize(stream), "gemm");
            for (size_t i = 0; i < batch.size(); ++i) {
                const int b = batch[i];
                out[b] = download(dc[i], static_cast<size_t>(s[b].m) * s[b].n);
                cudaFree(da[i]);
                cudaFree(db[i]);
                cudaFree(dc[i]);
            }
        }
        return out;
    };
    const auto all = groupings(boxes, g);
    const auto ref = run(all[0]);
    for (size_t k = 1; k < all.size(); ++k) {
        const auto got = run(all[k]);
        Report rep;
        for (int b = 0; b < boxes; ++b) compare(ref[b], got[b], b, rep, "C");
        print(label, Scalar<T>::name(), all[k].name, rep, boxes);
    }
}


template<typename T>
void test_type(int boxes, uint64_t seed) {
    std::mt19937_64 g(seed);
    test_qrcp<T>(boxes, g, 1e-4, false);
    test_qrcp<T>(boxes, g, 1e-4, true);
    test_getrf<T>(boxes, g, false);
    test_getrf<T>(boxes, g, true);
    test_trsm<T>(boxes, g, MagmaUpper, MagmaNonUnit, "trsm (right, upper)");
    test_trsm<T>(boxes, g, MagmaLower, MagmaUnit, "trsm (right, lower, unit)");
    test_trsm<T>(boxes, g, MagmaUpper, MagmaNonUnit, "trsm (right, upper) by class", true);
    test_trsm<T>(boxes, g, MagmaLower, MagmaUnit, "trsm (right, lower, unit) by class", true);
    tensor_core_gemm() = false;
    test_gemm<T>(boxes, g, MagmaNoTrans, MagmaNoTrans, "gemm NN (MAGMA)");
    test_gemm<T>(boxes, g, MagmaTrans, MagmaNoTrans, "gemm TN (MAGMA)");
    test_gemm<T>(boxes, g, MagmaNoTrans, MagmaTrans, "gemm NT (MAGMA)");
    if constexpr (std::is_same_v<T, double>) {
        tensor_core_gemm() = true;
        test_gemm<T>(boxes, g, MagmaNoTrans, MagmaNoTrans, "gemm NN (tensor cores)");
        test_gemm<T>(boxes, g, MagmaTrans, MagmaNoTrans, "gemm TN (tensor cores)");
        test_gemm<T>(boxes, g, MagmaNoTrans, MagmaTrans, "gemm NT (tensor cores)");
        tensor_core_gemm() = false;
    }
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    const int boxes = argc > 1 ? std::atoi(argv[1]) : 400;
    Context::instance().activate();
    std::printf("batch determinism: %d boxes per kernel; every grouping against each box alone\n", boxes);
    test_type<double>(boxes, 12345);
    test_type<dcomplex>(boxes, 67890);
    MPI_Finalize();
    return 0;
}
