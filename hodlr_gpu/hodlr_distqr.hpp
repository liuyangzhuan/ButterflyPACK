#pragma once
// A truncated SVD by TSQR of a tall matrix whose rows are spread over the
// ranks of a group, on the GPUs, for the merges of the pieces of a HODLR
// block compressed by several ranks (LR_HMerge on the GPU, tsqr): each
// rank's QR of its rows, up a binary tree over the group the QR of pairs of
// stacked R factors, on the root the SVD of the last R and the rank by the
// rule of PSVD_Truncate, down the tree the Q of each pair applied to the
// root's left singular vectors, then each rank's rows of them.  The merges'
// factors stay in device memory (dm_*: device matrices by id for the Fortran
// side, their redistribution between the ranks, the products of the other
// side); messages go between the ranks' device memories when MPI is
// CUDA-aware (in the exchange arena), else through pinned host buffers.
// Column major throughout.

#ifdef H2_HAVE_GPU

#include "color_gpu/device_heap.hpp"
#include "color_gpu/device_scalar.hpp"
#include "color_gpu/gpu_runtime.hpp"
#include "hodlr_svd.hpp"

#include <cusolverDn.h>
#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace bpack {
namespace gpu {

template<typename T>
class DistQr {
public:
    DistQr() = default;
    DistQr(const DistQr&) = delete;
    DistQr& operator=(const DistQr&) = delete;
    ~DistQr() {
        for (Factor& f : tree_) release(f);
        if (xnew_ != nullptr) fmm::gpu::DeviceHeap::instance().free(xnew_);
        for (auto& kv : mats_) {
            if (kv.second.owned) fmm::gpu::DeviceHeap::instance().free(kv.second.d);
        }
        if (params_ != nullptr) cusolverDnDestroyParams(params_);
        if (handle_ != nullptr) cusolverDnDestroy(handle_);
    }

    // TSQR of X (r columns) whose rows are spread over the ranks of comm, x
    // this rank's l rows (device, leading dimension l; overwritten by its
    // Householder vectors): returns the rank rn by the rule of PSVD_Truncate
    // (relative tolerance tol, values up to underflow as zero, then zero
    // factors of rank 1), and keeps this rank's rows xnew (l x rn) of the left
    // singular vectors kept and sw = diag(s) Vs^H (rn x r) for dm_tsqr.
    // tsec (5): seconds of the local QR, up the tree, the root's SVD, down
    // the tree, the local product
    int tsqr(MPI_Comm comm, int l, int r, T* x, double tol, double underflow, double* tsec) {
        using clock = std::chrono::steady_clock;
        using fmm::gpu::check_cuda;
        constexpr int tag_up = 7501, tag_down = 7502;
        fmm::gpu::Context& ctx = fmm::gpu::Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        if (r <= 0) throw std::invalid_argument("DistQr::tsqr: no columns");
        if (xnew_ != nullptr) heap.free(xnew_);
        xnew_ = nullptr;
        int p = 1, crank = 0;
        MPI_Comm_size(comm, &p);
        MPI_Comm_rank(comm, &crank);
        auto t0 = clock::now();
        auto lap = [&](int i) {
            const auto t = clock::now();
            tsec[i] = std::chrono::duration<double>(t - t0).count();
            t0 = t;
        };
        const size_t rr = static_cast<size_t>(r) * r;

        // this rank's QR, and its R (r x r, zero below the diagonal and from row l on)
        Factor loc;
        loc.m = l;
        loc.n = r;
        loc.a = x;
        loc.owned = false;
        Msg R = msg_alloc(rr);
        laset(r, r, R.d, r);
        if (l > 0) {
            loc.tau = heap.alloc<T>(static_cast<size_t>(std::min(l, r)) * sizeof(T));
            geqrf(loc);
            lacpy_upper(std::min(l, r), r, x, l, R.d, r);
        }
        check_cuda(cudaStreamSynchronize(stream), "TSQR local R");
        lap(0);

        // up the tree: rank c takes c + stride's R at the steps where c is a
        // multiple of 2 stride; stride ends as the step this rank sends at
        int stride = 1;
        while (stride < p) {
            if (crank % (2 * stride) != 0) break;
            if (crank + stride < p) {
                Msg P = msg_alloc(rr);
                recv(comm, crank + stride, tag_up, P);
                Factor f;
                f.m = 2 * r;
                f.n = r;
                f.a = heap.alloc<T>(2 * rr * sizeof(T));
                f.tau = heap.alloc<T>(static_cast<size_t>(r) * sizeof(T));
                const size_t col = static_cast<size_t>(r) * sizeof(T);
                check_cuda(cudaMemcpy2DAsync(f.a, 2 * col, R.d, col, col, r, cudaMemcpyDeviceToDevice, stream), "TSQR stack");
                check_cuda(cudaMemcpy2DAsync(f.a + r, 2 * col, P.d, col, col, r, cudaMemcpyDeviceToDevice, stream),
                           "TSQR stack");
                geqrf(f);
                msg_free(P);
                laset(r, r, R.d, r);
                lacpy_upper(r, r, f.a, 2 * r, R.d, r);
                tree_.push_back(f);
            }
            stride *= 2;
        }
        check_cuda(cudaStreamSynchronize(stream), "TSQR up");
        if (crank != 0) send(comm, crank - stride, tag_up, R);
        lap(1);

        // the root: the SVD of the last R, the rank, b = Us(:, 0 .. rn - 1)
        // (r x rn) and w = Vs(:, 0 .. rn - 1)^H (rn x r)
        int rn = 1;
        std::vector<double> s(static_cast<size_t>(r), 0.0);
        std::vector<T> w;
        Msg B{};
        if (crank == 0) {
            T* du = heap.alloc<T>(rr * sizeof(T));
            T* dv = heap.alloc<T>(rr * sizeof(T));
            double* ds = heap.alloc<double>(static_cast<size_t>(r) * sizeof(double));
            svd_.svd_device(r, r, R.d, ds, du, dv);
            flops_ += kFlopScale * (8.0 + 4.0 / 3.0) * r * static_cast<double>(r) * r;  // (as the CPU counts gesdd)
            check_cuda(cudaMemcpy(s.data(), ds, static_cast<size_t>(r) * sizeof(double), cudaMemcpyDeviceToHost), "TSQR s");
            const bool zero = !(s[0] > underflow);
            if (!zero) {
                rn = r;
                for (int i = 0; i < r; ++i) {
                    if (s[i] / s[0] <= tol || s[i] <= underflow) {
                        rn = i + 1;
                        if (s[i] <= underflow) break;
                        if (s[i] < s[0] * tol / 10) rn = i;
                        break;
                    }
                }
            }
            const size_t nb = static_cast<size_t>(r) * rn;
            B = msg_alloc(nb);
            std::vector<T> vt(nb);
            if (zero) {  // (zero factors of rank 1)
                std::fill(s.begin(), s.end(), 0.0);
                laset(r, rn, B.d, r);
                std::fill(vt.begin(), vt.end(), T(0.0));
            } else {
                check_cuda(cudaMemcpyAsync(B.d, du, nb * sizeof(T), cudaMemcpyDeviceToDevice, stream), "TSQR b");
                check_cuda(cudaMemcpyAsync(vt.data(), dv, nb * sizeof(T), cudaMemcpyDeviceToHost, stream), "TSQR v");
            }
            check_cuda(cudaStreamSynchronize(stream), "TSQR root");
            w.resize(nb);
            for (int i = 0; i < rn; ++i) {
                for (int j = 0; j < r; ++j) {
                    const T v = vt[static_cast<size_t>(j) + static_cast<size_t>(i) * r];
                    if constexpr (fmm::gpu::is_complex_scalar<T>) {
                        w[static_cast<size_t>(i) + static_cast<size_t>(j) * rn] = T{v.re, -v.im};
                    } else {
                        w[static_cast<size_t>(i) + static_cast<size_t>(j) * rn] = v;
                    }
                }
            }
            heap.free(ds);
            heap.free(dv);
            heap.free(du);
        }
        msg_free(R);
        lap(2);

        // down the tree, the stacked pairs' Q in reverse order
        // (the collectives first: a rank waiting for its parent's message
        // must not hold up a broadcast the parent is still in)
        MPI_Bcast(&rn, 1, MPI_INT, 0, comm);
        MPI_Bcast(s.data(), r, MPI_DOUBLE, 0, comm);
        const size_t nb = static_cast<size_t>(r) * rn;
        if (crank != 0) w.resize(nb);
        bcast_host(comm, w.data(), nb);
        if (crank != 0) {
            B = msg_alloc(nb);
            recv(comm, crank - stride, tag_down, B);
        }
        for (stride /= 2; stride >= 1; stride /= 2) {
            if (crank + stride >= p) continue;
            if (tree_.empty()) throw std::logic_error("DistQr::tsqr: tree out of step");
            Factor f = tree_.back();
            tree_.pop_back();
            T* c = apply_q_device(f, rn, B.d, r);  // (2r x rn)
            const size_t col = static_cast<size_t>(r) * sizeof(T);
            Msg C = msg_alloc(nb);
            check_cuda(cudaMemcpy2DAsync(B.d, col, c, 2 * col, col, rn, cudaMemcpyDeviceToDevice, stream), "TSQR top");
            check_cuda(cudaMemcpy2DAsync(C.d, col, c + r, 2 * col, col, rn, cudaMemcpyDeviceToDevice, stream), "TSQR bottom");
            check_cuda(cudaStreamSynchronize(stream), "TSQR down");
            heap.free(c);
            release(f);
            send(comm, crank + stride, tag_down, C);
            msg_free(C);
        }
        if (!tree_.empty()) throw std::logic_error("DistQr::tsqr: tree out of step");
        lap(3);

        // this rank's rows: xnew = Q_loc [b; 0], and sw = diag(s) w
        if (l > 0) {
            xnew_ = apply_q_device(loc, rn, B.d, r);
            heap.free(loc.tau);
        }
        msg_free(B);
        sw_.resize(nb);
        for (int j = 0; j < r; ++j) {
            for (int i = 0; i < rn; ++i) {
                const size_t k = static_cast<size_t>(i) + static_cast<size_t>(j) * rn;
                sw_[k] = s[static_cast<size_t>(i)] * w[k];
            }
        }
        l_ = l;
        rn_ = rn;
        r_ = r;
        lap(4);
        return rn;
    }

    // ---- the device matrices of the merges (ids for the Fortran side) ----

    // a rows x cols device matrix (leading dimension rows): zero, or host's
    // (column major); returns its id
    int dm_alloc(int rows, int cols, const T* host) {
        const int id = dm_alloc_raw(rows, cols);
        const cudaStream_t stream = fmm::gpu::Context::instance().stream();
        DMat& m = mat(id);
        const size_t count = static_cast<size_t>(m.rows) * m.cols;
        if (count > 0) {
            if (host != nullptr) {
                fmm::gpu::check_cuda(cudaMemcpyAsync(m.d, host, count * sizeof(T), cudaMemcpyHostToDevice, stream), "merge upload");
            } else {
                fmm::gpu::check_cuda(cudaMemsetAsync(m.d, 0, count * sizeof(T), stream), "merge zero");
            }
        }
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "merge matrix");
        return id;
    }
    // (uninitialized)
    int dm_alloc_raw(int rows, int cols) {
        fmm::gpu::Context::instance().activate();
        DMat m;
        m.rows = std::max(rows, 0);
        m.cols = std::max(cols, 0);
        m.d = fmm::gpu::DeviceHeap::instance().alloc<T>(std::max<size_t>(static_cast<size_t>(m.rows) * m.cols, 1) * sizeof(T));
        const int id = next_id_++;
        mats_[id] = m;
        return id;
    }
    // a rows x cols device matrix copied from device memory src (leading
    // dimension rows); returns its id
    int dm_copy(int rows, int cols, const T* src) {
        const int id = dm_alloc_raw(rows, cols);
        const DMat& m = mat(id);
        const size_t count = static_cast<size_t>(m.rows) * m.cols;
        if (count > 0) fmm::gpu::check_cuda(cudaMemcpy(m.d, src, count * sizeof(T), cudaMemcpyDeviceToDevice), "merge copy");
        return id;
    }
    void dm_download(int id, T* host) {
        const DMat& m = mat(id);
        fmm::gpu::Context::instance().activate();
        const size_t count = static_cast<size_t>(m.rows) * m.cols;
        if (count > 0) fmm::gpu::check_cuda(cudaMemcpy(host, m.d, count * sizeof(T), cudaMemcpyDeviceToHost), "merge download");
    }
    // Device matrix id, whose copy on the host is host: to keep as the device
    // HODLR's mirror of that host array (its memory handed over)
    template<typename Keeper>
    void dm_keep(int id, const T* host, Keeper& keep, bool stale = false) {
        auto it = mats_.find(id);
        if (it == mats_.end()) return;
        const size_t bytes = static_cast<size_t>(it->second.rows) * it->second.cols * sizeof(T);
        if (bytes > 0) keep.add_mirror(host, it->second.d, bytes, stale);
        keep.own_mirror_memory(it->second.d);
        mats_.erase(it);
    }
    void dm_free(int id) {
        auto it = mats_.find(id);
        if (it == mats_.end()) return;
        fmm::gpu::Context::instance().activate();
        if (it->second.owned) fmm::gpu::DeviceHeap::instance().free(it->second.d);
        mats_.erase(it);
    }
    // a device matrix id viewing d (rows x cols, leading dimension rows; not freed by dm_free)
    int dm_wrap(const T* d, int rows, int cols) {
        DMat m;
        m.d = const_cast<T*>(d);
        m.rows = rows;
        m.cols = cols;
        m.owned = false;
        const int id = next_id_++;
        mats_[id] = m;
        return id;
    }
    // device matrix id's memory, handed over to the caller (the id forgotten)
    T* dm_give(int id) {
        auto it = mats_.find(id);
        if (it == mats_.end() || !it->second.owned) throw std::logic_error("merge: no device matrix to give");
        T* d = it->second.d;
        mats_.erase(it);
        return d;
    }

    // Rows of device matrix src (ncols columns; 0: none on this rank) into
    // columns c0 .. c0 + ncols - 1 of dst (0: none) between the ranks of
    // comm: send i takes rows s_off[i] .. s_off[i] + s_rows[i] - 1 to rank
    // s_rank[i], receive j puts rows from rank r_rank[j] at row r_off[j]
    // (the plan of Redistribute1Dto1D; this rank's own part copied directly)
    void dm_redistribute(MPI_Comm comm, int src, int ncols, int dst, int c0, int nsend, const int* s_rank, const int* s_off,
                         const int* s_rows, int nrecv, const int* r_rank, const int* r_off, const int* r_rows, int tag) {
        using fmm::gpu::check_cuda;
        fmm::gpu::Context::instance().activate();
        const cudaStream_t stream = fmm::gpu::Context::instance().stream();
        int me = 0;
        MPI_Comm_rank(comm, &me);
        const DMat* S = src > 0 ? &mat(src) : nullptr;
        const DMat* D = dst > 0 ? &mat(dst) : nullptr;
        const size_t sz = sizeof(T);
        if (ncols <= 0) return;
        std::vector<Msg> smsg(static_cast<size_t>(std::max(nsend, 0))), rmsg(static_cast<size_t>(std::max(nrecv, 0)));
        std::vector<size_t> sstage(smsg.size(), 0), rstage(rmsg.size(), 0);
        size_t staged = 0;  // (elements, each segment rounded to 64)
        auto stage = [&](const Msg& m) {
            const size_t at = staged;
            staged += (m.count + 63) / 64 * 64;
            return at;
        };
        for (int j = 0; j < nrecv; ++j) {
            if (r_rank[j] == me || r_rows[j] <= 0) continue;
            if (D == nullptr) throw std::logic_error("merge redistribution: rows to a rank without the matrix");
            rmsg[j] = msg_alloc(static_cast<size_t>(r_rows[j]) * ncols);
            if (!rmsg[j].arena) rstage[j] = stage(rmsg[j]);
        }
        for (int i = 0; i < nsend; ++i) {
            if (s_rows[i] <= 0) continue;
            if (S == nullptr) throw std::logic_error("merge redistribution: rows from a rank without the matrix");
            if (s_rank[i] == me) {  // (this rank's own rows)
                int j = 0;
                while (j < nrecv && r_rank[j] != me) ++j;
                if (j == nrecv || r_rows[j] != s_rows[i] || D == nullptr) {
                    throw std::logic_error("merge redistribution: own rows without their place");
                }
                check_cuda(cudaMemcpy2DAsync(D->d + r_off[j] + static_cast<size_t>(c0) * D->rows, D->rows * sz, S->d + s_off[i],
                                             S->rows * sz, s_rows[i] * sz, ncols, cudaMemcpyDeviceToDevice, stream),
                           "merge own rows");
                continue;
            }
            smsg[i] = msg_alloc(static_cast<size_t>(s_rows[i]) * ncols);
            if (!smsg[i].arena) sstage[i] = stage(smsg[i]);
            check_cuda(cudaMemcpy2DAsync(smsg[i].d, s_rows[i] * sz, S->d + s_off[i], S->rows * sz, s_rows[i] * sz, ncols,
                                         cudaMemcpyDeviceToDevice, stream),
                       "merge pack");
        }
        check_cuda(cudaStreamSynchronize(stream), "merge pack");
        T* host = staged > 0 ? reinterpret_cast<T*>(staging(staged)) : nullptr;
        std::vector<MPI_Request> req;
        req.reserve(smsg.size() + rmsg.size());
        for (size_t j = 0; j < rmsg.size(); ++j) {
            if (rmsg[j].d == nullptr) continue;
            void* buf = rmsg[j].arena ? static_cast<void*>(rmsg[j].d) : static_cast<void*>(host + rstage[j]);
            req.emplace_back();
            MPI_Irecv(buf, doubles(rmsg[j].count), MPI_DOUBLE, r_rank[j], tag, comm, &req.back());
        }
        for (size_t i = 0; i < smsg.size(); ++i) {
            if (smsg[i].d == nullptr) continue;
            void* buf = smsg[i].d;
            if (!smsg[i].arena) {
                buf = host + sstage[i];
                check_cuda(cudaMemcpy(buf, smsg[i].d, smsg[i].count * sz, cudaMemcpyDeviceToHost), "merge send staging");
            }
            req.emplace_back();
            MPI_Isend(buf, doubles(smsg[i].count), MPI_DOUBLE, s_rank[i], tag, comm, &req.back());
        }
        if (!req.empty()) MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
        for (size_t j = 0; j < rmsg.size(); ++j) {
            if (rmsg[j].d == nullptr) continue;
            if (!rmsg[j].arena) {
                check_cuda(cudaMemcpyAsync(rmsg[j].d, host + rstage[j], rmsg[j].count * sz, cudaMemcpyHostToDevice, stream),
                           "merge receive staging");
            }
            check_cuda(cudaMemcpy2DAsync(D->d + r_off[j] + static_cast<size_t>(c0) * D->rows, D->rows * sz, rmsg[j].d,
                                         r_rows[j] * sz, r_rows[j] * sz, ncols, cudaMemcpyDeviceToDevice, stream),
                       "merge unpack");
        }
        check_cuda(cudaStreamSynchronize(stream), "merge unpack");
        for (Msg& m : smsg) msg_free(m);
        for (Msg& m : rmsg) msg_free(m);
    }

    // tsqr() of device matrix x (overwritten; its rows this rank's): the
    // rank rn; *xnew: the id of this rank's rows of the left singular vectors
    // kept (rows of x x rn); sw for dm_fetch_sw
    int dm_tsqr(MPI_Comm comm, int x, double tol, double underflow, double* tsec, int* xnew) {
        DMat& X = mat(x);
        const int rn = tsqr(comm, X.rows, X.cols, X.d, tol, underflow, tsec);
        DMat m;
        m.rows = X.rows;
        m.cols = rn;
        m.d = xnew_ != nullptr ? xnew_ : fmm::gpu::DeviceHeap::instance().alloc<T>(sizeof(T));
        xnew_ = nullptr;
        *xnew = next_id_++;
        mats_[*xnew] = m;
        return rn;
    }
    // sw of the last tsqr (rn x r, host)
    void dm_fetch_sw(T* sw) {
        std::copy(sw_.begin(), sw_.end(), sw);
        sw_.clear();
    }

    // the flops since the last call (TSQR and dm_gemm_nt, counted as the CPU
    // counts them: 8 m n k for a complex gemm)
    double take_flops() {
        const double f = flops_;
        flops_ = 0.0;
        return f;
    }

    // rows c_row0 .. c_row0 + m - 1 of device matrix c = (device matrix a's
    // first m rows, k columns) b^T, b (n x k) on the host
    void dm_gemm_nt(int a, int m, int k, const T* b, int n, int c, int c_row0) {
        if (m <= 0 || n <= 0 || k <= 0) return;
        flops_ += kFlopScale * 2.0 * m * static_cast<double>(n) * k;
        const DMat& A = mat(a);
        DMat& C = mat(c);
        fmm::gpu::Context& ctx = fmm::gpu::Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        T* db = heap.alloc<T>(static_cast<size_t>(n) * k * sizeof(T));
        fmm::gpu::check_cuda(cudaMemcpyAsync(db, b, static_cast<size_t>(n) * k * sizeof(T), cudaMemcpyHostToDevice, stream),
                             "merge product B");
        T* dc = C.d + c_row0;
        if constexpr (std::is_same_v<T, double>) {
            magma_dgemm(MagmaNoTrans, MagmaTrans, m, n, k, 1.0, A.d, A.rows, db, n, 0.0, dc, C.rows, ctx.queue());
        } else {
            magma_zgemm(MagmaNoTrans, MagmaTrans, m, n, k, MAGMA_Z_MAKE(1.0, 0.0),
                        reinterpret_cast<const magmaDoubleComplex*>(A.d), A.rows, reinterpret_cast<magmaDoubleComplex*>(db), n,
                        MAGMA_Z_MAKE(0.0, 0.0), reinterpret_cast<magmaDoubleComplex*>(dc), C.rows, ctx.queue());
        }
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "merge product");
        heap.free(db);
    }

private:
    struct Factor {
        int m = 0, n = 0;     // rows, columns
        T* a = nullptr;       // m x n: R and the Householder vectors
        T* tau = nullptr;     // min(m, n)
        bool owned = true;    // (a from the heap, freed by release)
    };

    // a device buffer MPI reads or writes: in the exchange arena when MPI is
    // CUDA-aware and the arena has room (then sent from and received into
    // directly), else from the heap (staged through pinned host memory)
    struct Msg {
        T* d = nullptr;
        size_t count = 0;
        bool arena = false;
    };

    static constexpr int kPer = fmm::gpu::is_complex_scalar<T> ? 2 : 1;  // doubles per scalar
    static constexpr double kFlopScale = fmm::gpu::is_complex_scalar<T> ? 4.0 : 1.0;
    double flops_ = 0.0;

    struct DMat {
        T* d = nullptr;
        int rows = 0, cols = 0;
        bool owned = true;  // (false: a view of memory kept elsewhere)
    };
    DMat& mat(int id) {
        auto it = mats_.find(id);
        if (it == mats_.end()) throw std::logic_error("merge: no device matrix " + std::to_string(id));
        return it->second;
    }

    static Msg msg_alloc(size_t count) {
        Msg m;
        m.count = count;
        const size_t bytes = std::max<size_t>(count, 1) * sizeof(T);
        if (fmm::gpu::device_exchange_enabled()) {
            m.d = reinterpret_cast<T*>(fmm::gpu::DeviceHeap::exchange_arena().try_alloc(bytes));
        }
        m.arena = m.d != nullptr;
        if (m.d == nullptr) m.d = fmm::gpu::DeviceHeap::instance().alloc<T>(bytes);
        return m;
    }
    static void msg_free(Msg& m) {
        if (m.d == nullptr) return;
        if (m.arena) {
            fmm::gpu::DeviceHeap::exchange_arena().free(m.d);
        } else {
            fmm::gpu::DeviceHeap::instance().free(m.d);
        }
        m = Msg{};
    }
    static int doubles(size_t count) {
        const size_t n = count * kPer;
        if (n > static_cast<size_t>(std::numeric_limits<int>::max())) throw std::runtime_error("TSQR: message too long for MPI");
        return static_cast<int>(n);
    }
    static double* staging(size_t count) {
        static fmm::gpu::PinnedBuffer* buffer = new fmm::gpu::PinnedBuffer;  // never freed (CUDA shutdown order)
        return static_cast<double*>(buffer->reserve(std::max<size_t>(count, 1) * sizeof(T)));
    }
    // (the stream's work on m done)
    static void send(MPI_Comm comm, int dest, int tag, const Msg& m) {
        if (m.arena) {
            MPI_Send(m.d, doubles(m.count), MPI_DOUBLE, dest, tag, comm);
            return;
        }
        double* h = staging(m.count);
        fmm::gpu::check_cuda(cudaMemcpy(h, m.d, m.count * sizeof(T), cudaMemcpyDeviceToHost), "TSQR send");
        MPI_Send(h, doubles(m.count), MPI_DOUBLE, dest, tag, comm);
    }
    static void recv(MPI_Comm comm, int src, int tag, Msg& m) {
        if (m.arena) {
            MPI_Recv(m.d, doubles(m.count), MPI_DOUBLE, src, tag, comm, MPI_STATUS_IGNORE);
            return;
        }
        double* h = staging(m.count);
        MPI_Recv(h, doubles(m.count), MPI_DOUBLE, src, tag, comm, MPI_STATUS_IGNORE);
        fmm::gpu::check_cuda(cudaMemcpy(m.d, h, m.count * sizeof(T), cudaMemcpyHostToDevice), "TSQR receive");
    }
    static void bcast_host(MPI_Comm comm, T* a, size_t count) {
        if (count > 0) MPI_Bcast(a, doubles(count), MPI_DOUBLE, 0, comm);
    }

    static void check(cusolverStatus_t st, const char* what) {
        if (st != CUSOLVER_STATUS_SUCCESS) {
            throw std::runtime_error(std::string(what) + " failed with status " + std::to_string(static_cast<int>(st)));
        }
    }

    cusolverDnHandle_t solver() {
        fmm::gpu::Context& ctx = fmm::gpu::Context::instance();
        ctx.activate();
        if (handle_ == nullptr) {
            check(cusolverDnCreate(&handle_), "cusolverDnCreate");
            check(cusolverDnCreateParams(&params_), "cusolverDnCreateParams");
        }
        check(cusolverDnSetStream(handle_, ctx.stream()), "cusolverDnSetStream");
        return handle_;
    }

    void release(Factor& f) {
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        if (f.a != nullptr && f.owned) heap.free(f.a);
        if (f.tau != nullptr) heap.free(f.tau);
        f = Factor{};
    }

    // a (m x n, leading dimension ld, device) = 0
    static void laset(int m, int n, T* a, int ld) {
        if (m <= 0 || n <= 0) return;
        const magma_queue_t q = fmm::gpu::Context::instance().queue();
        if constexpr (std::is_same_v<T, double>) {
            magmablas_dlaset(MagmaFull, m, n, 0.0, 0.0, a, ld, q);
        } else {
            magmablas_zlaset(MagmaFull, m, n, MAGMA_Z_ZERO, MAGMA_Z_ZERO, reinterpret_cast<magmaDoubleComplex*>(a), ld, q);
        }
    }
    // the upper trapezoid of a (m x n, leading dimension lda) into b (leading dimension ldb)
    static void lacpy_upper(int m, int n, const T* a, int lda, T* b, int ldb) {
        if (m <= 0 || n <= 0) return;
        const magma_queue_t q = fmm::gpu::Context::instance().queue();
        if constexpr (std::is_same_v<T, double>) {
            magmablas_dlacpy(MagmaUpper, m, n, a, lda, b, ldb, q);
        } else {
            magmablas_zlacpy(MagmaUpper, m, n, reinterpret_cast<const magmaDoubleComplex*>(a), lda,
                             reinterpret_cast<magmaDoubleComplex*>(b), ldb, q);
        }
    }

    // QR of f.a (on the device) in place
    void geqrf(Factor& f) {
        {
            const double k = std::min(f.m, f.n), mx = std::max(f.m, f.n);
            flops_ += kFlopScale * (2.0 * mx * k * k - 2.0 / 3.0 * k * k * k);
        }
        cusolverDnHandle_t h = solver();
        const cudaStream_t stream = fmm::gpu::Context::instance().stream();
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        const cudaDataType type = fmm::gpu::is_complex_scalar<T> ? CUDA_C_64F : CUDA_R_64F;
        size_t dbytes = 0, hbytes = 0;
        check(cusolverDnXgeqrf_bufferSize(h, params_, f.m, f.n, type, f.a, f.m, type, f.tau, type, &dbytes, &hbytes),
              "cusolverDnXgeqrf_bufferSize");
        const size_t wbytes = fmm::gpu::align_up(std::max<size_t>(dbytes, 1));
        char* work = heap.alloc<char>(wbytes + sizeof(int));
        std::vector<char> hwork(std::max<size_t>(hbytes, 1));
        check(cusolverDnXgeqrf(h, params_, f.m, f.n, type, f.a, f.m, type, f.tau, type, work, dbytes, hwork.data(), hbytes,
                               reinterpret_cast<int*>(work + wbytes)),
              "cusolverDnXgeqrf");
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "TSQR factor");
        heap.free(work);
    }

    // c (f.m x rn, from the heap) = Q_f [b(0 .. min(f.m, f.n) - 1, :); 0], b
    // (device or host) with leading dimension ldb
    T* apply_q_device(const Factor& f, int rn, const T* b, int ldb) {
        {
            const double k = std::min(f.m, f.n);
            flops_ += kFlopScale * (4.0 * f.m * k * rn - 2.0 * k * k * rn);
        }
        cusolverDnHandle_t h = solver();
        const cudaStream_t stream = fmm::gpu::Context::instance().stream();
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        const int k = std::min(f.m, f.n);
        const size_t count = static_cast<size_t>(f.m) * rn;
        T* c = heap.alloc<T>(std::max<size_t>(count, 1) * sizeof(T));
        fmm::gpu::check_cuda(cudaMemsetAsync(c, 0, count * sizeof(T), stream), "TSQR apply");
        fmm::gpu::check_cuda(cudaMemcpy2DAsync(c, static_cast<size_t>(f.m) * sizeof(T), b, static_cast<size_t>(ldb) * sizeof(T),
                                               static_cast<size_t>(k) * sizeof(T), rn, cudaMemcpyDefault, stream),
                             "TSQR apply input");
        int lwork = 0;
        if constexpr (std::is_same_v<T, double>) {
            check(cusolverDnDormqr_bufferSize(h, CUBLAS_SIDE_LEFT, CUBLAS_OP_N, f.m, rn, k, f.a, f.m, f.tau, c, f.m, &lwork),
                  "cusolverDnDormqr_bufferSize");
        } else {
            check(cusolverDnZunmqr_bufferSize(h, CUBLAS_SIDE_LEFT, CUBLAS_OP_N, f.m, rn, k,
                                              reinterpret_cast<const cuDoubleComplex*>(f.a), f.m,
                                              reinterpret_cast<const cuDoubleComplex*>(f.tau),
                                              reinterpret_cast<const cuDoubleComplex*>(c), f.m, &lwork),
                  "cusolverDnZunmqr_bufferSize");
        }
        T* work = heap.alloc<T>((static_cast<size_t>(std::max(lwork, 1)) + 2) * sizeof(T));
        int* info = reinterpret_cast<int*>(work + std::max(lwork, 1));
        if constexpr (std::is_same_v<T, double>) {
            check(cusolverDnDormqr(h, CUBLAS_SIDE_LEFT, CUBLAS_OP_N, f.m, rn, k, f.a, f.m, f.tau, c, f.m, work, lwork, info),
                  "cusolverDnDormqr");
        } else {
            check(cusolverDnZunmqr(h, CUBLAS_SIDE_LEFT, CUBLAS_OP_N, f.m, rn, k, reinterpret_cast<const cuDoubleComplex*>(f.a),
                                   f.m, reinterpret_cast<const cuDoubleComplex*>(f.tau),
                                   reinterpret_cast<cuDoubleComplex*>(c), f.m, reinterpret_cast<cuDoubleComplex*>(work),
                                   lwork, info),
                  "cusolverDnZunmqr");
        }
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "TSQR apply");
        heap.free(work);
        return c;
    }

    cusolverDnHandle_t handle_ = nullptr;
    cusolverDnParams_t params_ = nullptr;
    std::vector<Factor> tree_;  // the stacked pairs' factors up the tree, the last on top
    T* xnew_ = nullptr;         // tsqr's result: this rank's rows of the left singular vectors (l_ x rn_)
    std::vector<T> sw_;         // tsqr's result: diag(s) Vs^H (rn_ x r_)
    int l_ = 0, rn_ = 0, r_ = 0;
    DeviceSvd<T> svd_;
    std::unordered_map<int, DMat> mats_;  // the merges' device matrices
    int next_id_ = 1;
};

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
