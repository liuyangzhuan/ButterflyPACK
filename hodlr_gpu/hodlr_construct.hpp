#pragma once
// Construction of a HODLR (format 1, LRlevel 0) on the GPU from the device
// kernel the application registered: the entries of the dense leaves, and
// the low-rank blocks of a level by BACA without overlap (LR_BACA_noOverlap,
// RecLR_leaf=5), all blocks of the level in lockstep.  The Fortran side
// (HODLR_gpu_construct_level, BPACK_constr.f90) keeps the host part of the
// CPU algorithm, step for step: the random first columns (rperm, in the
// CPU's block order), the rank-revealing QR of the r x r cores
// (geqp3modf90) and W = R11^{-1} Q1^H, the norm estimates and the stopping
// test, and the SVD of the recompression.  The device holds U (m x rank) and
// V^T (n x rank) of every block and does, per iteration:
//   panels: C = A(:, J) - U V(:, J); I = the first r pivots of the
//     column-pivoted QR of C^T with the chosen rows masked (geqp3 on
//     column_RT); R = A(I, :) - U(I, :) V; the core C(I, :) to the host;
//   append: U(:, new) = C(:, jpvt), V(new, :) = W R; the Gram matrices of
//     the norm update (LR_Fnorm, LR_FnormUp) to the host; the next columns J
//     = the pivots of R with the chosen columns masked;
// and at the end the recompression (LR_ReCompression): QR of U and V^T
// (MAGMA batched geqrf, zero padded), with the products by Q in compact WY
// form, Q [X; 0] = [X; 0] - Y (T (Y1^H X)), T and Y1^H X from the host.
// The optional k-nearest-neighbour step (option%knn) seeds the blocks from
// the panels of the listed rows and columns first.
// Point slots are tree indices - 1; rows and columns of a block are local
// and 0-based here (1-based in the Fortran arrays).

#ifdef H2_HAVE_GPU

#include "color_gpu/device_kernels.hpp"
#include "color_gpu/evaluator.hpp"
#include "hodlr_batch.hpp"
#include "hodlr_device.hpp"
#include "hodlr_kernels.hpp"
#include "hodlr_svd.hpp"

#include <cusolverDn.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace bpack {
namespace gpu {

template<typename T>
class HodlrConstruct {
public:
    struct Stats {
        double t_eval = 0.0;     // dense leaves (s)
        double t_panels = 0.0;   // panel evaluation, pivots, cores
        double t_append = 0.0;   // factor updates, Gram matrices, next columns
        double t_recomp = 0.0;   // QR and products of the recompression
        double flops = 0.0;      // device flops (estimate)
        double t_entry = 0.0;    // device time of the entry evaluations of the panels (part of t_panels)
    };

    HodlrConstruct() = default;
    HodlrConstruct(const HodlrConstruct&) = delete;
    HodlrConstruct& operator=(const HodlrConstruct&) = delete;
    ~HodlrConstruct() {
        end();
        release_points();
        if (solver_params_ != nullptr) cusolverDnDestroyParams(solver_params_);
        if (solver_ != nullptr) cusolverDnDestroy(solver_);
    }

    const Stats& stats() const { return stats_; }
    void reset_stats() { stats_ = Stats{}; }

    // The points of the matrix in tree order: coordinates (dim per point; dim
    // 0: none) and the 0-based original indices (the global ids the
    // application's evaluator receives).
    void set_points(int64_t n, int dim, const double* xyz, const int64_t* ids) {
        Context::instance().activate();
        release_points();
        DeviceHeap& heap = DeviceHeap::instance();
        const size_t count = static_cast<size_t>(std::max<int64_t>(n, 1));
        dim_ = std::max(dim, 0);
        xyz_ = heap.alloc_resident<double>(std::max<size_t>(static_cast<size_t>(dim_) * count, 1) * sizeof(double));
        ids_ = heap.alloc_resident<int64_t>(count * sizeof(int64_t));
        if (n > 0) {
            if (dim_ > 0) {
                check_cuda(cudaMemcpy(xyz_, xyz, static_cast<size_t>(dim_) * static_cast<size_t>(n) * sizeof(double),
                                      cudaMemcpyHostToDevice),
                           "HODLR GPU points");
            }
            check_cuda(cudaMemcpy(ids_, ids, static_cast<size_t>(n) * sizeof(int64_t), cudaMemcpyHostToDevice),
                       "HODLR GPU point ids");
        }
        n_points_ = n;
    }
    bool has_points() const { return ids_ != nullptr; }
    int dim() const { return dim_; }
    int64_t n_points() const { return n_points_; }

    // out[b] (m[b] x n[b], leading dimension m[b], host) = scale * A(rows
    // r0[b] .., columns c0[b] ..), for the dense leaves
    void eval_dense(const fmm::gpu::Evaluator& evaluator, double scale, int count, const int64_t* r0, const int* m,
                    const int64_t* c0, const int* n, T* const* out) {
        auto t0 = std::chrono::steady_clock::now();
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        DeviceHeap& heap = DeviceHeap::instance();
        int i = 0;
        while (i < count) {
            // a chunk of at most kChunkBytes (at least one block)
            int j = i;
            size_t bytes = 0;
            while (j < count) {
                const size_t b = aligned(static_cast<size_t>(m[j]) * n[j] * sizeof(T));
                if (j > i && bytes + b > kChunkBytes) break;
                bytes += b;
                ++j;
            }
            char* buf = heap.alloc<char>(bytes);
            std::vector<T*> dst;
            ItemList<fmm::gpu::EvalItemT<T>> evals;
            ItemList<AxpbyItem<T>> scales;
            size_t at = 0;
            for (int b = i; b < j; ++b) {
                T* d = reinterpret_cast<T*>(buf + at);
                at += aligned(static_cast<size_t>(m[b]) * n[b] * sizeof(T));
                dst.push_back(d);
                check_point_range(r0[b], m[b]);
                check_point_range(c0[b], n[b]);
                evals.add({d, m[b], m[b], n[b], contiguous(r0[b]), contiguous(c0[b])}, m[b], n[b]);
                if (scale != 1.0) scales.add({d, m[b], d, m[b], nullptr, 0, m[b], n[b], T(scale), T(0.0)}, m[b], n[b]);
            }
            meta_.clear();
            evals.stage(meta_);
            scales.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            evaluator.eval<T>(evals.items, evals.device(md), evals.max_m, evals.max_n, md, points(), stream);
            launch_axpby(scales.device(md), scales.count(), scales.max_m, scales.max_n, stream);
            check_cuda(cudaStreamSynchronize(stream), "HODLR GPU dense leaves");
            for (int b = i; b < j; ++b) {
                copy_down(out[b], dst[static_cast<size_t>(b - i)], static_cast<size_t>(m[b]) * n[b]);
            }
            heap.free(buf);
            i = j;
        }
        stats_.t_eval += seconds_since(t0);
    }

    // ---- BACA without overlap on the blocks of one level ----

    // Blocks b = 0 .. nb-1: rows r0[b] .. r0[b]+m[b]-1, columns c0[b] ..
    // c0[b]+n[b]-1 (point slots), panel width r_est[b] (min(BACA_Batch, m, n))
    void begin(const fmm::gpu::Evaluator& evaluator, double scale, int variant, int nb, const int64_t* r0, const int* m,
               const int64_t* c0, const int* n, const int* r_est) {
        end();
        Context::instance().activate();
        if (variant != 4 && variant != 5) throw std::invalid_argument("HODLR GPU BACA: variant must be 4 or 5");
        evaluator_ = &evaluator;
        scale_ = scale;
        variant_ = variant;
        blocks_.assign(static_cast<size_t>(nb), Block{});
        size_t sel_count = 0;
        for (int b = 0; b < nb; ++b) {
            Block& B = blocks_[static_cast<size_t>(b)];
            B.r0 = r0[b];
            B.c0 = c0[b];
            B.m = m[b];
            B.n = n[b];
            B.r = r_est[b];
            if (B.m <= 0 || B.n <= 0 || B.r <= 0 || B.r > std::min(B.m, B.n)) {
                throw std::invalid_argument("HODLR GPU BACA: invalid block or panel size");
            }
            check_point_range(B.r0, B.m);
            check_point_range(B.c0, B.n);
            B.sel_offset = sel_count;
            sel_count += static_cast<size_t>(B.r);
        }
        DeviceHeap& heap = DeviceHeap::instance();
        sel_rows_d_ = heap.alloc_resident<int>(std::max<size_t>(sel_count, 1) * sizeof(int));
        sel_cols_d_ = heap.alloc_resident<int>(std::max<size_t>(sel_count, 1) * sizeof(int));
        sel_j2_d_ = heap.alloc_resident<int>(std::max<size_t>(sel_count, 1) * sizeof(int));
        sel_host_.assign(sel_count, 0);
        for (Block& B : blocks_) {
            const size_t mn = static_cast<size_t>(std::max(B.m, B.n));
            const size_t bytes = aligned(static_cast<size_t>(B.m) * B.r * sizeof(T)) +
                                 aligned(static_cast<size_t>(B.r) * B.n * sizeof(T)) +
                                 aligned(static_cast<size_t>(B.r) * mn * sizeof(T)) + aligned(3 * mn * sizeof(double)) +
                                 aligned(mn * sizeof(int));
            B.store = heap.alloc_resident<char>(bytes);
            char* at = B.store;
            auto take = [&](size_t b) {
                char* p = at;
                at += aligned(b);
                return p;
            };
            B.c = reinterpret_cast<T*>(take(static_cast<size_t>(B.m) * B.r * sizeof(T)));
            B.rp = reinterpret_cast<T*>(take(static_cast<size_t>(B.r) * B.n * sizeof(T)));
            B.work = reinterpret_cast<T*>(take(static_cast<size_t>(B.r) * mn * sizeof(T)));
            B.norms = reinterpret_cast<double*>(take(3 * mn * sizeof(double)));
            B.perm = reinterpret_cast<int*>(take(mn * sizeof(int)));
            if (at != B.store + bytes) throw std::logic_error("HODLR GPU BACA: panel storage overflow");
            B.sel_rows_d = sel_rows_d_ + B.sel_offset;
            B.sel_cols_d = sel_cols_d_ + B.sel_offset;
            B.j2_d = sel_j2_d_ + B.sel_offset;
            B.sel_rows.assign(static_cast<size_t>(B.r), 0);
            B.sel_cols.assign(static_cast<size_t>(B.r), 0);
            grow(B, B.r);
        }
    }

    int rank(int b) const { return block(b).rank; }

    // The panel columns of block b for the next panels() (1-based; r_est of them)
    void set_columns(int b, const int* cols1) {
        Block& B = block(b);
        for (int j = 0; j < B.r; ++j) {
            const int c = cols1[j] - 1;
            if (c < 0 || c >= B.n) throw std::out_of_range("HODLR GPU BACA: column index out of range");
            B.sel_cols[static_cast<size_t>(j)] = c;
        }
    }

    // Panels of the listed blocks from their current columns; core_out gets
    // the r x r cores C(I, :), packed in list order.  BACA (variant 4, LR_BACA)
    // takes I without masking the chosen rows, then J2 = the pivots of R,
    // the core from C = A(:, J2) - U V(:, J2), and the next columns from R
    // with J2 masked.
    void panels(int na, const int* bl, T* core_out) {
        const bool baca = variant_ == 4;
        auto t0 = std::chrono::steady_clock::now();
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        DeviceHeap& heap = DeviceHeap::instance();
        size_t gbytes = 0, cores = 0;
        for (int a = 0; a < na; ++a) {
            const Block& B = block(bl[a]);
            gbytes += 3 * aligned(static_cast<size_t>(B.r) * std::max(B.rank, 1) * sizeof(T));
            cores += static_cast<size_t>(B.r) * B.r;
        }
        char* gbuf = heap.alloc<char>(std::max<size_t>(gbytes, 1));
        T* core = heap.alloc<T>(std::max<size_t>(cores, 1) * sizeof(T));
        meta_.clear();
        ItemList<fmm::gpu::EvalItemT<T>> evc, evr;
        ItemList<AxpbyItem<T>> scc, scr;
        ItemList<fmm::gpu::GatherItemT<T>> gc, gr, gcore, gc2;
        PivotList piv, piv_j2, piv_next;
        ItemList<fmm::gpu::EvalItemT<T>> evc2;
        VBatch<T> bc, br, bc2;
        size_t at = 0, core_at = 0;
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            const int r = B.r, k = B.rank;
            const size_t off_cols = meta_.append(B.sel_cols);
            const size_t off_rows = meta_.append(B.rows);
            evc.add({B.c, B.m, B.m, r, contiguous(B.r0), listed(off_cols, B.c0)}, B.m, r);
            if (scale_ != 1.0) scc.add({B.c, B.m, B.c, B.m, nullptr, 0, B.m, r, T(scale_), T(0.0)}, B.m, r);
            T* gcb = reinterpret_cast<T*>(gbuf + at);
            at += aligned(static_cast<size_t>(r) * std::max(k, 1) * sizeof(T));
            T* grb = reinterpret_cast<T*>(gbuf + at);
            at += aligned(static_cast<size_t>(r) * std::max(k, 1) * sizeof(T));
            if (k > 0) {
                // V(:, J)^T = V^T(J, :) (r x k); C -= U V(:, J)
                gc.add({gcb, r, r, k, B.vt, 1, B.n, listed(off_cols, 0), contiguous(0)}, r, k);
                bc.entries.push_back({B.u, gcb, B.c, B.m, r, k, B.m, r, B.m});
                // U(I, :) (r x k); R -= U(I, :) V
                gr.add({grb, r, r, k, B.u, 1, B.m, device_list(B.sel_rows_d), contiguous(0)}, r, k);
                br.entries.push_back({grb, B.vt, B.rp, r, B.n, k, r, B.n, r});
                stats_.flops += kFlopScale * 2.0 * (static_cast<double>(B.m) + B.n) * r * k;
            }
            piv.add({B.c, B.m, 1, r, B.m, r, baca ? 0 : k, static_cast<int64_t>(off_rows), nullptr, B.work, B.norms,
                     B.perm, B.sel_rows_d});
            evr.add({B.rp, r, r, B.n, fmm::gpu::IndexList{-1, static_cast<int>(B.r0), B.sel_rows_d},
                     contiguous(B.c0)},
                    r, B.n);
            if (scale_ != 1.0) scr.add({B.rp, r, B.rp, r, nullptr, 0, r, B.n, T(scale_), T(0.0)}, r, B.n);
            if (baca) {
                // J2 = the pivots of R; C = A(:, J2) - U V(:, J2); the next columns: pivots of R with J2 masked
                piv_j2.add({B.rp, r, 0, r, B.n, r, 0, 0, nullptr, B.work, B.norms, B.perm, B.j2_d});
                evc2.add({B.c, B.m, B.m, r, contiguous(B.r0), fmm::gpu::IndexList{-1, static_cast<int>(B.c0), B.j2_d}},
                         B.m, r);
                T* g2 = reinterpret_cast<T*>(gbuf + at);
                at += aligned(static_cast<size_t>(r) * std::max(k, 1) * sizeof(T));
                if (k > 0) {
                    gc2.add({g2, r, r, k, B.vt, 1, B.n, device_list(B.j2_d), contiguous(0)}, r, k);
                    bc2.entries.push_back({B.u, g2, B.c, B.m, r, k, B.m, r, B.m});
                    stats_.flops += kFlopScale * 2.0 * static_cast<double>(B.m) * r * k;
                }
                piv_next.add({B.rp, r, 0, r, B.n, r, r, 0, B.j2_d, B.work, B.norms, B.perm, B.sel_cols_d});
                stats_.flops += 2.0 * qrcp_flops(r, B.n);  // (the pivots of R, twice)
            }
            gcore.add({core + core_at, r, r, r, B.c, 1, B.m, device_list(B.sel_rows_d), contiguous(0)}, r, r);
            core_at += static_cast<size_t>(r) * r;
            stats_.flops += qrcp_flops(r, B.m);  // (the pivots of C)
        }
        evc.stage(meta_);
        scc.stage(meta_);
        gc.stage(meta_);
        bc.stage(meta_);
        piv.stage(meta_);
        evr.stage(meta_);
        scr.stage(meta_);
        gr.stage(meta_);
        br.stage(meta_);
        piv_j2.stage(meta_);
        evc2.stage(meta_);
        gc2.stage(meta_);
        bc2.stage(meta_);
        piv_next.stage(meta_);
        gcore.stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        const fmm::gpu::PointTable pts = points();
        entry_begin(stream);
        evaluator_->eval<T>(evc.items, evc.device(md), evc.max_m, evc.max_n, md, pts, stream);
        launch_axpby(scc.device(md), scc.count(), scc.max_m, scc.max_n, stream);
        entry_end(stream);
        fmm::gpu::launch_gather(gc.device(md), gc.count(), gc.max_m, gc.max_n, md, stream);
        bc.gemm(md, MagmaNoTrans, MagmaTrans, T(-1.0), T(1.0), queue);
        piv.launch(md, stream);
        entry_begin(stream);
        evaluator_->eval<T>(evr.items, evr.device(md), evr.max_m, evr.max_n, md, pts, stream);
        launch_axpby(scr.device(md), scr.count(), scr.max_m, scr.max_n, stream);
        entry_end(stream);
        fmm::gpu::launch_gather(gr.device(md), gr.count(), gr.max_m, gr.max_n, md, stream);
        br.gemm(md, MagmaNoTrans, MagmaTrans, T(-1.0), T(1.0), queue);
        if (baca) {
            piv_j2.launch(md, stream);
            entry_begin(stream);
            evaluator_->eval<T>(evc2.items, evc2.device(md), evc2.max_m, evc2.max_n, md, pts, stream);
            launch_axpby(scc.device(md), scc.count(), scc.max_m, scc.max_n, stream);
            entry_end(stream);
            fmm::gpu::launch_gather(gc2.device(md), gc2.count(), gc2.max_m, gc2.max_n, md, stream);
            bc2.gemm(md, MagmaNoTrans, MagmaTrans, T(-1.0), T(1.0), queue);
            piv_next.launch(md, stream);
        }
        fmm::gpu::launch_gather(gcore.device(md), gcore.count(), gcore.max_m, gcore.max_n, md, stream);
        if (cores > 0) {
            check_cuda(cudaMemcpyAsync(core_out, core, cores * sizeof(T), cudaMemcpyDeviceToHost, stream),
                       "HODLR GPU BACA cores");
        }
        download_selection(baca ? sel_cols_d_ : sel_rows_d_, stream);
        check_cuda(cudaStreamSynchronize(stream), "HODLR GPU BACA panels");
        entry_collect();
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            std::copy_n(sel_host_.begin() + static_cast<std::ptrdiff_t>(B.sel_offset), B.r,
                        baca ? B.sel_cols.begin() : B.sel_rows.begin());
        }
        heap.free(core);
        heap.free(gbuf);
        stats_.t_panels += seconds_since(t0);
    }

    // The panels of the k-nearest-neighbour step: C = A(:, cols) (m x nc),
    // R = A(rows, :) (nr x n), core = C(rows, :) (nr x nc) to core_out,
    // packed in list order; cols, rows 1-based, packed.
    void knn_panels(int na, const int* bl, const int* nc, const int* cols1, const int* nr, const int* rows1,
                    T* core_out) {
        auto t0 = std::chrono::steady_clock::now();
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        DeviceHeap& heap = DeviceHeap::instance();
        size_t cores = 0, ci = 0, ri = 0;
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            if (B.rank != 0 || B.knn_store != nullptr) throw std::logic_error("HODLR GPU BACA: knn step out of order");
            if (nc[a] <= 0 || nr[a] <= 0 || nc[a] > B.n || nr[a] > B.m) {
                throw std::invalid_argument("HODLR GPU BACA: invalid knn panel size");
            }
            B.knn_cols.resize(static_cast<size_t>(nc[a]));
            B.knn_rows.resize(static_cast<size_t>(nr[a]));
            for (int j = 0; j < nc[a]; ++j) B.knn_cols[static_cast<size_t>(j)] = checked(cols1[ci + j] - 1, B.n);
            for (int i = 0; i < nr[a]; ++i) B.knn_rows[static_cast<size_t>(i)] = checked(rows1[ri + i] - 1, B.m);
            ci += static_cast<size_t>(nc[a]);
            ri += static_cast<size_t>(nr[a]);
            cores += static_cast<size_t>(nr[a]) * nc[a];
            const size_t kc = static_cast<size_t>(nc[a]), kr = static_cast<size_t>(nr[a]), n = static_cast<size_t>(B.n);
            const size_t bytes = aligned(B.m * kc * sizeof(T)) + aligned(kr * n * sizeof(T)) + aligned(kr * n * sizeof(T)) +
                                 aligned(3 * n * sizeof(double)) + aligned(n * sizeof(int));
            B.knn_store = heap.alloc<char>(bytes);
            char* at = B.knn_store;
            auto take = [&](size_t b) {
                char* p = at;
                at += aligned(b);
                return p;
            };
            B.ck = reinterpret_cast<T*>(take(B.m * kc * sizeof(T)));
            B.rk = reinterpret_cast<T*>(take(kr * n * sizeof(T)));
            B.kwork = reinterpret_cast<T*>(take(kr * n * sizeof(T)));
            B.knorms = reinterpret_cast<double*>(take(3 * n * sizeof(double)));
            B.kperm = reinterpret_cast<int*>(take(n * sizeof(int)));
        }
        T* core = heap.alloc<T>(std::max<size_t>(cores, 1) * sizeof(T));
        meta_.clear();
        ItemList<fmm::gpu::EvalItemT<T>> ev;
        ItemList<AxpbyItem<T>> sc;
        ItemList<fmm::gpu::GatherItemT<T>> gcore;
        size_t core_at = 0;
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            const int kc = static_cast<int>(B.knn_cols.size()), kr = static_cast<int>(B.knn_rows.size());
            const size_t off_c = meta_.append(B.knn_cols);
            const size_t off_r = meta_.append(B.knn_rows);
            ev.add({B.ck, B.m, B.m, kc, contiguous(B.r0), listed(off_c, B.c0)}, B.m, kc);
            ev.add({B.rk, kr, kr, B.n, listed(off_r, B.r0), contiguous(B.c0)}, kr, B.n);
            if (scale_ != 1.0) {
                sc.add({B.ck, B.m, B.ck, B.m, nullptr, 0, B.m, kc, T(scale_), T(0.0)}, B.m, kc);
                sc.add({B.rk, kr, B.rk, kr, nullptr, 0, kr, B.n, T(scale_), T(0.0)}, kr, B.n);
            }
            gcore.add({core + core_at, kr, kr, kc, B.ck, 1, B.m, listed(off_r, 0), contiguous(0)}, kr, kc);
            core_at += static_cast<size_t>(kr) * kc;
        }
        ev.stage(meta_);
        sc.stage(meta_);
        gcore.stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        entry_begin(stream);
        evaluator_->eval<T>(ev.items, ev.device(md), ev.max_m, ev.max_n, md, points(), stream);
        launch_axpby(sc.device(md), sc.count(), sc.max_m, sc.max_n, stream);
        entry_end(stream);
        fmm::gpu::launch_gather(gcore.device(md), gcore.count(), gcore.max_m, gcore.max_n, md, stream);
        if (cores > 0) {
            check_cuda(cudaMemcpyAsync(core_out, core, cores * sizeof(T), cudaMemcpyDeviceToHost, stream),
                       "HODLR GPU BACA knn cores");
        }
        check_cuda(cudaStreamSynchronize(stream), "HODLR GPU BACA knn panels");
        entry_collect();
        heap.free(core);
        stats_.t_panels += seconds_since(t0);
    }

    // Append ru[a] new rank-one terms to each listed block: U(:, new) =
    // C(:, jpvt) and V(new, :) = W R, with jpvt (1-based, packed) and W
    // (ru x p, packed, p = r_est, or the knn row count when knn != 0) from
    // the host.  The chosen rows and columns grow by I(1:ru) and J(jpvt).
    // Without knn, grams_out gets per block (packed) the Gram matrices of the
    // norm update (LR_Fnorm, LR_FnormUp): U_new^H U_new, V_new^* V_new^T
    // (ru x ru each), and U_old^H U_new, V_old^* V_new^T with U_old = U(:,
    // rskip+1:rank), V_old = V(rskip+1:rank, :) ((rank - rskip) x ru each,
    // rank before the append).  Then the next panel columns: the pivots of R
    // with the chosen columns masked.
    void append(int knn, int na, const int* bl, const int* ru, const int* jpvt1, const T* w, const int* rskip,
                T* grams_out) {
        const bool baca = variant_ == 4;  // (the next columns come from panels())
        if (baca && knn) throw std::logic_error("HODLR GPU BACA: no knn step in BACA");
        auto t0 = std::chrono::steady_clock::now();
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        DeviceHeap& heap = DeviceHeap::instance();
        // sizes, the W upload and the capacity of U and V^T
        size_t wcount = 0, gcount = 0;
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            const int p = knn ? static_cast<int>(B.knn_rows.size()) : B.r;
            if (knn && B.knn_store == nullptr) throw std::logic_error("HODLR GPU BACA: knn append without knn panels");
            if (ru[a] <= 0 || ru[a] > p || B.rank + ru[a] > std::min(B.m, B.n)) {
                throw std::invalid_argument("HODLR GPU BACA: invalid rank update");
            }
            wcount += static_cast<size_t>(ru[a]) * p;
            if (!knn) {
                const size_t nold = static_cast<size_t>(B.rank - rskip[a]);
                gcount += 2 * static_cast<size_t>(ru[a]) * ru[a] + 2 * nold * ru[a];
            }
            grow(B, B.rank + ru[a]);
        }
        T* wd = heap.alloc<T>(std::max<size_t>(wcount, 1) * sizeof(T));
        T* gd = heap.alloc<T>(std::max<size_t>(gcount, 1) * sizeof(T));
        check_cuda(cudaMemsetAsync(gd, 0, std::max<size_t>(gcount, 1) * sizeof(T), stream), "HODLR GPU BACA Grams");
        if (wcount > 0) {
            check_cuda(cudaMemcpyAsync(wd, w, wcount * sizeof(T), cudaMemcpyHostToDevice, stream), "HODLR GPU BACA W");
        }
        meta_.clear();
        ItemList<fmm::gpu::GatherItemT<T>> gu;
        PivotList piv;
        VBatch<T> vnew;
        std::vector<GramJob> gconj;
        size_t w_at = 0, g_at = 0, j_at = 0;
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            const int k = B.rank, u = ru[a];
            const int p = knn ? static_cast<int>(B.knn_rows.size()) : B.r;
            const int pc = knn ? static_cast<int>(B.knn_cols.size()) : B.r;
            T* cpanel = knn ? B.ck : B.c;
            T* rpanel = knn ? B.rk : B.rp;
            std::vector<int> jp(static_cast<size_t>(u));
            for (int j = 0; j < u; ++j) jp[static_cast<size_t>(j)] = checked(jpvt1[j_at + j] - 1, pc);
            j_at += static_cast<size_t>(u);
            const size_t off_jp = meta_.append(jp);
            T* unew = B.u + static_cast<size_t>(k) * B.m;
            T* vnew_t = B.vt + static_cast<size_t>(k) * B.n;
            gu.add({unew, B.m, B.m, u, cpanel, 1, B.m, contiguous(0), listed(off_jp, 0)}, B.m, u);
            // V_new^T = R^T W^T (n x u)
            vnew.entries.push_back({rpanel, wd + w_at, vnew_t, B.n, u, p, p, u, B.n});
            w_at += static_cast<size_t>(u) * p;
            stats_.flops += kFlopScale * 2.0 * B.n * static_cast<double>(u) * p;
            if (!knn) {
                const int nold = k - rskip[a];
                gconj.push_back({unew, B.m, unew, B.m, u, u, B.m, gd + g_at});
                g_at += static_cast<size_t>(u) * u;
                gconj.push_back({vnew_t, B.n, vnew_t, B.n, u, u, B.n, gd + g_at});
                g_at += static_cast<size_t>(u) * u;
                if (nold > 0) {
                    gconj.push_back({B.u + static_cast<size_t>(rskip[a]) * B.m, B.m, unew, B.m, nold, u, B.m, gd + g_at});
                    g_at += static_cast<size_t>(nold) * u;
                    gconj.push_back(
                        {B.vt + static_cast<size_t>(rskip[a]) * B.n, B.n, vnew_t, B.n, nold, u, B.n, gd + g_at});
                    g_at += static_cast<size_t>(nold) * u;
                }
                stats_.flops += kFlopScale * 2.0 * (B.m + B.n) * static_cast<double>(u) * (u + nold);
            }
            // the chosen rows and columns (masks of BACA without overlap)
            if (!baca) {
                const std::vector<int>& rows_from = knn ? B.knn_rows : B.sel_rows;
                const std::vector<int>& cols_from = knn ? B.knn_cols : B.sel_cols;
                for (int j = 0; j < u; ++j) {
                    B.rows.push_back(rows_from[static_cast<size_t>(j)]);
                    B.cols.push_back(cols_from[static_cast<size_t>(jp[static_cast<size_t>(j)])]);
                }
            }
            B.rank = k + u;
            // next panel columns
            const size_t off_mask = meta_.append(B.cols);
            if (!baca) {
                piv.add({rpanel, p, 0, p, B.n, B.r, B.rank, static_cast<int64_t>(off_mask), nullptr,
                         knn ? B.kwork : B.work, knn ? B.knorms : B.norms, knn ? B.kperm : B.perm, B.sel_cols_d});
            }
            if (!baca) stats_.flops += kFlopScale * 4.0 * static_cast<double>(p) * std::min(p, B.r) * B.n;
        }
        gu.stage(meta_);
        vnew.stage(meta_);
        piv.stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        fmm::gpu::launch_gather(gu.device(md), gu.count(), gu.max_m, gu.max_n, md, stream);
        vnew.gemm(md, MagmaTrans, MagmaTrans, T(1.0), T(0.0), queue);
        piv.launch(md, stream);
        run_grams(gconj, MagmaConjTrans);
        if (gcount > 0) {
            check_cuda(cudaMemcpyAsync(grams_out, gd, gcount * sizeof(T), cudaMemcpyDeviceToHost, stream),
                       "HODLR GPU BACA Gram matrices");
        }
        if (!baca) download_selection(sel_cols_d_, stream);
        check_cuda(cudaStreamSynchronize(stream), "HODLR GPU BACA append");
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            if (!baca) {
                std::copy_n(sel_host_.begin() + static_cast<std::ptrdiff_t>(B.sel_offset), B.r, B.sel_cols.begin());
            }
            if (knn) {
                heap.free(B.knn_store);
                B.knn_store = nullptr;
                B.ck = B.rk = B.kwork = nullptr;
                B.knorms = nullptr;
                B.kperm = nullptr;
            }
        }
        heap.free(gd);
        heap.free(wd);
        stats_.t_append += seconds_since(t0);
    }

    // The QR factorizations of U (m x rank) and V^T (n x rank) of the listed
    // blocks, in place (Householder vectors below the diagonal, unit lower
    // trapezoidal after this call).  Per block, packed in list order, U's
    // part before V's: top the leading rank x rank part of the ?geqrf output
    // (R and the top of the vectors), tau the rank scalar factors, gram
    // Y^H Y (rank x rank) of the unit lower vectors Y, left on the device
    struct QrOut {
        T* top = nullptr;   // per block (list order) at top_at: U's rank x rank top, then V's
        T* tau = nullptr;   // per block at tau_at: U's rank taus, then V's
        T* gram = nullptr;  // per block at top_at: U's Y^H Y, then V's
        std::vector<size_t> top_at, tau_at;
        size_t tops = 0, taus = 0;
        void release() {
            DeviceHeap& heap = DeviceHeap::instance();
            if (gram != nullptr) heap.free(gram);
            if (tau != nullptr) heap.free(tau);
            if (top != nullptr) heap.free(top);
            top = tau = gram = nullptr;
        }
    };
    QrOut qr_device(int na, const int* bl) {
        auto t0 = std::chrono::steady_clock::now();
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        DeviceHeap& heap = DeviceHeap::instance();
        size_t tops = 0, taus = 0;
        // BPACK_TRACE=hodlr-qr prints the blocks and chunks; BPACK_TRACE=sync
        // synchronizes after each step (to name the kernel of an asynchronous fault)
        static const bool debug = fmm::env::trace("hodlr-qr");
        static const bool sync = fmm::env::trace("sync");
        auto step = [&](const char* what) {
            if (sync) check_cuda(cudaStreamSynchronize(stream), what);
        };
        step("HODLR GPU BACA QR: before (a fault of an earlier step)");
        for (int a = 0; a < na; ++a) {
            const Block& B = block(bl[a]);
            if (debug) {
                std::printf(" HODLR GPU BACA qr: block %d of %d: m %d n %d rank %d cap %d r %d u %p vt %p\n", bl[a], na, B.m,
                            B.n, B.rank, B.cap, B.r, static_cast<const void*>(B.u), static_cast<const void*>(B.vt));
                std::fflush(stdout);
            }
            if (B.rank <= 0) throw std::logic_error("HODLR GPU BACA: recompression of a block of rank 0");
            if (B.rank > std::min(B.m, B.n) || B.rank > B.cap || B.u == nullptr || B.vt == nullptr) {
                throw std::logic_error("HODLR GPU BACA: a block's rank exceeds its size or its storage");
            }
            tops += 2 * static_cast<size_t>(B.rank) * B.rank;
            taus += 2 * static_cast<size_t>(B.rank);
        }
        T* topd = heap.alloc<T>(std::max<size_t>(tops, 1) * sizeof(T));
        T* taud = heap.alloc<T>(std::max<size_t>(taus, 1) * sizeof(T));
        T* gramd = heap.alloc<T>(std::max<size_t>(tops, 1) * sizeof(T));
        std::vector<size_t> tau_at(static_cast<size_t>(na)), top_at(static_cast<size_t>(na));
        {
            size_t t = 0, q = 0;
            for (int a = 0; a < na; ++a) {
                const Block& B = block(bl[a]);
                tau_at[static_cast<size_t>(a)] = t;
                top_at[static_cast<size_t>(a)] = q;
                t += 2 * static_cast<size_t>(B.rank);
                q += 2 * static_cast<size_t>(B.rank) * B.rank;
            }
        }
        // the blocks of large rank one by one (cuSOLVER geqrf in place), the
        // others (tall ones included) by MAGMA's batched geqrf of zero-padded
        // copies, in chunks
        std::vector<int> small;
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            if (B.rank < kBigQrRank) {
                small.push_back(a);
                continue;
            }
            T* td = taud + tau_at[static_cast<size_t>(a)];
            geqrf_one(B.m, B.rank, B.u, B.m, td);
            geqrf_one(B.n, B.rank, B.vt, B.n, td + B.rank);
            step("HODLR GPU BACA QR: cuSOLVER geqrf");
            stats_.flops += kFlopScale * 2.0 * (static_cast<double>(B.m) + B.n) * B.rank * B.rank;
        }
        const int ns = static_cast<int>(small.size());
        int a0 = 0;
        while (a0 < ns) {
            int a1 = a0, mp = 0, rp = 0;
            while (a1 < ns) {
                const Block& B = block(bl[small[static_cast<size_t>(a1)]]);
                const int mp2 = std::max(mp, std::max(B.m, B.n)), rp2 = std::max(rp, B.rank);
                const size_t bytes = 2 * static_cast<size_t>(a1 - a0 + 1) * mp2 * rp2 * sizeof(T);
                if (a1 > a0 && bytes > kChunkBytes) break;
                mp = mp2;
                rp = rp2;
                ++a1;
            }
            const int count = 2 * (a1 - a0);
            const size_t each = static_cast<size_t>(mp) * rp;
            T* pad = heap.alloc<T>(each * count * sizeof(T));
            T* ptau = heap.alloc<T>(static_cast<size_t>(rp) * count * sizeof(T));
            magma_int_t* info = heap.alloc<magma_int_t>(static_cast<size_t>(count) * sizeof(magma_int_t));
            check_cuda(cudaMemsetAsync(pad, 0, each * count * sizeof(T), stream), "HODLR GPU BACA QR padding");
            meta_.clear();
            ItemList<AxpbyItem<T>> in, out, tau_copy;
            std::vector<T*> aptr, tptr;
            for (int q = a0; q < a1; ++q) {
                const int a = small[static_cast<size_t>(q)];
                Block& B = block(bl[a]);
                const int k = B.rank;
                const size_t i = 2 * static_cast<size_t>(q - a0);
                T* pu = pad + i * each;
                T* pv = pad + (i + 1) * each;
                in.add({pu, mp, B.u, B.m, nullptr, 0, B.m, k, T(1.0), T(0.0)}, B.m, k);
                in.add({pv, mp, B.vt, B.n, nullptr, 0, B.n, k, T(1.0), T(0.0)}, B.n, k);
                out.add({B.u, B.m, pu, mp, nullptr, 0, B.m, k, T(1.0), T(0.0)}, B.m, k);
                out.add({B.vt, B.n, pv, mp, nullptr, 0, B.n, k, T(1.0), T(0.0)}, B.n, k);
                T* td = taud + tau_at[static_cast<size_t>(a)];
                tau_copy.add({td, k, ptau + i * rp, rp, nullptr, 0, k, 1, T(1.0), T(0.0)}, k, 1);
                tau_copy.add({td + k, k, ptau + (i + 1) * rp, rp, nullptr, 0, k, 1, T(1.0), T(0.0)}, k, 1);
                aptr.push_back(pu);
                aptr.push_back(pv);
                tptr.push_back(ptau + i * rp);
                tptr.push_back(ptau + (i + 1) * rp);
                stats_.flops += kFlopScale * 2.0 * (static_cast<double>(B.m) + B.n) * k * k;
            }
            in.stage(meta_);
            out.stage(meta_);
            tau_copy.stage(meta_);
            const size_t off_a = meta_.append(aptr);
            const size_t off_t = meta_.append(tptr);
            char* md = meta_.upload(meta_device_, stream);
            if (debug) {
                std::printf(" HODLR GPU BACA qr chunk: %d matrices of %d x %d (lda %d), %zu bytes padded:", count, mp, rp,
                            mp, each * count * sizeof(T));
                for (int q = a0; q < a1; ++q) {
                    const Block& B = block(bl[small[static_cast<size_t>(q)]]);
                    std::printf(" (%d, %d, %d)", B.m, B.n, B.rank);
                }
                std::printf("\n");
                std::fflush(stdout);
            }
            step("HODLR GPU BACA QR: metadata upload");
            launch_axpby(in.device(md), in.count(), in.max_m, in.max_n, stream);
            step("HODLR GPU BACA QR: padded copies in");
            geqrf_batched(mp, rp, reinterpret_cast<T**>(md + off_a), mp, reinterpret_cast<T**>(md + off_t), info,
                          count, queue);
            step("HODLR GPU BACA QR: MAGMA geqrf_batched");
            launch_axpby(out.device(md), out.count(), out.max_m, out.max_n, stream);
            launch_axpby(tau_copy.device(md), tau_copy.count(), tau_copy.max_m, tau_copy.max_n, stream);
            step("HODLR GPU BACA QR: copies out");
            std::vector<magma_int_t> info_h(static_cast<size_t>(count));
            check_cuda(cudaMemcpyAsync(info_h.data(), info, count * sizeof(magma_int_t), cudaMemcpyDeviceToHost,
                                       stream),
                       "HODLR GPU BACA QR info");
            check_cuda(cudaStreamSynchronize(stream), "HODLR GPU BACA QR");
            for (magma_int_t v : info_h) {
                if (v != 0) throw std::runtime_error("HODLR GPU BACA: batched geqrf failed, info " + std::to_string(v));
            }
            heap.free(info);
            heap.free(ptau);
            heap.free(pad);
            a0 = a1;
        }
        // the leading parts, the unit lower vectors and their Gram matrices
        meta_.clear();
        ItemList<AxpbyItem<T>> top;
        ItemList<UnitLowerItem<T>> unit;
        std::vector<GramJob> gram;
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            const int k = B.rank;
            T* tu = topd + top_at[static_cast<size_t>(a)];
            T* tv = tu + static_cast<size_t>(k) * k;
            T* gu = gramd + top_at[static_cast<size_t>(a)];
            T* gv = gu + static_cast<size_t>(k) * k;
            top.add({tu, k, B.u, B.m, nullptr, 0, k, k, T(1.0), T(0.0)}, k, k);
            top.add({tv, k, B.vt, B.n, nullptr, 0, k, k, T(1.0), T(0.0)}, k, k);
            unit.add({B.u, B.m, k}, k, k);
            unit.add({B.vt, B.n, k}, k, k);
            gram.push_back({B.u, B.m, B.u, B.m, k, k, B.m, gu});
            gram.push_back({B.vt, B.n, B.vt, B.n, k, k, B.n, gv});
            stats_.flops += kFlopScale * (static_cast<double>(B.m) + B.n) * k * k;
        }
        top.stage(meta_);
        unit.stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        launch_axpby(top.device(md), top.count(), top.max_m, top.max_n, stream);
        launch_unit_lower(unit.device(md), unit.count(), unit.max_n, stream);
        check_cuda(cudaMemsetAsync(gramd, 0, std::max<size_t>(tops, 1) * sizeof(T), stream), "HODLR GPU BACA QR Grams");
        run_grams(gram, MagmaConjTrans);
        check_cuda(cudaStreamSynchronize(stream), "HODLR GPU BACA QR factors");
        stats_.t_recomp += seconds_since(t0);
        QrOut out;
        out.top = topd;
        out.tau = taud;
        out.gram = gramd;
        out.top_at = std::move(top_at);
        out.tau_at = std::move(tau_at);
        out.tops = tops;
        out.taus = taus;
        return out;
    }

    // After qr_device(): the new factors [X_U; 0] - Y_U Z_U (= Q_U [X_U; 0], m x rn)
    // and [X_V; 0] - Y_V Z_V (n x rn), X and Z (rank x rn each; U's then V's per
    // block, packed) on the device; they stay on the device for download()
    void finish_device(int na, const int* bl, const int* rn, const T* xd, const T* zd) {
        auto t0 = std::chrono::steady_clock::now();
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        DeviceHeap& heap = DeviceHeap::instance();
        release_out();
        size_t bytes = 0;
        for (int a = 0; a < na; ++a) {
            const Block& B = block(bl[a]);
            if (rn[a] <= 0 || rn[a] > B.rank) throw std::invalid_argument("HODLR GPU BACA: invalid recompressed rank");
            bytes += aligned(static_cast<size_t>(B.m) * rn[a] * sizeof(T)) +
                     aligned(static_cast<size_t>(B.n) * rn[a] * sizeof(T));
        }
        out_ = heap.alloc<char>(std::max<size_t>(bytes, 1));
        check_cuda(cudaMemsetAsync(out_, 0, std::max<size_t>(bytes, 1), stream), "HODLR GPU BACA factors");
        meta_.clear();
        ItemList<AxpbyItem<T>> put;
        VBatch<T> prod;
        out_u_.clear();
        out_v_.clear();
        size_t x_at = 0, o_at = 0;
        for (int a = 0; a < na; ++a) {
            Block& B = block(bl[a]);
            const int k = B.rank, q = rn[a];
            T* uo = reinterpret_cast<T*>(out_ + o_at);
            o_at += aligned(static_cast<size_t>(B.m) * q * sizeof(T));
            T* vo = reinterpret_cast<T*>(out_ + o_at);
            o_at += aligned(static_cast<size_t>(B.n) * q * sizeof(T));
            const T* xu = xd + x_at;
            const T* zu = zd + x_at;
            x_at += static_cast<size_t>(k) * q;
            const T* xv = xd + x_at;
            const T* zv = zd + x_at;
            x_at += static_cast<size_t>(k) * q;
            put.add({uo, B.m, xu, k, nullptr, 0, k, q, T(1.0), T(0.0)}, k, q);
            put.add({vo, B.n, xv, k, nullptr, 0, k, q, T(1.0), T(0.0)}, k, q);
            prod.entries.push_back({B.u, zu, uo, B.m, q, k, B.m, k, B.m});
            prod.entries.push_back({B.vt, zv, vo, B.n, q, k, B.n, k, B.n});
            out_u_.push_back(uo);
            out_v_.push_back(vo);
            stats_.flops += kFlopScale * 2.0 * (static_cast<double>(B.m) + B.n) * k * q;
        }
        put.stage(meta_);
        prod.stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        launch_axpby(put.device(md), put.count(), put.max_m, put.max_n, stream);
        prod.gemm(md, MagmaNoTrans, MagmaNoTrans, T(-1.0), T(1.0), queue);
        check_cuda(cudaStreamSynchronize(stream), "HODLR GPU BACA recompression");
        stats_.t_recomp += seconds_since(t0);
    }

    // The new factors of the last recompress() to the host
    // (u_out[a]: m x rn, v_out[a]: n x rn); with keep, their device copies go
    // to keep (the device HODLR's mirrors of these host arrays), else freed;
    // !fill (with keep): the host arrays are left unfilled, stale mirrors
    void download(int na, const int* bl, const int* rn, T* const* u_out, T* const* v_out, HodlrDevice<T>* keep = nullptr,
                  bool fill = true) {
        if (static_cast<size_t>(na) != out_u_.size()) throw std::logic_error("HODLR GPU BACA: download out of step");
        if (!fill && keep == nullptr) throw std::invalid_argument("HODLR GPU BACA: factors neither downloaded nor kept");
        Context::instance().activate();
        for (int a = 0; a < na; ++a) {
            const Block& B = block(bl[a]);
            const size_t cu = static_cast<size_t>(B.m) * rn[a], cv = static_cast<size_t>(B.n) * rn[a];
            if (fill) {
                copy_down(u_out[a], out_u_[static_cast<size_t>(a)], cu);
                copy_down(v_out[a], out_v_[static_cast<size_t>(a)], cv);
            }
            if (keep != nullptr) {
                keep->add_mirror(u_out[a], out_u_[static_cast<size_t>(a)], cu * sizeof(T), !fill);
                keep->add_mirror(v_out[a], out_v_[static_cast<size_t>(a)], cv * sizeof(T), !fill);
            }
        }
        if (keep != nullptr && out_ != nullptr) {
            keep->own_mirror_memory(out_);
            out_ = nullptr;
        }
        release_out();
    }

    // The new factors of the last recompress() as device matrices of q
    // (DistQr::dm_copy): ids[2a], ids[2a + 1] those of block a's U (m x rn)
    // and V (n x rn); the copies here are freed
    template<typename Q>
    void to_dm(int na, const int* bl, const int* rn, Q& q, int* ids) {
        if (static_cast<size_t>(na) != out_u_.size()) throw std::logic_error("HODLR GPU BACA: to_dm out of step");
        Context::instance().activate();
        for (int a = 0; a < na; ++a) {
            const Block& B = block(bl[a]);
            ids[2 * a] = q.dm_copy(B.m, rn[a], out_u_[static_cast<size_t>(a)]);
            ids[2 * a + 1] = q.dm_copy(B.n, rn[a], out_v_[static_cast<size_t>(a)]);
        }
        release_out();
    }

    // The recompression (LR_ReCompression) of the listed blocks after their
    // BACA, all on the device: the QR of U and V^T (qr_device), the SVD of
    // R1 R2^T (MAGMA's batched one-sided Jacobi for ranks below kBigSvdRank,
    // zero padded by size class; cuSOLVER gesvdp one by one for the others),
    // the new rank rn by the rule of SVD_Truncate (tol relative, underflow
    // absolute), and the new factors [X; 0] - Y Z (= Q [X; 0]) with X_U =
    // Us(:, 1:rn) diag(s), X_V = conj(Vs(:, 1:rn)) and Z = T Y1^H X, T^-1 =
    // striu(Y^H Y) + diag(1/tau) (the larft recurrence on the host for a
    // block with a zero tau).  rn_out: the new ranks; the factors stay for
    // download()
    void recompress(int na, const int* bl, double tol, double underflow, DeviceSvd<T>& svd, int* rn_out) {
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        DeviceHeap& heap = DeviceHeap::instance();
        QrOut q = qr_device(na, bl);
        auto t0 = std::chrono::steady_clock::now();
        const size_t n_a = static_cast<size_t>(std::max(na, 0));
        std::vector<int> ks(n_a), rn(n_a, 1);
        std::vector<size_t> sq_at(n_a), s_at(n_a);
        size_t sq = 0, sc = 0;
        for (size_t a = 0; a < n_a; ++a) {
            const int k = block(bl[a]).rank;
            ks[a] = k;
            sq_at[a] = sq;
            sq += static_cast<size_t>(k) * k;
            s_at[a] = sc;
            sc += static_cast<size_t>(k);
        }
        T* rb = heap.alloc<T>(std::max<size_t>(2 * sq, 1) * sizeof(T));  // R1, R2; then T^-1 of U, V
        T* mb = heap.alloc<T>(std::max<size_t>(sq, 1) * sizeof(T));      // R1 R2^T
        T* ub = heap.alloc<T>(std::max<size_t>(sq, 1) * sizeof(T));      // its U
        T* vb = heap.alloc<T>(std::max<size_t>(sq, 1) * sizeof(T));      // its V
        double* sb = heap.alloc<double>(std::max<size_t>(sc, 1) * sizeof(double));
        check_cuda(cudaMemsetAsync(sb, 0, std::max<size_t>(sc, 1) * sizeof(double), stream), "HODLR GPU recompression s");

        // R1 R2^T (plain transpose, as the host's gemm 'N', 'T')
        {
            meta_.clear();
            ItemList<AxpbyItem<T>> rc;
            ItemList<TriuItem<T>> tz;
            VBatch<T> mm;
            for (size_t a = 0; a < n_a; ++a) {
                const int k = ks[a];
                T* r1 = rb + 2 * sq_at[a];
                T* r2 = r1 + static_cast<size_t>(k) * k;
                const T* tu = q.top + q.top_at[a];
                const T* tv = tu + static_cast<size_t>(k) * k;
                rc.add({r1, k, tu, k, nullptr, 0, k, k, T(1.0), T(0.0)}, k, k);
                rc.add({r2, k, tv, k, nullptr, 0, k, k, T(1.0), T(0.0)}, k, k);
                tz.add({r1, k, k}, k, k);
                tz.add({r2, k, k}, k, k);
                mm.entries.push_back({r1, r2, mb + sq_at[a], k, k, k, k, k, k});
                stats_.flops += kFlopScale * 2.0 * static_cast<double>(k) * k * k;
            }
            rc.stage(meta_);
            tz.stage(meta_);
            mm.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            launch_axpby(rc.device(md), rc.count(), rc.max_m, rc.max_n, stream);
            launch_zero_lower(tz.device(md), tz.count(), tz.max_n, stream);
            mm.gemm(md, MagmaNoTrans, MagmaTrans, T(1.0), T(0.0), queue);
            check_cuda(cudaStreamSynchronize(stream), "HODLR GPU recompression R1 R2^T");
        }

        // the SVDs: the small ones batched by size class, the others one by one
        std::map<int, std::vector<size_t>> classes;
        std::vector<size_t> large;
        for (size_t a = 0; a < n_a; ++a) {
            if (ks[a] < kBigSvdRank) {
                classes[(ks[a] + 31) / 32 * 32].push_back(a);
            } else {
                large.push_back(a);
            }
        }
        for (auto& cls : classes) {
            const int K = cls.first;
            const std::vector<size_t>& list = cls.second;
            const size_t each = static_cast<size_t>(K) * K;
            const size_t per = std::max<size_t>(1, kChunkBytes / (3 * each * sizeof(T) + K * sizeof(double)));
            for (size_t c0 = 0; c0 < list.size(); c0 += per) {
                const size_t cnt = std::min(per, list.size() - c0);
                T* pa = heap.alloc<T>(each * cnt * sizeof(T));
                T* pu = heap.alloc<T>(each * cnt * sizeof(T));
                T* pv = heap.alloc<T>(each * cnt * sizeof(T));
                double* ps = heap.alloc<double>(static_cast<size_t>(K) * cnt * sizeof(double));
                magma_int_t* info = heap.alloc<magma_int_t>(cnt * sizeof(magma_int_t));
                check_cuda(cudaMemsetAsync(pa, 0, each * cnt * sizeof(T), stream), "HODLR GPU SVD padding");
                meta_.clear();
                ItemList<AxpbyItem<T>> in, out;
                std::vector<T*> aptr, uptr, vptr;
                std::vector<double*> sptr;
                for (size_t i = 0; i < cnt; ++i) {
                    const size_t a = list[c0 + i];
                    const int k = ks[a];
                    in.add({pa + i * each, K, mb + sq_at[a], k, nullptr, 0, k, k, T(1.0), T(0.0)}, k, k);
                    out.add({ub + sq_at[a], k, pu + i * each, K, nullptr, 0, k, k, T(1.0), T(0.0)}, k, k);
                    out.add({vb + sq_at[a], k, pv + i * each, K, nullptr, 0, k, k, T(1.0), T(0.0)}, k, k);
                    aptr.push_back(pa + i * each);
                    uptr.push_back(pu + i * each);
                    vptr.push_back(pv + i * each);
                    sptr.push_back(ps + i * K);
                    stats_.flops += svd_flops(k);
                }
                in.stage(meta_);
                out.stage(meta_);
                const size_t off_a = meta_.append(aptr), off_u = meta_.append(uptr), off_v = meta_.append(vptr);
                const size_t off_s = meta_.append(sptr);
                char* md = meta_.upload(meta_device_, stream);
                launch_axpby(in.device(md), in.count(), in.max_m, in.max_n, stream);
                gesvj_batched(K, reinterpret_cast<T**>(md + off_a), reinterpret_cast<double**>(md + off_s),
                              reinterpret_cast<T**>(md + off_u), reinterpret_cast<T**>(md + off_v), info,
                              static_cast<magma_int_t>(cnt), queue);
                launch_axpby(out.device(md), out.count(), out.max_m, out.max_n, stream);
                for (size_t i = 0; i < cnt; ++i) {
                    const size_t a = list[c0 + i];
                    check_cuda(cudaMemcpyAsync(sb + s_at[a], ps + i * K, ks[a] * sizeof(double), cudaMemcpyDeviceToDevice, stream),
                               "HODLR GPU SVD values");
                }
                std::vector<magma_int_t> info_h(cnt);
                check_cuda(cudaMemcpyAsync(info_h.data(), info, cnt * sizeof(magma_int_t), cudaMemcpyDeviceToHost, stream),
                           "HODLR GPU SVD info");
                check_cuda(cudaStreamSynchronize(stream), "HODLR GPU batched SVD");
                for (magma_int_t v : info_h) {
                    if (v < 0) throw std::runtime_error("HODLR GPU: batched SVD failed, info " + std::to_string(v));
                }
                heap.free(info);
                heap.free(ps);
                heap.free(pv);
                heap.free(pu);
                heap.free(pa);
            }
        }
        for (size_t a : large) {
            const int k = ks[a];
            svd.svd_device(k, k, mb + sq_at[a], sb + s_at[a], ub + sq_at[a], vb + sq_at[a]);
            stats_.flops += svd_flops(k);
        }

        // the new ranks (SVD_Truncate)
        std::vector<double> hs(std::max<size_t>(sc, 1));
        if (sc > 0) check_cuda(cudaMemcpy(hs.data(), sb, sc * sizeof(double), cudaMemcpyDeviceToHost), "HODLR GPU SVD values");
        std::vector<char> zero(n_a, 0);
        size_t xc = 0;
        for (size_t a = 0; a < n_a; ++a) {
            const double* s = hs.data() + s_at[a];
            const int k = ks[a];
            int r = k;
            if (s[0] < underflow) {
                r = 0;
            } else {
                const double thr = std::max(s[0] * tol, underflow);
                for (int i = 0; i < k; ++i) {
                    if (s[i] <= thr) {
                        r = i + 1;
                        if (s[i] < thr / 10) r = i;
                        break;
                    }
                }
            }
            zero[a] = r == 0;
            rn[a] = std::max(r, 1);
            rn_out[a] = rn[a];
            xc += 2 * static_cast<size_t>(k) * rn[a];
        }

        // X, then Y1^H X in place of Z, then Z = T (Y1^H X) by the triangular solve with T^-1
        T* xb = heap.alloc<T>(std::max<size_t>(2 * xc, 1) * sizeof(T));
        T* zb = xb + xc;
        check_cuda(cudaMemsetAsync(xb, 0, std::max<size_t>(xc, 1) * sizeof(T), stream), "HODLR GPU recompression X");
        std::vector<size_t> x_at(n_a);
        std::vector<char> host_t(n_a, 0);  // (a zero tau: T by the larft recurrence on the host)
        {
            std::vector<T> ht(std::max<size_t>(q.taus, 1));
            if (q.taus > 0) check_cuda(cudaMemcpy(ht.data(), q.tau, q.taus * sizeof(T), cudaMemcpyDeviceToHost), "HODLR GPU tau");
            for (size_t a = 0; a < n_a; ++a) {
                for (int i = 0; i < 2 * ks[a]; ++i) {
                    if (ht[q.tau_at[a] + i] == T(0.0)) host_t[a] = 1;
                }
            }
        }
        meta_.clear();
        ItemList<ColScaleItem<T>> xs;
        ItemList<TinvItem<T>> ti;
        VBatch<T> yp;
        std::vector<magma_int_t> tm, tn, tlda, tldb;
        std::vector<T*> tA, tB;
        size_t at = 0;
        for (size_t a = 0; a < n_a; ++a) {
            const Block& B = block(bl[a]);
            const int k = ks[a], r = rn[a];
            x_at[a] = at;
            T* xu = xb + at;
            T* xv = xu + static_cast<size_t>(k) * r;
            T* yu = zb + at;
            T* yv = yu + static_cast<size_t>(k) * r;
            at += 2 * static_cast<size_t>(k) * r;
            if (!zero[a]) {
                xs.add({xu, k, ub + sq_at[a], k, k, r, sb + s_at[a], 0}, k, r);
                xs.add({xv, k, vb + sq_at[a], k, k, r, nullptr, 1}, k, r);
            }
            yp.entries.push_back({B.u, xu, yu, k, r, k, B.m, k, k});
            yp.entries.push_back({B.vt, xv, yv, k, r, k, B.n, k, k});
            stats_.flops += kFlopScale * 4.0 * static_cast<double>(k) * k * r;
            if (host_t[a]) continue;
            T* t1 = rb + 2 * sq_at[a];
            T* t2 = t1 + static_cast<size_t>(k) * k;
            const T* gu = q.gram + q.top_at[a];
            ti.add({t1, gu, k, q.tau + q.tau_at[a], k}, k, k);
            ti.add({t2, gu + static_cast<size_t>(k) * k, k, q.tau + q.tau_at[a] + k, k}, k, k);
            for (int side = 0; side < 2; ++side) {
                tm.push_back(k);
                tn.push_back(r);
                tlda.push_back(k);
                tldb.push_back(k);
                tA.push_back(side == 0 ? t1 : t2);
                tB.push_back(side == 0 ? yu : yv);
            }
            stats_.flops += kFlopScale * 2.0 * static_cast<double>(k) * k * r;
        }
        const magma_int_t ntr = static_cast<magma_int_t>(tA.size());
        magma_int_t max_tm = 0, max_tn = 0;
        for (magma_int_t i = 0; i < ntr; ++i) {
            max_tm = std::max(max_tm, tm[static_cast<size_t>(i)]);
            max_tn = std::max(max_tn, tn[static_cast<size_t>(i)]);
        }
        tm.push_back(0);
        tn.push_back(0);
        tlda.push_back(0);
        tldb.push_back(0);
        xs.stage(meta_);
        ti.stage(meta_);
        yp.stage(meta_);
        const size_t off_tm = meta_.append(tm), off_tn = meta_.append(tn), off_tlda = meta_.append(tlda);
        const size_t off_tldb = meta_.append(tldb), off_tA = meta_.append(tA), off_tB = meta_.append(tB);
        char* md = meta_.upload(meta_device_, stream);
        launch_col_scale(xs.device(md), xs.count(), xs.max_m, xs.max_n, stream);
        launch_tinv(ti.device(md), ti.count(), ti.max_n, stream);
        yp.gemm(md, MagmaConjTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
        if (ntr > 0) {
            trsm_vbatched(max_tm, max_tn, reinterpret_cast<magma_int_t*>(md + off_tm), reinterpret_cast<magma_int_t*>(md + off_tn),
                          reinterpret_cast<T**>(md + off_tA), reinterpret_cast<magma_int_t*>(md + off_tlda),
                          reinterpret_cast<T**>(md + off_tB), reinterpret_cast<magma_int_t*>(md + off_tldb), ntr, queue);
        }
        check_cuda(cudaStreamSynchronize(stream), "HODLR GPU recompression Z");
        for (size_t a = 0; a < n_a; ++a) {
            if (host_t[a]) host_wy(ks[a], rn[a], q.gram + q.top_at[a], q.tau + q.tau_at[a], zb + x_at[a]);
        }
        heap.free(sb);
        heap.free(vb);
        heap.free(ub);
        heap.free(mb);
        heap.free(rb);
        q.release();
        stats_.t_recomp += seconds_since(t0);
        finish_device(na, bl, rn.data(), xb, zb);
        heap.free(xb);
    }

    // Free the blocks of the level
    void end() {
        release_out();
        if (blocks_.empty() && sel_rows_d_ == nullptr) return;
        Context::instance().activate();
        check_cuda(cudaStreamSynchronize(Context::instance().stream()), "HODLR GPU BACA end");
        DeviceHeap& heap = DeviceHeap::instance();
        for (Block& B : blocks_) {
            if (B.u != nullptr) heap.free(B.u);
            if (B.vt != nullptr) heap.free(B.vt);
            if (B.store != nullptr) heap.free(B.store);
            if (B.knn_store != nullptr) heap.free(B.knn_store);
        }
        blocks_.clear();
        if (sel_rows_d_ != nullptr) heap.free(sel_rows_d_);
        if (sel_cols_d_ != nullptr) heap.free(sel_cols_d_);
        if (sel_j2_d_ != nullptr) heap.free(sel_j2_d_);
        sel_rows_d_ = sel_cols_d_ = sel_j2_d_ = nullptr;
    }

    // Before the first GPU construction (GpuState::warm_up): BACA's QRs
    // (cuSOLVER's, with its handle, and MAGMA's batched one) and the batched
    // Jacobi SVD of the recompression on small matrices of a few sizes (MAGMA
    // picks its kernels by size), so that no level pays their first use.  The
    // statistics stay as they were.
    void warm_up() {
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        DeviceHeap& heap = DeviceHeap::instance();
        const Stats saved = stats_;
        magma_int_t* info = heap.alloc<magma_int_t>(sizeof(magma_int_t));
        void** ptrs = heap.alloc<void*>(4 * sizeof(void*));
        for (const int qr : {0, 1}) {
            const int m = qr == 0 ? 64 : 1024, n = qr == 0 ? 8 : 64;
            T* a = heap_identity<T>(m, n);
            T* tau = heap.alloc<T>(static_cast<size_t>(n) * sizeof(T));
            geqrf_one(m, n, a, m, tau);
            const void* host[2] = {a, tau};
            check_cuda(cudaMemcpyAsync(ptrs, host, sizeof(host), cudaMemcpyHostToDevice, stream), "warm-up");
            geqrf_batched(m, n, reinterpret_cast<T**>(ptrs), m, reinterpret_cast<T**>(ptrs + 1), info, 1, queue);
            check_cuda(cudaStreamSynchronize(stream), "warm-up");
            heap.free(tau);
            heap.free(a);
        }
        for (const int k : {8, 32, 128}) {
            T* a = heap_identity<T>(k, k);
            double* s = heap.alloc<double>(static_cast<size_t>(k) * sizeof(double));
            T* u = heap.alloc<T>(static_cast<size_t>(k) * k * sizeof(T));
            T* v = heap.alloc<T>(static_cast<size_t>(k) * k * sizeof(T));
            const void* host[4] = {a, s, u, v};
            check_cuda(cudaMemcpyAsync(ptrs, host, sizeof(host), cudaMemcpyHostToDevice, stream), "warm-up");
            gesvj_batched(k, reinterpret_cast<T**>(ptrs), reinterpret_cast<double**>(ptrs + 1),
                          reinterpret_cast<T**>(ptrs + 2), reinterpret_cast<T**>(ptrs + 3), info, 1, queue);
            check_cuda(cudaStreamSynchronize(stream), "warm-up");
            for (void* p : {static_cast<void*>(v), static_cast<void*>(u), static_cast<void*>(s), static_cast<void*>(a)}) {
                heap.free(p);
            }
        }
        heap.free(ptrs);
        heap.free(info);
        stats_ = saved;
    }

    // ... and the application's evaluator on up to 32 of the points
    // (fmm::gpu::Evaluator::warm_up, once per evaluator); its seconds
    double warm_up_evaluator(const fmm::gpu::Evaluator& evaluator) {
        if (ids_ == nullptr || n_points_ <= 0) return 0.0;
        Context& ctx = Context::instance();
        ctx.activate();
        return evaluator.warm_up<T>(points(), static_cast<int>(std::min<int64_t>(n_points_, 32)), ctx.stream());
    }

private:
    struct Block {
        int64_t r0 = 0, c0 = 0;  // first point slots of the rows and the columns
        int m = 0, n = 0, r = 0;  // sizes, panel width
        int cap = 0, rank = 0;
        T* u = nullptr;   // m x cap
        T* vt = nullptr;  // n x cap: V^T
        char* store = nullptr;
        T* c = nullptr;    // m x r: column panel
        T* rp = nullptr;   // r x n: row panel
        T* work = nullptr;  // pivot scratch: r x max(m, n), 2 max(m, n) norms, max(m, n) ints
        double* norms = nullptr;  // 3 max(m, n)
        int* perm = nullptr;
        size_t sel_offset = 0;
        int* sel_rows_d = nullptr;  // r panel rows (device)
        int* sel_cols_d = nullptr;  // r panel columns (device)
        int* j2_d = nullptr;        // BACA: the columns of the second column panel (device)
        std::vector<int> rows, cols;          // chosen rows and columns, in order
        std::vector<int> sel_rows, sel_cols;  // the current panel rows and columns
        // the knn step: panels C (m x nc), R (nr x n) and the pivot scratch
        char* knn_store = nullptr;
        T* ck = nullptr;
        T* rk = nullptr;
        T* kwork = nullptr;
        double* knorms = nullptr;
        int* kperm = nullptr;
        std::vector<int> knn_rows, knn_cols;
    };

    static constexpr double kFlopScale = fmm::gpu::is_complex_scalar<T> ? 4.0 : 1.0;
    static constexpr size_t kChunkBytes = size_t{1} << 30;
    static constexpr int kBigQrRank = 256;  // QR of blocks of this rank and more one by one (cuSOLVER)
    // the flops of the SVD of a k x k matrix, as the CPU counts its gesdd
    // (flops_dgesdd: 8 k^3 + 4/3 k^3), whichever SVD runs
    static double svd_flops(int k) { return kFlopScale * (8.0 + 4.0 / 3.0) * k * static_cast<double>(k) * k; }
    // the flops of the column-pivoted QR of an r x n panel, as the CPU counts
    // its geqp3 (flops_dgeqpfmod with k = min(r, n))
    static double qrcp_flops(int r, int n) {
        const double k = std::min(r, n), mx = std::max(r, n);
        return kFlopScale * (2.0 * mx * k * k - 2.0 / 3.0 * k * k * k);
    }
    // SVD of blocks of this rank and more one by one (cuSOLVER gesvdp), the
    // others by MAGMA's batched one-sided Jacobi: on an A100 the batched SVD
    // takes 0.1 ms a matrix of 100 x 100, 9 ms of 448 x 448 and 31 ms of
    // 640 x 640, gesvdp 5, 15 and 27 ms
    static constexpr int kBigSvdRank = 512;

    // Pivot items of one launch: the ones with many columns spread over
    // thread blocks (launch_qrcp_pivots_wide), the others one block each
    struct PivotList {
        static constexpr int kWideMinL = 8192;     // always wide from this many columns
        static constexpr int kWideMinLFew = 1024;  // ... and from this many when the launch has few items
        static constexpr int kFewItems = 128;
        std::vector<PivotItem<T>> all;
        ItemList<PivotItem<T>> narrow, wide;
        int max_r = 0, max_steps = 0, max_p = 0;
        void add(const PivotItem<T>& it) { all.push_back(it); }
        int count() const { return static_cast<int>(all.size()); }
        void stage(MetaBuilder& meta) {
            const bool few = all.size() < static_cast<size_t>(kFewItems);
            for (const PivotItem<T>& it : all) {
                if (it.L >= kWideMinL || (few && it.L >= kWideMinLFew)) {
                    wide.add(it, it.r, it.L);
                    max_r = std::max(max_r, it.r);
                    max_steps = std::max(max_steps, std::min(it.p, std::min(it.r, it.L)));
                    max_p = std::max(max_p, it.p);
                } else {
                    narrow.add(it, it.r, it.L);
                }
            }
            narrow.stage(meta);
            wide.stage(meta);
        }
        void launch(char* md, cudaStream_t stream) const {
            launch_qrcp_pivots(narrow.device(md), narrow.count(), md, stream);
            launch_qrcp_pivots_wide(wide.device(md), wide.count(), max_r, wide.max_n, max_steps, max_p, md, stream);
        }
    };

    // out (m x n, leading dimension m, zero on entry) = op(a)^T b, a (K x m)
    // and b (K x n) with leading dimensions lda, ldb
    struct GramJob {
        const T* a;
        int lda;
        const T* b;
        int ldb;
        int m, n, K;
        T* out;
    };

    // The Gram products of jobs (op = MagmaTrans or MagmaConjTrans).  Their
    // outputs are small and K long, so each product with few output tiles is
    // split along K into chunks computed side by side, and the chunks are
    // summed in order.
    void run_grams(const std::vector<GramJob>& jobs, magma_trans_t op) {
        if (jobs.empty()) return;
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        DeviceHeap& heap = DeviceHeap::instance();
        constexpr int kMinChunk = 2048;   // rows per chunk, at least
        constexpr int kTargetTiles = 512;  // 32 x 32 output tiles per product, at least (if K allows)
        std::vector<int> chunks(jobs.size(), 1), rows(jobs.size(), 0);
        size_t scratch = 0;
        for (size_t i = 0; i < jobs.size(); ++i) {
            const GramJob& j = jobs[i];
            const int tiles = std::max(1, ((j.m + 31) / 32) * ((j.n + 31) / 32));
            int nch = std::min((j.K + kMinChunk - 1) / kMinChunk, std::max(1, kTargetTiles / tiles));
            nch = std::max(nch, 1);
            const int kc = (j.K + nch - 1) / nch;
            nch = (j.K + kc - 1) / kc;
            chunks[i] = nch;
            rows[i] = kc;
            if (nch > 1) scratch += static_cast<size_t>(nch) * aligned(static_cast<size_t>(j.m) * j.n * sizeof(T));
        }
        char* sbuf = scratch > 0 ? heap.alloc<char>(scratch) : nullptr;
        meta_.clear();
        VBatch<T> g;
        ItemList<fmm::gpu::SumAddItemT<T>> sums;
        size_t at = 0;
        for (size_t i = 0; i < jobs.size(); ++i) {
            const GramJob& j = jobs[i];
            if (chunks[i] == 1) {
                g.entries.push_back({j.a, j.b, j.out, j.m, j.n, j.K, j.lda, j.ldb, j.m});
                continue;
            }
            std::vector<T*> parts;
            for (int c = 0; c < chunks[i]; ++c) {
                const int k0 = c * rows[i], kk = std::min(rows[i], j.K - k0);
                T* p = reinterpret_cast<T*>(sbuf + at);
                at += aligned(static_cast<size_t>(j.m) * j.n * sizeof(T));
                g.entries.push_back({j.a + k0, j.b + k0, p, j.m, j.n, kk, j.lda, j.ldb, j.m});
                parts.push_back(p);
            }
            const size_t off = meta_.append(parts);
            sums.add({j.out, j.m, j.m, j.n, static_cast<int64_t>(off), chunks[i]}, j.m, j.n);
        }
        g.stage(meta_);
        sums.stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        g.gemm(md, op, MagmaNoTrans, T(1.0), T(0.0), ctx.queue());
        fmm::gpu::launch_sum_add(sums.device(md), sums.count(), sums.max_m, sums.max_n, md, stream);
        if (sbuf != nullptr) heap.free(sbuf);
    }

    static size_t aligned(size_t bytes) { return fmm::gpu::align_up(std::max<size_t>(bytes, 1)); }

    // the new factors of finish_device() (for download())
    char* out_ = nullptr;
    std::vector<T*> out_u_, out_v_;
    void release_out() {
        if (out_ != nullptr) DeviceHeap::instance().free(out_);
        out_ = nullptr;
        out_u_.clear();
        out_v_.clear();
    }

    // MAGMA's batched one-sided Jacobi SVD of count k x k matrices (U, V and
    // the singular values, descending; a destroyed)
    static void gesvj_batched(magma_int_t k, T** a, double** s, T** u, T** v, magma_int_t* info, magma_int_t count,
                              magma_queue_t queue) {
        magma_int_t status = 0;
        if constexpr (std::is_same_v<T, double>) {
            status = magma_dgesvj_batched(MagmaSomeVec, MagmaSomeVec, k, k, a, k, s, u, k, v, k, info, count, queue);
        } else {
            status = magma_zgesvj_batched(MagmaSomeVec, MagmaSomeVec, k, k, reinterpret_cast<magmaDoubleComplex**>(a), k, s,
                                          reinterpret_cast<magmaDoubleComplex**>(u), k,
                                          reinterpret_cast<magmaDoubleComplex**>(v), k, info, count, queue);
        }
        if (status != 0) throw std::runtime_error("magma gesvj_batched failed with status " + std::to_string(status));
    }

    // B := A^-1 B for count upper triangular A (m x m) and B (m x n), sizes per entry (device arrays)
    static void trsm_vbatched(magma_int_t max_m, magma_int_t max_n, magma_int_t* m, magma_int_t* n, T** a, magma_int_t* lda,
                              T** b, magma_int_t* ldb, magma_int_t count, magma_queue_t queue) {
        if constexpr (std::is_same_v<T, double>) {
            magmablas_dtrsm_vbatched_max_nocheck(MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit, max_m, max_n, m, n, 1.0, a,
                                                 lda, b, ldb, count, queue);
        } else {
            magmablas_ztrsm_vbatched_max_nocheck(MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit, max_m, max_n, m, n,
                                                 MAGMA_Z_MAKE(1.0, 0.0), reinterpret_cast<magmaDoubleComplex**>(a), lda,
                                                 reinterpret_cast<magmaDoubleComplex**>(b), ldb, count, queue);
        }
    }

    // Z = T Y1^H X for a block with a zero tau (U's then V's, each k x rn,
    // Y1^H X in z on entry), T by the recurrence of ?larft ('F', 'C') on the host
    static void host_wy(int k, int rn, const T* gram, const T* tau, T* z) {
        const size_t kk = static_cast<size_t>(k) * k, kr = static_cast<size_t>(k) * rn;
        std::vector<T> g(2 * kk), tv(2 * static_cast<size_t>(k)), y(2 * kr), zz(2 * kr), t(kk), v(static_cast<size_t>(k));
        check_cuda(cudaMemcpy(g.data(), gram, 2 * kk * sizeof(T), cudaMemcpyDeviceToHost), "HODLR GPU WY gram");
        check_cuda(cudaMemcpy(tv.data(), tau, 2 * static_cast<size_t>(k) * sizeof(T), cudaMemcpyDeviceToHost), "HODLR GPU WY tau");
        check_cuda(cudaMemcpy(y.data(), z, 2 * kr * sizeof(T), cudaMemcpyDeviceToHost), "HODLR GPU WY rhs");
        for (int side = 0; side < 2; ++side) {
            const T* gs = g.data() + side * kk;
            const T* ts = tv.data() + side * k;
            std::fill(t.begin(), t.end(), T(0.0));
            for (int i = 0; i < k; ++i) {
                if (!(ts[i] == T(0.0))) {
                    for (int j = 0; j < i; ++j) v[static_cast<size_t>(j)] = T(0.0) - ts[i] * gs[j + static_cast<size_t>(i) * k];
                    for (int j = 0; j < i; ++j) {
                        T acc = T(0.0);
                        for (int l = j; l < i; ++l) acc += t[j + static_cast<size_t>(l) * k] * v[static_cast<size_t>(l)];
                        t[j + static_cast<size_t>(i) * k] = acc;
                    }
                }
                t[i + static_cast<size_t>(i) * k] = ts[i];
            }
            for (int c = 0; c < rn; ++c) {
                for (int i = 0; i < k; ++i) {
                    T acc = T(0.0);
                    for (int l = i; l < k; ++l) acc += t[i + static_cast<size_t>(l) * k] * y[side * kr + l + static_cast<size_t>(c) * k];
                    zz[side * kr + i + static_cast<size_t>(c) * k] = acc;
                }
            }
        }
        check_cuda(cudaMemcpy(z, zz.data(), 2 * kr * sizeof(T), cudaMemcpyHostToDevice), "HODLR GPU WY Z");
    }
    // device time of the entry evaluations (Stats::t_entry): an event pair
    // around each launch, read after the synchronization that follows
    void entry_begin(cudaStream_t stream) {
        cudaEvent_t a = nullptr, b = nullptr;
        check_cuda(cudaEventCreate(&a), "HODLR GPU entry timer");
        check_cuda(cudaEventCreate(&b), "HODLR GPU entry timer");
        check_cuda(cudaEventRecord(a, stream), "HODLR GPU entry timer");
        entry_events_.push_back({a, b});
    }
    void entry_end(cudaStream_t stream) {
        check_cuda(cudaEventRecord(entry_events_.back().second, stream), "HODLR GPU entry timer");
    }
    void entry_collect() {
        for (auto& e : entry_events_) {
            float ms = 0.0f;
            check_cuda(cudaEventElapsedTime(&ms, e.first, e.second), "HODLR GPU entry timer");
            stats_.t_entry += 1.0e-3 * ms;
            cudaEventDestroy(e.first);
            cudaEventDestroy(e.second);
        }
        entry_events_.clear();
    }
    std::vector<std::pair<cudaEvent_t, cudaEvent_t>> entry_events_;

    static double seconds_since(std::chrono::steady_clock::time_point t0) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    static fmm::gpu::IndexList contiguous(int64_t base) { return fmm::gpu::IndexList{-1, static_cast<int>(base), nullptr}; }
    static fmm::gpu::IndexList listed(size_t offset, int64_t base) {
        return fmm::gpu::IndexList{static_cast<int64_t>(offset), static_cast<int>(base), nullptr};
    }
    static fmm::gpu::IndexList device_list(const int* ptr) { return fmm::gpu::IndexList{-1, 0, ptr}; }
    static int checked(int i, int n) {
        if (i < 0 || i >= n) throw std::out_of_range("HODLR GPU BACA: index out of range");
        return i;
    }
    static void copy_down(T* host, const T* dev, size_t count) {
        if (count > 0) check_cuda(cudaMemcpy(host, dev, count * sizeof(T), cudaMemcpyDeviceToHost), "download");
    }

    Block& block(int b) {
        if (b < 0 || static_cast<size_t>(b) >= blocks_.size()) throw std::out_of_range("HODLR GPU BACA: no such block");
        return blocks_[static_cast<size_t>(b)];
    }
    const Block& block(int b) const {
        if (b < 0 || static_cast<size_t>(b) >= blocks_.size()) throw std::out_of_range("HODLR GPU BACA: no such block");
        return blocks_[static_cast<size_t>(b)];
    }

    void check_point_range(int64_t first, int count) const {
        if (ids_ == nullptr) throw std::logic_error("HODLR GPU construction: no points (c_bpack_hodlr_gpu_set_points)");
        if (first < 0 || first + count > n_points_ || n_points_ > std::numeric_limits<int>::max()) {
            throw std::out_of_range("HODLR GPU construction: point range outside the point table");
        }
    }

    fmm::gpu::PointTable points() const {
        fmm::gpu::PointTable p;
        p.xyz = dim_ > 0 ? xyz_ : nullptr;
        p.ids = ids_;
        p.dim = dim_;
        return p;
    }

    // capacity of U and V^T for `need` columns (the CPU's rule: double, at most min(m, n))
    void grow(Block& B, int need) {
        if (need <= B.cap) return;
        const int cap = std::min(std::max(need, 2 * B.cap), std::min(B.m, B.n));
        DeviceHeap& heap = DeviceHeap::instance();
        const cudaStream_t stream = Context::instance().stream();
        T* u = heap.alloc_resident<T>(static_cast<size_t>(B.m) * cap * sizeof(T));
        T* vt = heap.alloc_resident<T>(static_cast<size_t>(B.n) * cap * sizeof(T));
        if (B.rank > 0) {
            check_cuda(cudaMemcpyAsync(u, B.u, static_cast<size_t>(B.m) * B.rank * sizeof(T), cudaMemcpyDeviceToDevice,
                                       stream),
                       "HODLR GPU BACA grow U");
            check_cuda(cudaMemcpyAsync(vt, B.vt, static_cast<size_t>(B.n) * B.rank * sizeof(T),
                                       cudaMemcpyDeviceToDevice, stream),
                       "HODLR GPU BACA grow V");
        }
        if (B.u != nullptr) heap.free(B.u);
        if (B.vt != nullptr) heap.free(B.vt);
        B.u = u;
        B.vt = vt;
        B.cap = cap;
    }

    // all panel rows or columns (every block) into sel_host_
    void download_selection(const int* d, cudaStream_t stream) {
        if (!sel_host_.empty()) {
            check_cuda(cudaMemcpyAsync(sel_host_.data(), d, sel_host_.size() * sizeof(int), cudaMemcpyDeviceToHost,
                                       stream),
                       "HODLR GPU BACA panel indices");
        }
    }

    // QR of one m x n matrix in place (LAPACK layout: R and the Householder
    // vectors, tau), with cuSOLVER on the stream
    void geqrf_one(int m, int n, T* a, int lda, T* tau) {
        auto check = [](cusolverStatus_t st, const char* what) {
            if (st != CUSOLVER_STATUS_SUCCESS) {
                throw std::runtime_error(std::string(what) + " failed with status " + std::to_string(static_cast<int>(st)));
            }
        };
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        if (solver_ == nullptr) {
            check(cusolverDnCreate(&solver_), "cusolverDnCreate");
            check(cusolverDnCreateParams(&solver_params_), "cusolverDnCreateParams");
        }
        check(cusolverDnSetStream(solver_, stream), "cusolverDnSetStream");
        const cudaDataType type = fmm::gpu::is_complex_scalar<T> ? CUDA_C_64F : CUDA_R_64F;
        size_t dbytes = 0, hbytes = 0;
        check(cusolverDnXgeqrf_bufferSize(solver_, solver_params_, m, n, type, a, lda, type, tau, type, &dbytes, &hbytes),
              "cusolverDnXgeqrf_bufferSize");
        DeviceHeap& heap = DeviceHeap::instance();
        const size_t wbytes = fmm::gpu::align_up(std::max<size_t>(dbytes, 1));
        char* dwork = heap.alloc<char>(wbytes + sizeof(int));
        int* dinfo = reinterpret_cast<int*>(dwork + wbytes);
        std::vector<char> hwork(std::max<size_t>(hbytes, 1));
        check(cusolverDnXgeqrf(solver_, solver_params_, m, n, type, a, lda, type, tau, type, dwork, dbytes, hwork.data(),
                               hbytes, dinfo),
              "cusolverDnXgeqrf");
        heap.free(dwork);
    }

    // MAGMA's batched geqrf through its _work variant with the workspace of
    // its blocked path: magma_[dz]geqrf_batched (MAGMA 2.9.0) sizes the
    // workspace for the fused panel path when its tuning picks that one (0
    // bytes, so none is allocated), but falls back to the blocked path when
    // the fused kernel declines the size (tall matrices: n <= 96 and m of
    // about 12k and more in double, 6k in complex), which then writes its
    // pointer arrays through the null workspace (illegal memory access; also
    // square matrices of order 31 and 32)
    static void geqrf_batched(magma_int_t m, magma_int_t n, T** a, magma_int_t lda, T** tau, magma_int_t* info,
                              magma_int_t count, magma_queue_t queue) {
        if (count <= 0 || m <= 0 || n <= 0) return;
        constexpr bool real = std::is_same_v<T, double>;
        const magma_int_t nb = real ? magma_get_dgeqrf_batched_nb(m) : magma_get_zgeqrf_batched_nb(m);
        const size_t ld = static_cast<size_t>(std::min(nb, std::min(m, n)));
        const size_t ptrs = (4 * static_cast<size_t>(count) + 15) / 16 * 16;  // (dR, dT, 2 dW arrays; 128-byte rounding)
        const size_t bytes =
            (2 * ld * ld * count + 2 * static_cast<size_t>(nb) * n * count) * sizeof(T) + ptrs * sizeof(void*);
        if (bytes > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("HODLR GPU BACA: batched geqrf workspace too large");
        }
        DeviceHeap& heap = DeviceHeap::instance();
        void* work = heap.alloc<char>(bytes);
        magma_int_t lwork = static_cast<magma_int_t>(bytes);
        magma_int_t status = 0;
        if constexpr (real) {
            status = magma_dgeqrf_batched_work(m, n, a, lda, tau, info, work, &lwork, count, queue);
        } else {
            status = magma_zgeqrf_batched_work(m, n, reinterpret_cast<magmaDoubleComplex**>(a), lda,
                                               reinterpret_cast<magmaDoubleComplex**>(tau), info, work, &lwork, count,
                                               queue);
        }
        heap.free(work);  // (stream ordered: the next user of the memory runs after MAGMA's kernels)
        if (status != 0) throw std::runtime_error("magma geqrf_batched failed with status " + std::to_string(status));
    }

    void release_points() {
        if (ids_ == nullptr) return;
        Context::instance().activate();
        DeviceHeap& heap = DeviceHeap::instance();
        heap.free(xyz_);
        heap.free(ids_);
        xyz_ = nullptr;
        ids_ = nullptr;
        n_points_ = 0;
    }

    double* xyz_ = nullptr;
    int64_t* ids_ = nullptr;
    int64_t n_points_ = 0;
    int dim_ = 0;
    const fmm::gpu::Evaluator* evaluator_ = nullptr;  // the application's (GpuState::evaluator)
    double scale_ = 1.0;
    std::vector<Block> blocks_;
    int* sel_rows_d_ = nullptr;
    int* sel_cols_d_ = nullptr;
    int* sel_j2_d_ = nullptr;
    int variant_ = 5;  // RecLR_leaf: 5 BACA without overlap, 4 BACA
    std::vector<int> sel_host_;
    MetaBuilder meta_;
    DeviceBuffer meta_device_;
    Stats stats_;
    cusolverDnHandle_t solver_ = nullptr;
    cusolverDnParams_t solver_params_ = nullptr;
};

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
