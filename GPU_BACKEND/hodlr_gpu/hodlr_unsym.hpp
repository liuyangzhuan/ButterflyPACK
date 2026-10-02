#pragma once
// Unsymmetric HODLR factorization and solve on the GPU (option%sym=0), the
// device form of HODLR_factorization with LRlevel 0 (SRC/BPACK_factor.f90)
// and of HODLR_Inv_Apply (SRC/BPACK_solve_mul.f90), operation for operation:
//   - dense leaves: LU with the pivot rule of getrfmodf90, explicit inverse;
//   - Sblock (LR_Sblock): U' = A_ii^{-1} U for every off-diagonal block,
//     applying the leaf inverses and then, level by level from the bottom,
//     each lower level's 2 x 2 inverse (the CPU applies the same operators
//     in the same order to each block; here a level applies its inverse to
//     the blocks of all levels above it at once);
//   - the 2 x 2 inverse of a node (LR_minusBC, LR_SMW): with the updated
//     blocks B12 = U1' V1^T (rows of child 0) and B21 = U2' V2^T (rows of
//     child 1), Uo = -(-U1' (V1^T U2')) (I + V2^T (-U1' (V1^T U2')))^{-1},
//     so that (I - B12 B21)^{-1} = I + Uo V2^T.
// A node's inverse applied to X = [X1; X2] ('N', BF_block_MVP_inverse_dat):
//   X1 -= U1' V1^T X2;  X1 += Uo V2^T X1;  X2 -= U2' V2^T X1.
// The factors stay on the device for the solve.
//
// With several ranks, a rank holds of the node of a shared level (one per
// level) the rows of one child: in child 0, U1', V2 and Uo; in child 1, U2'
// and V1 (HodlrDevice moves each V to the rows it multiplies).  Each k x k
// or k x ncols product V^T U is then a sum over the node's ranks, and the LU
// and inverse of I + V2^T (-U1' V1^T U2') are formed on the node's head
// (rank 0 of its communicator) and broadcast.  The log-determinant here is
// this rank's part (shared nodes on the head).

#ifdef H2_HAVE_GPU

#include "hodlr_batch.hpp"
#include "hodlr_device.hpp"
#include "hodlr_lu.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bpack {
namespace gpu {

template<typename T>
class HodlrUnsymFactor {
public:
    struct Stats {
        double logabsdet = 0.0;
        std::complex<double> phase{1.0, 0.0};
        double flops = 0.0;
        double mbytes = 0.0;  // resident factor data (beyond the forward blocks)
        double t_leaf = 0.0, t_sblock = 0.0, t_nodes = 0.0;
        int raised_leaves = 0, raised_nodes = 0;  // LUs with pivots raised to the threshold
        int maxrank = 0;
    };

    HodlrUnsymFactor() = default;
    HodlrUnsymFactor(const HodlrUnsymFactor&) = delete;
    HodlrUnsymFactor& operator=(const HodlrUnsymFactor&) = delete;
    ~HodlrUnsymFactor() { release(); }

    bool ready() const { return ready_; }
    const Stats& stats() const { return stats_; }

    void factor(const HodlrDevice<T>& dev, double jitter_option) {
        using clock = std::chrono::steady_clock;
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        release();
        stats_ = Stats{};
        build_topology(dev);
        allocate();

        // dense leaves: LU with the CPU's pivot rule, log-determinant, inverse
        auto t0 = clock::now();
        {
            std::vector<const T*> src;
            std::vector<T*> lu, inv;
            std::vector<int*> piv;
            std::vector<int> sizes;
            for (Leaf& lf : leaves_) {
                src.push_back(lf.d);
                lu.push_back(lf.lu);
                inv.push_back(lf.inv);
                piv.push_back(lf.ipiv);
                sizes.push_back(lf.m);
            }
            stats_.raised_leaves = factor_lu(src, lu, piv, sizes, jitter_option, "dense leaf");
            lu_inverses(std::vector<const T*>(lu.begin(), lu.end()), std::vector<const int*>(piv.begin(), piv.end()),
                        sizes, inv, meta_, meta_device_);
            for (int m : sizes) stats_.flops += kFlopScale * (2.0 / 3.0 + 2.0) * std::pow(static_cast<double>(m), 3);
        }
        // Sblock, leaf part: U' = D^{-1} U on the rows of each leaf of each block
        {
            VBatch<T> b;
            for (auto& level : blocks_) {
                for (Block& bl : level) {
                    if (bl.k == 0 || bl.m == 0) continue;
                    for (const Leaf* lf : leaves_in(bl.row0, bl.m)) {
                        const int64_t off = lf->row0 - bl.row0;
                        b.entries.push_back({lf->inv, bl.uf + off, bl.up + off, lf->m, bl.k, lf->m, lf->m, bl.m, bl.m});
                        stats_.flops += kFlopScale * 2.0 * lf->m * static_cast<double>(lf->m) * bl.k;
                    }
                }
            }
            meta_.clear();
            b.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            b.gemm(md, MagmaNoTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
            check_cuda(cudaStreamSynchronize(stream), "HODLR unsym leaf Sblock");
        }
        auto t1 = clock::now();

        // levels, bottom up: the node inverses, then their application to
        // the updated blocks of the levels above
        double t_apply = 0.0;
        for (int l = maxlevel_; l >= 1; --l) {
            Node* shared = shared_node(l);
            if (shared != nullptr) {
                factor_shared(*shared, jitter_option);
            } else {
                factor_level(l, jitter_option);
            }
            if (l > 1) {
                auto ta = clock::now();
                std::vector<Target> targets;
                for (const Node& nd : nodes_[static_cast<size_t>(l)]) {
                    // (the rows of the children this rank holds: one child on a shared level)
                    int64_t lo = nd.c0, hi = nd.c1 + nd.n1;
                    if (nd.n0 == 0) {
                        lo = nd.c1;
                    } else if (nd.n1 == 0) {
                        hi = nd.c0 + nd.n0;
                    }
                    for (int al = 1; al < l; ++al) {
                        const Block* anc = block_holding(al, lo, hi);
                        if (anc->k == 0) continue;
                        targets.push_back({&nd, anc->up + (nd.c0 - anc->row0), anc->up + (nd.c1 - anc->row0), anc->m, anc->k});
                    }
                }
                if (shared != nullptr) {
                    apply_shared('N', *shared, targets);
                } else {
                    apply_nodes('N', targets);
                }
                check_cuda(cudaStreamSynchronize(stream), "HODLR unsym Sblock");
                t_apply += std::chrono::duration<double>(clock::now() - ta).count();
            }
        }
        check_cuda(cudaStreamSynchronize(stream), "HODLR unsym factorization");
        auto t2 = clock::now();
        stats_.t_leaf = std::chrono::duration<double>(t1 - t0).count();
        stats_.t_sblock = t_apply;
        stats_.t_nodes = std::chrono::duration<double>(t2 - t1).count() - t_apply;
        ready_ = true;
    }

    // y = op(A)^{-1} x ('N' or 'T') for the local n_loc x nrhs host arrays x and y
    double solve_flops() const { return solve_flops_; }  // of the last solve()
    void solve(char trans, int nrhs, const T* x, T* y) {
        if (!ready_) throw std::logic_error("HODLR GPU: unsymmetric solve before the factorization");
        if (trans != 'N' && trans != 'T') throw std::invalid_argument("HODLR GPU: solve takes 'N' or 'T'");
        const double flops0 = stats_.flops;  // (the level applications count there; the factorization's count is kept)
        for (const Leaf& lf : leaves_) stats_.flops += kFlopScale * 2.0 * lf.m * static_cast<double>(lf.m) * nrhs;
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const size_t bytes = static_cast<size_t>(n_loc_) * nrhs * sizeof(T);
        DeviceHeap& heap = DeviceHeap::instance();
        T* dx = heap.alloc<T>(bytes);
        T* dy = heap.alloc<T>(bytes);
        check_cuda(cudaMemcpyAsync(dx, x, bytes, cudaMemcpyHostToDevice, stream), "HODLR solve input");
        const int ld = static_cast<int>(n_loc_);
        auto leaf_apply = [&](const T* in, T* out) {
            VBatch<T> b;
            for (const Leaf& lf : leaves_) b.entries.push_back({lf.inv, in + lf.row0, out + lf.row0, lf.m, nrhs, lf.m, lf.m, ld, ld});
            meta_.clear();
            b.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            b.gemm(md, trans == 'N' ? MagmaNoTrans : MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), ctx.queue());
        };
        auto level_apply = [&](int l, T* v) {
            std::vector<Target> targets;
            for (const Node& nd : nodes_[static_cast<size_t>(l)]) targets.push_back({&nd, v + nd.c0, v + nd.c1, ld, nrhs});
            Node* shared = shared_node(l);
            if (shared != nullptr) {
                apply_shared(trans, *shared, targets);
            } else {
                apply_nodes(trans, targets);
            }
        };
        T* result = dy;
        if (trans == 'N') {
            leaf_apply(dx, dy);
            for (int l = maxlevel_; l >= 1; --l) level_apply(l, dy);
        } else {
            for (int l = 1; l <= maxlevel_; ++l) level_apply(l, dx);
            leaf_apply(dx, dy);
        }
        check_cuda(cudaMemcpyAsync(y, result, bytes, cudaMemcpyDeviceToHost, stream), "HODLR solve output");
        check_cuda(cudaStreamSynchronize(stream), "HODLR unsym solve");
        heap.free(dy);
        heap.free(dx);
        solve_flops_ = stats_.flops - flops0;
        stats_.flops = flops0;
    }

    // Copies for the comparison with the CPU (false if there is no such block):
    // U' of the block of `level` starting at local rows row0 / columns col0,
    // Uo of the node whose child 0 starts at c0, the inverse of the leaf at row0
    bool download_block(int level, int64_t row0, int64_t col0, int m, int k, T* up) const {
        if (!ready_ || level < 1 || level > maxlevel_) return false;
        for (const Block& bl : blocks_[static_cast<size_t>(level)]) {
            if (bl.row0 != row0 || bl.col0 != col0) continue;
            check_sizes(bl.m, bl.k, m, k, "block");
            copy_down(up, bl.up, static_cast<size_t>(m) * k);
            return true;
        }
        return false;
    }
    bool download_schur(int level, int64_t c0, int n0, int k, T* uo) const {
        if (!ready_ || level < 1 || level > maxlevel_) return false;
        for (const Node& nd : nodes_[static_cast<size_t>(level)]) {
            if (nd.c0 != c0) continue;
            check_sizes(nd.n0, nd.off2->k, n0, k, "Schur correction");
            copy_down(uo, nd.uo, static_cast<size_t>(n0) * k);
            return true;
        }
        return false;
    }
    bool download_leaf_inverse(int64_t row0, int m, T* inv) const {
        if (!ready_) return false;
        for (const Leaf& lf : leaves_) {
            if (lf.row0 != row0) continue;
            check_sizes(lf.m, lf.m, m, m, "leaf");
            copy_down(inv, lf.inv, static_cast<size_t>(m) * m);
            return true;
        }
        return false;
    }

private:
    double solve_flops_ = 0.0;
    struct Block {
        int level;
        int64_t row0, col0;
        int m, n, k;
        const T* uf;  // forward U (m x k)
        const T* v;   // forward V (n x k)
        T* up;        // updated U' = A_ii^{-1} U (m x k)
    };
    struct Node {
        int64_t c0, c1;  // first rows of child 0 and child 1 (of those this rank holds)
        int n0, n1;      // (0 for a child whose rows this rank does not hold)
        const Block* off1;  // rows of child 0, columns of child 1
        const Block* off2;  // rows of child 1, columns of child 0
        T* uo;              // n0 x k2 (k2 = off2->k)
        T* mlu;             // LU of I + V2^T (-U1' V1^T U2'), k2 x k2
        T* minv;            // its inverse
        int* ipiv;
        MPI_Comm comm = MPI_COMM_NULL;  // a shared node's ranks (MPI_COMM_NULL: this rank's alone)
        bool head = true;               // this rank factors I + V2^T Uo (rank 0 of comm, or a local node)
    };
    struct Leaf {
        int64_t row0;
        int m;
        const T* d;  // forward block
        T* lu;
        int* ipiv;
        T* inv;
    };
    struct Target {
        const Node* node;
        T* x1;  // n0 x ncols (rows of child 0)
        T* x2;  // n1 x ncols
        int ld, ncols;
    };

    static constexpr double kFlopScale = fmm::gpu::is_complex_scalar<T> ? 4.0 : 1.0;

    static size_t aligned(size_t bytes) { return fmm::gpu::align_up(std::max<size_t>(bytes, 1)); }

    static void check_sizes(int m, int k, int m_host, int k_host, const char* what) {
        if (m != m_host || k != k_host) {
            std::ostringstream oss;
            oss << "HODLR GPU: " << what << " of size " << m << " x " << k << " on the device but " << m_host << " x "
                << k_host << " on the host";
            throw std::logic_error(oss.str());
        }
    }
    static void copy_down(T* host, const T* dev, size_t count) {
        if (count > 0) check_cuda(cudaMemcpy(host, dev, count * sizeof(T), cudaMemcpyDeviceToHost), "download");
    }

    void build_topology(const HodlrDevice<T>& dev) {
        n_loc_ = dev.n_loc();
        maxlevel_ = dev.maxlevel();
        blocks_.assign(static_cast<size_t>(maxlevel_ + 1), {});
        for (const auto& b : dev.blocks()) {
            if (b.sym) throw std::logic_error("HODLR GPU: the unsymmetric factorization needs both blocks of a node");
            blocks_[static_cast<size_t>(b.level)].push_back(Block{b.level, b.row0, b.col0, b.m, b.n, b.k, b.u, b.v, nullptr});
        }
        nodes_.assign(static_cast<size_t>(maxlevel_ + 1), {});
        for (int l = 1; l <= maxlevel_; ++l) {
            auto& v = blocks_[static_cast<size_t>(l)];
            if (dev.shared_level(l)) {  // (A12, A21 of this rank's node, as uploaded)
                if (v.size() != 2) throw std::logic_error("HODLR GPU: a shared level must hold one node on each rank");
                for (const Block& b : v) stats_.maxrank = std::max(stats_.maxrank, b.k);
                continue;
            }
            std::sort(v.begin(), v.end(), [](const Block& a, const Block& b) { return a.row0 < b.row0; });
            if (!dev.distributed() && v.size() != (size_t{1} << l)) {
                throw std::logic_error("HODLR GPU: level " + std::to_string(l) + " has " + std::to_string(v.size()) +
                                       " blocks; the one-rank backend needs the whole tree");
            }
            for (const Block& b : v) stats_.maxrank = std::max(stats_.maxrank, b.k);
        }
        // pair the blocks of each node: off1 (rows before columns) and its transpose position off2
        for (int l = 1; l <= maxlevel_; ++l) {
            const auto& v = blocks_[static_cast<size_t>(l)];
            if (dev.shared_level(l)) {
                const Block& a12 = v[0];
                const Block& a21 = v[1];
                if (a12.n != a21.m || a21.n != a12.m) {
                    throw std::logic_error("HODLR GPU: the rows of a shared node's blocks on this rank do not match");
                }
                Node nd{a12.row0, a21.row0, a12.m, a21.m, &a12, &a21, nullptr, nullptr, nullptr, nullptr};
                nd.comm = dev.level_comm(l);
                int rank = 0;
                MPI_Comm_rank(nd.comm, &rank);
                nd.head = rank == 0;
                nodes_[static_cast<size_t>(l)].push_back(nd);
                continue;
            }
            for (const Block& b1 : v) {
                if (b1.row0 > b1.col0) continue;
                const Block* b2 = nullptr;
                for (const Block& c : v) {
                    if (c.row0 == b1.col0 && c.col0 == b1.row0) b2 = &c;
                }
                if (b2 == nullptr || b2->m != b1.n || b2->n != b1.m) {
                    throw std::logic_error("HODLR GPU: an off-diagonal block has no transposed partner");
                }
                nodes_[static_cast<size_t>(l)].push_back(Node{b1.row0, b1.col0, b1.m, b1.n, &b1, b2, nullptr, nullptr,
                                                              nullptr, nullptr});
            }
        }
        leaves_.clear();
        for (const auto& lf : dev.leaves()) leaves_.push_back(Leaf{lf.row0, lf.m, lf.d, nullptr, nullptr, nullptr});
        std::sort(leaves_.begin(), leaves_.end(), [](const Leaf& a, const Leaf& b) { return a.row0 < b.row0; });
    }

    void allocate() {
        size_t bytes = 0;
        for (const auto& level : blocks_)
            for (const Block& b : level) bytes += aligned(static_cast<size_t>(b.m) * b.k * sizeof(T));
        for (const auto& level : nodes_) {
            for (const Node& nd : level) {
                const size_t k = static_cast<size_t>(nd.off2->k);
                bytes += aligned(nd.n0 * k * sizeof(T)) + 2 * aligned(k * k * sizeof(T)) + aligned(k * sizeof(int));
            }
        }
        for (const Leaf& lf : leaves_) {
            const size_t m = static_cast<size_t>(lf.m);
            bytes += 2 * aligned(m * m * sizeof(T)) + aligned(m * sizeof(int));
        }
        store_ = DeviceHeap::instance().alloc_resident<char>(bytes);
        stats_.mbytes = static_cast<double>(bytes) / 1.0e6;
        char* at = store_;
        auto take = [&](size_t b) {
            char* p = at;
            at += aligned(b);
            return p;
        };
        for (auto& level : blocks_)
            for (Block& b : level) b.up = reinterpret_cast<T*>(take(static_cast<size_t>(b.m) * b.k * sizeof(T)));
        for (auto& level : nodes_) {
            for (Node& nd : level) {
                const size_t k = static_cast<size_t>(nd.off2->k);
                nd.uo = reinterpret_cast<T*>(take(nd.n0 * k * sizeof(T)));
                nd.mlu = reinterpret_cast<T*>(take(k * k * sizeof(T)));
                nd.minv = reinterpret_cast<T*>(take(k * k * sizeof(T)));
                nd.ipiv = reinterpret_cast<int*>(take(k * sizeof(int)));
            }
        }
        for (Leaf& lf : leaves_) {
            const size_t m = static_cast<size_t>(lf.m);
            lf.lu = reinterpret_cast<T*>(take(m * m * sizeof(T)));
            lf.inv = reinterpret_cast<T*>(take(m * m * sizeof(T)));
            lf.ipiv = reinterpret_cast<int*>(take(m * sizeof(int)));
        }
        if (at != store_ + bytes) throw std::logic_error("HODLR GPU: unsymmetric factors placed outside their allocation");
    }

    void release() {
        if (store_ != nullptr) {
            Context::instance().activate();
            DeviceHeap::instance().free(store_);
            store_ = nullptr;
        }
        blocks_.clear();
        nodes_.clear();
        leaves_.clear();
        ready_ = false;
    }

    // the leaves covering rows [row0, row0 + m)
    std::vector<const Leaf*> leaves_in(int64_t row0, int m) const {
        std::vector<const Leaf*> out;
        auto it = std::lower_bound(leaves_.begin(), leaves_.end(), row0,
                                   [](const Leaf& lf, int64_t r) { return lf.row0 < r; });
        int64_t covered = 0;
        for (; it != leaves_.end() && it->row0 < row0 + m; ++it) {
            out.push_back(&*it);
            covered += it->m;
        }
        if (covered != m) throw std::logic_error("HODLR GPU: the leaves do not tile a block's rows");
        return out;
    }

    // the block of level al whose rows hold [lo, hi)
    const Block* block_holding(int al, int64_t lo, int64_t hi) const {
        const auto& v = blocks_[static_cast<size_t>(al)];
        for (const Block& b : v) {  // a shared level: the block whose U rows this rank holds
            if (is_shared(al) && b.m > 0 && lo >= b.row0 && hi <= b.row0 + b.m) return &b;
        }
        auto it = std::upper_bound(v.begin(), v.end(), lo, [](int64_t r, const Block& b) { return r < b.row0; });
        if (it != v.begin()) {
            const Block& b = *std::prev(it);
            if (lo >= b.row0 && hi <= b.row0 + b.m) return &b;
        }
        std::ostringstream oss;
        oss << "HODLR GPU: no block of level " << al << " holds rows [" << lo << ", " << hi << ")";
        throw std::logic_error(oss.str());
    }

    // LU of the matrices src[i] into lu[i] (which may be src[i]) with the
    // pivot rule of getrfmodf90 (pivots below jitter * |A|_F raised to that
    // size), and their log-determinants; the matrices whose GPU LU has such
    // a pivot (or an exact zero) are redone on the host.  Returns their count.
    int factor_lu(const std::vector<const T*>& src, const std::vector<T*>& lu, const std::vector<int*>& piv,
                  const std::vector<int>& sizes, double jitter_option, const char* what) {
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const size_t count = src.size();
        if (count == 0) return 0;
        DeviceHeap& heap = DeviceHeap::instance();
        double* d_norms = heap.alloc<double>(sizeof(double) * count);
        ItemList<NormItem<T>> norms;
        ItemList<AxpbyItem<T>> copies;
        GetrfBatch<T> getrf;
        for (size_t i = 0; i < count; ++i) {
            const int n = sizes[i];
            norms.add({src[i], n, n, n, d_norms + i}, n, n);
            if (lu[i] != src[i]) copies.add({lu[i], n, src[i], n, nullptr, 0, n, n, T(1.0), T(0.0)}, n, n);
            getrf.add(lu[i], n, n, piv[i]);
        }
        meta_.clear();
        norms.stage(meta_);
        copies.stage(meta_);
        getrf.stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        launch_fnorm(norms.device(md), norms.count(), stream);
        launch_axpby(copies.device(md), copies.count(), copies.max_m, copies.max_n, stream);
        std::vector<double> norm(count);
        check_cuda(cudaMemcpyAsync(norm.data(), d_norms, sizeof(double) * count, cudaMemcpyDeviceToHost, stream), "norms");
        std::vector<int> failed = run_getrf(getrf, md, getrf_work_);  // (synchronizes)
        heap.free(d_norms);

        std::vector<T> diag;
        std::vector<int> pv;
        std::vector<char> redo(count, 0);
        for (int f : failed) redo[static_cast<size_t>(f)] = 1;
        {
            std::vector<const T*> m(lu.begin(), lu.end());
            std::vector<const int*> p(piv.begin(), piv.end());
            lu_diagonals(m, p, sizes, meta_, meta_device_, diag, pv, what);
        }
        size_t off = 0;
        for (size_t i = 0; i < count; ++i) {
            const double thresh = norm[i] * jitter_option;
            for (int r = 0; r < sizes[i] && !redo[i]; ++r) {
                if (std::abs(to_std(diag[off + static_cast<size_t>(r)])) < thresh) redo[i] = 1;
            }
            off += static_cast<size_t>(sizes[i]);
        }
        int raised = 0;
        for (size_t i = 0; i < count; ++i) {  // the host LU of getrfmodf90
            if (!redo[i]) continue;
            ++raised;
            const int n = sizes[i];
            std::vector<T> a(static_cast<size_t>(n) * n);
            std::vector<int> ip;
            check_cuda(cudaMemcpy(a.data(), src[i] == lu[i] ? lu_source_[i] : src[i], a.size() * sizeof(T),
                                  cudaMemcpyDeviceToHost),
                       "host LU");
            host_lu_threshold(a, n, norm[i] * jitter_option, ip);
            check_cuda(cudaMemcpy(lu[i], a.data(), a.size() * sizeof(T), cudaMemcpyHostToDevice), "host LU");
            check_cuda(cudaMemcpy(piv[i], ip.data(), ip.size() * sizeof(int), cudaMemcpyHostToDevice), "host LU");
        }
        if (raised > 0) {
            std::vector<const T*> m(lu.begin(), lu.end());
            std::vector<const int*> p(piv.begin(), piv.end());
            lu_diagonals(m, p, sizes, meta_, meta_device_, diag, pv, what);
        }
        off = 0;
        for (size_t i = 0; i < count; ++i) {  // log-determinants, in order
            std::complex<double> phase;
            double logabs = 0.0;
            lu_slogdet(diag.data() + off, pv.data() + off, sizes[i], phase, logabs);
            stats_.phase *= phase;
            stats_.logabsdet += logabs;
            off += static_cast<size_t>(sizes[i]);
        }
        return raised;
    }

    // the 2 x 2 inverses of the nodes of level l (LR_minusBC, LR_SMW)
    void factor_level(int l, double jitter_option) {
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        std::vector<const Node*> active;
        for (const Node& nd : nodes_[static_cast<size_t>(l)]) {
            if (nd.off2->k > 0) active.push_back(&nd);
        }
        if (active.empty()) return;
        size_t ws = 0;
        for (const Node* nd : active) {
            const size_t k1 = static_cast<size_t>(nd->off1->k), k2 = static_cast<size_t>(nd->off2->k);
            ws += k1 * k2 + static_cast<size_t>(nd->n0) * k2;
        }
        DeviceHeap& heap = DeviceHeap::instance();
        T* work = heap.alloc<T>(sizeof(T) * std::max<size_t>(ws, 1));
        std::vector<T*> uo_tmp;
        {
            VBatch<T> b_t, b_u, b_m;
            ItemList<DiagItem<T>> ones;
            size_t off = 0;
            for (const Node* nd : active) {
                const Block& o1 = *nd->off1;
                const Block& o2 = *nd->off2;
                const int k1 = o1.k, k2 = o2.k;
                T* t = work + off;
                T* u = t + static_cast<size_t>(k1) * k2;
                off += static_cast<size_t>(k1) * k2 + static_cast<size_t>(nd->n0) * k2;
                uo_tmp.push_back(u);
                if (k1 > 0) {
                    b_t.entries.push_back({o1.v, o2.up, t, k1, k2, nd->n1, nd->n1, nd->n1, k1});
                    b_u.entries.push_back({o1.up, t, u, nd->n0, k2, k1, nd->n0, k1, nd->n0});
                } else {
                    check_cuda(cudaMemsetAsync(u, 0, sizeof(T) * static_cast<size_t>(nd->n0) * k2, stream), "Uo");
                }
                b_m.entries.push_back({o2.v, u, nd->mlu, k2, k2, nd->n0, nd->n0, nd->n0, k2});
                ones.add({nd->mlu, k2, k2, T(1.0)}, k2, k2);
                stats_.flops += kFlopScale * 2.0 *
                                (static_cast<double>(k1) * k2 * (nd->n1 + nd->n0) + static_cast<double>(k2) * k2 * nd->n0);
            }
            meta_.clear();
            b_t.stage(meta_);
            b_u.stage(meta_);
            b_m.stage(meta_);
            ones.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            b_t.gemm(md, MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
            b_u.gemm(md, MagmaNoTrans, MagmaNoTrans, T(-1.0), T(0.0), queue);
            b_m.gemm(md, MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
            launch_add_diagonal(ones.device(md), ones.count(), ones.max_n, stream);
        }
        // LU of I + V2^T Uo (kept as the source of a host redo in minv), inverse
        {
            std::vector<const T*> src;
            std::vector<T*> lu, inv;
            std::vector<int*> piv;
            std::vector<int> sizes;
            ItemList<AxpbyItem<T>> keep;
            for (const Node* nd : active) {
                const int k2 = nd->off2->k;
                keep.add({nd->minv, k2, nd->mlu, k2, nullptr, 0, k2, k2, T(1.0), T(0.0)}, k2, k2);
                src.push_back(nd->mlu);
                lu.push_back(nd->mlu);
                inv.push_back(nd->minv);
                piv.push_back(nd->ipiv);
                sizes.push_back(k2);
                stats_.flops += kFlopScale * (2.0 / 3.0 + 2.0) * std::pow(static_cast<double>(k2), 3);
            }
            meta_.clear();
            keep.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            launch_axpby(keep.device(md), keep.count(), keep.max_m, keep.max_n, stream);
            lu_source_.assign(inv.begin(), inv.end());  // the unfactored matrices, for a host redo
            stats_.raised_nodes += factor_lu(src, lu, piv, sizes, jitter_option, "Schur correction");
            lu_inverses(std::vector<const T*>(lu.begin(), lu.end()), std::vector<const int*>(piv.begin(), piv.end()),
                        sizes, inv, meta_, meta_device_);
        }
        // Uo = -Uo_tmp M^{-1}
        {
            VBatch<T> b;
            for (size_t i = 0; i < active.size(); ++i) {
                const Node* nd = active[i];
                const int k2 = nd->off2->k;
                b.entries.push_back({uo_tmp[i], nd->minv, nd->uo, nd->n0, k2, k2, nd->n0, k2, nd->n0});
                stats_.flops += kFlopScale * 2.0 * nd->n0 * static_cast<double>(k2) * k2;
            }
            meta_.clear();
            b.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            b.gemm(md, MagmaNoTrans, MagmaNoTrans, T(-1.0), T(0.0), queue);
        }
        heap.free(work);
    }

    bool is_shared(int l) const {
        const auto& v = nodes_[static_cast<size_t>(l)];
        return v.size() == 1 && v[0].comm != MPI_COMM_NULL;
    }
    // the node of level l if the level is shared (its only node on this rank)
    Node* shared_node(int l) { return is_shared(l) ? &nodes_[static_cast<size_t>(l)][0] : nullptr; }

    // device buffer for a message: in the exchange arena when MPI is
    // CUDA-aware (then *in_arena), else from the heap (staged by sum_over)
    static T* message_buffer(size_t count, bool* in_arena) {
        T* p = nullptr;
        if (fmm::gpu::device_exchange_enabled()) {
            p = reinterpret_cast<T*>(DeviceHeap::exchange_arena().try_alloc(std::max<size_t>(count, 1) * sizeof(T)));
        }
        *in_arena = p != nullptr;
        return p != nullptr ? p : DeviceHeap::instance().alloc<T>(std::max<size_t>(count, 1) * sizeof(T));
    }
    static void free_message(T* p, bool in_arena) {
        if (in_arena) {
            DeviceHeap::exchange_arena().free(p);
        } else {
            DeviceHeap::instance().free(p);
        }
    }

    // The 2 x 2 inverse of a shared node: T = V1^T U2' (summed over the
    // node's ranks), W = -U1' T, M = I + V2^T W (summed), its LU and
    // inverse on the head (broadcast), Uo = W M^{-1} with the sign of
    // factor_level (Uo = -W_tmp M^{-1}, W_tmp = -U1' T)
    void factor_shared(Node& nd, double jitter_option) {
        const Block& o1 = *nd.off1;
        const Block& o2 = *nd.off2;
        const int k1 = o1.k, k2 = o2.k;
        if (k2 == 0) return;
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        DeviceHeap& heap = DeviceHeap::instance();
        bool t_arena = false, m_arena = false;
        T* t = message_buffer(static_cast<size_t>(std::max(k1, 1)) * k2, &t_arena);
        T* w = heap.alloc<T>(sizeof(T) * std::max<size_t>(static_cast<size_t>(nd.n0) * k2, 1));
        check_cuda(cudaMemsetAsync(t, 0, sizeof(T) * static_cast<size_t>(std::max(k1, 1)) * k2, stream), "HODLR T");
        check_cuda(cudaMemsetAsync(w, 0, sizeof(T) * std::max<size_t>(static_cast<size_t>(nd.n0) * k2, 1), stream), "HODLR W");
        if (k1 > 0) {
            // T = V1^T U2' on the ranks of child 1
            VBatch<T> b_t;
            if (nd.n1 > 0) b_t.entries.push_back({o1.v, o2.up, t, k1, k2, nd.n1, nd.n1, nd.n1, k1});
            meta_.clear();
            b_t.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            if (b_t.count() > 0) b_t.gemm(md, MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
            HodlrDevice<T>::sum_over(nd.comm, t, static_cast<int64_t>(k1) * k2, t_arena, stream);
            // W = -U1' T on the ranks of child 0
            VBatch<T> b_u;
            if (nd.n0 > 0) b_u.entries.push_back({o1.up, t, w, nd.n0, k2, k1, nd.n0, k1, nd.n0});
            meta_.clear();
            b_u.stage(meta_);
            md = meta_.upload(meta_device_, stream);
            if (b_u.count() > 0) b_u.gemm(md, MagmaNoTrans, MagmaNoTrans, T(-1.0), T(0.0), queue);
            stats_.flops += kFlopScale * 2.0 * static_cast<double>(k1) * k2 * (nd.n1 + nd.n0);
        }
        // M = I + V2^T W (V2 on the ranks of child 0)
        T* m = message_buffer(static_cast<size_t>(k2) * k2, &m_arena);
        check_cuda(cudaMemsetAsync(m, 0, sizeof(T) * static_cast<size_t>(k2) * k2, stream), "HODLR M");
        {
            VBatch<T> b_m;
            if (nd.n0 > 0) b_m.entries.push_back({o2.v, w, m, k2, k2, nd.n0, nd.n0, nd.n0, k2});
            meta_.clear();
            b_m.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            if (b_m.count() > 0) b_m.gemm(md, MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
            stats_.flops += kFlopScale * 2.0 * static_cast<double>(k2) * k2 * nd.n0;
        }
        HodlrDevice<T>::sum_over(nd.comm, m, static_cast<int64_t>(k2) * k2, m_arena, stream);
        int ok = 1;
        if (nd.head) {
            try {
                ItemList<DiagItem<T>> ones;
                ItemList<AxpbyItem<T>> keep;
                ones.add({m, k2, k2, T(1.0)}, k2, k2);
                keep.add({nd.mlu, k2, m, k2, nullptr, 0, k2, k2, T(1.0), T(0.0)}, k2, k2);
                keep.add({nd.minv, k2, m, k2, nullptr, 0, k2, k2, T(1.0), T(0.0)}, k2, k2);
                meta_.clear();
                ones.stage(meta_);
                keep.stage(meta_);
                char* md = meta_.upload(meta_device_, stream);
                launch_add_diagonal(ones.device(md), ones.count(), ones.max_n, stream);
                launch_axpby(keep.device(md), keep.count(), keep.max_m, keep.max_n, stream);
                lu_source_.assign(1, nd.minv);  // the unfactored matrix, for a host redo
                stats_.raised_nodes += factor_lu({nd.mlu}, {nd.mlu}, {nd.ipiv}, {k2}, jitter_option, "Schur correction");
                lu_inverses(std::vector<const T*>{nd.mlu}, std::vector<const int*>{nd.ipiv}, std::vector<int>{k2},
                            std::vector<T*>{nd.minv}, meta_, meta_device_);
                stats_.flops += kFlopScale * (2.0 / 3.0 + 2.0) * std::pow(static_cast<double>(k2), 3);
                // (M^{-1} goes out through m)
                ItemList<AxpbyItem<T>> out;
                out.add({m, k2, nd.minv, k2, nullptr, 0, k2, k2, T(1.0), T(0.0)}, k2, k2);
                meta_.clear();
                out.stage(meta_);
                md = meta_.upload(meta_device_, stream);
                launch_axpby(out.device(md), out.count(), out.max_m, out.max_n, stream);
            } catch (const std::runtime_error&) {
                ok = 0;
            }
        }
        MPI_Bcast(&ok, 1, MPI_INT, 0, nd.comm);
        if (!ok) throw std::runtime_error("unsymmetric HODLR LU failed for a Schur correction (GPU)");
        HodlrDevice<T>::bcast_from_root(nd.comm, m, static_cast<int64_t>(k2) * k2, m_arena, stream);
        {
            // M^{-1} on every rank of the node; Uo = -W M^{-1} on the ranks of child 0
            ItemList<AxpbyItem<T>> keep;
            if (!nd.head) keep.add({nd.minv, k2, m, k2, nullptr, 0, k2, k2, T(1.0), T(0.0)}, k2, k2);
            VBatch<T> b;
            if (nd.n0 > 0) b.entries.push_back({w, m, nd.uo, nd.n0, k2, k2, nd.n0, k2, nd.n0});
            meta_.clear();
            keep.stage(meta_);
            b.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            launch_axpby(keep.device(md), keep.count(), keep.max_m, keep.max_n, stream);
            if (b.count() > 0) b.gemm(md, MagmaNoTrans, MagmaNoTrans, T(-1.0), T(0.0), queue);
            stats_.flops += kFlopScale * 2.0 * nd.n0 * static_cast<double>(k2) * k2;
        }
        check_cuda(cudaStreamSynchronize(stream), "HODLR unsym shared node");
        free_message(m, m_arena);
        heap.free(w);
        free_message(t, t_arena);
    }

    // The node inverse of a shared node on its targets (apply_nodes with
    // several ranks): each of the three updates sums its coefficients over
    // the node's ranks first.
    void apply_shared(char trans, Node& nd, const std::vector<Target>& targets) {
        if (targets.empty()) return;
        const Block& o1 = *nd.off1;
        const Block& o2 = *nd.off2;
        const int k1 = o1.k, k2 = o2.k, n0 = nd.n0, n1 = nd.n1;
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        size_t cols = 0;
        for (const Target& t : targets) cols += static_cast<size_t>(t.ncols);
        const int kmax = std::max(std::max(k1, k2), 1);
        bool in_arena = false;
        T* c = message_buffer(static_cast<size_t>(kmax) * cols, &in_arena);
        // one update: coefficients c = A_part^T X_part (k rows) summed over the
        // node's ranks, then Y_part += alpha B_part c
        auto phase = [&](int k, const T* a, int na, bool from_x1, const T* b, int nb, bool to_x1, T alpha) {
            if (k == 0) return;
            check_cuda(cudaMemsetAsync(c, 0, sizeof(T) * static_cast<size_t>(k) * cols, stream), "HODLR coefficients");
            VBatch<T> bc, bu;
            size_t off = 0;
            for (const Target& t : targets) {
                T* ct = c + off;
                off += static_cast<size_t>(k) * t.ncols;
                if (na > 0) bc.entries.push_back({a, from_x1 ? t.x1 : t.x2, ct, k, t.ncols, na, na, t.ld, k});
                if (nb > 0) bu.entries.push_back({b, ct, to_x1 ? t.x1 : t.x2, nb, t.ncols, k, nb, k, t.ld});
                stats_.flops += kFlopScale * 2.0 * k * static_cast<double>(na + nb) * t.ncols;
            }
            meta_.clear();
            bc.stage(meta_);
            bu.stage(meta_);
            char* md = meta_.upload(meta_device_, stream);
            if (bc.count() > 0) bc.gemm(md, MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
            HodlrDevice<T>::sum_over(nd.comm, c, static_cast<int64_t>(k) * static_cast<int64_t>(cols), in_arena, stream);
            if (bu.count() > 0) bu.gemm(md, MagmaNoTrans, MagmaNoTrans, alpha, T(1.0), queue);
        };
        if (trans == 'N') {
            // X1 -= U1' (V1^T X2); X1 += Uo (V2^T X1); X2 -= U2' (V2^T X1)
            phase(k1, o1.v, n1, false, o1.up, n0, true, T(-1.0));
            phase(k2, o2.v, n0, true, nd.uo, n0, true, T(1.0));
            phase(k2, o2.v, n0, true, o2.up, n1, false, T(-1.0));
        } else {
            // X1 -= V2 (U2'^T X2); X1 += V2 (Uo^T X1); X2 -= V1 (U1'^T X1)
            phase(k2, o2.up, n1, false, o2.v, n0, true, T(-1.0));
            phase(k2, nd.uo, n0, true, o2.v, n0, true, T(1.0));
            phase(k1, o1.up, n0, true, o1.v, n1, false, T(-1.0));
        }
        check_cuda(cudaStreamSynchronize(stream), "HODLR unsym shared node apply");
        free_message(c, in_arena);
    }

    // the 2 x 2 node inverse on each target (BF_block_MVP_inverse_dat)
    void apply_nodes(char trans, const std::vector<Target>& targets) {
        if (targets.empty()) return;
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        size_t ws = 0;
        for (const Target& t : targets) {
            const size_t kmax = static_cast<size_t>(std::max(t.node->off1->k, t.node->off2->k));
            ws += 3 * kmax * t.ncols;
        }
        DeviceHeap& heap = DeviceHeap::instance();
        T* work = heap.alloc<T>(sizeof(T) * std::max<size_t>(ws, 1));
        VBatch<T> c1, u1, c2, u2, c3, u3;
        size_t off = 0;
        for (const Target& t : targets) {
            const Node& nd = *t.node;
            const Block& o1 = *nd.off1;
            const Block& o2 = *nd.off2;
            const int k1 = o1.k, k2 = o2.k, nc = t.ncols, n0 = nd.n0, n1 = nd.n1;
            const int kmax = std::max(k1, k2);
            T* w1 = work + off;
            T* w2 = w1 + static_cast<size_t>(kmax) * nc;
            T* w3 = w2 + static_cast<size_t>(kmax) * nc;
            off += 3 * static_cast<size_t>(kmax) * nc;
            if (trans == 'N') {
                // X1 -= U1' (V1^T X2); X1 += Uo (V2^T X1); X2 -= U2' (V2^T X1)
                if (k1 > 0) {
                    c1.entries.push_back({o1.v, t.x2, w1, k1, nc, n1, n1, t.ld, k1});
                    u1.entries.push_back({o1.up, w1, t.x1, n0, nc, k1, n0, k1, t.ld});
                }
                if (k2 > 0) {
                    c2.entries.push_back({o2.v, t.x1, w2, k2, nc, n0, n0, t.ld, k2});
                    u2.entries.push_back({nd.uo, w2, t.x1, n0, nc, k2, n0, k2, t.ld});
                    c3.entries.push_back({o2.v, t.x1, w3, k2, nc, n0, n0, t.ld, k2});
                    u3.entries.push_back({o2.up, w3, t.x2, n1, nc, k2, n1, k2, t.ld});
                }
            } else {
                // X1 -= V2 (U2'^T X2); X1 += V2 (Uo^T X1); X2 -= V1 (U1'^T X1)
                if (k2 > 0) {
                    c1.entries.push_back({o2.up, t.x2, w1, k2, nc, n1, n1, t.ld, k2});
                    u1.entries.push_back({o2.v, w1, t.x1, n0, nc, k2, n0, k2, t.ld});
                    c2.entries.push_back({nd.uo, t.x1, w2, k2, nc, n0, n0, t.ld, k2});
                    u2.entries.push_back({o2.v, w2, t.x1, n0, nc, k2, n0, k2, t.ld});
                }
                if (k1 > 0) {
                    c3.entries.push_back({o1.up, t.x1, w3, k1, nc, n0, n0, t.ld, k1});
                    u3.entries.push_back({o1.v, w3, t.x2, n1, nc, k1, n1, k1, t.ld});
                }
            }
            stats_.flops += kFlopScale * 2.0 * nc *
                            (static_cast<double>(k1) * (n0 + n1) + static_cast<double>(k2) * (3.0 * n0 + n1));
        }
        meta_.clear();
        for (VBatch<T>* b : {&c1, &u1, &c2, &u2, &c3, &u3}) b->stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        const T one(1.0), zero(0.0), minus(-1.0);
        c1.gemm(md, MagmaTrans, MagmaNoTrans, one, zero, queue);
        u1.gemm(md, MagmaNoTrans, MagmaNoTrans, minus, one, queue);
        c2.gemm(md, MagmaTrans, MagmaNoTrans, one, zero, queue);
        u2.gemm(md, MagmaNoTrans, MagmaNoTrans, one, one, queue);
        c3.gemm(md, MagmaTrans, MagmaNoTrans, one, zero, queue);
        u3.gemm(md, MagmaNoTrans, MagmaNoTrans, minus, one, queue);
        heap.free(work);  // (reused by later launches only after these: one stream)
    }

    int64_t n_loc_ = 0;
    int maxlevel_ = 0;
    std::vector<std::vector<Block>> blocks_;  // [level][block], by first row
    std::vector<std::vector<Node>> nodes_;    // [level][node], by first row
    std::vector<Leaf> leaves_;
    std::vector<const T*> lu_source_;         // factor_lu: unfactored copies of in-place LUs
    char* store_ = nullptr;
    bool ready_ = false;
    Stats stats_;
    MetaBuilder meta_;
    DeviceBuffer meta_device_;
    DeviceBuffer getrf_work_;
};

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
