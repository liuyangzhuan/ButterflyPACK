#pragma once
// Symmetric HODLR factorization and solve on the GPU (option%sym=1), the
// device form of HODLR_factorization_sym and HODLR_Sym_Inv_Apply
// (SRC/BPACK_factor.f90, SRC/BPACK_solve_mul.f90), operation for operation,
// with the dense leaves factored by LU with partial pivoting instead of
// Bunch-Kaufman (MAGMA has no batched Bunch-Kaufman).
//
// A node of level l stores A21 = U1 V0^T (HodlrDevice, flagged sym); its
// bases are Z0 = V0 (rows of child 0) and Z1 = U1 (rows of child 1), and
// Q = A_child^{-1} Z is built by applying every leaf solve and every
// descendant node's correction to all ancestor bases.  A node then forms
// G0 = Z0^T Q0, G1 = Z1^T Q1 (symmetrized), S = I - G0 G1 and its LU.  The
// work of a level runs as one set of batched launches.  The factors stay on
// the device for the solve.
//
// With several ranks, a rank holds of the node of a shared level (one per
// level) the rows of Z0 and Q0 or of Z1 and Q1 it owns.  As the CPU does,
// G0 and G1 are summed over the node's ranks and S is factored on the node's
// head (rank 0 of its communicator); a node correction sums c0 and c1 over
// the ranks, the head solves, and [delta; gamma] goes back to all of them.
// The log-determinant here is this rank's part (shared nodes on the head).

#ifdef H2_HAVE_GPU

#include "hodlr_batch.hpp"
#include "hodlr_device.hpp"
#include "hodlr_lu.hpp"

#include <chrono>
#include <cmath>
#include <string>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace bpack {
namespace gpu {

template<typename T>
class HodlrSymFactor {
public:
    struct Stats {
        double logabsdet = 0.0;
        std::complex<double> phase{1.0, 0.0};
        double flops = 0.0;
        double mbytes = 0.0;  // resident factor data (beyond the forward blocks)
        double t_leaf = 0.0, t_prop = 0.0, t_nodes = 0.0;
        int jitter_leaves = 0, jitter_nodes = 0;
        double max_jitter = 0.0;
        int maxrank = 0;
    };

    HodlrSymFactor() = default;
    HodlrSymFactor(const HodlrSymFactor&) = delete;
    HodlrSymFactor& operator=(const HodlrSymFactor&) = delete;
    ~HodlrSymFactor() { release(); }

    bool ready() const { return ready_; }
    const Stats& stats() const { return stats_; }

    void factor(const HodlrDevice<T>& dev, double jitter_option) {
        using clock = std::chrono::steady_clock;
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        release();
        stats_ = Stats{};
        build_topology(dev);
        allocate();

        // Q0 = Z0, Q1 = Z1, LU = D
        {
            ItemList<AxpbyItem<T>> copies;
            for (const auto& level : nodes_) {
                for (const Node& nd : level) {
                    if (nd.k == 0) continue;
                    copies.add({nd.q0, nd.n0, nd.z0, nd.n0, nullptr, 0, nd.n0, nd.k, T(1.0), T(0.0)}, nd.n0, nd.k);
                    copies.add({nd.q1, nd.n1, nd.z1, nd.n1, nullptr, 0, nd.n1, nd.k, T(1.0), T(0.0)}, nd.n1, nd.k);
                }
            }
            for (const Leaf& lf : leaves_) copies.add({lf.lu, lf.m, lf.d, lf.m, nullptr, 0, lf.m, lf.m, T(1.0), T(0.0)}, lf.m, lf.m);
            meta_.clear();
            copies.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            launch_axpby(copies.device(meta), copies.count(), copies.max_m, copies.max_n, stream);
        }

        // dense leaves: LU (with the CPU's jitter retries), then their solves on all ancestor bases
        auto t0 = clock::now();
        factor_leaves(jitter_option);
        auto t1 = clock::now();
        {
            std::vector<LeafTarget> targets;
            for (const Leaf& lf : leaves_) {
                for (int al = 1; al <= maxlevel_; ++al) {
                    int side = 0;
                    const Node* anc = find_node(al, lf.row0, lf.row0 + lf.m, side);
                    if (anc->k == 0) continue;
                    targets.push_back(slice(lf, *anc, side));
                }
            }
            apply_leaves(targets);
            check_cuda(cudaStreamSynchronize(stream), "HODLR sym leaf propagation");
        }
        auto t2 = clock::now();

        // nodes, bottom up
        for (int l = maxlevel_; l >= 1; --l) {
            Node* shared = shared_node(l);
            if (shared != nullptr) {
                factor_shared(*shared, jitter_option);
            } else {
                factor_level(l, jitter_option);
            }
            if (l > 1) {
                std::vector<NodeTarget> targets;
                for (const Node& nd : nodes_[l]) {
                    if (nd.k == 0) continue;
                    // (the rows of the children this rank holds: one child on a shared level)
                    int64_t lo = std::min(nd.c0, nd.c1), hi = std::max(nd.c0 + nd.n0, nd.c1 + nd.n1);
                    if (nd.n0 == 0) {
                        lo = nd.c1;
                        hi = nd.c1 + nd.n1;
                    } else if (nd.n1 == 0) {
                        lo = nd.c0;
                        hi = nd.c0 + nd.n0;
                    }
                    for (int al = 1; al < l; ++al) {
                        int side = 0;
                        const Node* anc = find_node(al, lo, hi, side);
                        if (anc->k == 0) continue;
                        T* q = side == 0 ? anc->q0 : anc->q1;
                        const int64_t start = side == 0 ? anc->c0 : anc->c1;
                        const int ld = side == 0 ? anc->n0 : anc->n1;
                        targets.push_back({&nd, q + (nd.c0 - start), q + (nd.c1 - start), ld, anc->k});
                    }
                }
                if (shared != nullptr) {
                    apply_shared(*shared, targets);
                } else {
                    apply_nodes(targets);
                }
            }
        }
        check_cuda(cudaStreamSynchronize(stream), "HODLR sym factorization");
        auto t3 = clock::now();
        stats_.t_leaf = std::chrono::duration<double>(t1 - t0).count();
        stats_.t_prop = std::chrono::duration<double>(t2 - t1).count();
        stats_.t_nodes = std::chrono::duration<double>(t3 - t2).count();
        ready_ = true;
    }

    // y = A^{-1} x for the local n_loc x nrhs host arrays x and y
    void solve(int nrhs, const T* x, T* y) {
        if (!ready_) throw std::logic_error("HODLR GPU: symmetric solve before the factorization");
        const double flops0 = stats_.flops;  // (the applications count there; the factorization's count is kept)
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        const size_t bytes = static_cast<size_t>(n_loc_) * nrhs * sizeof(T);
        DeviceHeap& heap = DeviceHeap::instance();
        T* dv = heap.alloc<T>(bytes);
        check_cuda(cudaMemcpyAsync(dv, x, bytes, cudaMemcpyHostToDevice, stream), "HODLR solve input");
        const int ld = static_cast<int>(n_loc_);
        std::vector<LeafTarget> leaf_targets;
        for (const Leaf& lf : leaves_) leaf_targets.push_back({&lf, dv + lf.row0, ld, nrhs});
        apply_leaves(leaf_targets);
        for (int l = maxlevel_; l >= 1; --l) {
            std::vector<NodeTarget> targets;
            for (const Node& nd : nodes_[l]) {
                if (nd.k > 0) targets.push_back({&nd, dv + nd.c0, dv + nd.c1, ld, nrhs});
            }
            Node* shared = shared_node(l);
            if (shared != nullptr) {
                apply_shared(*shared, targets);
            } else {
                apply_nodes(targets);
            }
        }
        check_cuda(cudaMemcpyAsync(y, dv, bytes, cudaMemcpyDeviceToHost, stream), "HODLR solve output");
        check_cuda(cudaStreamSynchronize(stream), "HODLR sym solve");
        heap.free(dv);
        solve_flops_ = stats_.flops - flops0;
        stats_.flops = flops0;
    }
    double solve_flops() const { return solve_flops_; }  // of the last solve()

    // Copies of the factor of the node of level `level` whose child 0 starts
    // at local row c0 (for the comparison with the CPU); false if none.
    // n0, n1, k: the caller's sizes, checked before copying
    // (on a shared level c0 is the first local row of whichever child this
    // rank holds rows of; g0, g1, s, ipiv are the head's)
    bool download_node(int level, int64_t c0, int n0, int n1, T* q0, T* q1, T* g0, T* g1, T* s, int* ipiv,
                       int* k, T* z0 = nullptr, T* z1 = nullptr) const {
        if (!ready_ || level < 1 || level > maxlevel_) return false;
        for (const Node& nd : nodes_[level]) {
            const int64_t start = nd.comm != MPI_COMM_NULL && nd.n0 == 0 ? nd.c1 : nd.c0;
            if (start != c0) continue;
            if (nd.n0 != n0 || nd.n1 != n1 || nd.k != *k) {
                std::ostringstream oss;
                oss << "HODLR GPU: node of level " << level << " at row " << c0 << " has sizes " << nd.n0 << ", "
                    << nd.n1 << ", rank " << nd.k << " on the device but " << n0 << ", " << n1 << ", " << *k
                    << " on the host";
                throw std::logic_error(oss.str());
            }
            if (nd.k == 0) return true;
            const size_t kk = static_cast<size_t>(nd.k) * nd.k * sizeof(T);
            check_cuda(cudaMemcpy(q0, nd.q0, static_cast<size_t>(nd.n0) * nd.k * sizeof(T), cudaMemcpyDeviceToHost), "download");
            check_cuda(cudaMemcpy(q1, nd.q1, static_cast<size_t>(nd.n1) * nd.k * sizeof(T), cudaMemcpyDeviceToHost), "download");
            check_cuda(cudaMemcpy(g0, nd.g0, kk, cudaMemcpyDeviceToHost), "download");
            check_cuda(cudaMemcpy(g1, nd.g1, kk, cudaMemcpyDeviceToHost), "download");
            check_cuda(cudaMemcpy(s, nd.s, kk, cudaMemcpyDeviceToHost), "download");
            check_cuda(cudaMemcpy(ipiv, nd.ipiv, static_cast<size_t>(nd.k) * sizeof(int), cudaMemcpyDeviceToHost), "download");
            if (z0 != nullptr)
                check_cuda(cudaMemcpy(z0, nd.z0, static_cast<size_t>(nd.n0) * nd.k * sizeof(T), cudaMemcpyDeviceToHost), "download");
            if (z1 != nullptr)
                check_cuda(cudaMemcpy(z1, nd.z1, static_cast<size_t>(nd.n1) * nd.k * sizeof(T), cudaMemcpyDeviceToHost), "download");
            return true;
        }
        return false;
    }

private:
    double solve_flops_ = 0.0;
    struct Node {
        int64_t c0, c1;  // first rows of child 0 and child 1
        int n0, n1, k;
        const T* z0;     // n0 x k (forward V0)
        const T* z1;     // n1 x k (forward U1)
        MPI_Comm comm;   // a shared node's ranks (MPI_COMM_NULL: this rank's alone)
        bool head;       // this rank factors S (rank 0 of comm, or the node is local)
        T* q0;
        T* q1;
        T *g0, *g1, *s;  // k x k
        int* ipiv;       // k
    };
    struct Leaf {
        int64_t row0;
        int m;
        const T* d;  // forward block
        T* lu;
        int* ipiv;
    };
    struct LeafTarget {
        const Leaf* leaf;
        T* x;  // m x ncols
        int ld, ncols;
    };
    struct NodeTarget {
        const Node* node;
        T* x0;  // n0 x ncols (rows of child 0)
        T* x1;  // n1 x ncols
        int ld, ncols;
    };

    static constexpr double kFlopScale = fmm::gpu::is_complex_scalar<T> ? 4.0 : 1.0;

    static size_t aligned(size_t bytes) { return fmm::gpu::align_up(std::max<size_t>(bytes, 1)); }

    void build_topology(const HodlrDevice<T>& dev) {
        n_loc_ = dev.n_loc();
        maxlevel_ = dev.maxlevel();
        nodes_.assign(static_cast<size_t>(maxlevel_ + 1), {});
        for (const auto& b : dev.blocks()) {
            if (!b.sym) throw std::logic_error("HODLR GPU: the symmetric factorization needs symmetric blocks");
            MPI_Comm comm = dev.level_comm(b.level);
            bool head = true;
            if (comm != MPI_COMM_NULL) {
                int rank = 0;
                MPI_Comm_rank(comm, &rank);
                head = rank == 0;
            }
            nodes_[static_cast<size_t>(b.level)].push_back(Node{b.col0, b.row0, b.n, b.m, b.k, b.v, b.u, comm, head,
                                                                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});
        }
        for (int l = 1; l <= maxlevel_; ++l) {
            auto& v = nodes_[static_cast<size_t>(l)];
            if (dev.shared_level(l) && v.size() != 1) {
                throw std::logic_error("HODLR GPU: a shared level must hold one node on each rank");
            }
            if (!dev.distributed() && v.size() != (size_t{1} << (l - 1))) {
                throw std::logic_error("HODLR GPU: level " + std::to_string(l) + " has " + std::to_string(v.size()) +
                                       " nodes; the one-rank backend needs the whole tree");
            }
            std::sort(v.begin(), v.end(), [](const Node& a, const Node& b) { return a.c0 < b.c0; });
            for (const Node& nd : v) stats_.maxrank = std::max(stats_.maxrank, nd.k);
        }
        leaves_.clear();
        for (const auto& lf : dev.leaves()) leaves_.push_back(Leaf{lf.row0, lf.m, lf.d, nullptr, nullptr});
        std::sort(leaves_.begin(), leaves_.end(), [](const Leaf& a, const Leaf& b) { return a.row0 < b.row0; });
    }

    void allocate() {
        size_t bytes = 0;
        for (const auto& level : nodes_) {
            for (const Node& nd : level) {
                const size_t k = static_cast<size_t>(nd.k);
                bytes += aligned(nd.n0 * k * sizeof(T)) + aligned(nd.n1 * k * sizeof(T)) + 3 * aligned(k * k * sizeof(T)) +
                         aligned(k * sizeof(int));  // as carved below
            }
        }
        for (const Leaf& lf : leaves_) {
            const size_t m = static_cast<size_t>(lf.m);
            bytes += aligned(m * m * sizeof(T)) + aligned(m * sizeof(int));
        }
        DeviceHeap& heap = DeviceHeap::instance();
        store_ = heap.alloc_resident<char>(bytes);
        stats_.mbytes = static_cast<double>(bytes) / 1.0e6;
        char* at = store_;
        auto take = [&](size_t b) {
            char* p = at;
            at += aligned(b);
            return p;
        };
        for (auto& level : nodes_) {
            for (Node& nd : level) {
                const size_t k = static_cast<size_t>(nd.k);
                nd.q0 = reinterpret_cast<T*>(take(nd.n0 * k * sizeof(T)));
                nd.q1 = reinterpret_cast<T*>(take(nd.n1 * k * sizeof(T)));
                nd.g0 = reinterpret_cast<T*>(take(k * k * sizeof(T)));
                nd.g1 = reinterpret_cast<T*>(take(k * k * sizeof(T)));
                nd.s = reinterpret_cast<T*>(take(k * k * sizeof(T)));
                nd.ipiv = reinterpret_cast<int*>(take(k * sizeof(int)));
            }
        }
        for (Leaf& lf : leaves_) {
            const size_t m = static_cast<size_t>(lf.m);
            lf.lu = reinterpret_cast<T*>(take(m * m * sizeof(T)));
            lf.ipiv = reinterpret_cast<int*>(take(m * sizeof(int)));
        }
        if (at != store_ + bytes) throw std::logic_error("HODLR GPU: symmetric factors placed outside their allocation");
    }

    void release() {
        if (store_ != nullptr) {
            Context::instance().activate();
            DeviceHeap::instance().free(store_);
            store_ = nullptr;
        }
        nodes_.clear();
        leaves_.clear();
        ready_ = false;
    }

    // the node of level al whose child 0 or child 1 holds rows [lo, hi)
    const Node* find_node(int al, int64_t lo, int64_t hi, int& side) const {
        const auto& v = nodes_[static_cast<size_t>(al)];
        auto it = std::upper_bound(v.begin(), v.end(), lo, [](int64_t r, const Node& nd) { return r < nd.c0; });
        if (it != v.begin()) {
            const Node& nd = *std::prev(it);
            if (lo >= nd.c0 && hi <= nd.c0 + nd.n0) {
                side = 0;
                return &nd;
            }
        }
        for (const Node& nd : v) {  // child 1 precedes child 0 only in unusual trees
            if (lo >= nd.c0 && hi <= nd.c0 + nd.n0) {
                side = 0;
                return &nd;
            }
            if (lo >= nd.c1 && hi <= nd.c1 + nd.n1) {
                side = 1;
                return &nd;
            }
        }
        std::ostringstream oss;
        oss << "HODLR GPU: no node of level " << al << " holds rows [" << lo << ", " << hi << ")";
        throw std::logic_error(oss.str());
    }

    LeafTarget slice(const Leaf& lf, const Node& anc, int side) const {
        T* q = side == 0 ? anc.q0 : anc.q1;
        const int64_t start = side == 0 ? anc.c0 : anc.c1;
        const int ld = side == 0 ? anc.n0 : anc.n1;
        return {&lf, q + (lf.row0 - start), ld, anc.k};
    }

    // LU of the leaves; a failed leaf is refactored with a growing diagonal
    // jitter, as the CPU does, and the log-determinant is accumulated.
    void factor_leaves(double jitter_option) {
        const cudaStream_t stream = Context::instance().stream();
        DeviceHeap& heap = DeviceHeap::instance();
        const int count = static_cast<int>(leaves_.size());
        std::vector<int> todo(static_cast<size_t>(count));
        for (int i = 0; i < count; ++i) todo[static_cast<size_t>(i)] = i;
        std::vector<double> norms;  // of the leaves in todo, from the first retry on
        for (int attempt = 0; attempt <= 8 && !todo.empty(); ++attempt) {
            if (attempt == 1) {  // Frobenius norms of the failed leaves, once
                ItemList<NormItem<T>> norm_items;
                double* d_norms = heap.alloc<double>(sizeof(double) * todo.size());
                for (size_t t = 0; t < todo.size(); ++t) {
                    const Leaf& lf = leaves_[static_cast<size_t>(todo[t])];
                    norm_items.add({lf.d, lf.m, lf.m, lf.m, d_norms + t}, lf.m, lf.m);
                }
                meta_.clear();
                norm_items.stage(meta_);
                char* meta = meta_.upload(meta_device_, stream);
                launch_fnorm(norm_items.device(meta), norm_items.count(), stream);
                norms.resize(todo.size());
                check_cuda(cudaMemcpyAsync(norms.data(), d_norms, sizeof(double) * todo.size(), cudaMemcpyDeviceToHost,
                                           stream),
                           "leaf norms");
                check_cuda(cudaStreamSynchronize(stream), "leaf norms");
                heap.free(d_norms);
            }
            GetrfBatch<T> lu;
            ItemList<AxpbyItem<T>> restore;
            ItemList<DiagItem<T>> jitter_items;
            for (size_t t = 0; t < todo.size(); ++t) {
                const Leaf& lf = leaves_[static_cast<size_t>(todo[t])];
                if (attempt > 0) {
                    const double scale = std::max(norms[t], 1.0);
                    const double base = std::max(jitter_option * scale, std::numeric_limits<double>::epsilon() * scale);
                    const double jitter = base * std::pow(10.0, attempt - 1);
                    restore.add({lf.lu, lf.m, lf.d, lf.m, nullptr, 0, lf.m, lf.m, T(1.0), T(0.0)}, lf.m, lf.m);
                    jitter_items.add({lf.lu, lf.m, lf.m, T(jitter)}, lf.m, lf.m);
                    stats_.max_jitter = std::max(stats_.max_jitter, jitter);
                }
                lu.add(lf.lu, lf.m, lf.m, lf.ipiv);
                stats_.flops += kFlopScale * 2.0 / 3.0 * std::pow(static_cast<double>(lf.m), 3);
            }
            meta_.clear();
            lu.stage(meta_);
            restore.stage(meta_);
            jitter_items.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            if (attempt > 0) {
                launch_axpby(restore.device(meta), restore.count(), restore.max_m, restore.max_n, stream);
                launch_add_diagonal(jitter_items.device(meta), jitter_items.count(), jitter_items.max_n, stream);
            }
            std::vector<int> failed = run_getrf(lu, meta);
            std::vector<int> next;
            std::vector<double> kept;
            for (int f : failed) {
                next.push_back(todo[static_cast<size_t>(f)]);
                if (attempt > 0) kept.push_back(norms[static_cast<size_t>(f)]);
            }
            if (attempt > 0) {
                stats_.jitter_leaves += static_cast<int>(todo.size() - next.size());
                norms.swap(kept);
            }
            todo.swap(next);
        }
        if (!todo.empty()) throw std::runtime_error("symmetric HODLR LU failed on a dense leaf (GPU)");
        // log-determinant, in leaf order as the CPU accumulates it
        std::vector<const T*> mats;
        std::vector<const int*> pivs;
        std::vector<int> sizes;
        for (const Leaf& lf : leaves_) {
            mats.push_back(lf.lu);
            pivs.push_back(lf.ipiv);
            sizes.push_back(lf.m);
        }
        accumulate_logdet(mats, pivs, sizes, "dense leaf");
    }

    std::vector<int> run_getrf(const GetrfBatch<T>& lu, char* meta) { return bpack::gpu::run_getrf(lu, meta, getrf_work_); }

    void accumulate_logdet(const std::vector<const T*>& mats, const std::vector<const int*>& pivs,
                           const std::vector<int>& sizes, const char* what) {
        std::vector<T> diag;
        std::vector<int> piv;
        lu_diagonals(mats, pivs, sizes, meta_, meta_device_, diag, piv, what);
        size_t off = 0;
        for (size_t i = 0; i < mats.size(); ++i) {
            std::complex<double> phase;
            double logabs = 0.0;
            if (!lu_slogdet(diag.data() + off, piv.data() + off, sizes[i], phase, logabs)) {
                throw std::runtime_error(std::string("symmetric HODLR: zero pivot in a ") + what + " (GPU)");
            }
            stats_.phase *= phase;
            stats_.logabsdet += logabs;
            off += static_cast<size_t>(sizes[i]);
        }
    }

    // G0, G1, S = I - G0 G1 and the LU of S for the nodes of level l
    void factor_level(int l, double jitter_option) {
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        std::vector<const Node*> active;
        for (const Node& nd : nodes_[static_cast<size_t>(l)]) {
            if (nd.k > 0) active.push_back(&nd);
        }
        if (active.empty()) return;
        {
            VBatch<T> g;
            ItemList<SquareItem<T>> sym;
            for (const Node* nd : active) {
                g.entries.push_back({nd->z0, nd->q0, nd->g0, nd->k, nd->k, nd->n0, nd->n0, nd->n0, nd->k});
                g.entries.push_back({nd->z1, nd->q1, nd->g1, nd->k, nd->k, nd->n1, nd->n1, nd->n1, nd->k});
                sym.add({nd->g0, nd->k, nd->k}, nd->k, nd->k);
                sym.add({nd->g1, nd->k, nd->k}, nd->k, nd->k);
                stats_.flops += kFlopScale * 2.0 * nd->k * nd->k * static_cast<double>(nd->n0 + nd->n1);
            }
            meta_.clear();
            g.stage(meta_);
            sym.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            g.gemm(meta, MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
            launch_symmetrize(sym.device(meta), sym.count(), sym.max_n, stream);
        }
        factor_s(active, jitter_option);
    }

    // S = I - G0 G1 and its LU for the nodes (with the CPU's jitter
    // retries), and their log-determinants
    void factor_s(const std::vector<const Node*>& active, double jitter_option) {
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        const double jitter_base = std::max(jitter_option, std::numeric_limits<double>::epsilon());
        std::vector<const Node*> todo = active;
        for (int attempt = 0; attempt <= 8 && !todo.empty(); ++attempt) {
            VBatch<T> s;
            ItemList<DiagItem<T>> one, jit;
            GetrfBatch<T> lu;
            const double jitter = attempt > 0 ? jitter_base * std::pow(10.0, attempt - 1) : 0.0;
            for (const Node* nd : todo) {
                s.entries.push_back({nd->g0, nd->g1, nd->s, nd->k, nd->k, nd->k, nd->k, nd->k, nd->k});
                one.add({nd->s, nd->k, nd->k, T(1.0)}, nd->k, nd->k);
                if (attempt > 0) jit.add({nd->s, nd->k, nd->k, T(jitter)}, nd->k, nd->k);
                lu.add(nd->s, nd->k, nd->k, nd->ipiv);
                stats_.flops += kFlopScale * (2.0 + 2.0 / 3.0) * std::pow(static_cast<double>(nd->k), 3);
            }
            meta_.clear();
            s.stage(meta_);
            one.stage(meta_);
            jit.stage(meta_);
            lu.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            s.gemm(meta, MagmaNoTrans, MagmaNoTrans, T(-1.0), T(0.0), queue);
            launch_add_diagonal(one.device(meta), one.count(), one.max_n, stream);
            if (attempt > 0) launch_add_diagonal(jit.device(meta), jit.count(), jit.max_n, stream);
            std::vector<int> failed = run_getrf(lu, meta);
            std::vector<const Node*> next;
            for (int f : failed) next.push_back(todo[static_cast<size_t>(f)]);
            if (attempt > 0) {
                stats_.jitter_nodes += static_cast<int>(todo.size() - next.size());
                stats_.max_jitter = std::max(stats_.max_jitter, jitter);
            }
            todo.swap(next);
        }
        if (!todo.empty()) throw std::runtime_error("symmetric HODLR LU failed for a rank-r Schur correction (GPU)");
        std::vector<const T*> mats;
        std::vector<const int*> pivs;
        std::vector<int> sizes;
        for (const Node* nd : active) {
            mats.push_back(nd->s);
            pivs.push_back(nd->ipiv);
            sizes.push_back(nd->k);
        }
        accumulate_logdet(mats, pivs, sizes, "Schur correction");
    }

    // the node of level l if the level is shared (its only node on this rank)
    Node* shared_node(int l) {
        auto& v = nodes_[static_cast<size_t>(l)];
        return v.size() == 1 && v[0].comm != MPI_COMM_NULL ? &v[0] : nullptr;
    }

    // G0, G1 of a shared node, summed over its ranks (one message of both,
    // in the exchange arena when MPI is CUDA-aware, else staged through the
    // host), then S and its LU on the head (whose success all ranks learn)
    void factor_shared(Node& nd, double jitter_option) {
        if (nd.k == 0) return;
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const size_t kk = static_cast<size_t>(nd.k) * nd.k;
        DeviceHeap& heap = DeviceHeap::instance();
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        T* buf = nullptr;
        if (fmm::gpu::device_exchange_enabled()) buf = reinterpret_cast<T*>(arena.try_alloc(2 * kk * sizeof(T)));
        const bool in_arena = buf != nullptr;
        if (buf == nullptr) buf = heap.alloc<T>(2 * kk * sizeof(T));
        check_cuda(cudaMemsetAsync(buf, 0, 2 * kk * sizeof(T), stream), "HODLR sym G0, G1");
        {
            VBatch<T> g;
            if (nd.n0 > 0) g.entries.push_back({nd.z0, nd.q0, buf, nd.k, nd.k, nd.n0, nd.n0, nd.n0, nd.k});
            if (nd.n1 > 0) g.entries.push_back({nd.z1, nd.q1, buf + kk, nd.k, nd.k, nd.n1, nd.n1, nd.n1, nd.k});
            stats_.flops += kFlopScale * 2.0 * nd.k * nd.k * static_cast<double>(nd.n0 + nd.n1);
            meta_.clear();
            g.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            if (g.count() > 0) g.gemm(meta, MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), ctx.queue());
        }
        HodlrDevice<T>::sum_over(nd.comm, buf, static_cast<int64_t>(2 * kk), in_arena, stream);
        check_cuda(cudaMemcpyAsync(nd.g0, buf, kk * sizeof(T), cudaMemcpyDeviceToDevice, stream), "HODLR sym G0");
        check_cuda(cudaMemcpyAsync(nd.g1, buf + kk, kk * sizeof(T), cudaMemcpyDeviceToDevice, stream), "HODLR sym G1");
        check_cuda(cudaStreamSynchronize(stream), "HODLR sym G0, G1");
        if (in_arena) {
            arena.free(buf);
        } else {
            heap.free(buf);
        }
        {
            ItemList<SquareItem<T>> sym;
            sym.add({nd.g0, nd.k, nd.k}, nd.k, nd.k);
            sym.add({nd.g1, nd.k, nd.k}, nd.k, nd.k);
            meta_.clear();
            sym.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            launch_symmetrize(sym.device(meta), sym.count(), sym.max_n, stream);
        }
        int ok = 1;
        if (nd.head) {
            try {
                factor_s({&nd}, jitter_option);
            } catch (const std::runtime_error&) {
                ok = 0;
            }
        }
        MPI_Bcast(&ok, 1, MPI_INT, 0, nd.comm);
        if (!ok) throw std::runtime_error("symmetric HODLR LU failed for a rank-r Schur correction (GPU)");
    }

    // the node corrections of a shared node on its targets
    // (HODLR_Sym_Node_Apply with several ranks): the partial c0 = Z0^T X0,
    // c1 = Z1^T X1 summed over the node's ranks, delta and gamma on the head
    // (as apply_nodes), [c; delta] broadcast, then the local rows updated
    void apply_shared(Node& nd, const std::vector<NodeTarget>& targets) {
        if (targets.empty() || nd.k == 0) return;
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        const int k = nd.k;
        size_t cols = 0;
        for (const NodeTarget& t : targets) cols += static_cast<size_t>(t.ncols);
        // [c0 c1 per target | delta per target] (exchanged), then w per target
        const size_t ex = 3 * static_cast<size_t>(k) * cols, all = ex + static_cast<size_t>(k) * cols;
        DeviceHeap& heap = DeviceHeap::instance();
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        T* buf = nullptr;
        if (fmm::gpu::device_exchange_enabled()) buf = reinterpret_cast<T*>(arena.try_alloc(ex * sizeof(T)));
        const bool in_arena = buf != nullptr;
        if (buf == nullptr) buf = heap.alloc<T>(ex * sizeof(T));
        T* wbuf = heap.alloc<T>((all - ex) * sizeof(T));
        check_cuda(cudaMemsetAsync(buf, 0, ex * sizeof(T), stream), "HODLR sym node coefficients");
        std::vector<T*> c0s, c1s, ds, ws;
        {
            size_t cat = 0, dat = 2 * static_cast<size_t>(k) * cols, wat = 0;
            for (const NodeTarget& t : targets) {
                const size_t kc = static_cast<size_t>(k) * t.ncols;
                c0s.push_back(buf + cat);
                c1s.push_back(buf + cat + kc);
                cat += 2 * kc;
                ds.push_back(buf + dat);
                dat += kc;
                ws.push_back(wbuf + wat);
                wat += kc;
            }
        }
        {
            VBatch<T> b_c;
            for (size_t i = 0; i < targets.size(); ++i) {
                const NodeTarget& t = targets[i];
                if (nd.n0 > 0) b_c.entries.push_back({nd.z0, t.x0, c0s[i], k, t.ncols, nd.n0, nd.n0, t.ld, k});
                if (nd.n1 > 0) b_c.entries.push_back({nd.z1, t.x1, c1s[i], k, t.ncols, nd.n1, nd.n1, t.ld, k});
                stats_.flops += kFlopScale * 2.0 * k * static_cast<double>(nd.n0 + nd.n1) * t.ncols;
            }
            meta_.clear();
            b_c.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            if (b_c.count() > 0) b_c.gemm(meta, MagmaTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
        }
        HodlrDevice<T>::sum_over(nd.comm, buf, static_cast<int64_t>(2 * static_cast<size_t>(k) * cols), in_arena, stream);
        if (nd.head) {
            VBatch<T> b_w, b_d, b_gamma;
            ItemList<AxpbyItem<T>> copy_w, sub_c0;
            ItemList<RowSwapItem<T>> swaps;
            TrsmBatch<T> lower, upper;
            for (size_t i = 0; i < targets.size(); ++i) {
                const int nc = targets[i].ncols;
                copy_w.add({ws[i], k, c1s[i], k, nullptr, 0, k, nc, T(1.0), T(0.0)}, k, nc);
                b_w.entries.push_back({nd.g1, c0s[i], ws[i], k, nc, k, k, k, k});
                b_d.entries.push_back({nd.g0, ws[i], ds[i], k, nc, k, k, k, k});
                swaps.add({ds[i], k, k, nc, nd.ipiv}, k, nc);
                lower.add(nd.s, k, ds[i], k, k, nc);
                upper.add(nd.s, k, ds[i], k, k, nc);
                sub_c0.add({ds[i], k, ds[i], k, c0s[i], k, k, nc, T(1.0), T(-1.0)}, k, nc);
                b_gamma.entries.push_back({nd.g1, ds[i], c1s[i], k, nc, k, k, k, k});
                stats_.flops += kFlopScale * 6.0 * k * static_cast<double>(k) * nc;
            }
            meta_.clear();
            copy_w.stage(meta_);
            b_w.stage(meta_);
            b_d.stage(meta_);
            swaps.stage(meta_);
            lower.stage(meta_);
            upper.stage(meta_);
            sub_c0.stage(meta_);
            b_gamma.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            const T one(1.0), zero(0.0), minus(-1.0);
            launch_axpby(copy_w.device(meta), copy_w.count(), copy_w.max_m, copy_w.max_n, stream);
            b_w.gemm(meta, MagmaNoTrans, MagmaNoTrans, minus, one, queue);
            b_d.gemm(meta, MagmaNoTrans, MagmaNoTrans, one, zero, queue);
            launch_row_swaps(swaps.device(meta), swaps.count(), swaps.max_n, stream);
            lower.solve(meta, MagmaLower, MagmaUnit, queue);
            upper.solve(meta, MagmaUpper, MagmaNonUnit, queue);
            launch_axpby(sub_c0.device(meta), sub_c0.count(), sub_c0.max_m, sub_c0.max_n, stream);
            b_gamma.gemm(meta, MagmaNoTrans, MagmaNoTrans, one, one, queue);
        }
        HodlrDevice<T>::bcast_from_root(nd.comm, buf, static_cast<int64_t>(ex), in_arena, stream);
        {
            VBatch<T> b_x0, b_x1;
            for (size_t i = 0; i < targets.size(); ++i) {
                const NodeTarget& t = targets[i];
                if (nd.n0 > 0) b_x0.entries.push_back({nd.q0, c1s[i], t.x0, nd.n0, t.ncols, k, nd.n0, k, t.ld});
                if (nd.n1 > 0) b_x1.entries.push_back({nd.q1, ds[i], t.x1, nd.n1, t.ncols, k, nd.n1, k, t.ld});
                stats_.flops += kFlopScale * 2.0 * k * static_cast<double>(nd.n0 + nd.n1) * t.ncols;
            }
            meta_.clear();
            b_x0.stage(meta_);
            b_x1.stage(meta_);
            char* meta = meta_.upload(meta_device_, stream);
            const T one(1.0), minus(-1.0);
            if (b_x0.count() > 0) b_x0.gemm(meta, MagmaNoTrans, MagmaNoTrans, minus, one, queue);
            if (b_x1.count() > 0) b_x1.gemm(meta, MagmaNoTrans, MagmaNoTrans, one, one, queue);
        }
        check_cuda(cudaStreamSynchronize(stream), "HODLR sym shared node");
        heap.free(wbuf);
        if (in_arena) {
            arena.free(buf);
        } else {
            heap.free(buf);
        }
    }

    // X = D^{-1} X for each target (HODLR_Sym_Leaf_Apply)
    void apply_leaves(const std::vector<LeafTarget>& targets) {
        if (targets.empty()) return;
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        ItemList<RowSwapItem<T>> swaps;
        TrsmBatch<T> lower, upper;
        for (const LeafTarget& t : targets) {
            const Leaf& lf = *t.leaf;
            swaps.add({t.x, t.ld, lf.m, t.ncols, lf.ipiv}, lf.m, t.ncols);
            lower.add(lf.lu, lf.m, t.x, t.ld, lf.m, t.ncols);
            upper.add(lf.lu, lf.m, t.x, t.ld, lf.m, t.ncols);
            stats_.flops += kFlopScale * 2.0 * lf.m * static_cast<double>(lf.m) * t.ncols;
        }
        meta_.clear();
        swaps.stage(meta_);
        lower.stage(meta_);
        upper.stage(meta_);
        char* meta = meta_.upload(meta_device_, stream);
        launch_row_swaps(swaps.device(meta), swaps.count(), swaps.max_n, stream);
        lower.solve(meta, MagmaLower, MagmaUnit, ctx.queue());
        upper.solve(meta, MagmaUpper, MagmaNonUnit, ctx.queue());
    }

    // the node correction on each target (HODLR_Sym_Node_Apply):
    //   c0 = Z0^T X0, c1 = Z1^T X1, delta = S^{-1} G0 (c1 - G1 c0) - c0,
    //   gamma = c1 + G1 delta, X0 -= Q0 gamma, X1 += Q1 delta
    void apply_nodes(const std::vector<NodeTarget>& targets) {
        if (targets.empty()) return;
        Context& ctx = Context::instance();
        const cudaStream_t stream = ctx.stream();
        const magma_queue_t queue = ctx.queue();
        size_t ws = 0;
        for (const NodeTarget& t : targets) ws += 4 * static_cast<size_t>(t.node->k) * t.ncols;
        DeviceHeap& heap = DeviceHeap::instance();
        T* work = heap.alloc<T>(sizeof(T) * std::max<size_t>(ws, 1));
        VBatch<T> b_c, b_w, b_d, b_gamma, b_x0, b_x1;
        ItemList<AxpbyItem<T>> copy_w, sub_c0;
        ItemList<RowSwapItem<T>> swaps;
        TrsmBatch<T> lower, upper;
        size_t off = 0;
        for (const NodeTarget& t : targets) {
            const Node& nd = *t.node;
            const int k = nd.k, nc = t.ncols;
            T* c0 = work + off;
            T* c1 = c0 + static_cast<size_t>(k) * nc;
            T* w = c1 + static_cast<size_t>(k) * nc;
            T* d = w + static_cast<size_t>(k) * nc;
            off += 4 * static_cast<size_t>(k) * nc;
            b_c.entries.push_back({nd.z0, t.x0, c0, k, nc, nd.n0, nd.n0, t.ld, k});
            b_c.entries.push_back({nd.z1, t.x1, c1, k, nc, nd.n1, nd.n1, t.ld, k});
            copy_w.add({w, k, c1, k, nullptr, 0, k, nc, T(1.0), T(0.0)}, k, nc);
            b_w.entries.push_back({nd.g1, c0, w, k, nc, k, k, k, k});
            b_d.entries.push_back({nd.g0, w, d, k, nc, k, k, k, k});
            swaps.add({d, k, k, nc, nd.ipiv}, k, nc);
            lower.add(nd.s, k, d, k, k, nc);
            upper.add(nd.s, k, d, k, k, nc);
            sub_c0.add({d, k, d, k, c0, k, k, nc, T(1.0), T(-1.0)}, k, nc);
            b_gamma.entries.push_back({nd.g1, d, c1, k, nc, k, k, k, k});
            b_x0.entries.push_back({nd.q0, c1, t.x0, nd.n0, nc, k, nd.n0, k, t.ld});
            b_x1.entries.push_back({nd.q1, d, t.x1, nd.n1, nc, k, nd.n1, k, t.ld});
            stats_.flops += kFlopScale * static_cast<double>(nc) *
                            (4.0 * k * static_cast<double>(nd.n0 + nd.n1) + 6.0 * k * static_cast<double>(k));
        }
        meta_.clear();
        b_c.stage(meta_);
        copy_w.stage(meta_);
        b_w.stage(meta_);
        b_d.stage(meta_);
        swaps.stage(meta_);
        lower.stage(meta_);
        upper.stage(meta_);
        sub_c0.stage(meta_);
        b_gamma.stage(meta_);
        b_x0.stage(meta_);
        b_x1.stage(meta_);
        char* meta = meta_.upload(meta_device_, stream);
        const T one(1.0), zero(0.0), minus(-1.0);
        b_c.gemm(meta, MagmaTrans, MagmaNoTrans, one, zero, queue);
        launch_axpby(copy_w.device(meta), copy_w.count(), copy_w.max_m, copy_w.max_n, stream);
        b_w.gemm(meta, MagmaNoTrans, MagmaNoTrans, minus, one, queue);
        b_d.gemm(meta, MagmaNoTrans, MagmaNoTrans, one, zero, queue);
        launch_row_swaps(swaps.device(meta), swaps.count(), swaps.max_n, stream);
        lower.solve(meta, MagmaLower, MagmaUnit, queue);
        upper.solve(meta, MagmaUpper, MagmaNonUnit, queue);
        launch_axpby(sub_c0.device(meta), sub_c0.count(), sub_c0.max_m, sub_c0.max_n, stream);
        b_gamma.gemm(meta, MagmaNoTrans, MagmaNoTrans, one, one, queue);
        b_x0.gemm(meta, MagmaNoTrans, MagmaNoTrans, minus, one, queue);
        b_x1.gemm(meta, MagmaNoTrans, MagmaNoTrans, one, one, queue);
        // the workspace is reused by the next launch only after these (one stream)
        heap.free(work);
    }

    int64_t n_loc_ = 0;
    int maxlevel_ = 0;
    std::vector<std::vector<Node>> nodes_;  // [level][node], by first row
    std::vector<Leaf> leaves_;
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
