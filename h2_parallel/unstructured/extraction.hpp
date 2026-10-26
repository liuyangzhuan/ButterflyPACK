#pragma once
// Batched entry extraction from a compression-only H2 (unstructured Color;
// GPU_BACKEND/h2_gpu/GPU_CA_PLAN.md, "Entries of a compressed H2"): the
// blocks A(I_n, J_n) of a list of intersections, from the stored
// representation alone, never the kernel.  The representation
// (compression_impl.hpp) holds, on the owner of each box t of a level:
//   - the interpolation T_t (k x r): x[S] += T x[R] going up;
//   - K(S_t, S_s) for each s in t's interaction list (levels 2..leaf);
//   - at the leaf, K(B_t, B_s) for s = t and each neighbour s.
// An entry A(a, b), a in leaf box t and b in leaf box s, is
//   - K_near(t, s)(a, b) if s is t or a neighbour of t;
//   - else u_a^T K(S_tl, S_sl) u_b at the one level l at which the
//     ancestors t_l and s_l are in each other's interaction lists (the
//     coarsest level, at least 2, at which they are not neighbours), u_x
//     being x's interpolation vector at level l: at the leaf e_i for
//     skeleton slot i, T(:, j) for redundant slot j; one level up, that
//     vector placed at its box's offset in the parent's point list (the
//     children's skeletons in child order), then v[S] + T v[R].
// The extraction therefore runs like the compressed matvec on the requested
// points only: the u vectors go up the tree (to the parents' owners at a
// process reduction), the column side's go to the row side's owner of each
// level, which computes the block, and the blocks go to their intersections'
// owners.
//
// Collective over the caller's communicator, every rank passing the same
// list, and inside it over the source's, a subset of the caller's ranks (the
// solver null elsewhere).  Indices are the user's, as the kernel callbacks
// take them, 1-based.  A nonsymmetric matrix's blocks are each box's own
// rows' (K(S_t, S_s) on t, K(B_t, B_s) on t), as the representation keeps
// them for any matrix.

#include "core/butterfly_types.hpp"
#include "core/blas_declare.hpp"
#include "core/butterfly_matvec.hpp"
#include "core/extraction_profile.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace butterfly {
namespace color_unstructured {

namespace extraction_detail {

inline bool boxes_adjacent(int dimension, int64_t a, int64_t b) {
    uint32_t ax = 0, ay = 0, az = 0, bx = 0, by = 0, bz = 0;
    morton::decode_nd(dimension, static_cast<uint64_t>(a), ax, ay, az);
    morton::decode_nd(dimension, static_cast<uint64_t>(b), bx, by, bz);
    auto close = [](uint32_t p, uint32_t q) { return (p > q ? p - q : q - p) <= 1; };
    return close(ax, bx) && close(ay, by) && close(az, bz);
}

// The level at which leaf boxes t and s (neither equal nor neighbours) meet
// in each other's interaction lists.
inline int meeting_level(int dimension, int leaf, int64_t t, int64_t s) {
    int level = leaf;
    while (level > 2) {
        const int64_t tp = t >> dimension, sp = s >> dimension;
        if (tp == sp || boxes_adjacent(dimension, tp, sp)) break;
        t = tp;
        s = sp;
        --level;
    }
    return level;
}

// The requested points of each intersection's rows (or columns), grouped by
// leaf box: perm holds offsets into the stacked index array, by
// intersection, then leaf box, then position.
struct Groups {
    std::vector<int64_t> perm;
    std::vector<int64_t> start;   // group g: perm[start[g] .. start[g + 1])
    std::vector<int64_t> first;   // intersection n: groups first[n] .. first[n + 1]
    std::vector<int64_t> morton;  // group g's leaf box
};

inline Groups make_groups(int64_t ninter, const std::vector<int64_t>& offset,
                          const std::vector<int64_t>& point_of, const std::vector<int64_t>& leaf_of) {
    Groups g;
    g.perm.resize(point_of.size());
    g.first.assign(static_cast<size_t>(ninter) + 1, 0);
    g.start.push_back(0);
    for (int64_t n = 0; n < ninter; ++n) {
        const int64_t a = offset[static_cast<size_t>(n)], b = offset[static_cast<size_t>(n) + 1];
        for (int64_t i = a; i < b; ++i) g.perm[static_cast<size_t>(i)] = i;
        std::stable_sort(g.perm.begin() + a, g.perm.begin() + b, [&](int64_t x, int64_t y) {
            return leaf_of[static_cast<size_t>(point_of[static_cast<size_t>(x)])] <
                   leaf_of[static_cast<size_t>(point_of[static_cast<size_t>(y)])];
        });
        g.first[static_cast<size_t>(n)] = static_cast<int64_t>(g.morton.size());
        for (int64_t i = a; i < b; ++i) {
            const int64_t m = leaf_of[static_cast<size_t>(point_of[static_cast<size_t>(g.perm[static_cast<size_t>(i)])])];
            if (i == a || m != g.morton.back()) {
                if (i != a) g.start.push_back(i);
                g.morton.push_back(m);
            }
        }
        if (b > a) g.start.push_back(b);
    }
    g.first[static_cast<size_t>(ninter)] = static_cast<int64_t>(g.morton.size());
    if (g.start.size() != g.morton.size() + 1) {
        throw std::logic_error("extract_compressed_blocks: inconsistent groups");
    }
    return g;
}

// The interpolation vectors of the requested points of one box of a level,
// on its owner: column j is points[j]'s, k entries each.
template<typename DataType>
struct BoxVectors {
    int64_t k = 0;
    std::vector<int64_t> points;
    std::vector<DataType> u;
    std::unordered_map<int64_t, int64_t> column;  // point -> column, built on first use
    const DataType* of(int64_t point) {
        if (column.empty()) {
            for (size_t j = 0; j < points.size(); ++j) column.emplace(points[j], static_cast<int64_t>(j));
        }
        const auto it = column.find(point);
        if (it == column.end()) throw std::logic_error("extract_compressed_blocks: a missing interpolation vector");
        return u.data() + it->second * k;
    }
};

inline int checked_count(size_t count, const char* what) {
    if (count > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(std::string("extract_compressed_blocks: ") + what + " exceeds INT_MAX");
    }
    return static_cast<int>(count);
}

// One all-to-all of per-rank buffers of T over comm (mpi_seconds: plus the
// time of the all-to-alls themselves).
template<typename T>
std::vector<std::vector<T>> exchange(MPI_Comm comm, const std::vector<std::vector<T>>& send, MPI_Datatype type,
                                     double* mpi_seconds = nullptr) {
    int size = 1;
    MPI_Comm_size(comm, &size);
    std::vector<int> sc(static_cast<size_t>(size)), rc(static_cast<size_t>(size)), sd(static_cast<size_t>(size)),
        rd(static_cast<size_t>(size));
    size_t total = 0;
    for (int p = 0; p < size; ++p) {
        sc[static_cast<size_t>(p)] = checked_count(send[static_cast<size_t>(p)].size(), "a message");
        sd[static_cast<size_t>(p)] = checked_count(total, "a send buffer");
        total += send[static_cast<size_t>(p)].size();
    }
    std::vector<T> sbuf;
    sbuf.reserve(total);
    for (const auto& v : send) sbuf.insert(sbuf.end(), v.begin(), v.end());
    const double t_mpi = MPI_Wtime();
    MPI_Alltoall(sc.data(), 1, MPI_INT, rc.data(), 1, MPI_INT, comm);
    size_t rtotal = 0;
    for (int p = 0; p < size; ++p) {
        rd[static_cast<size_t>(p)] = checked_count(rtotal, "a receive buffer");
        rtotal += static_cast<size_t>(rc[static_cast<size_t>(p)]);
    }
    std::vector<T> rbuf(rtotal);
    MPI_Alltoallv(sbuf.data(), sc.data(), sd.data(), type, rbuf.data(), rc.data(), rd.data(), type, comm);
    if (mpi_seconds != nullptr) *mpi_seconds += MPI_Wtime() - t_mpi;
    std::vector<std::vector<T>> recv(static_cast<size_t>(size));
    for (int p = 0; p < size; ++p) {
        recv[static_cast<size_t>(p)].assign(rbuf.begin() + rd[static_cast<size_t>(p)],
                                            rbuf.begin() + rd[static_cast<size_t>(p)] + rc[static_cast<size_t>(p)]);
    }
    return recv;
}

}  // namespace extraction_detail

namespace extraction_detail {

// One task: an intersection's row group by column group, on the owner of its
// row side's box of `level` (the leaf for a near task, else the meeting level).
struct Task {
    int64_t inter, rgroup, cgroup;
    int level;
    bool near;
    int64_t t, s;  // the row and column sides' boxes of that level
    int computer;  // its source rank
};

// What both extractions share: the requested points resolved to their leaf
// boxes and grouped, the tasks, and the points' interpolation vectors carried
// up the tree on their owners (on the source's ranks; tree null elsewhere).
template<typename CoordType, typename DataType>
struct Plan {
    ParallelTree<CoordType, DataType>* tree = nullptr;
    int crank = 0, csize = 1, srank = -1, ssize = 0, dimension = 0, leaf = 0, coarsest = 0;
    MPI_Datatype value_type = MPI_DATATYPE_NULL;
    int64_t ninter = 0;
    std::vector<int64_t> row_offset, col_offset, distinct, row_point, col_point, leaf_box, leaf_position;
    Groups rg, cg;
    std::vector<Task> tasks;
    std::vector<std::unordered_map<int64_t, BoxVectors<DataType>>> vectors;  // level -> box -> vectors

    // f(point, stacked index) for each point of a row (column) group
    template<typename F>
    void row_points(int64_t group, F&& f) const { points_of(rg, row_point, group, f); }
    template<typename F>
    void col_points(int64_t group, F&& f) const { points_of(cg, col_point, group, f); }
    int64_t row_count(int64_t group) const {
        return rg.start[static_cast<size_t>(group) + 1] - rg.start[static_cast<size_t>(group)];
    }
    int64_t col_count(int64_t group) const {
        return cg.start[static_cast<size_t>(group) + 1] - cg.start[static_cast<size_t>(group)];
    }

private:
    template<typename F>
    static void points_of(const Groups& g, const std::vector<int64_t>& point_of, int64_t group, F& f) {
        for (int64_t i = g.start[static_cast<size_t>(group)]; i < g.start[static_cast<size_t>(group) + 1]; ++i) {
            const int64_t stacked = g.perm[static_cast<size_t>(i)];
            f(point_of[static_cast<size_t>(stacked)], stacked);
        }
    }
};

template<typename CoordType, typename DataType, typename Lap>
Plan<CoordType, DataType> plan_extraction(
    H2<CoordType, DataType>* source, MPI_Comm caller_comm, int64_t ninter,
    const int* allrows, const int* allcols, const int* rowidx, const int* colidx, Lap&& lap) {
    Plan<CoordType, DataType> p;
    p.ninter = ninter;
    MPI_Comm_rank(caller_comm, &p.crank);
    MPI_Comm_size(caller_comm, &p.csize);
    p.value_type = fmm::mpi_datatype_for<DataType>();

    // the source's ranks among the caller's
    p.tree = source != nullptr ? source->tree.get() : nullptr;
    auto* tree = p.tree;
    if (tree != nullptr) MPI_Comm_rank(tree->comm, &p.srank);
    std::vector<int> source_rank_of(static_cast<size_t>(p.csize));
    MPI_Allgather(&p.srank, 1, MPI_INT, source_rank_of.data(), 1, MPI_INT, caller_comm);
    int nsource = 0;
    for (int r : source_rank_of) nsource += r >= 0 ? 1 : 0;
    std::vector<int> caller_rank_of(static_cast<size_t>(std::max(nsource, 1)), -1);
    for (int c = 0; c < p.csize; ++c) {
        const int r = source_rank_of[static_cast<size_t>(c)];
        if (r < 0) continue;
        if (r >= nsource || caller_rank_of[static_cast<size_t>(r)] != -1) {
            throw std::invalid_argument("extract_compressed_blocks: the source's ranks are not a subset of the caller's");
        }
        caller_rank_of[static_cast<size_t>(r)] = c;
    }
    if (nsource == 0) throw std::invalid_argument("extract_compressed_blocks: no rank holds the source");
    if (tree != nullptr) {
        if (source->build_state != H2BuildState::H2_COMPRESSED || source->options.unstructured != 1) {
            throw std::invalid_argument(
                "extract_compressed_blocks: the source must be a compression-only unstructured H2 (precon 2)");
        }
        MPI_Comm_size(tree->comm, &p.ssize);
        if (p.ssize != nsource) {
            throw std::invalid_argument("extract_compressed_blocks: the source's ranks are not a subset of the caller's");
        }
    }

    // the intersections' points: the distinct indices, and each stacked
    // index's place among them
    p.row_offset.assign(static_cast<size_t>(ninter) + 1, 0);
    p.col_offset.assign(static_cast<size_t>(ninter) + 1, 0);
    for (int64_t n = 0; n < ninter; ++n) {
        if (rowidx[n] < 0 || colidx[n] < 0) throw std::invalid_argument("extract_compressed_blocks: a negative size");
        p.row_offset[static_cast<size_t>(n) + 1] = p.row_offset[static_cast<size_t>(n)] + rowidx[n];
        p.col_offset[static_cast<size_t>(n) + 1] = p.col_offset[static_cast<size_t>(n)] + colidx[n];
    }
    const int64_t nrows = p.row_offset.back(), ncols = p.col_offset.back();
    auto& distinct = p.distinct;
    distinct.reserve(static_cast<size_t>(nrows + ncols));
    for (int64_t i = 0; i < nrows; ++i) distinct.push_back(static_cast<int64_t>(allrows[i]) - 1);
    for (int64_t j = 0; j < ncols; ++j) distinct.push_back(static_cast<int64_t>(allcols[j]) - 1);
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
    const int64_t npoints = static_cast<int64_t>(distinct.size());
    auto place = [&](int64_t id) {
        const auto it = std::lower_bound(distinct.begin(), distinct.end(), id);
        return static_cast<int64_t>(it - distinct.begin());
    };
    p.row_point.resize(static_cast<size_t>(nrows));
    p.col_point.resize(static_cast<size_t>(ncols));
    for (int64_t i = 0; i < nrows; ++i) p.row_point[static_cast<size_t>(i)] = place(static_cast<int64_t>(allrows[i]) - 1);
    for (int64_t j = 0; j < ncols; ++j) p.col_point[static_cast<size_t>(j)] = place(static_cast<int64_t>(allcols[j]) - 1);

    // each point's leaf box and its position there, from the leaf owners
    p.leaf_box.assign(static_cast<size_t>(npoints), -1);
    p.leaf_position.assign(static_cast<size_t>(npoints), -1);
    if (tree != nullptr) {
        p.dimension = tree->dimension;
        p.leaf = tree->num_levels - 1;
        auto& level = tree->levels[static_cast<size_t>(p.leaf)];
        if (level.is_process_active) {
            for (auto& box : level.local_boxes) {
                for (size_t q = 0; q < box.point_indices.size(); ++q) {
                    const auto it = std::lower_bound(distinct.begin(), distinct.end(), box.point_indices[q]);
                    if (it == distinct.end() || *it != box.point_indices[q]) continue;
                    const size_t d = static_cast<size_t>(it - distinct.begin());
                    p.leaf_box[d] = box.morton_index;
                    p.leaf_position[d] = static_cast<int64_t>(q);
                }
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, p.leaf_box.data(), static_cast<int>(npoints), MPI_INT64_T, MPI_MAX, caller_comm);
    MPI_Allreduce(MPI_IN_PLACE, p.leaf_position.data(), static_cast<int>(npoints), MPI_INT64_T, MPI_MAX, caller_comm);
    for (int64_t d = 0; d < npoints; ++d) {
        if (p.leaf_box[static_cast<size_t>(d)] < 0) {
            throw std::invalid_argument("extract_compressed_blocks: index " +
                                        std::to_string(distinct[static_cast<size_t>(d)] + 1) +
                                        " is not a point of the source");
        }
    }
    p.rg = make_groups(ninter, p.row_offset, p.row_point, p.leaf_box);
    p.cg = make_groups(ninter, p.col_offset, p.col_point, p.leaf_box);
    lap(fmm::ExtractionProfile::kResolve);
    if (tree == nullptr) return p;

    // ---- the tasks (every source rank: the same list), the coarsest level
    // each point's interpolation vector is needed at
    const int dimension = p.dimension, leaf = p.leaf, srank = p.srank;
    for (int64_t n = 0; n < ninter; ++n) {
        for (int64_t a = p.rg.first[static_cast<size_t>(n)]; a < p.rg.first[static_cast<size_t>(n) + 1]; ++a) {
            for (int64_t b = p.cg.first[static_cast<size_t>(n)]; b < p.cg.first[static_cast<size_t>(n) + 1]; ++b) {
                Task task{n, a, b, leaf, true, p.rg.morton[static_cast<size_t>(a)], p.cg.morton[static_cast<size_t>(b)], -1};
                if (task.t != task.s && !boxes_adjacent(dimension, task.t, task.s)) {
                    task.near = false;
                    task.level = meeting_level(dimension, leaf, task.t, task.s);
                    task.t >>= dimension * (leaf - task.level);
                    task.s >>= dimension * (leaf - task.level);
                }
                task.computer = tree->levels[static_cast<size_t>(task.level)].owner_rank_of(task.t);
                p.tasks.push_back(task);
            }
        }
    }
    std::vector<int> need(static_cast<size_t>(npoints), leaf + 1);
    p.coarsest = leaf + 1;
    for (const Task& task : p.tasks) {
        if (task.near) continue;
        p.coarsest = std::min(p.coarsest, task.level);
        p.row_points(task.rgroup, [&](int64_t d, int64_t) {
            need[static_cast<size_t>(d)] = std::min(need[static_cast<size_t>(d)], task.level);
        });
        p.col_points(task.cgroup, [&](int64_t d, int64_t) {
            need[static_cast<size_t>(d)] = std::min(need[static_cast<size_t>(d)], task.level);
        });
    }
    lap(fmm::ExtractionProfile::kTasks);

    // ---- the interpolation vectors, leaf up
    auto& vectors = p.vectors;
    vectors.resize(static_cast<size_t>(leaf) + 1);
    {
        auto& level = tree->levels[static_cast<size_t>(leaf)];
        for (int64_t d = 0; d < npoints; ++d) {
            if (need[static_cast<size_t>(d)] > leaf || level.owner_rank_of(p.leaf_box[static_cast<size_t>(d)]) != srank) {
                continue;
            }
            auto* box = level.find_local_box(p.leaf_box[static_cast<size_t>(d)]);
            if (box == nullptr) throw std::logic_error("extract_compressed_blocks: a leaf box is not local");
            auto& bv = vectors[static_cast<size_t>(leaf)][box->morton_index];
            bv.k = static_cast<int64_t>(box->skeleton_indices.size());
            bv.points.push_back(d);
            const int64_t q = p.leaf_position[static_cast<size_t>(d)];
            std::vector<DataType> col(static_cast<size_t>(bv.k), DataType{0});
            const auto sit = std::find(box->skeleton_indices.begin(), box->skeleton_indices.end(), q);
            if (sit != box->skeleton_indices.end()) {
                col[static_cast<size_t>(sit - box->skeleton_indices.begin())] = DataType{1};
            } else {
                const auto rit = std::find(box->redundant_indices.begin(), box->redundant_indices.end(), q);
                if (rit == box->redundant_indices.end()) {
                    throw std::logic_error("extract_compressed_blocks: a point neither skeleton nor redundant");
                }
                const int64_t j = rit - box->redundant_indices.begin();
                const auto& T = box->interpolation_matrix;
                for (int64_t i = 0; i < bv.k; ++i) col[static_cast<size_t>(i)] = T.data[static_cast<size_t>(i + j * T.lda)];
            }
            bv.u.insert(bv.u.end(), col.begin(), col.end());
        }
    }
    const int ssize = p.ssize;
    for (int l = leaf; l > std::max(p.coarsest, 2); --l) {
        // to the parents (their owners at a process reduction): parent,
        // offset in its point list, k, m; the points; the vectors
        auto& child_level = tree->levels[static_cast<size_t>(l)];
        auto& parent_level = tree->levels[static_cast<size_t>(l) - 1];
        std::vector<std::vector<int64_t>> meta(static_cast<size_t>(ssize));
        std::vector<std::vector<DataType>> data(static_cast<size_t>(ssize));
        for (auto& [c, bv] : vectors[static_cast<size_t>(l)]) {
            std::vector<int64_t> up;
            for (size_t j = 0; j < bv.points.size(); ++j) {
                if (need[static_cast<size_t>(bv.points[j])] < l) up.push_back(static_cast<int64_t>(j));
            }
            if (up.empty()) continue;
            const int64_t parent = c >> dimension;
            int64_t offset = 0;
            for (int64_t sib = parent << dimension; sib < c; ++sib) {
                if (const auto* sb = child_level.find_local_box(sib)) offset += static_cast<int64_t>(sb->skeleton_indices.size());
            }
            const int to = parent_level.owner_rank_of(parent);
            auto& mv = meta[static_cast<size_t>(to)];
            mv.insert(mv.end(), {parent, offset, bv.k, static_cast<int64_t>(up.size())});
            for (int64_t j : up) mv.push_back(bv.points[static_cast<size_t>(j)]);
            auto& dv = data[static_cast<size_t>(to)];
            for (int64_t j : up) dv.insert(dv.end(), bv.u.begin() + j * bv.k, bv.u.begin() + (j + 1) * bv.k);
        }
        const auto rmeta = exchange(tree->comm, meta, MPI_INT64_T);
        const auto rdata = exchange(tree->comm, data, p.value_type);
        for (int from = 0; from < ssize; ++from) {
            const auto& mv = rmeta[static_cast<size_t>(from)];
            const auto& dv = rdata[static_cast<size_t>(from)];
            size_t mi = 0, di = 0;
            while (mi < mv.size()) {
                const int64_t parent = mv[mi], offset = mv[mi + 1], kc = mv[mi + 2], m = mv[mi + 3];
                mi += 4;
                auto* pb = parent_level.find_local_box(parent);
                if (pb == nullptr) throw std::logic_error("extract_compressed_blocks: a parent box is not local");
                const int64_t np = static_cast<int64_t>(pb->point_indices.size());
                const int64_t kp = static_cast<int64_t>(pb->skeleton_indices.size());
                if (offset + kc > np) throw std::logic_error("extract_compressed_blocks: a child outside its parent");
                auto& bv = vectors[static_cast<size_t>(l) - 1][parent];
                bv.k = kp;
                std::vector<DataType> v(static_cast<size_t>(np));
                for (int64_t j = 0; j < m; ++j) {
                    std::fill(v.begin(), v.end(), DataType{0});
                    std::copy_n(dv.begin() + static_cast<std::ptrdiff_t>(di + j * kc), kc, v.begin() + offset);
                    // u = v[S] + T v[R]
                    std::vector<DataType> u(static_cast<size_t>(kp));
                    for (int64_t i = 0; i < kp; ++i) u[static_cast<size_t>(i)] = v[static_cast<size_t>(pb->skeleton_indices[static_cast<size_t>(i)])];
                    const auto& T = pb->interpolation_matrix;
                    for (size_t r = 0; r < pb->redundant_indices.size(); ++r) {
                        const DataType x = v[static_cast<size_t>(pb->redundant_indices[r])];
                        if (x == DataType{0}) continue;
                        for (int64_t i = 0; i < kp; ++i) u[static_cast<size_t>(i)] += T.data[static_cast<size_t>(i + static_cast<int64_t>(r) * T.lda)] * x;
                    }
                    bv.points.push_back(mv[mi + static_cast<size_t>(j)]);
                    bv.u.insert(bv.u.end(), u.begin(), u.end());
                }
                mi += static_cast<size_t>(m);
                di += static_cast<size_t>(m * kc);
            }
        }
    }
    // (the vectors' point maps, before threads read them)
    for (auto& per_level : vectors) {
        for (auto& [m, bv] : per_level) {
            if (!bv.points.empty()) (void)bv.of(bv.points.front());
        }
    }
    lap(fmm::ExtractionProfile::kUp);
    return p;
}

// The near block K(B_t, B_s) of leaf box t (local) and its column map (the
// global index of each column to its position).
template<typename CoordType, typename DataType>
const H2Block<DataType>& near_block(TreeLevel<CoordType, DataType>& level, BoxData<CoordType, DataType>* tb, int64_t s,
                                    std::unordered_map<int64_t, int64_t>& column_of) {
    const H2Block<DataType>* block = nullptr;
    for (const auto& b : tb->h2_near_blocks) {
        if (b.source_morton == s) block = &b;
    }
    if (block == nullptr) throw std::logic_error("extract_compressed_blocks: a missing near block");
    const auto columns = h2_box_global_indices(level, s, false);
    column_of.clear();
    for (size_t j = 0; j < columns.size(); ++j) column_of.emplace(columns[j], static_cast<int64_t>(j));
    return *block;
}

}  // namespace extraction_detail

// The blocks A(I_n, J_n), n < ninter, each to the caller's rank owner[n]
// (its blocks in alldat_loc, in intersection order, each column major).
// allrows / allcols: the intersections' 1-based indices, stacked (rowidx[n]
// and colidx[n] of them).  source: the compression-only H2, null on the
// caller's ranks outside it; its ranks must be a subset of caller_comm's.
template<typename CoordType, typename DataType>
void extract_compressed_blocks(
    H2<CoordType, DataType>* source, MPI_Comm caller_comm, int64_t ninter,
    const int* allrows, const int* allcols, const int* rowidx, const int* colidx,
    const std::vector<int>& owner, DataType* alldat_loc) {
    using namespace extraction_detail;
    // (the time by step: core/extraction_profile.hpp)
    fmm::ExtractionProfile& profile = fmm::extraction_profile();
    ++profile.calls;
    double t_step = MPI_Wtime();
    auto lap = [&](int step) {
        const double now = MPI_Wtime();
        profile.seconds[step] += now - t_step;
        t_step = now;
    };
    Plan<CoordType, DataType> p = plan_extraction(source, caller_comm, ninter, allrows, allcols, rowidx, colidx, lap);
    auto* tree = p.tree;
    const int srank = p.srank, csize = p.csize, crank = p.crank;
    const auto& tasks = p.tasks;
    const auto& rg = p.rg;
    const auto& cg = p.cg;

    // ---- the source: each task on the owner of its row side's box, the
    // values per task
    std::vector<std::vector<DataType>> values(tasks.size());
    if (tree != nullptr) {
        auto& vectors = p.vectors;
        // ---- the column sides' vectors to the row sides' owners, by level
        std::vector<std::unordered_map<int64_t, std::unordered_map<int64_t, std::vector<DataType>>>> remote(
            static_cast<size_t>(p.leaf) + 1);  // level -> box -> point -> vector
        for (int l = p.leaf; l >= std::max(p.coarsest, 2); --l) {
            auto& level = tree->levels[static_cast<size_t>(l)];
            std::vector<std::vector<int64_t>> meta(static_cast<size_t>(p.ssize));
            std::vector<std::vector<DataType>> data(static_cast<size_t>(p.ssize));
            std::vector<std::unordered_map<int64_t, std::vector<int64_t>>> sent(static_cast<size_t>(p.ssize));
            for (const Task& task : tasks) {
                if (task.near || task.level != l || task.computer == srank || level.owner_rank_of(task.s) != srank) continue;
                auto& bv = vectors[static_cast<size_t>(l)].at(task.s);
                auto& done = sent[static_cast<size_t>(task.computer)][task.s];
                p.col_points(task.cgroup, [&](int64_t d, int64_t) {
                    if (std::find(done.begin(), done.end(), d) != done.end()) return;
                    done.push_back(d);
                    auto& mv = meta[static_cast<size_t>(task.computer)];
                    mv.insert(mv.end(), {task.s, d, bv.k});
                    const DataType* u = bv.of(d);
                    data[static_cast<size_t>(task.computer)].insert(data[static_cast<size_t>(task.computer)].end(), u, u + bv.k);
                });
            }
            const auto rmeta = exchange(tree->comm, meta, MPI_INT64_T);
            const auto rdata = exchange(tree->comm, data, p.value_type);
            for (int from = 0; from < p.ssize; ++from) {
                const auto& mv = rmeta[static_cast<size_t>(from)];
                const auto& dv = rdata[static_cast<size_t>(from)];
                size_t di = 0;
                for (size_t mi = 0; mi + 3 <= mv.size(); mi += 3) {
                    const int64_t k = mv[mi + 2];
                    remote[static_cast<size_t>(l)][mv[mi]][mv[mi + 1]].assign(dv.begin() + static_cast<std::ptrdiff_t>(di),
                                                                              dv.begin() + static_cast<std::ptrdiff_t>(di + k));
                    di += static_cast<size_t>(k);
                }
            }
        }
        lap(fmm::ExtractionProfile::kColumns);

        // ---- this rank's tasks
        std::vector<size_t> mine;
        for (size_t i = 0; i < tasks.size(); ++i) {
            if (tasks[i].computer == srank) mine.push_back(i);
        }
        std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic)
        for (size_t q = 0; q < mine.size(); ++q) {
            if (failure) continue;
            try {
                const Task& task = tasks[mine[q]];
                const int64_t nr = p.row_count(task.rgroup), nc = p.col_count(task.cgroup);
                std::vector<DataType>& out = values[mine[q]];
                out.assign(static_cast<size_t>(nr * nc), DataType{0});
                auto& level = tree->levels[static_cast<size_t>(task.level)];
                auto* tb = level.find_local_box(task.t);
                if (tb == nullptr) throw std::logic_error("extract_compressed_blocks: a task's box is not local");
                if (task.near) {
                    std::unordered_map<int64_t, int64_t> column_of;
                    const auto& K = near_block(level, tb, task.s, column_of).matrix;
                    std::vector<int64_t> rpos, cpos;
                    p.row_points(task.rgroup, [&](int64_t d, int64_t) { rpos.push_back(p.leaf_position[static_cast<size_t>(d)]); });
                    p.col_points(task.cgroup, [&](int64_t d, int64_t) {
                        cpos.push_back(column_of.at(p.distinct[static_cast<size_t>(d)]));
                    });
                    for (int64_t j = 0; j < nc; ++j) {
                        for (int64_t i = 0; i < nr; ++i) {
                            out[static_cast<size_t>(i + j * nr)] =
                                K.data[static_cast<size_t>(rpos[static_cast<size_t>(i)] + cpos[static_cast<size_t>(j)] * K.lda)];
                        }
                    }
                    continue;
                }
                // far: U_t^T K(S_t, S_s) U_s
                const H2Block<DataType>* block = nullptr;
                for (const auto& b : tb->h2_interaction_blocks) {
                    if (b.source_morton == task.s) block = &b;
                }
                auto& tv = vectors[static_cast<size_t>(task.level)];
                const auto tit = tv.find(task.t);
                const int64_t kt = static_cast<int64_t>(tb->skeleton_indices.size());
                if (block == nullptr || kt == 0) continue;  // an empty skeleton: no coupling, zeros
                const auto& K = block->matrix;
                const int64_t ks = K.cols;
                if (tit == tv.end() || K.rows != kt) throw std::logic_error("extract_compressed_blocks: inconsistent far data");
                std::vector<DataType> ut(static_cast<size_t>(kt * nr)), us(static_cast<size_t>(ks * nc));
                int64_t i = 0;
                p.row_points(task.rgroup, [&](int64_t d, int64_t) {
                    std::copy_n(tit->second.of(d), kt, ut.begin() + i * kt);
                    ++i;
                });
                int64_t j = 0;
                const bool local_s = level.owner_rank_of(task.s) == srank;
                p.col_points(task.cgroup, [&](int64_t d, int64_t) {
                    const DataType* u = local_s ? tv.at(task.s).of(d)
                                                : remote[static_cast<size_t>(task.level)].at(task.s).at(d).data();
                    std::copy_n(u, ks, us.begin() + j * ks);
                    ++j;
                });
                std::vector<DataType> w(static_cast<size_t>(kt * nc));
                {
                    const int m = static_cast<int>(kt), n = static_cast<int>(nc), k = static_cast<int>(ks);
                    const int lda = static_cast<int>(K.lda), ldb = static_cast<int>(ks), ldc = static_cast<int>(kt);
                    const DataType one{1}, zero{0};
                    const char no = 'N';
                    if (m > 0 && n > 0 && k > 0) {
                        gemm_(&no, &no, &m, &n, &k, &one, K.data.data(), &lda, us.data(), &ldb, &zero, w.data(), &ldc);
                    }
                }
                {
                    const int m = static_cast<int>(nr), n = static_cast<int>(nc), k = static_cast<int>(kt);
                    const int lda = static_cast<int>(kt), ldb = static_cast<int>(kt), ldc = static_cast<int>(nr);
                    const DataType one{1}, zero{0};
                    const char tr = 'T', no = 'N';
                    if (m > 0 && n > 0 && k > 0) {
                        gemm_(&tr, &no, &m, &n, &k, &one, ut.data(), &lda, w.data(), &ldb, &zero, out.data(), &ldc);
                    }
                }
            } catch (...) {
#pragma omp critical(extract_compressed_blocks_failure)
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
        for (size_t i : mine) profile.entries_computed += static_cast<double>(values[i].size());
        lap(fmm::ExtractionProfile::kCompute);
    }

    // ---- the blocks to their owners: per task, its intersection, row and
    // column groups, then its values
    std::vector<std::vector<int64_t>> meta(static_cast<size_t>(csize));
    std::vector<std::vector<DataType>> data(static_cast<size_t>(csize));
    for (size_t i = 0; i < tasks.size(); ++i) {
        const Task& task = tasks[i];
        if (task.computer != srank) continue;
        const int to = owner[static_cast<size_t>(task.inter)];
        meta[static_cast<size_t>(to)].insert(meta[static_cast<size_t>(to)].end(), {task.inter, task.rgroup, task.cgroup});
        data[static_cast<size_t>(to)].insert(data[static_cast<size_t>(to)].end(), values[i].begin(), values[i].end());
        profile.bytes_routed += static_cast<double>(values[i].size() * sizeof(DataType));
    }
    lap(fmm::ExtractionProfile::kPack);
    const auto rmeta = exchange(caller_comm, meta, MPI_INT64_T, &profile.route_mpi);
    const auto rdata = exchange(caller_comm, data, p.value_type, &profile.route_mpi);
    lap(fmm::ExtractionProfile::kRoute);
    std::vector<int64_t> base(static_cast<size_t>(ninter), -1);
    int64_t expected = 0;
    for (int64_t n = 0; n < ninter; ++n) {
        if (owner[static_cast<size_t>(n)] != crank) continue;
        base[static_cast<size_t>(n)] = expected;
        expected += static_cast<int64_t>(rowidx[n]) * colidx[n];
    }
    int64_t written = 0;
    for (int from = 0; from < csize; ++from) {
        const auto& mv = rmeta[static_cast<size_t>(from)];
        const auto& dv = rdata[static_cast<size_t>(from)];
        size_t di = 0;
        for (size_t mi = 0; mi + 3 <= mv.size(); mi += 3) {
            const int64_t n = mv[mi], a = mv[mi + 1], b = mv[mi + 2];
            if (base[static_cast<size_t>(n)] < 0) throw std::logic_error("extract_compressed_blocks: a block for another rank");
            const int64_t nr = p.row_count(a), nc = p.col_count(b);
            for (int64_t j = 0; j < nc; ++j) {
                const int64_t col = cg.perm[static_cast<size_t>(cg.start[static_cast<size_t>(b)] + j)] - p.col_offset[static_cast<size_t>(n)];
                for (int64_t i = 0; i < nr; ++i) {
                    const int64_t row = rg.perm[static_cast<size_t>(rg.start[static_cast<size_t>(a)] + i)] - p.row_offset[static_cast<size_t>(n)];
                    alldat_loc[base[static_cast<size_t>(n)] + row + col * rowidx[n]] = dv[di + static_cast<size_t>(i + j * nr)];
                }
            }
            di += static_cast<size_t>(nr * nc);
            written += nr * nc;
        }
    }
    if (written != expected) {
        throw std::logic_error("extract_compressed_blocks: " + std::to_string(written) + " of " +
                               std::to_string(expected) + " entries received");
    }
    lap(fmm::ExtractionProfile::kUnpack);
}

// The sketched blocks Y_n = Omega_n A(I_n, J_n), n < nsketch, each to the
// caller's rank owner[n] (in alldat_loc, in request order, each dims[n] x
// |J_n|, column major).  Omega_n (dims[n] x |I_n|) is given by its entries
// omega_off[n] .. omega_off[n + 1]: omega_row, the position in I_n, and
// omega_dest, the sketch row (both 0-based), and the value.  Indices and the
// source as for extract_compressed_blocks.  A is never formed far from the
// diagonal: the owner of a row box t forms (Omega U_t^T) K(S_t, S_s), d x k_s,
// summed on that rank over its boxes t, then over the ranks on the owner of
// s, which multiplies by U_s of the columns; near blocks are sketched where
// they are stored.
//
// transposed: Y_n = Omega_n A(J_n, I_n)^T (a nonsymmetric source; the
// rows I_n still Omega's).  The blocks of A(J, I) are stored on the column
// side's boxes s, so the owner of t forms M = U_t Omega^T (k_t x d) and sends
// it to the owner of s, which sums (K(S_s, S_t) M)^T over the boxes t and
// multiplies by U_s; the near blocks K(B_s, B_t) are sketched on s's owner.
template<typename CoordType, typename DataType>
void extract_compressed_sketches(
    H2<CoordType, DataType>* source, MPI_Comm caller_comm, int64_t nsketch,
    const int* allrows, const int* allcols, const int* rowidx, const int* colidx, const int* dims,
    const int64_t* omega_off, const int* omega_row, const int* omega_dest, const double* omega_val,
    const std::vector<int>& owner, DataType* alldat_loc, bool transposed = false) {
    using namespace extraction_detail;
    fmm::ExtractionProfile& profile = fmm::extraction_profile();
    ++profile.calls;
    double t_step = MPI_Wtime();
    auto lap = [&](int step) {
        const double now = MPI_Wtime();
        profile.seconds[step] += now - t_step;
        t_step = now;
    };
    for (int64_t n = 0; n < nsketch; ++n) {
        if (dims[n] < 0) throw std::invalid_argument("extract_compressed_sketches: a negative sketch size");
        for (int64_t e = omega_off[n]; e < omega_off[n + 1]; ++e) {
            if (omega_row[e] < 0 || omega_row[e] >= rowidx[n] || omega_dest[e] < 0 || omega_dest[e] >= dims[n]) {
                throw std::invalid_argument("extract_compressed_sketches: a sketch entry out of range");
            }
        }
    }
    Plan<CoordType, DataType> p = plan_extraction(source, caller_comm, nsketch, allrows, allcols, rowidx, colidx, lap);
    auto* tree = p.tree;
    const int srank = p.srank, csize = p.csize, crank = p.crank;
    const auto& tasks = p.tasks;

    // the contributions to the owners: per piece, its request, column group,
    // sketch rows and columns, then its values (d x columns)
    struct Piece {
        int64_t n, cgroup, d, nc;
        std::vector<DataType> y;
    };
    std::vector<Piece> pieces;
    if (tree != nullptr) {
        auto& vectors = p.vectors;
        // each stacked row's sketch entries
        const int64_t nrows = p.row_offset.back();
        std::vector<int64_t> row_omega(static_cast<size_t>(nrows) + 1, 0), omega_of(static_cast<size_t>(omega_off[nsketch]));
        for (int64_t n = 0; n < nsketch; ++n) {
            for (int64_t e = omega_off[n]; e < omega_off[n + 1]; ++e) {
                ++row_omega[static_cast<size_t>(p.row_offset[static_cast<size_t>(n)] + omega_row[e]) + 1];
            }
        }
        for (int64_t i = 0; i < nrows; ++i) row_omega[static_cast<size_t>(i) + 1] += row_omega[static_cast<size_t>(i)];
        {
            std::vector<int64_t> at(row_omega.begin(), row_omega.end() - 1);
            for (int64_t n = 0; n < nsketch; ++n) {
                for (int64_t e = omega_off[n]; e < omega_off[n + 1]; ++e) {
                    omega_of[static_cast<size_t>(at[static_cast<size_t>(p.row_offset[static_cast<size_t>(n)] + omega_row[e])]++)] = e;
                }
            }
        }
        // M (k x d) += Omega's columns of the group's rows times their vectors
        auto sketch_vectors = [&](int64_t n, int64_t group, BoxVectors<DataType>& bv, int64_t d, DataType* M) {
            const int64_t k = bv.k;
            p.row_points(group, [&](int64_t pt, int64_t stacked) {
                const DataType* u = bv.of(pt);
                for (int64_t q = row_omega[static_cast<size_t>(stacked)]; q < row_omega[static_cast<size_t>(stacked) + 1]; ++q) {
                    const int64_t e = omega_of[static_cast<size_t>(q)];
                    const DataType w = static_cast<DataType>(omega_val[e]);
                    DataType* col = M + static_cast<int64_t>(omega_dest[e]) * k;
                    for (int64_t i = 0; i < k; ++i) col[i] += w * u[i];
                }
            });
            (void)n;
            (void)d;
        };

        // the row groups of each (n, l, s, t) on t's owner (the task's computer)
        std::map<std::tuple<int64_t, int, int64_t, int64_t>, std::vector<int64_t>> row_groups;  // (n, l, s, t) -> groups
        for (const Task& task : tasks) {
            if (task.near || task.computer != srank) continue;
            auto& g = row_groups[{task.inter, task.level, task.s, task.t}];
            if (g.empty() || g.back() != task.rgroup) g.push_back(task.rgroup);
        }
        std::map<std::tuple<int64_t, int, int64_t>, std::pair<int64_t, std::vector<DataType>>> sums;  // (n, l, s) -> ks, W
        std::exception_ptr failure;
        if (!transposed) {
            // ---- far, on the row sides' owners: W(n, l, s) = sum over this
            // rank's boxes t of (Omega U_t^T) K(S_t, S_s), d x k_s
            struct Partial {
                int64_t n;
                int l;
                int64_t s, ks;
                std::vector<std::pair<int64_t, const std::vector<int64_t>*>> boxes;  // t, its row groups
                std::vector<DataType> w;
            };
            std::vector<Partial> partials;
            for (const auto& [key, groups] : row_groups) {
                const auto& [n, l, s, t] = key;
                if (partials.empty() || partials.back().n != n || partials.back().l != l || partials.back().s != s) {
                    partials.push_back(Partial{n, l, s, 0, {}, {}});
                }
                partials.back().boxes.emplace_back(t, &groups);
            }
#pragma omp parallel for schedule(dynamic)
            for (size_t q = 0; q < partials.size(); ++q) {
                if (failure) continue;
                try {
                    Partial& part = partials[q];
                    auto& level = tree->levels[static_cast<size_t>(part.l)];
                    const int64_t d = dims[part.n];
                    std::vector<DataType> M;
                    for (const auto& [t, groups] : part.boxes) {
                        auto* tb = level.find_local_box(t);
                        if (tb == nullptr) throw std::logic_error("extract_compressed_sketches: a task's box is not local");
                        const H2Block<DataType>* block = nullptr;
                        for (const auto& b : tb->h2_interaction_blocks) {
                            if (b.source_morton == part.s) block = &b;
                        }
                        const int64_t kt = static_cast<int64_t>(tb->skeleton_indices.size());
                        if (block == nullptr || kt == 0 || d == 0) continue;  // no coupling
                        const auto& K = block->matrix;
                        if (K.rows != kt) throw std::logic_error("extract_compressed_sketches: inconsistent far data");
                        if (part.w.empty()) {
                            part.ks = K.cols;
                            part.w.assign(static_cast<size_t>(d * part.ks), DataType{0});
                        } else if (K.cols != part.ks) {
                            throw std::logic_error("extract_compressed_sketches: inconsistent column skeletons");
                        }
                        auto& bv = vectors[static_cast<size_t>(part.l)].at(t);
                        M.assign(static_cast<size_t>(kt * d), DataType{0});
                        for (int64_t g : *groups) sketch_vectors(part.n, g, bv, d, M.data());
                        // W += M^T K
                        const int m = static_cast<int>(d), nn = static_cast<int>(part.ks), k = static_cast<int>(kt);
                        const int lda = static_cast<int>(kt), ldb = static_cast<int>(K.lda), ldc = static_cast<int>(d);
                        const DataType one{1};
                        const char tr = 'T', no = 'N';
                        if (m > 0 && nn > 0) {
                            gemm_(&tr, &no, &m, &nn, &k, &one, M.data(), &lda, K.data.data(), &ldb, &one, part.w.data(), &ldc);
                        }
                    }
                } catch (...) {
#pragma omp critical(extract_compressed_sketches_failure)
                    if (!failure) failure = std::current_exception();
                }
            }
            if (failure) std::rethrow_exception(failure);
            lap(fmm::ExtractionProfile::kCompute);

            // ---- the partials to the column sides' owners, summed there
            std::vector<std::vector<int64_t>> meta(static_cast<size_t>(p.ssize));
            std::vector<std::vector<DataType>> data(static_cast<size_t>(p.ssize));
            for (const Partial& part : partials) {
                if (part.w.empty()) continue;
                const int to = tree->levels[static_cast<size_t>(part.l)].owner_rank_of(part.s);
                meta[static_cast<size_t>(to)].insert(meta[static_cast<size_t>(to)].end(), {part.n, part.l, part.s, part.ks});
                data[static_cast<size_t>(to)].insert(data[static_cast<size_t>(to)].end(), part.w.begin(), part.w.end());
            }
            partials.clear();
            const auto rmeta = exchange(tree->comm, meta, MPI_INT64_T);
            const auto rdata = exchange(tree->comm, data, p.value_type);
            for (int from = 0; from < p.ssize; ++from) {
                const auto& mv = rmeta[static_cast<size_t>(from)];
                const auto& dv = rdata[static_cast<size_t>(from)];
                size_t di = 0;
                for (size_t mi = 0; mi + 4 <= mv.size(); mi += 4) {
                    const int64_t n = mv[mi], s = mv[mi + 2], ks = mv[mi + 3];
                    const int l = static_cast<int>(mv[mi + 1]);
                    const size_t len = static_cast<size_t>(dims[n] * ks);
                    auto& [k, w] = sums[{n, l, s}];
                    if (w.empty()) {
                        k = ks;
                        w.assign(dv.begin() + static_cast<std::ptrdiff_t>(di), dv.begin() + static_cast<std::ptrdiff_t>(di + len));
                    } else {
                        if (k != ks) throw std::logic_error("extract_compressed_sketches: inconsistent partials");
                        for (size_t i = 0; i < len; ++i) w[i] += dv[di + i];
                    }
                    di += len;
                }
            }
            lap(fmm::ExtractionProfile::kColumns);
        } else {
            // ---- far, transposed: M = U_t Omega^T (k_t x d) on t's owner, to
            // s's owner
            struct Rows {
                std::tuple<int64_t, int, int64_t, int64_t> key;  // (n, l, s, t)
                const std::vector<int64_t>* groups;
                int64_t kt = 0;
                std::vector<DataType> m;
            };
            std::vector<Rows> rows;
            for (const auto& [key, groups] : row_groups) rows.push_back(Rows{key, &groups, 0, {}});
#pragma omp parallel for schedule(dynamic)
            for (size_t q = 0; q < rows.size(); ++q) {
                if (failure) continue;
                try {
                    Rows& row = rows[q];
                    const auto& [n, l, s, t] = row.key;
                    auto* tb = tree->levels[static_cast<size_t>(l)].find_local_box(t);
                    if (tb == nullptr) throw std::logic_error("extract_compressed_sketches: a task's box is not local");
                    const int64_t kt = static_cast<int64_t>(tb->skeleton_indices.size()), d = dims[n];
                    if (kt == 0 || d == 0) continue;  // no coupling
                    auto& bv = vectors[static_cast<size_t>(l)].at(t);
                    row.kt = kt;
                    row.m.assign(static_cast<size_t>(kt * d), DataType{0});
                    for (int64_t g : *row.groups) sketch_vectors(n, g, bv, d, row.m.data());
                } catch (...) {
#pragma omp critical(extract_compressed_sketches_failure)
                    if (!failure) failure = std::current_exception();
                }
            }
            if (failure) std::rethrow_exception(failure);
            std::vector<std::vector<int64_t>> meta(static_cast<size_t>(p.ssize));
            std::vector<std::vector<DataType>> data(static_cast<size_t>(p.ssize));
            for (const Rows& row : rows) {
                if (row.m.empty()) continue;
                const auto& [n, l, s, t] = row.key;
                const int to = tree->levels[static_cast<size_t>(l)].owner_rank_of(s);
                meta[static_cast<size_t>(to)].insert(meta[static_cast<size_t>(to)].end(), {n, l, s, t, row.kt});
                data[static_cast<size_t>(to)].insert(data[static_cast<size_t>(to)].end(), row.m.begin(), row.m.end());
            }
            rows.clear();
            lap(fmm::ExtractionProfile::kCompute);
            const auto rmeta = exchange(tree->comm, meta, MPI_INT64_T);
            const auto rdata = exchange(tree->comm, data, p.value_type);
            lap(fmm::ExtractionProfile::kColumns);
            // ---- on s's owner: W(n, l, s) = sum over the boxes t of
            // (K(S_s, S_t) M)^T, d x k_s
            struct Received {
                int64_t t, kt;
                const DataType* m;
            };
            std::map<std::tuple<int64_t, int, int64_t>, std::vector<Received>> received;  // (n, l, s) -> its M
            for (int from = 0; from < p.ssize; ++from) {
                const auto& mv = rmeta[static_cast<size_t>(from)];
                const auto& dv = rdata[static_cast<size_t>(from)];
                size_t di = 0;
                for (size_t mi = 0; mi + 5 <= mv.size(); mi += 5) {
                    const int64_t n = mv[mi], s = mv[mi + 2], t = mv[mi + 3], kt = mv[mi + 4];
                    const int l = static_cast<int>(mv[mi + 1]);
                    received[{n, l, s}].push_back(Received{t, kt, dv.data() + di});
                    di += static_cast<size_t>(kt * dims[n]);
                }
            }
            std::vector<std::pair<std::tuple<int64_t, int, int64_t>, const std::vector<Received>*>> targets;
            for (const auto& [key, list] : received) targets.emplace_back(key, &list);
            std::vector<std::pair<int64_t, std::vector<DataType>>> formed(targets.size());
#pragma omp parallel for schedule(dynamic)
            for (size_t q = 0; q < targets.size(); ++q) {
                if (failure) continue;
                try {
                    const auto& [n, l, s] = targets[q].first;
                    auto* sb = tree->levels[static_cast<size_t>(l)].find_local_box(s);
                    if (sb == nullptr) throw std::logic_error("extract_compressed_sketches: a column box is not local");
                    const int64_t ks = static_cast<int64_t>(sb->skeleton_indices.size()), d = dims[n];
                    if (ks == 0) continue;
                    auto& [k, w] = formed[q];
                    for (const Received& r : *targets[q].second) {
                        const H2Block<DataType>* block = nullptr;
                        for (const auto& b : sb->h2_interaction_blocks) {
                            if (b.source_morton == r.t) block = &b;
                        }
                        if (block == nullptr) continue;  // no coupling
                        const auto& K = block->matrix;
                        if (K.rows != ks || K.cols != r.kt) {
                            throw std::logic_error("extract_compressed_sketches: inconsistent transposed far data");
                        }
                        if (w.empty()) {
                            k = ks;
                            w.assign(static_cast<size_t>(d * ks), DataType{0});
                        }
                        // W += M^T K^T
                        const int m = static_cast<int>(d), nn = static_cast<int>(ks), kk = static_cast<int>(r.kt);
                        const int lda = static_cast<int>(r.kt), ldb = static_cast<int>(K.lda), ldc = static_cast<int>(d);
                        const DataType one{1};
                        const char tr = 'T';
                        if (m > 0 && nn > 0 && kk > 0) {
                            gemm_(&tr, &tr, &m, &nn, &kk, &one, r.m, &lda, K.data.data(), &ldb, &one, w.data(), &ldc);
                        }
                    }
                } catch (...) {
#pragma omp critical(extract_compressed_sketches_failure)
                    if (!failure) failure = std::current_exception();
                }
            }
            if (failure) std::rethrow_exception(failure);
            for (size_t q = 0; q < targets.size(); ++q) {
                if (!formed[q].second.empty()) sums[targets[q].first] = std::move(formed[q]);
            }
            lap(fmm::ExtractionProfile::kCompute);
        }

        // ---- the columns: Y = W U_s for each column group of n under s at l
        std::map<std::tuple<int64_t, int, int64_t>, std::vector<int64_t>> col_groups;  // (n, l, s) -> groups
        for (const Task& task : tasks) {
            if (task.near) continue;
            const auto it = sums.find({task.inter, task.level, task.s});
            if (it == sums.end()) continue;
            auto& g = col_groups[{task.inter, task.level, task.s}];
            if (std::find(g.begin(), g.end(), task.cgroup) == g.end()) g.push_back(task.cgroup);
        }
        struct Columns {
            int64_t n;
            int l;
            int64_t s, cgroup;
            const std::pair<int64_t, std::vector<DataType>>* w;
        };
        std::vector<Columns> work;
        for (const auto& [key, groups] : col_groups) {
            const auto& w = sums.at(key);
            for (int64_t g : groups) work.push_back(Columns{std::get<0>(key), std::get<1>(key), std::get<2>(key), g, &w});
        }
        // (near: on the owner of the stored block's row box, t, or with
        // transposed s)
        std::vector<size_t> near_tasks;
        auto& leaf_level = tree->levels[static_cast<size_t>(p.leaf)];
        for (size_t i = 0; i < tasks.size(); ++i) {
            if (!tasks[i].near) continue;
            if ((transposed ? leaf_level.owner_rank_of(tasks[i].s) : tasks[i].computer) == srank) near_tasks.push_back(i);
        }
        pieces.resize(work.size() + near_tasks.size());
#pragma omp parallel for schedule(dynamic)
        for (size_t q = 0; q < pieces.size(); ++q) {
            if (failure) continue;
            try {
                Piece& piece = pieces[q];
                if (q < work.size()) {
                    const Columns& c = work[q];
                    const int64_t d = dims[c.n], ks = c.w->first, nc = p.col_count(c.cgroup);
                    piece = Piece{c.n, c.cgroup, d, nc, std::vector<DataType>(static_cast<size_t>(d * nc), DataType{0})};
                    auto& bv = vectors[static_cast<size_t>(c.l)].at(c.s);
                    std::vector<DataType> us(static_cast<size_t>(ks * nc));
                    int64_t j = 0;
                    p.col_points(c.cgroup, [&](int64_t pt, int64_t) {
                        std::copy_n(bv.of(pt), ks, us.begin() + j * ks);
                        ++j;
                    });
                    const int m = static_cast<int>(d), nn = static_cast<int>(nc), k = static_cast<int>(ks);
                    const int lda = static_cast<int>(d), ldb = static_cast<int>(ks), ldc = static_cast<int>(d);
                    const DataType one{1}, zero{0};
                    const char no = 'N';
                    if (m > 0 && nn > 0 && k > 0) {
                        gemm_(&no, &no, &m, &nn, &k, &one, c.w->second.data(), &lda, us.data(), &ldb, &zero, piece.y.data(), &ldc);
                    }
                    continue;
                }
                // near: Omega's columns of the rows times the stored block's rows
                const Task& task = tasks[near_tasks[q - work.size()]];
                const int64_t d = dims[task.inter], nc = p.col_count(task.cgroup);
                piece = Piece{task.inter, task.cgroup, d, nc, std::vector<DataType>(static_cast<size_t>(d * nc), DataType{0})};
                auto& level = tree->levels[static_cast<size_t>(task.level)];
                if (transposed) {
                    // K(B_s, B_t), stored on s: Omega's columns of the rows
                    // (t's points) times its columns
                    auto* sb = level.find_local_box(task.s);
                    if (sb == nullptr) throw std::logic_error("extract_compressed_sketches: a task's box is not local");
                    std::unordered_map<int64_t, int64_t> column_of;
                    const auto& K = near_block(level, sb, task.t, column_of).matrix;
                    std::vector<int64_t> spos;
                    p.col_points(task.cgroup, [&](int64_t pt, int64_t) {
                        spos.push_back(p.leaf_position[static_cast<size_t>(pt)]);
                    });
                    p.row_points(task.rgroup, [&](int64_t pt, int64_t stacked) {
                        const int64_t c = column_of.at(p.distinct[static_cast<size_t>(pt)]);
                        for (int64_t qq = row_omega[static_cast<size_t>(stacked)]; qq < row_omega[static_cast<size_t>(stacked) + 1]; ++qq) {
                            const int64_t e = omega_of[static_cast<size_t>(qq)];
                            const DataType w = static_cast<DataType>(omega_val[e]);
                            DataType* y = piece.y.data() + omega_dest[e];
                            for (int64_t j = 0; j < nc; ++j) {
                                y[j * d] += w * K.data[static_cast<size_t>(spos[static_cast<size_t>(j)] + c * K.lda)];
                            }
                        }
                    });
                    continue;
                }
                auto* tb = level.find_local_box(task.t);
                if (tb == nullptr) throw std::logic_error("extract_compressed_sketches: a task's box is not local");
                std::unordered_map<int64_t, int64_t> column_of;
                const auto& K = near_block(level, tb, task.s, column_of).matrix;
                std::vector<int64_t> cpos;
                p.col_points(task.cgroup, [&](int64_t pt, int64_t) {
                    cpos.push_back(column_of.at(p.distinct[static_cast<size_t>(pt)]));
                });
                p.row_points(task.rgroup, [&](int64_t pt, int64_t stacked) {
                    const int64_t r = p.leaf_position[static_cast<size_t>(pt)];
                    for (int64_t qq = row_omega[static_cast<size_t>(stacked)]; qq < row_omega[static_cast<size_t>(stacked) + 1]; ++qq) {
                        const int64_t e = omega_of[static_cast<size_t>(qq)];
                        const DataType w = static_cast<DataType>(omega_val[e]);
                        DataType* y = piece.y.data() + omega_dest[e];
                        for (int64_t j = 0; j < nc; ++j) {
                            y[j * d] += w * K.data[static_cast<size_t>(r + cpos[static_cast<size_t>(j)] * K.lda)];
                        }
                    }
                });
            } catch (...) {
#pragma omp critical(extract_compressed_sketches_failure)
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
        for (const Piece& piece : pieces) profile.entries_computed += static_cast<double>(piece.y.size());
        lap(fmm::ExtractionProfile::kCompute);
    }

    // ---- the pieces to the requests' owners, summed into their columns
    std::vector<std::vector<int64_t>> meta(static_cast<size_t>(csize));
    std::vector<std::vector<DataType>> data(static_cast<size_t>(csize));
    for (const Piece& piece : pieces) {
        const int to = owner[static_cast<size_t>(piece.n)];
        meta[static_cast<size_t>(to)].insert(meta[static_cast<size_t>(to)].end(), {piece.n, piece.cgroup, piece.d, piece.nc});
        data[static_cast<size_t>(to)].insert(data[static_cast<size_t>(to)].end(), piece.y.begin(), piece.y.end());
        profile.bytes_routed += static_cast<double>(piece.y.size() * sizeof(DataType));
    }
    pieces.clear();
    lap(fmm::ExtractionProfile::kPack);
    const auto rmeta = exchange(caller_comm, meta, MPI_INT64_T, &profile.route_mpi);
    const auto rdata = exchange(caller_comm, data, p.value_type, &profile.route_mpi);
    lap(fmm::ExtractionProfile::kRoute);
    std::vector<int64_t> base(static_cast<size_t>(nsketch), -1);
    int64_t total = 0;
    for (int64_t n = 0; n < nsketch; ++n) {
        if (owner[static_cast<size_t>(n)] != crank) continue;
        base[static_cast<size_t>(n)] = total;
        total += static_cast<int64_t>(dims[n]) * colidx[n];
    }
    std::fill(alldat_loc, alldat_loc + total, DataType{0});
    const auto& cg = p.cg;
    for (int from = 0; from < csize; ++from) {
        const auto& mv = rmeta[static_cast<size_t>(from)];
        const auto& dv = rdata[static_cast<size_t>(from)];
        size_t di = 0;
        for (size_t mi = 0; mi + 4 <= mv.size(); mi += 4) {
            const int64_t n = mv[mi], b = mv[mi + 1], d = mv[mi + 2], nc = mv[mi + 3];
            if (base[static_cast<size_t>(n)] < 0) throw std::logic_error("extract_compressed_sketches: a piece for another rank");
            for (int64_t j = 0; j < nc; ++j) {
                const int64_t col = cg.perm[static_cast<size_t>(cg.start[static_cast<size_t>(b)] + j)] - p.col_offset[static_cast<size_t>(n)];
                DataType* y = alldat_loc + base[static_cast<size_t>(n)] + col * d;
                const DataType* v = dv.data() + di + static_cast<size_t>(j * d);
                for (int64_t i = 0; i < d; ++i) y[i] += v[i];
            }
            di += static_cast<size_t>(d * nc);
        }
    }
    lap(fmm::ExtractionProfile::kUnpack);
}

}  // namespace color_unstructured
}  // namespace butterfly
