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
// take them, 1-based.  The matrix is transpose-symmetric (option%sym=1).

#include "core/butterfly_types.hpp"
#include "core/blas_declare.hpp"
#include "core/butterfly_matvec.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
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

// One all-to-all of per-rank buffers of T over comm.
template<typename T>
std::vector<std::vector<T>> exchange(MPI_Comm comm, const std::vector<std::vector<T>>& send, MPI_Datatype type) {
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
    MPI_Alltoall(sc.data(), 1, MPI_INT, rc.data(), 1, MPI_INT, comm);
    size_t rtotal = 0;
    for (int p = 0; p < size; ++p) {
        rd[static_cast<size_t>(p)] = checked_count(rtotal, "a receive buffer");
        rtotal += static_cast<size_t>(rc[static_cast<size_t>(p)]);
    }
    std::vector<T> rbuf(rtotal);
    MPI_Alltoallv(sbuf.data(), sc.data(), sd.data(), type, rbuf.data(), rc.data(), rd.data(), type, comm);
    std::vector<std::vector<T>> recv(static_cast<size_t>(size));
    for (int p = 0; p < size; ++p) {
        recv[static_cast<size_t>(p)].assign(rbuf.begin() + rd[static_cast<size_t>(p)],
                                            rbuf.begin() + rd[static_cast<size_t>(p)] + rc[static_cast<size_t>(p)]);
    }
    return recv;
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
    int crank = 0, csize = 1;
    MPI_Comm_rank(caller_comm, &crank);
    MPI_Comm_size(caller_comm, &csize);
    const MPI_Datatype value_type = fmm::mpi_datatype_for<DataType>();

    // the source's ranks among the caller's
    ParallelTree<CoordType, DataType>* tree = source != nullptr ? source->tree.get() : nullptr;
    int srank = -1;
    if (tree != nullptr) MPI_Comm_rank(tree->comm, &srank);
    std::vector<int> source_rank_of(static_cast<size_t>(csize));
    MPI_Allgather(&srank, 1, MPI_INT, source_rank_of.data(), 1, MPI_INT, caller_comm);
    int nsource = 0;
    for (int r : source_rank_of) nsource += r >= 0 ? 1 : 0;
    std::vector<int> caller_rank_of(static_cast<size_t>(std::max(nsource, 1)), -1);
    for (int c = 0; c < csize; ++c) {
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
        int ssize = 0;
        MPI_Comm_size(tree->comm, &ssize);
        if (ssize != nsource) {
            throw std::invalid_argument("extract_compressed_blocks: the source's ranks are not a subset of the caller's");
        }
    }

    // the intersections' points: the distinct indices, and each stacked
    // index's place among them
    std::vector<int64_t> row_offset(static_cast<size_t>(ninter) + 1, 0), col_offset(static_cast<size_t>(ninter) + 1, 0);
    for (int64_t n = 0; n < ninter; ++n) {
        if (rowidx[n] < 0 || colidx[n] < 0) throw std::invalid_argument("extract_compressed_blocks: a negative size");
        row_offset[static_cast<size_t>(n) + 1] = row_offset[static_cast<size_t>(n)] + rowidx[n];
        col_offset[static_cast<size_t>(n) + 1] = col_offset[static_cast<size_t>(n)] + colidx[n];
    }
    const int64_t nrows = row_offset.back(), ncols = col_offset.back();
    std::vector<int64_t> distinct;
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
    std::vector<int64_t> row_point(static_cast<size_t>(nrows)), col_point(static_cast<size_t>(ncols));
    for (int64_t i = 0; i < nrows; ++i) row_point[static_cast<size_t>(i)] = place(static_cast<int64_t>(allrows[i]) - 1);
    for (int64_t j = 0; j < ncols; ++j) col_point[static_cast<size_t>(j)] = place(static_cast<int64_t>(allcols[j]) - 1);

    // each point's leaf box and its position there, from the leaf owners
    std::vector<int64_t> leaf_box(static_cast<size_t>(npoints), -1), leaf_position(static_cast<size_t>(npoints), -1);
    int dimension = 0, leaf = 0;
    if (tree != nullptr) {
        dimension = tree->dimension;
        leaf = tree->num_levels - 1;
        auto& level = tree->levels[static_cast<size_t>(leaf)];
        if (level.is_process_active) {
            for (auto& box : level.local_boxes) {
                for (size_t p = 0; p < box.point_indices.size(); ++p) {
                    const auto it = std::lower_bound(distinct.begin(), distinct.end(), box.point_indices[p]);
                    if (it == distinct.end() || *it != box.point_indices[p]) continue;
                    const size_t d = static_cast<size_t>(it - distinct.begin());
                    leaf_box[d] = box.morton_index;
                    leaf_position[d] = static_cast<int64_t>(p);
                }
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, leaf_box.data(), static_cast<int>(npoints), MPI_INT64_T, MPI_MAX, caller_comm);
    MPI_Allreduce(MPI_IN_PLACE, leaf_position.data(), static_cast<int>(npoints), MPI_INT64_T, MPI_MAX, caller_comm);
    for (int64_t d = 0; d < npoints; ++d) {
        if (leaf_box[static_cast<size_t>(d)] < 0) {
            throw std::invalid_argument("extract_compressed_blocks: index " +
                                        std::to_string(distinct[static_cast<size_t>(d)] + 1) +
                                        " is not a point of the source");
        }
    }
    const Groups rg = make_groups(ninter, row_offset, row_point, leaf_box);
    const Groups cg = make_groups(ninter, col_offset, col_point, leaf_box);

    // ---- the source: each task (an intersection's row group by column
    // group) on the owner of its row side's box, the values per task
    struct Task {
        int64_t inter, rgroup, cgroup;
        int level;     // the leaf for a near task, else the meeting level
        bool near;
        int64_t t, s;  // the row and column sides' boxes of that level
        int computer;  // its source rank
    };
    std::vector<Task> tasks;
    std::vector<std::vector<DataType>> values;  // of the tasks this rank computes
    if (tree != nullptr) {
        for (int64_t n = 0; n < ninter; ++n) {
            for (int64_t a = rg.first[static_cast<size_t>(n)]; a < rg.first[static_cast<size_t>(n) + 1]; ++a) {
                for (int64_t b = cg.first[static_cast<size_t>(n)]; b < cg.first[static_cast<size_t>(n) + 1]; ++b) {
                    Task task{n, a, b, leaf, true, rg.morton[static_cast<size_t>(a)], cg.morton[static_cast<size_t>(b)], -1};
                    if (task.t != task.s && !boxes_adjacent(dimension, task.t, task.s)) {
                        task.near = false;
                        task.level = meeting_level(dimension, leaf, task.t, task.s);
                        task.t >>= dimension * (leaf - task.level);
                        task.s >>= dimension * (leaf - task.level);
                    }
                    task.computer = tree->levels[static_cast<size_t>(task.level)].owner_rank_of(task.t);
                    tasks.push_back(task);
                }
            }
        }
        values.resize(tasks.size());

        // the coarsest level each point's interpolation vector is needed at
        std::vector<int> need(static_cast<size_t>(npoints), leaf + 1);
        auto group_points = [&](const Groups& g, const std::vector<int64_t>& point_of, int64_t group, auto&& f) {
            for (int64_t i = g.start[static_cast<size_t>(group)]; i < g.start[static_cast<size_t>(group) + 1]; ++i) {
                f(point_of[static_cast<size_t>(g.perm[static_cast<size_t>(i)])], i);
            }
        };
        int coarsest = leaf + 1;
        for (const Task& task : tasks) {
            if (task.near) continue;
            coarsest = std::min(coarsest, task.level);
            group_points(rg, row_point, task.rgroup, [&](int64_t d, int64_t) {
                need[static_cast<size_t>(d)] = std::min(need[static_cast<size_t>(d)], task.level);
            });
            group_points(cg, col_point, task.cgroup, [&](int64_t d, int64_t) {
                need[static_cast<size_t>(d)] = std::min(need[static_cast<size_t>(d)], task.level);
            });
        }

        // ---- the interpolation vectors, leaf up
        std::vector<std::unordered_map<int64_t, BoxVectors<DataType>>> vectors(static_cast<size_t>(leaf) + 1);
        {
            auto& level = tree->levels[static_cast<size_t>(leaf)];
            for (int64_t d = 0; d < npoints; ++d) {
                if (need[static_cast<size_t>(d)] > leaf || level.owner_rank_of(leaf_box[static_cast<size_t>(d)]) != srank) {
                    continue;
                }
                auto* box = level.find_local_box(leaf_box[static_cast<size_t>(d)]);
                if (box == nullptr) throw std::logic_error("extract_compressed_blocks: a leaf box is not local");
                auto& bv = vectors[static_cast<size_t>(leaf)][box->morton_index];
                bv.k = static_cast<int64_t>(box->skeleton_indices.size());
                bv.points.push_back(d);
                const int64_t p = leaf_position[static_cast<size_t>(d)];
                std::vector<DataType> col(static_cast<size_t>(bv.k), DataType{0});
                const auto sit = std::find(box->skeleton_indices.begin(), box->skeleton_indices.end(), p);
                if (sit != box->skeleton_indices.end()) {
                    col[static_cast<size_t>(sit - box->skeleton_indices.begin())] = DataType{1};
                } else {
                    const auto rit = std::find(box->redundant_indices.begin(), box->redundant_indices.end(), p);
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
        int ssize = 1;
        MPI_Comm_size(tree->comm, &ssize);
        for (int l = leaf; l > std::max(coarsest, 2); --l) {
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
            const auto rdata = exchange(tree->comm, data, value_type);
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

        // ---- the column sides' vectors to the row sides' owners, by level
        std::vector<std::unordered_map<int64_t, std::unordered_map<int64_t, std::vector<DataType>>>> remote(
            static_cast<size_t>(leaf) + 1);  // level -> box -> point -> vector
        for (int l = leaf; l >= std::max(coarsest, 2); --l) {
            auto& level = tree->levels[static_cast<size_t>(l)];
            std::vector<std::vector<int64_t>> meta(static_cast<size_t>(ssize));
            std::vector<std::vector<DataType>> data(static_cast<size_t>(ssize));
            std::vector<std::unordered_map<int64_t, std::vector<int64_t>>> sent(static_cast<size_t>(ssize));
            for (const Task& task : tasks) {
                if (task.near || task.level != l || task.computer == srank || level.owner_rank_of(task.s) != srank) continue;
                auto& bv = vectors[static_cast<size_t>(l)].at(task.s);
                auto& done = sent[static_cast<size_t>(task.computer)][task.s];
                group_points(cg, col_point, task.cgroup, [&](int64_t d, int64_t) {
                    if (std::find(done.begin(), done.end(), d) != done.end()) return;
                    done.push_back(d);
                    auto& mv = meta[static_cast<size_t>(task.computer)];
                    mv.insert(mv.end(), {task.s, d, bv.k});
                    const DataType* u = bv.of(d);
                    data[static_cast<size_t>(task.computer)].insert(data[static_cast<size_t>(task.computer)].end(), u, u + bv.k);
                });
            }
            const auto rmeta = exchange(tree->comm, meta, MPI_INT64_T);
            const auto rdata = exchange(tree->comm, data, value_type);
            for (int from = 0; from < ssize; ++from) {
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

        // ---- this rank's tasks
        std::vector<size_t> mine;
        for (size_t i = 0; i < tasks.size(); ++i) {
            if (tasks[i].computer == srank) mine.push_back(i);
        }
        // (the vectors' point maps, before the threads read them)
        for (auto& per_level : vectors) {
            for (auto& [m, bv] : per_level) {
                if (!bv.points.empty()) (void)bv.of(bv.points.front());
            }
        }
        std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic)
        for (size_t q = 0; q < mine.size(); ++q) {
            if (failure) continue;
            try {
                const Task& task = tasks[mine[q]];
                const int64_t nr = rg.start[static_cast<size_t>(task.rgroup) + 1] - rg.start[static_cast<size_t>(task.rgroup)];
                const int64_t nc = cg.start[static_cast<size_t>(task.cgroup) + 1] - cg.start[static_cast<size_t>(task.cgroup)];
                std::vector<DataType>& out = values[mine[q]];
                out.assign(static_cast<size_t>(nr * nc), DataType{0});
                auto& level = tree->levels[static_cast<size_t>(task.level)];
                auto* tb = level.find_local_box(task.t);
                if (tb == nullptr) throw std::logic_error("extract_compressed_blocks: a task's box is not local");
                if (task.near) {
                    const H2Block<DataType>* block = nullptr;
                    for (const auto& b : tb->h2_near_blocks) {
                        if (b.source_morton == task.s) block = &b;
                    }
                    if (block == nullptr) throw std::logic_error("extract_compressed_blocks: a missing near block");
                    const auto columns = h2_box_global_indices(level, task.s, false);
                    std::unordered_map<int64_t, int64_t> column_of;
                    for (size_t j = 0; j < columns.size(); ++j) column_of.emplace(columns[j], static_cast<int64_t>(j));
                    std::vector<int64_t> rpos, cpos;
                    group_points(rg, row_point, task.rgroup, [&](int64_t d, int64_t) { rpos.push_back(leaf_position[static_cast<size_t>(d)]); });
                    group_points(cg, col_point, task.cgroup, [&](int64_t d, int64_t) {
                        cpos.push_back(column_of.at(distinct[static_cast<size_t>(d)]));
                    });
                    const auto& K = block->matrix;
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
                group_points(rg, row_point, task.rgroup, [&](int64_t d, int64_t) {
                    std::copy_n(tit->second.of(d), kt, ut.begin() + i * kt);
                    ++i;
                });
                int64_t j = 0;
                const bool local_s = level.owner_rank_of(task.s) == srank;
                group_points(cg, col_point, task.cgroup, [&](int64_t d, int64_t) {
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
    }
    const auto rmeta = exchange(caller_comm, meta, MPI_INT64_T);
    const auto rdata = exchange(caller_comm, data, value_type);
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
            const int64_t nr = rg.start[static_cast<size_t>(a) + 1] - rg.start[static_cast<size_t>(a)];
            const int64_t nc = cg.start[static_cast<size_t>(b) + 1] - cg.start[static_cast<size_t>(b)];
            for (int64_t j = 0; j < nc; ++j) {
                const int64_t col = cg.perm[static_cast<size_t>(cg.start[static_cast<size_t>(b)] + j)] - col_offset[static_cast<size_t>(n)];
                for (int64_t i = 0; i < nr; ++i) {
                    const int64_t row = rg.perm[static_cast<size_t>(rg.start[static_cast<size_t>(a)] + i)] - row_offset[static_cast<size_t>(n)];
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
}

}  // namespace color_unstructured
}  // namespace butterfly
