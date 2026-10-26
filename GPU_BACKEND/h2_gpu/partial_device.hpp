#pragma once
// =============================================================================
// H2_partial on the device: the Schur complement's data and its product.
//
// The F22 records of the sparsified boxes (core/tree.hpp PartialRecord: the
// Bitcomp chunks under H2_compress_factors, else the raw entries) and the
// root's CB block go to the device at the first CB product after a
// factorization and stay for the later ones, when they fit: decided on every
// rank alike (each rank's records, root block and decompression scratch
// beside what the device holds), else the host product runs.
//
// partial_schur_multiply (unstructured/partial.hpp) keeps the up and down
// passes, the exchanges and the sums into the levels' accumulators on the
// host.  Per level it hands the device the raw vectors of its boxes (and the
// received ones of other ranks' neighbours) and each sparsified box's p and
// w, and takes back, per sparsified box,
//     gR = X_RR p + X_RS w + X_RN rn,   gS = X_SR p,   t_q = X_NR[q] p
// (rn: the neighbours' raw values at the box's X_RN columns, gathered on the
// device).  Compressed records are decompressed a segment of boxes at a time
// into one scratch (Bitcomp's batched decompression) before the segment's
// GEMM batches.  The root's block is one GEMM.
// =============================================================================

#ifdef H2_HAVE_GPU

#include "core/partial_record_codec.hpp"
#include "device_kernels.hpp"
#include "device_solve.hpp"

#ifdef H2_HAVE_NVCOMP
#include <nvcomp/native/bitcomp.h>
#endif

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// A record on the device: its entries (ld rows), or its compressed chunks
// end to end (chunk c at dev + chunk_off[c]), decompressed into the
// segment's scratch at scratch_off.
struct PartialDeviceRecord {
    const char* dev = nullptr;
    int64_t rows = 0, cols = 0;
    std::vector<size_t> chunk_off;
    size_t scratch_off = 0;
    bool compressed() const { return !chunk_off.empty(); }
};

// A sparsified box of a level: its records and its places in the level's
// layouts, in columns of the product (p then w; rn; gR, gS, then t_q in
// one-hop order)
struct PartialDeviceBox {
    int64_t bi = 0;                // its local index at the level
    int r = 0, kp = 0, total = 0;  // R22, S' (after the step), X_RN's columns
    PartialDeviceRecord rr, rs, sr, rn;
    std::vector<PartialDeviceRecord> nr;  // per one-hop slot (rows 0: none)
    int64_t pw_off = 0, rn_off = 0, out_off = 0;
    std::vector<int64_t> t_off;           // per one-hop slot, t_q's place after gR and gS
};

// Boxes whose records one decompression batch takes into the scratch
struct PartialDeviceSegment {
    int first = 0, count = 0;
    size_t nchunks = 0;
    const void** d_srcs = nullptr;
    void** d_dsts = nullptr;
#ifdef H2_HAVE_NVCOMP
    bitcompHandle_t plan{};
    bool have_plan = false;
#endif
};

struct PartialDeviceLevel {
    std::vector<PartialDeviceBox> boxes;  // the level's sparsified local boxes, in local order
    std::vector<int> box_of;              // local box index -> boxes (-1)
    int64_t pw_units = 0, rn_units = 0, out_units = 0;
    // rn's rows from a neighbour's raw vector: rows [col0, col0 + m) of box
    // q's rn from the neighbour's slots (a local box, or another rank's)
    struct Gather {
        int q = 0, col0 = 0, m = 0;
        bool remote = false;
        int64_t neighbour = 0;  // local index, or Morton index
        const int* d_slots = nullptr;
    };
    std::vector<Gather> gathers;
    int max_gather_m = 0;
    std::vector<PartialDeviceSegment> segments;
};

struct PartialDeviceStore {
    const void* tree = nullptr;
    bool decided = false;
    bool usable = false;
    std::vector<char*> blocks;  // resident: records, root block, lists, segment tables
    size_t bytes = 0;
    std::vector<PartialDeviceLevel> levels;
    const char* root_cb = nullptr;
    int64_t root_n = 0;
    char* scratch = nullptr;
    size_t scratch_bytes = 0;
    int64_t max_columns = 1;  // columns one pass takes (its vectors beside the resident data)
    MetaBuilder meta;
    DeviceBuffer meta_device;
    PinnedBuffer host_in, host_out;

    void release() {
        DeviceHeap& heap = DeviceHeap::instance();
#ifdef H2_HAVE_NVCOMP
        for (auto& lv : levels) {
            for (auto& s : lv.segments) {
                if (s.have_plan) bitcompDestroyPlan(s.plan);
                s.have_plan = false;
            }
        }
#endif
        for (char* b : blocks) heap.free(b);
        blocks.clear();
        if (scratch != nullptr) heap.free(scratch);
        scratch = nullptr;
        scratch_bytes = 0;
        levels.clear();
        root_cb = nullptr;
        root_n = 0;
        bytes = 0;
        max_columns = 1;
        tree = nullptr;
        decided = false;
        usable = false;
    }
};

// The store of the active operator (one per operator, as the solve's): its
// device data goes with the operator, and with a new factorization of it
// (invalidate_device_solve).
inline PartialDeviceStore& partial_device_store() {
    static std::unordered_map<const void*, PartialDeviceStore>* const stores = [] {
        auto* s = new std::unordered_map<const void*, PartialDeviceStore>;  // (never freed: see DeviceHeap)
        OperatorContext& c = operator_context();
        c.releasers.push_back([s](const void* tree) {
            auto it = s->find(tree);
            if (it == s->end()) return;
            it->second.release();
            s->erase(it);
        });
        c.holders.push_back([s](const void* tree) -> size_t {
            auto it = s->find(tree);
            return it == s->end() ? 0 : it->second.bytes;
        });
        device_store_invalidators().push_back([s] {
            auto it = s->find(operator_context().active);
            if (it != s->end()) it->second.release();
        });
        return s;
    }();
    return (*stores)[operator_context().active];
}

namespace partial_detail {

// Host pieces copied to device addresses through two pinned staging
// buffers, each filled by all threads (pieces split at kStage bytes).
struct UploadPiece {
    char* dst;
    const char* src;
    size_t bytes;
};
inline void staged_upload(const std::vector<UploadPiece>& pieces, cudaStream_t stream) {
    constexpr size_t kStage = size_t{128} << 20;
    struct Part { char* dst; const char* src; size_t bytes, at; };
    std::vector<std::vector<Part>> stages(1);
    size_t fill = 0;
    for (const UploadPiece& p : pieces) {
        for (size_t done = 0; done < p.bytes;) {
            if (fill == kStage) {
                stages.emplace_back();
                fill = 0;
            }
            const size_t n = std::min(p.bytes - done, kStage - fill);
            stages.back().push_back(Part{p.dst + done, p.src + done, n, fill});
            fill += n;
            done += n;
        }
    }
    cudaEvent_t done[2] = {nullptr, nullptr};
    for (cudaEvent_t& e : done) check_cuda(cudaEventCreateWithFlags(&e, cudaEventDisableTiming), "cudaEventCreate");
    bool used[2] = {false, false};
    for (size_t s = 0; s < stages.size(); ++s) {
        const int turn = static_cast<int>(s & 1);
        if (used[turn]) check_cuda(cudaEventSynchronize(done[turn]), "partial upload");
        char* h = static_cast<char*>(pinned_pool().staging[turn].reserve(kStage));
        const auto& parts = stages[s];
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t i = 0; i < static_cast<int64_t>(parts.size()); ++i) {
            const Part& p = parts[static_cast<size_t>(i)];
            std::memcpy(h + p.at, p.src, p.bytes);
        }
        for (const Part& p : parts) {
            check_cuda(cudaMemcpyAsync(p.dst, h + p.at, p.bytes, cudaMemcpyHostToDevice, stream), "partial upload");
        }
        check_cuda(cudaEventRecord(done[turn], stream), "cudaEventRecord");
        used[turn] = true;
    }
    check_cuda(cudaStreamSynchronize(stream), "partial upload");
    for (cudaEvent_t e : done) cudaEventDestroy(e);
}

}  // namespace partial_detail

// Whether the CB product of the partial factorization `tree` runs its record
// terms and root term on the device: the first call after a factorization
// decides (collectively) and uploads the records and the root block.
template<typename CoordType, typename DataType>
bool partial_device_ready(ParallelTree<CoordType, DataType>* tree, int verbosity) {
    activate_operator(tree);
    if (device_solve_suspended() || !color_gpu_enabled()) return false;
    PartialDeviceStore& st = partial_device_store();
    if (st.decided && st.tree == tree) return st.usable;
    st.release();
    st.tree = tree;
    st.decided = true;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    const size_t D = sizeof(DataType), I = sizeof(int);
    const int leaf = tree->num_levels - 1;
    int rank = 0;
    MPI_Comm_rank(tree->comm, &rank);
    const auto& pr = tree->partial_root;
    std::string reason;
    if (pr.cb_released) reason = "released";

    // ---- the plan: per level its sparsified boxes, their records, layouts and gathers
    st.levels.assign(static_cast<size_t>(leaf + 1), PartialDeviceLevel{});
    struct Placed { PartialDeviceRecord* rec; const PartialRecord<DataType>* host; };
    std::vector<Placed> placed;  // the records in device order
    size_t record_bytes = 0, list_ints = 0, largest_box_raw = 0, largest_level_raw = 0;
    bool compressed = false;
    std::vector<size_t> box_raw;  // per (level, box) in order: its records' entries' bytes
    double per_column = 0.0;      // the largest level's vectors a column
    for (int l = 2; l <= leaf && reason.empty(); ++l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        if (!level.is_process_active) continue;
        PartialDeviceLevel& L = st.levels[static_cast<size_t>(l)];
        L.box_of.assign(static_cast<size_t>(level.num_boxes_local), -1);
        int64_t points = 0;
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            points += box.num_points;
            const auto& f = box.f22;
            if (!f.sparsified || f.redundant.empty() || f.released) continue;
            PartialDeviceBox B;
            B.bi = bi;
            B.r = static_cast<int>(f.redundant.size());
            B.kp = static_cast<int>(box.skeleton_indices.size());
            B.total = static_cast<int>(f.X_RN.cols);
            const size_t nslots = box.one_hop.size();
            if (f.X_NR.size() != nslots || f.neighbor_nr_counts.size() != nslots || f.neighbor_f22_counts.size() != nslots) {
                reason = "a box's neighbour records do not match its neighbours";
                break;
            }
            auto shape = [&](const PartialRecord<DataType>& m, int64_t rows, int64_t cols) {
                if (rows * cols == 0) return true;
                return m.rows == rows && m.cols == cols && (m.compressed() || m.lda == rows);
            };
            bool ok = shape(f.X_RR, B.r, B.r) && shape(f.X_RS, B.r, B.kp) && shape(f.X_SR, B.kp, B.r) &&
                      shape(f.X_RN, B.r, B.total);
            int64_t tcols = 0;
            B.nr.resize(nslots);
            B.t_off.assign(nslots, 0);
            for (size_t q = 0; q < nslots; ++q) {
                const int64_t m = f.neighbor_nr_counts[q];
                ok = ok && shape(f.X_NR[q], m, B.r);
                B.t_off[q] = B.r + B.kp + tcols;
                tcols += m;
            }
            if (!ok) {
                reason = "a record's shape";
                break;
            }
            B.pw_off = L.pw_units;
            L.pw_units += B.r + B.kp;
            B.rn_off = L.rn_units;
            L.rn_units += B.total;
            B.out_off = L.out_units;
            L.out_units += B.r + B.kp + tcols;
            L.box_of[static_cast<size_t>(bi)] = static_cast<int>(L.boxes.size());
            L.boxes.push_back(std::move(B));
        }
        // (the boxes' addresses are final now)
        size_t level_raw = 0;
        for (PartialDeviceBox& B : L.boxes) {
            const auto& box = level.local_boxes[static_cast<size_t>(B.bi)];
            const auto& f = box.f22;
            size_t raw = 0;
            auto add = [&](PartialDeviceRecord& d, const PartialRecord<DataType>& h) {
                d.rows = h.rows;
                d.cols = h.cols;
                if (h.rows * h.cols == 0) return;
                if (h.compressed()) compressed = true;
                record_bytes += align_up(h.compressed() ? h.packed.size() : static_cast<size_t>(h.rows * h.cols) * D);
                raw += align_up(static_cast<size_t>(h.rows * h.cols) * D);
                placed.push_back(Placed{&d, &h});
            };
            add(B.rr, f.X_RR);
            add(B.rs, f.X_RS);
            add(B.sr, f.X_SR);
            add(B.rn, f.X_RN);
            for (size_t q = 0; q < B.nr.size(); ++q) add(B.nr[q], f.X_NR[q]);
            box_raw.push_back(raw);
            largest_box_raw = std::max(largest_box_raw, raw);
            level_raw += raw;
            for (size_t q = 0; q < box.one_hop.size(); ++q) list_ints += static_cast<size_t>(f.neighbor_f22_counts[q]);
        }
        largest_level_raw = std::max(largest_level_raw, level_raw);
        // vectors a column: the level's raw vectors (and as many again for the
        // received ones, at most), p|w, rn, and the outputs
        per_column = std::max(per_column, static_cast<double>(2 * points + L.pw_units + L.rn_units + L.out_units) * D);
    }
    size_t root_bytes = 0;
    if (pr.active && pr.n22 > 0) {
        if (pr.cb.rows != pr.n22 || pr.cb.cols != pr.n22) reason = "the root block's shape";
        root_bytes = align_up(static_cast<size_t>(pr.n22 * pr.n22) * D);
        per_column = std::max(per_column, static_cast<double>(2 * pr.n22) * D);
    }
    DeviceHeap& heap = DeviceHeap::instance();
    heap.ensure_initialized();
    const double room = static_cast<double>(heap.capacity() - heap.used()) - static_cast<double>(size_t{256} << 20);
    // the scratch: a segment of boxes' records uncompressed, an eighth of the
    // room left (at least the largest box; at most 2 GB, or a level's records)
    double scratch = 0.0;
    if (compressed) {
        const double left = room - static_cast<double>(record_bytes + root_bytes);
        scratch = std::max(static_cast<double>(largest_box_raw),
                           std::min({static_cast<double>(size_t{2} << 30), left / 8.0, static_cast<double>(largest_level_raw)}));
    }
    const double tables = static_cast<double>(list_ints * I) + 1e6;
    const double need = static_cast<double>(record_bytes + root_bytes) + scratch + tables + 4.0 * per_column;
    int bad = reason.empty() && need <= room ? 0 : 1;
    if (bad && reason.empty()) reason = "not enough device memory";
    MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, tree->comm);
    double sizes[3] = {static_cast<double>(record_bytes), static_cast<double>(root_bytes), need};
    MPI_Allreduce(MPI_IN_PLACE, sizes, 3, MPI_DOUBLE, MPI_MAX, tree->comm);
    if (bad) {
        if (rank == 0 && verbosity >= 0) {
            std::printf("GPU partial CB: host product (%s on some rank; records %.2f GB, root block %.2f GB, need %.2f GB "
                        "against %.2f GB on rank 0)\n", reason.empty() ? "another rank" : reason.c_str(), sizes[0] / 1e9,
                        sizes[1] / 1e9, sizes[2] / 1e9, room / 1e9);
            std::fflush(stdout);
        }
        st.levels.clear();
        return false;
    }

    // ---- the upload: the records end to end in blocks of up to 256 MB, then the root block
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    std::vector<partial_detail::UploadPiece> pieces;
    {
        constexpr size_t kBlock = size_t{256} << 20;
        std::vector<std::pair<size_t, size_t>> block_of(placed.size());  // (block, offset)
        std::vector<size_t> sizes_of;
        for (size_t i = 0; i < placed.size(); ++i) {
            const auto& h = *placed[i].host;
            const size_t b = align_up(h.compressed() ? h.packed.size() : static_cast<size_t>(h.rows * h.cols) * D);
            if (sizes_of.empty() || (sizes_of.back() > 0 && sizes_of.back() + b > kBlock)) sizes_of.push_back(0);
            block_of[i] = {sizes_of.size() - 1, sizes_of.back()};
            sizes_of.back() += b;
        }
        std::vector<char*> blocks(sizes_of.size());
        for (size_t b = 0; b < sizes_of.size(); ++b) {
            blocks[b] = heap.alloc_resident(std::max<size_t>(sizes_of[b], 1));
            st.blocks.push_back(blocks[b]);
            st.bytes += sizes_of[b];
        }
        for (size_t i = 0; i < placed.size(); ++i) {
            const auto& h = *placed[i].host;
            PartialDeviceRecord& d = *placed[i].rec;
            char* dst = blocks[block_of[i].first] + block_of[i].second;
            d.dev = dst;
            if (h.compressed()) {
                d.chunk_off.assign(h.strides.size(), 0);
                for (size_t c = 1; c < h.strides.size(); ++c) d.chunk_off[c] = d.chunk_off[c - 1] + h.strides[c - 1];
                pieces.push_back({dst, h.packed.data(), h.packed.size()});
            } else {
                pieces.push_back({dst, reinterpret_cast<const char*>(h.data.data()), static_cast<size_t>(h.rows * h.cols) * D});
            }
        }
    }
    if (root_bytes > 0) {
        char* rb = heap.alloc_resident(root_bytes);
        st.blocks.push_back(rb);
        st.bytes += root_bytes;
        const size_t col = static_cast<size_t>(pr.n22) * D;
        if (pr.cb.lda == pr.n22) {
            pieces.push_back({rb, reinterpret_cast<const char*>(pr.cb.data.data()), col * static_cast<size_t>(pr.n22)});
        } else {
            for (int64_t j = 0; j < pr.n22; ++j)
                pieces.push_back({rb + j * col, reinterpret_cast<const char*>(pr.cb.data.data() + j * pr.cb.lda), col});
        }
        st.root_cb = rb;
        st.root_n = pr.n22;
    }
    // the gathers' slot lists, and the gathers
    std::vector<int> lists;
    lists.reserve(list_ints);
    std::vector<std::vector<size_t>> list_at(static_cast<size_t>(leaf + 1));
    for (int l = 2; l <= leaf; ++l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        if (!level.is_process_active) continue;
        PartialDeviceLevel& L = st.levels[static_cast<size_t>(l)];
        for (size_t q = 0; q < L.boxes.size(); ++q) {
            const auto& box = level.local_boxes[static_cast<size_t>(L.boxes[q].bi)];
            const auto& f = box.f22;
            int col0 = 0;
            for (size_t s = 0; s < box.one_hop.size(); ++s) {
                const int m = static_cast<int>(f.neighbor_f22_counts[s]);
                if (m == 0) continue;
                PartialDeviceLevel::Gather g;
                g.q = static_cast<int>(q);
                g.col0 = col0;
                g.m = m;
                const int64_t nm = box.one_hop[s];
                g.remote = !level.is_box_on_process(nm);
                g.neighbour = g.remote ? nm : nm - level.local_morton_start;
                list_at[static_cast<size_t>(l)].push_back(lists.size());
                for (int i = 0; i < m; ++i) lists.push_back(static_cast<int>(f.neighbor_f22_slots[s][static_cast<size_t>(i)]));
                L.gathers.push_back(g);
                L.max_gather_m = std::max(L.max_gather_m, m);
                col0 += m;
            }
        }
    }
    int* d_lists = nullptr;
    if (!lists.empty()) {
        d_lists = reinterpret_cast<int*>(heap.alloc_resident(lists.size() * I));
        st.blocks.push_back(reinterpret_cast<char*>(d_lists));
        st.bytes += lists.size() * I;
        pieces.push_back({reinterpret_cast<char*>(d_lists), reinterpret_cast<const char*>(lists.data()), lists.size() * I});
    }
    const auto tu = clock::now();
    partial_detail::staged_upload(pieces, stream);
    const double upload_seconds = std::chrono::duration<double>(clock::now() - tu).count();
    for (int l = 2; l <= leaf; ++l) {
        PartialDeviceLevel& L = st.levels[static_cast<size_t>(l)];
        for (size_t g = 0; g < L.gathers.size(); ++g) L.gathers[g].d_slots = d_lists + list_at[static_cast<size_t>(l)][g];
    }

    // ---- the decompression segments: consecutive boxes whose records fit the scratch
    if (compressed) {
        st.scratch_bytes = static_cast<size_t>(scratch);
        st.scratch = heap.alloc_resident(std::max<size_t>(st.scratch_bytes, 1));
        st.bytes += st.scratch_bytes;
#ifdef H2_HAVE_NVCOMP
        const size_t chunk_elems = std::max<size_t>(1, kRecordChunkBytes / D);
        size_t bi_all = 0;
        for (int l = 2; l <= leaf; ++l) {
            PartialDeviceLevel& L = st.levels[static_cast<size_t>(l)];
            size_t used = st.scratch_bytes;  // (forces a new segment)
            for (size_t q = 0; q < L.boxes.size(); ++q, ++bi_all) {
                if (used + box_raw[bi_all] > st.scratch_bytes || L.segments.empty()) {
                    L.segments.emplace_back();
                    L.segments.back().first = static_cast<int>(q);
                    used = 0;
                }
                PartialDeviceSegment& s = L.segments.back();
                ++s.count;
                PartialDeviceBox& B = L.boxes[q];
                auto place = [&](PartialDeviceRecord& d) {
                    if (d.rows * d.cols == 0) return;
                    d.scratch_off = used;
                    used += align_up(static_cast<size_t>(d.rows * d.cols) * D);
                };
                place(B.rr);
                place(B.rs);
                place(B.sr);
                place(B.rn);
                for (auto& d : B.nr) place(d);
            }
            for (PartialDeviceSegment& s : L.segments) {
                std::vector<const void*> srcs;
                std::vector<void*> dsts;
                for (int q = s.first; q < s.first + s.count; ++q) {
                    PartialDeviceBox& B = L.boxes[static_cast<size_t>(q)];
                    auto chunks = [&](const PartialDeviceRecord& d) {
                        if (!d.compressed()) return;
                        for (size_t c = 0; c < d.chunk_off.size(); ++c) {
                            srcs.push_back(d.dev + d.chunk_off[c]);
                            dsts.push_back(st.scratch + d.scratch_off + c * chunk_elems * D);
                        }
                    };
                    chunks(B.rr);
                    chunks(B.rs);
                    chunks(B.sr);
                    chunks(B.rn);
                    for (const auto& d : B.nr) chunks(d);
                }
                s.nchunks = srcs.size();
                if (s.nchunks == 0) continue;
                char* t = heap.alloc_resident(2 * align_up(s.nchunks * sizeof(void*)));
                st.blocks.push_back(t);
                st.bytes += 2 * align_up(s.nchunks * sizeof(void*));
                s.d_srcs = reinterpret_cast<const void**>(t);
                s.d_dsts = reinterpret_cast<void**>(t + align_up(s.nchunks * sizeof(void*)));
                check_cuda(cudaMemcpyAsync(s.d_srcs, srcs.data(), s.nchunks * sizeof(void*), cudaMemcpyHostToDevice, stream),
                           "partial segment srcs");
                check_cuda(cudaMemcpyAsync(s.d_dsts, dsts.data(), s.nchunks * sizeof(void*), cudaMemcpyHostToDevice, stream),
                           "partial segment dsts");
                check_cuda(cudaStreamSynchronize(stream), "partial segment tables");
                if (bitcompCreateBatchPlanFromCompressedData(&s.plan, s.d_srcs, s.nchunks, nullptr, stream) != BITCOMP_SUCCESS) {
                    throw std::runtime_error("GPU partial CB: a decompression plan");
                }
                s.have_plan = true;
            }
        }
#else
        throw std::logic_error("GPU partial CB: compressed records without a compressor");
#endif
    } else {
        for (int l = 2; l <= leaf; ++l) {
            PartialDeviceLevel& L = st.levels[static_cast<size_t>(l)];
            if (L.boxes.empty()) continue;
            L.segments.emplace_back();
            L.segments.back().count = static_cast<int>(L.boxes.size());
        }
    }
    // columns a pass takes: its vectors in a quarter of what is left
    {
        double left = static_cast<double>(heap.capacity() - heap.used()) - static_cast<double>(size_t{256} << 20);
        double cols = per_column > 0.0 ? std::floor(left / 4.0 / per_column) : 1e9;
        MPI_Allreduce(MPI_IN_PLACE, &cols, 1, MPI_DOUBLE, MPI_MIN, tree->comm);
        st.max_columns = static_cast<int64_t>(std::max(1.0, std::min(cols, 1e9)));
    }
    st.usable = true;
    double seconds = std::chrono::duration<double>(clock::now() - t0).count();
    double held = static_cast<double>(st.bytes);
    MPI_Allreduce(MPI_IN_PLACE, &seconds, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
    MPI_Allreduce(MPI_IN_PLACE, &held, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
    if (rank == 0 && verbosity >= 0) {
        std::printf("GPU partial CB: the records%s and the root block on the device, %.2f GB per rank at most (records "
                    "%.2f, root %.2f), set up in %.2f s (upload %.2f s on rank 0); %lld columns a pass\n",
                    compressed ? " (compressed)" : "", held / 1e9, sizes[0] / 1e9, sizes[1] / 1e9, seconds, upload_seconds,
                    static_cast<long long>(st.max_columns));
        std::fflush(stdout);
    }
    return true;
}

// One level's record terms of the CB product (partial_device_ready true).
// x_local: per local box its raw vector (n x ncol); x_remote: the received
// raw vectors of other ranks' neighbours, Morton -> (n, vector); pw: per
// sparsified box (the level plan's order) p (r x ncol) then w (kp x ncol),
// at pw_off * ncol; out: per sparsified box gR (r x ncol), gS (kp x ncol),
// then t_q (m_q x ncol) per one-hop slot, at out_off * ncol.
template<typename CoordType, typename DataType>
void partial_device_terms(ParallelTree<CoordType, DataType>* tree, int level, int64_t ncol,
                          const std::vector<const DataType*>& x_local,
                          const std::unordered_map<int64_t, std::pair<int64_t, const DataType*>>& x_remote,
                          const DataType* pw, DataType* out) {
    using S = typename DeviceScalar<DataType>::type;
    PartialDeviceStore& st = partial_device_store();
    const PartialDeviceLevel& L = st.levels[static_cast<size_t>(level)];
    const auto& lvl = tree->levels[static_cast<size_t>(level)];
    if (L.boxes.empty()) return;
    const size_t D = sizeof(DataType);
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    DeviceHeap& heap = DeviceHeap::instance();
    // ---- the level vector: the local boxes, then the remote ones the gathers read
    std::vector<int64_t> xoff(x_local.size() + 1, 0);
    for (size_t b = 0; b < x_local.size(); ++b) xoff[b + 1] = xoff[b] + lvl.local_boxes[b].num_points;
    std::unordered_map<int64_t, int64_t> remote_off;
    int64_t xpoints = xoff.back();
    for (const auto& g : L.gathers) {
        if (!g.remote || remote_off.count(g.neighbour) != 0) continue;
        const auto it = x_remote.find(g.neighbour);
        if (it == x_remote.end()) throw std::runtime_error("GPU partial CB: no vector of neighbour " + std::to_string(g.neighbour));
        remote_off.emplace(g.neighbour, xpoints);
        xpoints += it->second.first;
    }
    const size_t x_bytes = static_cast<size_t>(xpoints * ncol) * D, pw_bytes = static_cast<size_t>(L.pw_units * ncol) * D;
    const size_t rn_bytes = static_cast<size_t>(L.rn_units * ncol) * D, out_bytes = static_cast<size_t>(L.out_units * ncol) * D;
    char* h = static_cast<char*>(st.host_in.reserve(std::max<size_t>(x_bytes + pw_bytes, 1)));
#pragma omp parallel for schedule(dynamic, 16)
    for (int64_t b = 0; b < static_cast<int64_t>(x_local.size()); ++b) {
        std::memcpy(h + static_cast<size_t>(xoff[static_cast<size_t>(b)] * ncol) * D, x_local[static_cast<size_t>(b)],
                    static_cast<size_t>((xoff[static_cast<size_t>(b) + 1] - xoff[static_cast<size_t>(b)]) * ncol) * D);
    }
    for (const auto& [m, off] : remote_off) {
        const auto& xr = x_remote.at(m);
        std::memcpy(h + static_cast<size_t>(off * ncol) * D, xr.second, static_cast<size_t>(xr.first * ncol) * D);
    }
    std::memcpy(h + x_bytes, pw, pw_bytes);
    char* d_x = heap.alloc(std::max<size_t>(align_up(x_bytes) + align_up(pw_bytes) + align_up(rn_bytes) + align_up(out_bytes), 1));
    char* d_pw = d_x + align_up(x_bytes);
    char* d_rn = d_pw + align_up(pw_bytes);
    char* d_out = d_rn + align_up(rn_bytes);
    check_cuda(cudaMemcpyAsync(d_x, h, x_bytes, cudaMemcpyHostToDevice, stream), "partial CB upload");
    check_cuda(cudaMemcpyAsync(d_pw, h + x_bytes, pw_bytes, cudaMemcpyHostToDevice, stream), "partial CB upload");
    check_cuda(cudaMemsetAsync(d_out, 0, out_bytes, stream), "partial CB outputs");
    // ---- the metadata: the gathers, then per segment its three GEMM batches
    MetaBuilder& meta = st.meta;
    meta.clear();
    std::vector<GatherItemT<S>> gathers;
    gathers.reserve(L.gathers.size());
    for (const auto& g : L.gathers) {
        const PartialDeviceBox& B = L.boxes[static_cast<size_t>(g.q)];
        GatherItemT<S> it{};
        it.out = reinterpret_cast<S*>(d_rn) + B.rn_off * ncol + g.col0;
        it.ld = B.total;
        it.m = g.m;
        it.n = static_cast<int>(ncol);
        int64_t off = 0, n = 0;
        if (g.remote) {
            off = remote_off.at(g.neighbour);
            n = x_remote.at(g.neighbour).first;
        } else {
            off = xoff[static_cast<size_t>(g.neighbour)];
            n = xoff[static_cast<size_t>(g.neighbour) + 1] - off;
        }
        it.src = reinterpret_cast<const S*>(d_x) + off * ncol;
        it.rs = 1;
        it.cs = n;
        it.rows.ptr = g.d_slots;
        gathers.push_back(it);
    }
    const size_t gathers_off = meta.append(gathers);
    auto rec = [&](const PartialDeviceRecord& d) -> const DataType* {
        return reinterpret_cast<const DataType*>(d.compressed() ? st.scratch + d.scratch_off : d.dev);
    };
    std::vector<VBatch<DataType>> batches(3 * L.segments.size());
    const int nc = static_cast<int>(ncol);
    for (size_t s = 0; s < L.segments.size(); ++s) {
        const PartialDeviceSegment& seg = L.segments[s];
        for (int q = seg.first; q < seg.first + seg.count; ++q) {
            const PartialDeviceBox& B = L.boxes[static_cast<size_t>(q)];
            const DataType* p = reinterpret_cast<const DataType*>(d_pw) + B.pw_off * ncol;
            const DataType* w = p + static_cast<int64_t>(B.r) * ncol;
            const DataType* rn = reinterpret_cast<const DataType*>(d_rn) + B.rn_off * ncol;
            DataType* gR = reinterpret_cast<DataType*>(d_out) + B.out_off * ncol;
            DataType* gS = gR + static_cast<int64_t>(B.r) * ncol;
            auto add = [&](VBatch<DataType>& b, const DataType* a, const DataType* x, DataType* c, int m, int k) {
                if (m == 0 || k == 0 || nc == 0) return;
                b.entries.push_back({a, x, c, m, nc, k, m, k, m});
            };
            // gR = X_RR p, gS = X_SR p, t_q = X_NR[q] p; then gR += X_RS w; gR += X_RN rn
            add(batches[3 * s], rec(B.rr), p, gR, B.r, B.r);
            add(batches[3 * s], rec(B.sr), p, gS, B.kp, B.r);
            for (size_t t = 0; t < B.nr.size(); ++t) {
                add(batches[3 * s], rec(B.nr[t]), p, gR + B.t_off[t] * ncol, static_cast<int>(B.nr[t].rows), B.r);
            }
            add(batches[3 * s + 1], rec(B.rs), w, gR, B.r, B.kp);
            add(batches[3 * s + 2], rec(B.rn), rn, gR, B.r, B.total);
        }
    }
    for (auto& b : batches) b.stage(meta);
    char* md = meta.upload(st.meta_device, stream);
    launch_gather(reinterpret_cast<const GatherItemT<S>*>(md + gathers_off), static_cast<int>(gathers.size()), L.max_gather_m,
                  nc, md, stream);
    for (size_t s = 0; s < L.segments.size(); ++s) {
#ifdef H2_HAVE_NVCOMP
        const PartialDeviceSegment& seg = L.segments[s];
        if (seg.have_plan) {
            bitcompSetStream(seg.plan, stream);
            if (bitcompBatchUncompress(seg.plan, seg.d_srcs, seg.d_dsts) != BITCOMP_SUCCESS) {
                throw std::runtime_error("GPU partial CB: decompressing a segment's records failed");
            }
        }
#endif
        for (int k = 0; k < 3; ++k) {
            const auto& b = batches[3 * s + static_cast<size_t>(k)];
            if (b.count() > 0) b.gemm(md, MagmaNoTrans, MagmaNoTrans, DataType(1), DataType(1), queue);
        }
    }
    char* ho = static_cast<char*>(st.host_out.reserve(std::max<size_t>(out_bytes, 1)));
    check_cuda(cudaMemcpyAsync(ho, d_out, out_bytes, cudaMemcpyDeviceToHost, stream), "partial CB download");
    check_cuda(cudaStreamSynchronize(stream), "partial CB terms");
    heap.free(d_x);
    std::memcpy(out, ho, out_bytes);
}

// The root's term: g = CB_root v (v, g: root_n x ncol, on the root's owner).
template<typename DataType>
void partial_device_root(const DataType* v, int64_t ncol, DataType* g) {
    PartialDeviceStore& st = partial_device_store();
    const int64_t n = st.root_n;
    if (n == 0 || ncol == 0) return;
    const size_t bytes = static_cast<size_t>(n * ncol) * sizeof(DataType);
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    DeviceHeap& heap = DeviceHeap::instance();
    char* d = heap.alloc(2 * align_up(bytes));
    char* h = static_cast<char*>(st.host_in.reserve(bytes));
    std::memcpy(h, v, bytes);
    check_cuda(cudaMemcpyAsync(d, h, bytes, cudaMemcpyHostToDevice, stream), "partial CB root upload");
    check_cuda(cudaMemsetAsync(d + align_up(bytes), 0, bytes, stream), "partial CB root");
    MetaBuilder& meta = st.meta;
    meta.clear();
    VBatch<DataType> b;
    b.entries.push_back({reinterpret_cast<const DataType*>(st.root_cb), reinterpret_cast<const DataType*>(d),
                         reinterpret_cast<DataType*>(d + align_up(bytes)), static_cast<magma_int_t>(n),
                         static_cast<magma_int_t>(ncol), static_cast<magma_int_t>(n), static_cast<magma_int_t>(n),
                         static_cast<magma_int_t>(n), static_cast<magma_int_t>(n)});
    b.stage(meta);
    char* md = meta.upload(st.meta_device, stream);
    b.gemm(md, MagmaNoTrans, MagmaNoTrans, DataType(1), DataType(1), ctx.queue());
    char* ho = static_cast<char*>(st.host_out.reserve(bytes));
    check_cuda(cudaMemcpyAsync(ho, d + align_up(bytes), bytes, cudaMemcpyDeviceToHost, stream), "partial CB root download");
    check_cuda(cudaStreamSynchronize(stream), "partial CB root");
    heap.free(d);
    std::memcpy(g, ho, bytes);
}

// c_bpack_partial_schur_release: the device copies go with the host's.
template<typename CoordType, typename DataType>
void partial_device_release(ParallelTree<CoordType, DataType>* tree) {
    activate_operator(tree);
    PartialDeviceStore& st = partial_device_store();
    st.release();
    st.tree = tree;
    st.decided = true;  // (and stays off: the records are gone)
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
