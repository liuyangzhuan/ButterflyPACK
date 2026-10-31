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
    PartialDeviceRecord t22;              // T22 (kp x r, raw: the CB entries' chain steps)
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
    // the records the factorization kept here (their blocks, owned), and how
    // to give them host copies before an eviction (partial_records_to_host)
    std::vector<char*> kept_blocks;
    size_t kept_bytes = 0;
    std::function<void()> spill;

    // what the first CB product or entries set up (the kept records stay)
    void release_prepared() {
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
    void release() {
        release_prepared();
        DeviceHeap& heap = DeviceHeap::instance();
        for (char* b : kept_blocks) heap.free(b);
        kept_blocks.clear();
        kept_bytes = 0;
        spill = nullptr;
    }
};

// The factorization keeps a chunk of records on the device while all it
// keeps stays under this fraction of the device pool (the higher levels'
// factorization runs beside them; decide_solve_keep sees what they take)
constexpr double kPartialKeepFraction = 0.2;

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
            return it == s->end() ? 0 : it->second.bytes + it->second.kept_bytes;
        });
        c.spillers.push_back([s](const void* tree) {
            auto it = s->find(tree);
            if (it != s->end() && !it->second.kept_blocks.empty() && it->second.spill) it->second.spill();
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
    st.release_prepared();
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
                      shape(f.X_RN, B.r, B.total) && shape(f.T, B.kp, B.r) && !f.T.compressed();
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
                raw += align_up(static_cast<size_t>(h.rows * h.cols) * D);
                if (h.device != nullptr) {  // kept by the factorization: read in place
                    d.dev = h.device;
                    if (h.compressed()) {
                        d.chunk_off.assign(h.strides.size(), 0);
                        for (size_t c = 1; c < h.strides.size(); ++c) d.chunk_off[c] = d.chunk_off[c - 1] + h.strides[c - 1];
                    }
                    return;
                }
                record_bytes += align_up(h.compressed() ? h.packed.size() : static_cast<size_t>(h.rows * h.cols) * D);
                placed.push_back(Placed{&d, &h});
            };
            add(B.rr, f.X_RR);
            add(B.rs, f.X_RS);
            add(B.sr, f.X_SR);
            add(B.rn, f.X_RN);
            for (size_t q = 0; q < B.nr.size(); ++q) add(B.nr[q], f.X_NR[q]);
            // (T22, raw, is never decompressed: outside the box's scratch)
            B.t22.rows = f.T.rows;
            B.t22.cols = f.T.cols;
            if (f.T.rows * f.T.cols > 0) {
                record_bytes += align_up(static_cast<size_t>(f.T.rows * f.T.cols) * D);
                placed.push_back(Placed{&B.t22, &f.T});
            }
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
        // the records kept here go to the host (the host paths read them) and
        // free the device for the other operators
        if (!st.kept_blocks.empty() && st.spill) st.spill();
        st.release();
        st.tree = tree;
        st.decided = true;
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

// ---------------------------------------------------------------------------
// The CB entries on the device (unstructured/partial.hpp
// partial_schur_blocks_device): transient device memory, transfers through
// the store's pinned buffers, and the pair terms' executor.
// ---------------------------------------------------------------------------
class PartialDeviceWork {
public:
    PartialDeviceWork() : st_(partial_device_store()), heap_(DeviceHeap::instance()) {
        Context& ctx = Context::instance();
        ctx.activate();
        stream_ = ctx.stream();
        queue_ = ctx.queue();
        for (cudaEvent_t& e : staged_) check_cuda(cudaEventCreateWithFlags(&e, cudaEventDisableTiming), "cudaEventCreate");
    }
    ~PartialDeviceWork() {
        cudaStreamSynchronize(stream_);
        for (char* p : held_) heap_.free(p);
        for (cudaEvent_t e : staged_) cudaEventDestroy(e);
    }
    PartialDeviceWork(const PartialDeviceWork&) = delete;
    PartialDeviceWork& operator=(const PartialDeviceWork&) = delete;

    // transient device memory, freed by release() or at the end (the stream
    // keeps a later user of the same bytes behind the earlier kernels)
    char* alloc(size_t bytes) {
        char* p = heap_.alloc(std::max<size_t>(align_up(bytes), 256));
        held_.push_back(p);
        return p;
    }
    void release(char* p) {
        if (p == nullptr) return;
        auto it = std::find(held_.begin(), held_.end(), p);
        if (it != held_.end()) held_.erase(it);
        heap_.free(p);
    }
    void zero(void* p, size_t bytes) {
        if (bytes > 0) check_cuda(cudaMemsetAsync(p, 0, bytes, stream_), "partial CB zero");
    }
    // (host copies through the pinned buffers by all threads)
    static void copy(void* dst, const void* src, size_t bytes) {
        constexpr size_t kPiece = size_t{4} << 20;
        const int64_t n = static_cast<int64_t>((bytes + kPiece - 1) / kPiece);
#pragma omp parallel for schedule(static) if (n > 1)
        for (int64_t i = 0; i < n; ++i) {
            const size_t a = static_cast<size_t>(i) * kPiece;
            std::memcpy(static_cast<char*>(dst) + a, static_cast<const char*>(src) + a, std::min(kPiece, bytes - a));
        }
    }
    // Transfers through the backend's two pinned staging buffers, a piece of
    // kStage bytes at a time (they stay at their size: growing a pinned
    // buffer costs about a second per GB)
    void upload(void* dst, const void* src, size_t bytes) {
        for (size_t at = 0, k = 0; at < bytes; at += kStage, ++k) {
            const size_t n = std::min(kStage, bytes - at);
            const int turn = static_cast<int>(k & 1);
            if (k >= 2) check_cuda(cudaEventSynchronize(staged_[turn]), "partial CB upload");
            char* h = static_cast<char*>(pinned_pool().staging[turn].reserve(kStage));
            copy(h, static_cast<const char*>(src) + at, n);
            check_cuda(cudaMemcpyAsync(static_cast<char*>(dst) + at, h, n, cudaMemcpyHostToDevice, stream_), "partial CB upload");
            check_cuda(cudaEventRecord(staged_[turn], stream_), "cudaEventRecord");
        }
        check_cuda(cudaStreamSynchronize(stream_), "partial CB upload");
    }
    // (a piece comes over while the host copies the one before out of the other buffer)
    void download(void* dst, const void* src, size_t bytes) {
        const size_t pieces = (bytes + kStage - 1) / kStage;
        char* h[2] = {nullptr, nullptr};
        auto start = [&](size_t k) {
            const size_t at = k * kStage, n = std::min(kStage, bytes - at);
            h[k & 1] = static_cast<char*>(pinned_pool().staging[k & 1].reserve(kStage));
            check_cuda(cudaMemcpyAsync(h[k & 1], static_cast<const char*>(src) + at, n, cudaMemcpyDeviceToHost, stream_),
                       "partial CB download");
            check_cuda(cudaEventRecord(staged_[k & 1], stream_), "cudaEventRecord");
        };
        if (pieces > 0) start(0);
        for (size_t k = 0; k < pieces; ++k) {
            if (k + 1 < pieces) start(k + 1);
            check_cuda(cudaEventSynchronize(staged_[k & 1]), "partial CB download");
            const size_t at = k * kStage;
            copy(static_cast<char*>(dst) + at, h[k & 1], std::min(kStage, bytes - at));
        }
    }
    MetaBuilder& meta() {
        st_.meta.clear();
        return st_.meta;
    }
    char* upload_meta() { return st_.meta.upload(st_.meta_device, stream_); }
    void sync() { check_cuda(cudaStreamSynchronize(stream_), "partial CB"); }
    cudaStream_t stream() const { return stream_; }
    magma_queue_t queue() const { return queue_; }
    PartialDeviceStore& store() { return st_; }
    // what a group of terms may take
    size_t budget() const {
        const size_t free = heap_.capacity() > heap_.used() ? heap_.capacity() - heap_.used() : 0;
        return std::max<size_t>(size_t{64} << 20, std::min<size_t>(size_t{2} << 30, free / 4));
    }

private:
    static constexpr size_t kStage = size_t{64} << 20;
    cudaEvent_t staged_[2] = {nullptr, nullptr};
    PartialDeviceStore& st_;
    DeviceHeap& heap_;
    cudaStream_t stream_ = nullptr;
    magma_queue_t queue_ = nullptr;
    std::vector<char*> held_;
};

// The values of an exchange of the CB entries (partial_schur_blocks_device,
// whose blocks' index lists go through the host): spans[p] the device values
// for rank p, in the order of its lists, and rcount[p] the elements from rank
// p (known from the lists received).  The values are gathered into one buffer
// on the device and, when MPI is GPU-aware (BPACK_GPU_AWARE_MPI), sent and
// received in device memory - in the exchange arena when they fit (MPI
// registers it once; device memory outside it makes the transfers much
// slower) - else through host copies.  Returns the received values in work
// memory, rank p's from element roff[p] (null when none arrive).  Collective
// over comm.
template<typename DataType>
char* partial_values_exchange(PartialDeviceWork& work, MPI_Comm comm,
                              const std::vector<std::vector<std::pair<const DataType*, size_t>>>& spans,
                              const std::vector<int64_t>& rcount, int tag, std::vector<int64_t>& roff) {
    int nproc = 1, me = 0;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &me);
    constexpr size_t D = sizeof(DataType);
    constexpr int64_t kWords = static_cast<int64_t>(D / sizeof(double));
    const size_t P = static_cast<size_t>(nproc), self = static_cast<size_t>(me);
    std::vector<int64_t> soff(P + 1, 0);
    roff.assign(P + 1, 0);
    for (size_t p = 0; p < P; ++p) {
        int64_t n = 0;
        for (const auto& s : spans[p]) n += static_cast<int64_t>(s.second);
        soff[p + 1] = soff[p] + n;
        roff[p + 1] = roff[p] + rcount[p];
        if (n * kWords > INT32_MAX || rcount[p] * kWords > INT32_MAX) {
            throw std::runtime_error("partial CB entries: a message over MPI's counts");
        }
    }
    if (rcount[self] != soff[self + 1] - soff[self]) throw std::logic_error("partial CB entries: a rank's own values miscounted");
    const size_t sbytes = static_cast<size_t>(soff[P]) * D, rbytes = static_cast<size_t>(roff[P]) * D;
    const bool direct = device_exchange_enabled();
    DeviceHeap& arena = DeviceHeap::exchange_arena();
    auto buffer = [&](size_t bytes, bool& in_arena) -> char* {
        in_arena = false;
        if (bytes == 0) return nullptr;
        if (direct && arena.initialized()) {
            if (char* p = arena.try_alloc(bytes)) {
                in_arena = true;
                return p;
            }
        }
        return work.alloc(bytes);
    };
    bool s_arena = false, r_arena = false;
    char* sbuf = buffer(sbytes, s_arena);
    char* rbuf = buffer(rbytes, r_arena);
    cudaStream_t stream = work.stream();
    if (sbytes > 0) {
        MetaBuilder& meta = work.meta();
        std::vector<CopySpan> items;
        int64_t max_words = 0;
        for (size_t p = 0; p < P; ++p) {
            int64_t at = soff[p];
            for (const auto& [src, e] : spans[p]) {
                const int64_t words = static_cast<int64_t>(e * D / 4);
                if (words > 0) items.push_back({src, sbuf + static_cast<size_t>(at) * D, words});
                max_words = std::max(max_words, words);
                at += static_cast<int64_t>(e);
            }
        }
        const size_t off = meta.append(items);
        char* md = work.upload_meta();
        launch_copy_spans(reinterpret_cast<const CopySpan*>(md + off), static_cast<int>(items.size()), max_words, stream);
    }
    if (rcount[self] > 0) {
        check_cuda(cudaMemcpyAsync(rbuf + static_cast<size_t>(roff[self]) * D, sbuf + static_cast<size_t>(soff[self]) * D,
                                   static_cast<size_t>(rcount[self]) * D, cudaMemcpyDeviceToDevice, stream),
                   "partial CB exchange");
    }
    // the messages: from and into device memory, or host copies of it
    std::vector<char> hs, hr;
    char* s_msg = sbuf;
    char* r_msg = rbuf;
    if (!direct) {
        hs.resize(sbytes);
        hr.resize(rbytes);
        if (sbytes > 0) work.download(hs.data(), sbuf, sbytes);
        s_msg = hs.data();
        r_msg = hr.data();
    }
    work.sync();
    std::vector<MPI_Request> reqs;
    for (int p = 0; p < nproc; ++p) {
        const size_t q = static_cast<size_t>(p);
        if (p == me || rcount[q] == 0) continue;
        reqs.emplace_back();
        MPI_Irecv(r_msg + static_cast<size_t>(roff[q]) * D, static_cast<int>(rcount[q] * kWords), MPI_DOUBLE, p, tag, comm,
                  &reqs.back());
    }
    for (int p = 0; p < nproc; ++p) {
        const size_t q = static_cast<size_t>(p);
        const int64_t n = soff[q + 1] - soff[q];
        if (p == me || n == 0) continue;
        reqs.emplace_back();
        MPI_Isend(s_msg + static_cast<size_t>(soff[q]) * D, static_cast<int>(n * kWords), MPI_DOUBLE, p, tag, comm, &reqs.back());
    }
    if (!reqs.empty()) MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
    if (!direct) {
        for (size_t p = 0; p < P; ++p) {
            if (p == self || rcount[p] == 0) continue;
            work.upload(rbuf + static_cast<size_t>(roff[p]) * D, hr.data() + static_cast<size_t>(roff[p]) * D,
                        static_cast<size_t>(rcount[p]) * D);
        }
    }
    char* out = rbuf;
    if (r_arena) {
        // (out of the arena: the values stay for the level)
        out = work.alloc(rbytes);
        check_cuda(cudaMemcpyAsync(out, rbuf, rbytes, cudaMemcpyDeviceToDevice, stream), "partial CB exchange");
        work.sync();
        arena.free(rbuf);
    }
    if (s_arena) {
        arena.free(sbuf);
    } else {
        work.release(sbuf);
    }
    return out;
}

// An operand of a term: a block's values on the device (ld rows), all of
// them, or rows of them gathered (a block's list, or the term's own)
template<typename DataType>
struct PartialTermOperand {
    const DataType* dev = nullptr;
    int64_t ld = 0;
    bool gathered = false;
    const std::vector<int64_t>* rows_ext = nullptr;
    std::vector<int64_t> rows_own;
    const std::vector<int64_t>& rows() const { return rows_ext != nullptr ? *rows_ext : rows_own; }
};

// A term G (ni x nj) += A^T X B of a pair, X (k1 x k2) = rec(rows, cols),
// A (k1 x ni), B (k2 x nj), in the cheaper order.  kind: its round (a pair
// has one term of a kind at most; the rounds run in kind order, which is the
// host's order of a pair's terms).  rec null: X from `raw` (the root block,
// leading dimension ld).
template<typename DataType>
struct PartialTerm {
    int kind = 0;
    const PartialDeviceRecord* rec = nullptr;
    const DataType* raw = nullptr;
    int64_t ld = 0;
    // X's rows and columns in the record: a block's lists, or the term's own
    const std::vector<int64_t>* rows_ext = nullptr;
    const std::vector<int64_t>* cols_ext = nullptr;
    std::vector<int64_t> rows_own, cols_own;
    PartialTermOperand<DataType> a, b;
    int k1 = 0, k2 = 0;
    const std::vector<int64_t>& rows() const { return rows_ext != nullptr ? *rows_ext : rows_own; }
    const std::vector<int64_t>& cols() const { return cols_ext != nullptr ? *cols_ext : cols_own; }
};

template<typename DataType>
struct PartialPair {
    int ni = 0, nj = 0;
    std::vector<PartialTerm<DataType>> terms;
};

constexpr int kPartialTermKinds = 6;  // RR, RS, SR, RN, NR, the root

// The chunks of compressed records a term reads (record -> sorted chunk indices)
template<typename DataType>
void partial_term_chunks(const PartialTerm<DataType>& t, std::vector<int64_t>& chunks) {
    chunks.clear();
    if (t.rec == nullptr || !t.rec->compressed() || t.rows().empty() || t.cols().empty()) return;
    const int64_t E = static_cast<int64_t>(std::max<size_t>(1, kRecordChunkBytes / sizeof(DataType)));
    const auto mm = std::minmax_element(t.rows().begin(), t.rows().end());
    for (int64_t c : t.cols()) {
        const int64_t lo = (*mm.first + c * t.ld) / E, hi = (*mm.second + c * t.ld) / E;
        for (int64_t k = lo; k <= hi; ++k) chunks.push_back(k);
    }
    std::sort(chunks.begin(), chunks.end());
    chunks.erase(std::unique(chunks.begin(), chunks.end()), chunks.end());
}

// Device bytes a pair takes in a group (X, gathered operands, products, G)
template<typename DataType>
size_t partial_pair_bytes(const PartialPair<DataType>& p) {
    size_t e = static_cast<size_t>(p.ni) * static_cast<size_t>(p.nj);
    for (const auto& t : p.terms) {
        e += static_cast<size_t>(t.k1) * static_cast<size_t>(t.k2);
        if (t.a.gathered) e += static_cast<size_t>(t.k1) * static_cast<size_t>(p.ni);
        if (t.b.gathered) e += static_cast<size_t>(t.k2) * static_cast<size_t>(p.nj);
        e += std::max(static_cast<size_t>(t.k1) * static_cast<size_t>(p.nj), static_cast<size_t>(p.ni) * static_cast<size_t>(t.k2));
    }
    return e * sizeof(DataType) + 256 * (4 + 3 * p.terms.size());
}

// One side of a pair's targets: entry a of the list takes row (or column)
// src[a] of G (src empty: a itself) to row dst[a] of the output, weighted by
// w[a] (w empty: 1).
struct PartialTargetList {
    std::vector<int> src, dst;
    std::vector<double> w;
};

// Where the G of a group go when they stay on the device: pair q's G is
// added into out at rows[q] x cols[q] (PartialTargetList), entry (a, b) of
// its block of out at base[q] + a ors[q] + b ocs[q].  With add 1 the pairs
// of a group write distinct entries: at a level, the row blocks partition
// the requested rows' points and the column blocks the columns'.  A
// sketched request's blocks share its sketch rows: add 2, atomic adds.
template<typename DataType>
struct PartialPairTargets {
    DataType* out = nullptr;
    std::vector<int64_t> base, ors, ocs;
    std::vector<const PartialTargetList*> rows, cols;
    int add = 1;
};

// The G (ni x nj, column major, end to end in pair order) of a group of
// pairs, into g (host); with targets, added into their blocks of the device
// output instead (g unused).
template<typename DataType>
void partial_device_pairs(PartialDeviceWork& work, const std::vector<const PartialPair<DataType>*>& pairs, DataType* g,
                          const PartialPairTargets<DataType>* targets = nullptr) {
    using S = typename DeviceScalar<DataType>::type;
    const size_t D = sizeof(DataType);
    if (pairs.empty()) return;
    // ---- the layout: G of every pair, then per term its X, gathered operands and product
    struct Place { size_t x = 0, a = 0, b = 0, tmp = 0; bool xb = true; };
    std::vector<size_t> g_off(pairs.size() + 1, 0);
    std::vector<std::vector<Place>> place(pairs.size());
    size_t at = 0;
    for (size_t q = 0; q < pairs.size(); ++q) {
        g_off[q] = at;
        at += align_up(static_cast<size_t>(pairs[q]->ni) * static_cast<size_t>(pairs[q]->nj) * D);
    }
    const size_t g_bytes = at;
    for (size_t q = 0; q < pairs.size(); ++q) {
        const auto& P = *pairs[q];
        place[q].resize(P.terms.size());
        for (size_t t = 0; t < P.terms.size(); ++t) {
            const auto& T = P.terms[t];
            Place& pl = place[q][t];
            const double k1 = T.k1, k2 = T.k2, ni = P.ni, nj = P.nj;
            pl.xb = k1 * nj * (k2 + ni) <= ni * k2 * (k1 + nj);
            pl.x = at;
            at += align_up(static_cast<size_t>(T.k1) * static_cast<size_t>(T.k2) * D);
            if (T.a.gathered) {
                pl.a = at;
                at += align_up(static_cast<size_t>(T.k1) * static_cast<size_t>(P.ni) * D);
            }
            if (T.b.gathered) {
                pl.b = at;
                at += align_up(static_cast<size_t>(T.k2) * static_cast<size_t>(P.nj) * D);
            }
            pl.tmp = at;
            at += align_up((pl.xb ? static_cast<size_t>(T.k1) * static_cast<size_t>(P.nj)
                                  : static_cast<size_t>(P.ni) * static_cast<size_t>(T.k2)) * D);
        }
    }
    char* d = work.alloc(at);
    work.zero(d, g_bytes);
    // ---- the compressed records' chunks the group reads, into the store's scratch
    PartialDeviceStore& st = work.store();
    const int64_t E = static_cast<int64_t>(std::max<size_t>(1, kRecordChunkBytes / D));
    std::unordered_map<const PartialDeviceRecord*, std::vector<int>> slots;  // record -> chunk -> slot
    std::vector<const void*> srcs;
    std::vector<void*> dsts;
    std::vector<int64_t> chunks;
    const size_t pool_slots = st.scratch_bytes / static_cast<size_t>(E * static_cast<int64_t>(D));
    for (const auto* P : pairs) {
        for (const auto& T : P->terms) {
            partial_term_chunks(T, chunks);
            if (chunks.empty()) continue;
            auto& sl = slots[T.rec];
            if (sl.empty()) sl.assign(T.rec->chunk_off.size(), -1);
            for (int64_t k : chunks) {
                if (sl[static_cast<size_t>(k)] >= 0) continue;
                if (srcs.size() >= pool_slots) throw std::logic_error("GPU partial CB entries: a group's chunks exceed the scratch");
                sl[static_cast<size_t>(k)] = static_cast<int>(srcs.size());
                srcs.push_back(T.rec->dev + T.rec->chunk_off[static_cast<size_t>(k)]);
                dsts.push_back(st.scratch + srcs.size() * static_cast<size_t>(E) * D - static_cast<size_t>(E) * D);
            }
        }
    }
    // ---- the metadata: lists, slot tables, gathers, the decompression's pointers, the GEMM batches
    MetaBuilder& meta = work.meta();
    std::unordered_map<const void*, size_t> list_at;
    auto list = [&](const std::vector<int64_t>* v) -> IndexList {
        auto it = list_at.find(v);
        if (it == list_at.end()) {
            std::vector<int> ints(v->begin(), v->end());
            it = list_at.emplace(v, meta.append(ints)).first;
        }
        IndexList l;
        l.offset = static_cast<int64_t>(it->second);
        return l;
    };
    std::unordered_map<const PartialDeviceRecord*, size_t> slot_at;
    for (const auto& [rec, sl] : slots) slot_at.emplace(rec, meta.append(sl));
    std::vector<GatherItemT<S>> gathers;
    std::vector<ChunkGatherItemT<S>> chunked;
    int gm = 0, gn = 0, cm = 0, cn = 0;
    std::vector<VBatch<DataType>> batches(4 * kPartialTermKinds);  // per kind: X B, A^T X; A^T tmp, tmp B
    for (size_t q = 0; q < pairs.size(); ++q) {
        const auto& P = *pairs[q];
        DataType* G = reinterpret_cast<DataType*>(d + g_off[q]);
        for (size_t t = 0; t < P.terms.size(); ++t) {
            const auto& T = P.terms[t];
            const Place& pl = place[q][t];
            DataType* X = reinterpret_cast<DataType*>(d + pl.x);
            if (T.rec != nullptr && T.rec->compressed()) {
                ChunkGatherItemT<S> it{};
                it.out = reinterpret_cast<S*>(X);
                it.ld = T.k1;
                it.m = T.k1;
                it.n = T.k2;
                it.pool = reinterpret_cast<const S*>(st.scratch);
                it.slot.offset = static_cast<int64_t>(slot_at.at(T.rec));
                it.chunk = E;
                it.lda = T.ld;
                it.rows = list(&T.rows());
                it.cols = list(&T.cols());
                chunked.push_back(it);
                cm = std::max(cm, T.k1);
                cn = std::max(cn, T.k2);
            } else {
                GatherItemT<S> it{};
                it.out = reinterpret_cast<S*>(X);
                it.ld = T.k1;
                it.m = T.k1;
                it.n = T.k2;
                it.src = reinterpret_cast<const S*>(T.rec != nullptr ? reinterpret_cast<const DataType*>(T.rec->dev) : T.raw);
                it.rs = 1;
                it.cs = T.ld;
                it.rows = list(&T.rows());
                it.cols = list(&T.cols());
                gathers.push_back(it);
                gm = std::max(gm, T.k1);
                gn = std::max(gn, T.k2);
            }
            const DataType* A = T.a.dev;
            int64_t lda = T.a.ld;
            if (T.a.gathered) {  // A's rows gathered
                GatherItemT<S> it{};
                it.out = reinterpret_cast<S*>(d + pl.a);
                it.ld = T.k1;
                it.m = T.k1;
                it.n = P.ni;
                it.src = reinterpret_cast<const S*>(T.a.dev);
                it.rs = 1;
                it.cs = T.a.ld;
                it.rows = list(&T.a.rows());
                gathers.push_back(it);
                gm = std::max(gm, T.k1);
                gn = std::max(gn, P.ni);
                A = reinterpret_cast<const DataType*>(d + pl.a);
                lda = T.k1;
            }
            const DataType* B = T.b.dev;
            int64_t ldb = T.b.ld;
            if (T.b.gathered) {
                GatherItemT<S> it{};
                it.out = reinterpret_cast<S*>(d + pl.b);
                it.ld = T.k2;
                it.m = T.k2;
                it.n = P.nj;
                it.src = reinterpret_cast<const S*>(T.b.dev);
                it.rs = 1;
                it.cs = T.b.ld;
                it.rows = list(&T.b.rows());
                gathers.push_back(it);
                gm = std::max(gm, T.k2);
                gn = std::max(gn, P.nj);
                B = reinterpret_cast<const DataType*>(d + pl.b);
                ldb = T.k2;
            }
            DataType* tmp = reinterpret_cast<DataType*>(d + pl.tmp);
            const magma_int_t k1 = T.k1, k2 = T.k2, ni = P.ni, nj = P.nj;
            auto& b = batches;
            const size_t base = static_cast<size_t>(4 * T.kind);
            if (pl.xb) {  // tmp = X B; G += A^T tmp
                b[base].entries.push_back({X, B, tmp, k1, nj, k2, k1, static_cast<magma_int_t>(ldb), k1});
                b[base + 2].entries.push_back({A, tmp, G, ni, nj, k1, static_cast<magma_int_t>(lda), k1, ni});
            } else {      // tmp = A^T X; G += tmp B
                b[base + 1].entries.push_back({A, X, tmp, ni, k2, k1, static_cast<magma_int_t>(lda), k1, ni});
                b[base + 3].entries.push_back({tmp, B, G, ni, nj, k2, ni, static_cast<magma_int_t>(ldb), ni});
            }
        }
    }
    // (with targets: each G into its block of out, in tiles - a launch's grid
    // covers its largest item)
    constexpr int kScatterTile = 64;
    std::vector<BlockCopyItemT<S>> tiles;
    if (targets != nullptr) {
        // (a list's dst, src and w in meta, once a list; -1: absent)
        struct ListAt { int64_t dst = -1, src = -1, w = -1; int n = 0; };
        std::unordered_map<const PartialTargetList*, ListAt> list_at;
        auto lists = [&](const PartialTargetList* v) {
            auto it = list_at.find(v);
            if (it == list_at.end()) {
                ListAt a;
                a.n = static_cast<int>(v->dst.size());
                a.dst = static_cast<int64_t>(meta.append(v->dst));
                if (!v->src.empty()) a.src = static_cast<int64_t>(meta.append(v->src));
                if (!v->w.empty()) a.w = static_cast<int64_t>(meta.append(v->w));
                it = list_at.emplace(v, a).first;
            }
            return it->second;
        };
        auto at = [](int64_t list, int i, size_t bytes) { return list < 0 ? list : list + static_cast<int64_t>(i) * static_cast<int64_t>(bytes); };
        for (size_t q = 0; q < pairs.size(); ++q) {
            const int ni = pairs[q]->ni;
            const ListAt ro = lists(targets->rows[q]), co = lists(targets->cols[q]);
            const S* G = reinterpret_cast<const S*>(d + g_off[q]);
            for (int j0 = 0; j0 < co.n; j0 += kScatterTile) {
                for (int i0 = 0; i0 < ro.n; i0 += kScatterTile) {
                    BlockCopyItemT<S> it{};
                    it.out = reinterpret_cast<S*>(targets->out + targets->base[q]);
                    it.ors = targets->ors[q];
                    it.ocs = targets->ocs[q];
                    it.out_rows.offset = at(ro.dst, i0, sizeof(int));
                    it.out_cols.offset = at(co.dst, j0, sizeof(int));
                    it.src = G;
                    it.rs = 1;
                    it.cs = ni;
                    // (an identity source list: the tile's own rows)
                    if (ro.src >= 0) it.src_rows.offset = at(ro.src, i0, sizeof(int));
                    else it.src_rows.base = i0;
                    if (co.src >= 0) it.src_cols.offset = at(co.src, j0, sizeof(int));
                    else it.src_cols.base = j0;
                    it.row_w = at(ro.w, i0, sizeof(double));
                    it.col_w = at(co.w, j0, sizeof(double));
                    it.m = std::min(kScatterTile, ro.n - i0);
                    it.n = std::min(kScatterTile, co.n - j0);
                    it.add = targets->add;
                    tiles.push_back(it);
                }
            }
        }
    }
    const size_t gathers_off = meta.append(gathers), chunked_off = meta.append(chunked);
    const size_t srcs_off = meta.append(srcs), dsts_off = meta.append(dsts), tiles_off = meta.append(tiles);
    for (auto& b : batches) b.stage(meta);
    char* md = work.upload_meta();
    cudaStream_t stream = work.stream();
    launch_gather(reinterpret_cast<const GatherItemT<S>*>(md + gathers_off), static_cast<int>(gathers.size()), gm, gn, md, stream);
    if (!srcs.empty()) {
#ifdef H2_HAVE_NVCOMP
        const void** d_srcs = reinterpret_cast<const void**>(md + srcs_off);
        void** d_dsts = reinterpret_cast<void**>(md + dsts_off);
        bitcompHandle_t plan{};
        if (bitcompCreateBatchPlanFromCompressedData(&plan, d_srcs, srcs.size(), nullptr, stream) != BITCOMP_SUCCESS) {
            throw std::runtime_error("GPU partial CB entries: a decompression plan");
        }
        bitcompSetStream(plan, stream);
        const bitcompResult_t r = bitcompBatchUncompress(plan, d_srcs, d_dsts);
        bitcompDestroyPlan(plan);
        if (r != BITCOMP_SUCCESS) throw std::runtime_error("GPU partial CB entries: decompressing records failed");
#else
        throw std::logic_error("GPU partial CB entries: compressed records without a compressor");
#endif
    }
    launch_chunk_gather(reinterpret_cast<const ChunkGatherItemT<S>*>(md + chunked_off), static_cast<int>(chunked.size()), cm, cn,
                        md, stream);
    magma_queue_t queue = work.queue();
    for (int k = 0; k < kPartialTermKinds; ++k) {
        const size_t base = static_cast<size_t>(4 * k);
        if (batches[base].count() > 0) batches[base].gemm(md, MagmaNoTrans, MagmaNoTrans, DataType(1), DataType(0), queue);
        if (batches[base + 1].count() > 0) batches[base + 1].gemm(md, MagmaTrans, MagmaNoTrans, DataType(1), DataType(0), queue);
        if (batches[base + 2].count() > 0) batches[base + 2].gemm(md, MagmaTrans, MagmaNoTrans, DataType(1), DataType(1), queue);
        if (batches[base + 3].count() > 0) batches[base + 3].gemm(md, MagmaNoTrans, MagmaNoTrans, DataType(1), DataType(1), queue);
    }
    if (targets != nullptr) {
        launch_block_copy(reinterpret_cast<const BlockCopyItemT<S>*>(md + tiles_off), static_cast<int>(tiles.size()), kScatterTile,
                          kScatterTile, md, stream);
        work.release(d);  // (the stream keeps a later user of the bytes behind the copy)
        return;
    }
    // G of the pairs, end to end
    std::vector<DataType> all(g_bytes / D);
    work.download(all.data(), d, g_bytes);
    size_t o = 0;
    for (size_t q = 0; q < pairs.size(); ++q) {
        const size_t n = static_cast<size_t>(pairs[q]->ni) * static_cast<size_t>(pairs[q]->nj);
        std::copy_n(all.data() + g_off[q] / D, n, g + o);
        o += n;
    }
    work.release(d);
}

// Bytes from the device to the host through the backend's pinned staging
// buffers (the kept records' host copies)
inline void partial_download(void* dst, const void* src, size_t bytes) {
    PartialDeviceWork work;
    work.download(dst, src, bytes);
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
