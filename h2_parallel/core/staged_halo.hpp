// Staged halo exchange (H2_CA_staged_halo=2)
//
// The pre-elimination halo is split into four stages sent in priority order:
//   stage 0: ghost box shells (geometry, points, index lists — no matrices)
//   stage 1: Blue(+orange) ghosts' entry blocks (schur + near-field list)
//   stage 2: Purple ghosts' entry blocks
//   stage 3: Green ghosts' entry blocks
// All receives are posted up front; senders pack and send stages in order, so
// later stages transfer while earlier groups compute. A box's stage is its
// geometric group in its OWNER's brick, which sender and receiver derive
// identically from their color lists (the receiver eliminates every halo box
// with the same global schedule).
//
// In mode 2, block stages are waited per
//     group inside the elimination loop, so later groups' payloads transfer
//     while earlier groups compute; writes into unarrived ghosts divert to
//     the level's pending state and merge at arrival (staged_halo_merge_stage
//     in factorization.hpp). Without lazy Schur updates, all stages complete
//     before elimination, which is equivalent to the unstaged schedule.
//
// For a symmetric matrix, staged pair pruning is enabled automatically. One
// orientation is transferred and the reciprocal view is reconstructed after
// its dependency has arrived.
#pragma once

#include "serialization.hpp"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <limits>
#include <type_traits>

namespace fmm {

// ============================================================================

// ----- staged wire format ----------------------------------------------------

template<typename CoordType, typename DataType>
size_t get_staged_shell_size(const BoxData<CoordType, DataType>& box) {
    size_t size = 0;
    size += sizeof(int64_t);          // morton_index
    size += sizeof(int32_t);          // level
    size += sizeof(int32_t) * 3;      // grid_coords
    size += sizeof(CoordType) * 6;    // bounds
    size += sizeof(CoordType) * 3;    // center
    size += sizeof(CoordType);        // size
    size += sizeof(int64_t);          // parent_morton
    size += sizeof(int64_t) * 8;      // children_morton
    size += sizeof(int32_t);          // num_children
    size += sizeof(int64_t);          // num_points
    size += sizeof(bool);             // on_boundary
    auto vec_size = [](const auto& v) {
        using V = typename std::decay_t<decltype(v)>::value_type;
        return sizeof(size_t) + v.size() * sizeof(V);
    };
    size += vec_size(box.point_indices);
    size += vec_size(box.point_coords);
    size += vec_size(box.redundant_indices);
    size += vec_size(box.skeleton_indices);
    size += vec_size(box.one_hop);
    size += vec_size(box.use_full_set);
    size += vec_size(box.two_hop);
    return size;
}

template<typename CoordType, typename DataType>
char* serialize_staged_shell(const BoxData<CoordType, DataType>& box, char* ptr) {
    std::memcpy(ptr, &box.morton_index, sizeof(int64_t)); ptr += sizeof(int64_t);
    std::memcpy(ptr, &box.level, sizeof(int32_t)); ptr += sizeof(int32_t);
    std::memcpy(ptr, box.grid_coords, sizeof(int32_t) * 3); ptr += sizeof(int32_t) * 3;
    std::memcpy(ptr, box.bounds, sizeof(CoordType) * 6); ptr += sizeof(CoordType) * 6;
    std::memcpy(ptr, box.center, sizeof(CoordType) * 3); ptr += sizeof(CoordType) * 3;
    std::memcpy(ptr, &box.size, sizeof(CoordType)); ptr += sizeof(CoordType);
    std::memcpy(ptr, &box.parent_morton, sizeof(int64_t)); ptr += sizeof(int64_t);
    std::memcpy(ptr, box.children_morton, sizeof(int64_t) * 8); ptr += sizeof(int64_t) * 8;
    std::memcpy(ptr, &box.num_children, sizeof(int32_t)); ptr += sizeof(int32_t);
    std::memcpy(ptr, &box.num_points, sizeof(int64_t)); ptr += sizeof(int64_t);
    std::memcpy(ptr, &box.on_boundary, sizeof(bool)); ptr += sizeof(bool);

    auto put_vec = [&ptr](const auto& v) {
        using V = typename std::decay_t<decltype(v)>::value_type;
        size_t count = v.size();
        std::memcpy(ptr, &count, sizeof(size_t)); ptr += sizeof(size_t);
        if (count > 0) {
            std::memcpy(ptr, v.data(), count * sizeof(V));
            ptr += count * sizeof(V);
        }
    };
    put_vec(box.point_indices);
    put_vec(box.point_coords);
    put_vec(box.redundant_indices);
    put_vec(box.skeleton_indices);
    put_vec(box.one_hop);
    put_vec(box.use_full_set);
    put_vec(box.two_hop);
    return ptr;
}

template<typename CoordType, typename DataType>
const char* deserialize_staged_shell(BoxData<CoordType, DataType>& box, const char* ptr) {
    std::memcpy(&box.morton_index, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    std::memcpy(&box.level, ptr, sizeof(int32_t)); ptr += sizeof(int32_t);
    std::memcpy(box.grid_coords, ptr, sizeof(int32_t) * 3); ptr += sizeof(int32_t) * 3;
    std::memcpy(box.bounds, ptr, sizeof(CoordType) * 6); ptr += sizeof(CoordType) * 6;
    std::memcpy(box.center, ptr, sizeof(CoordType) * 3); ptr += sizeof(CoordType) * 3;
    std::memcpy(&box.size, ptr, sizeof(CoordType)); ptr += sizeof(CoordType);
    std::memcpy(&box.parent_morton, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    std::memcpy(box.children_morton, ptr, sizeof(int64_t) * 8); ptr += sizeof(int64_t) * 8;
    std::memcpy(&box.num_children, ptr, sizeof(int32_t)); ptr += sizeof(int32_t);
    std::memcpy(&box.num_points, ptr, sizeof(int64_t)); ptr += sizeof(int64_t);
    std::memcpy(&box.on_boundary, ptr, sizeof(bool)); ptr += sizeof(bool);

    auto get_vec = [&ptr](auto& v) {
        using V = typename std::decay_t<decltype(v)>::value_type;
        size_t count = 0;
        std::memcpy(&count, ptr, sizeof(size_t)); ptr += sizeof(size_t);
        v.resize(count);
        if (count > 0) {
            std::memcpy(v.data(), ptr, count * sizeof(V));
            ptr += count * sizeof(V);
        }
    };
    get_vec(box.point_indices);
    get_vec(box.point_coords);
    get_vec(box.redundant_indices);
    get_vec(box.skeleton_indices);
    get_vec(box.one_hop);
    get_vec(box.use_full_set);
    get_vec(box.two_hop);

    box.num_near_field_interactions = 0;
    box.num_far_field_interactions = 0;
    return ptr;
}

template<typename CoordType, typename DataType>
size_t get_staged_blocks_size(const BoxData<CoordType, DataType>& box,
                              const std::vector<uint8_t>* omit = nullptr) {
    size_t size = get_serialized_size(box.schur_complement);
    size += sizeof(size_t);  // block count
    const auto& blocks = box.near_field_modified_interactions;
    for (size_t b = 0; b < blocks.size(); ++b) {
        size += get_serialized_size(blocks[b], omit != nullptr && (*omit)[b] != 0);
    }
    return size;
}

template<typename CoordType, typename DataType>
char* serialize_staged_blocks(const BoxData<CoordType, DataType>& box, char* ptr,
                              const std::vector<uint8_t>* omit = nullptr) {
    ptr = serialize(box.schur_complement, ptr);
    const auto& blocks = box.near_field_modified_interactions;
    size_t count = blocks.size();
    std::memcpy(ptr, &count, sizeof(size_t)); ptr += sizeof(size_t);
    for (size_t b = 0; b < blocks.size(); ++b) {
        ptr = serialize(blocks[b], ptr, omit != nullptr && (*omit)[b] != 0);
    }
    return ptr;
}

template<typename CoordType, typename DataType>
const char* deserialize_staged_blocks_into(BoxData<CoordType, DataType>& box, const char* ptr) {
    ptr = deserialize(box.schur_complement, ptr);
    size_t count = 0;
    std::memcpy(&count, ptr, sizeof(size_t)); ptr += sizeof(size_t);
    box.near_field_modified_interactions.clear();
    box.near_field_modified_interactions.resize(count);
    box.near_field_interaction_map.clear();
    for (size_t b = 0; b < count; ++b) {
        ptr = deserialize(box.near_field_modified_interactions[b], ptr);
        box.near_field_interaction_map
            [box.near_field_modified_interactions[b].neighbor_morton] =
            static_cast<int64_t>(b);
    }
    box.num_near_field_interactions = static_cast<int64_t>(count);
    return ptr;
}

// ----- staged state and protocol --------------------------------------------

constexpr int STAGED_HALO_STAGES = 4;

// Byte wrapper whose default constructor intentionally leaves the value
// uninitialized, so vector::resize on multi-GB staged payload buffers does
// not serially memset memory the (parallel) pack overwrites in full anyway.
struct StagedByte {
    char b;
    StagedByte() noexcept {}
};
static_assert(sizeof(StagedByte) == 1, "StagedByte must be byte-sized");

template<typename CoordType, typename DataType>
struct StagedHaloState {
    bool active = false;
    std::vector<int> peers;
    std::vector<std::vector<int64_t>> send_mortons[STAGED_HALO_STAGES];
    std::vector<std::vector<int64_t>> recv_mortons[STAGED_HALO_STAGES];
    std::vector<std::vector<StagedByte>> send_buffers[STAGED_HALO_STAGES];
    std::vector<std::vector<StagedByte>> recv_buffers[STAGED_HALO_STAGES];
    std::vector<size_t> send_sizes[STAGED_HALO_STAGES];
    std::vector<size_t> recv_sizes[STAGED_HALO_STAGES];
    std::vector<MPI_Request> recv_requests[STAGED_HALO_STAGES];
    std::vector<MPI_Request> send_requests[STAGED_HALO_STAGES];
    bool stage_arrived[STAGED_HALO_STAGES] = {false, false, false, false};

    // Phase timing (ms), printed via staged_halo_print_timing.
    double t_request_ms = 0.0;
    double t_sizes_ms = 0.0;
    double t_pack_ms = 0.0;
    double t_assist_ms = 0.0;
    double t_wait_ms[STAGED_HALO_STAGES] = {0.0, 0.0, 0.0, 0.0};
    double t_deser_ms[STAGED_HALO_STAGES] = {0.0, 0.0, 0.0, 0.0};
    double t_merge_ms[STAGED_HALO_STAGES] = {0.0, 0.0, 0.0, 0.0};
};

inline void staged_halo_print_timing_row(
    const double v[STAGED_HALO_STAGES], char* buf, size_t buflen) {
    snprintf(buf, buflen, "%.0f/%.0f/%.0f/%.0f", v[0], v[1], v[2], v[3]);
}

template<typename CoordType, typename DataType>
void staged_halo_print_timing(
    const StagedHaloState<CoordType, DataType>& state, int level) {
    char waits[64], desers[64], merges[64];
    staged_halo_print_timing_row(state.t_wait_ms, waits, sizeof(waits));
    staged_halo_print_timing_row(state.t_deser_ms, desers, sizeof(desers));
    staged_halo_print_timing_row(state.t_merge_ms, merges, sizeof(merges));
    printf("  [staged timing] level %d: request %.0f ms, sizes %.0f ms, "
           "pack+send %.0f ms, assisting %.0f ms, waits %s ms, "
           "deserialize %s ms, merge %s ms\n",
           level, state.t_request_ms, state.t_sizes_ms, state.t_pack_ms,
           state.t_assist_ms, waits, desers, merges);
    fflush(stdout);
}

/**
 * @brief Geometric stage of each boundary box (1 = blue/orange, 2 = purple,
 * 3 = green), derived from the level's color lists. Sender and receiver both
 * hold a box's group in their lists (the receiver eliminates every halo box),
 * so the classification agrees on both sides.
 */
template<typename CoordType, typename DataType>
void build_staged_stage_map(
    const TreeLevel<CoordType, DataType>& lvl,
    std::unordered_map<int64_t, int>& out) {
    out.clear();
    for (int64_t m : lvl.blue) out[m] = 1;
    for (int64_t m : lvl.orange) out[m] = 1;
    for (int64_t m : lvl.purple) out[m] = 2;
    for (int64_t m : lvl.green) out[m] = 3;
}

inline int staged_stage_of(
    const std::unordered_map<int64_t, int>& stage_map, int64_t morton) {
    auto it = stage_map.find(morton);
    return it == stage_map.end() ? 1 : it->second;  // unknown => earliest (safe)
}

/// Symmetric-pair pruning is enabled automatically for symmetric CA levels.
inline bool staged_halo_prune_enabled() {
    return ca_staged_prune_enabled();
}

/**
 * @brief Should this near block ship metadata only, letting the receiver
 * rebuild A_NS by transposing the pair's other view? True when the receiver
 * provably holds that other view no later than this box's own arrival: the
 * partner is the receiver's local box, or it is shipped to the same receiver
 * in an earlier stage (on stage ties the lower-morton box carries the pair).
 */
template<typename CoordType, typename DataType>
bool staged_should_omit_block(
    const BoxData<CoordType, DataType>& box,
    const ModifiedBlock<DataType>& block,
    int peer_rank,
    const std::unordered_set<int64_t>& peer_requests,
    const std::unordered_map<int64_t, int>& stage_map,
    int msg_stage,
    int partner_owner_rank) {
    if (block.neighbor_morton < 0) {
        return false;
    }
    // Only pure-symmetric storage is reconstructable by transpose.
    if (block.A_SN.is_allocated() && !block.A_SN.data.empty()) {
        return false;
    }
    if (!block.a_ns_is_allocated() || block.a_ns_data_empty()) {
        return false;
    }
    if (partner_owner_rank == peer_rank) {
        return true;
    }
    if (peer_requests.find(block.neighbor_morton) == peer_requests.end()) {
        return false;
    }
    const int partner_stage = staged_stage_of(stage_map, block.neighbor_morton);
    if (partner_stage != msg_stage) {
        return partner_stage < msg_stage;
    }
    return box.morton_index > block.neighbor_morton;
}

/**
 * @brief Initiate the staged halo exchange: request lists, per-stage sizes,
 * all receives posted, sends packed and issued in stage-priority order, then
 * the assisting point-data exchange (synchronous, small) while the staged
 * payloads are in flight. Completion happens in staged_halo_wait_stage.
 */
/// Request-list encoding: a ghost whose entry BLOCKS are wanted carries this
/// bit; every requested ghost gets its shell (stage 0).
constexpr int64_t STAGED_REQUEST_BLOCKS_BIT = int64_t(1) << 62;

template<typename CoordType, typename DataType>
void initiate_staged_halo_gather(
    ParallelTree<CoordType, DataType>* tree,
    int level,
    StagedHaloState<CoordType, DataType>& state,
    MPI_Comm level_comm,
    const std::unordered_set<int64_t>* block_ghosts = nullptr) {
    // block_ghosts (component-owner schedule): request entry blocks only for
    // these ghosts — the remote boxes of components this rank eliminates —
    // and shells for every ghost. nullptr: blocks for all ghosts (replicated).

    auto& lvl = tree->levels[level];
    if (!lvl.is_process_active) {
        return;
    }

    int rank = 0;
    MPI_Comm_rank(level_comm, &rank);
    int comm_size = 0;
    MPI_Comm_size(level_comm, &comm_size);
    std::vector<int> global_rank_by_comm(static_cast<size_t>(comm_size));
    MPI_Allgather(
        &tree->mpi_rank, 1, MPI_INT,
        global_rank_by_comm.data(), 1, MPI_INT,
        level_comm);
    std::unordered_map<int, int> comm_rank_by_global;
    for (int comm_rank = 0; comm_rank < comm_size; ++comm_rank) {
        comm_rank_by_global.emplace(
            global_rank_by_comm[static_cast<size_t>(comm_rank)], comm_rank);
    }
    const uint32_t grid_size = 1u << level;

    auto owner_rank_for_morton = [&](int64_t morton) -> int {
        std::vector<uint64_t> single_box = {static_cast<uint64_t>(morton)};
        std::vector<uint32_t> morton_ids;
        if (tree->dimension == 2) {
            morton_ids = morton::assign_to_processes_2d(
                single_box, lvl.num_active_processes, grid_size);
        } else {
            morton_ids = morton::assign_to_processes_3d(
                single_box, lvl.num_active_processes, grid_size);
        }
        const int owner_global =
            lvl.morton_to_rank.at(static_cast<int>(morton_ids[0]));
        auto it = comm_rank_by_global.find(owner_global);
        if (it == comm_rank_by_global.end()) {
            throw std::runtime_error(
                "initiate_staged_halo_gather: owner is absent from level communicator");
        }
        return it->second;
    };

    auto abort_on_mpi_error = [&](int err, const char* op, int peer, size_t bytes) {
        if (err == MPI_SUCCESS) return;
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        std::fprintf(stderr,
                     "[staged-halo-mpi-error] rank=%d level=%d op=%s peer=%d bytes=%zu err=%.*s\n",
                     rank, level, op, peer, bytes, errlen, errbuf);
        std::fflush(stderr);
        MPI_Abort(level_comm, err);
    };

    using staged_clock = std::chrono::high_resolution_clock;
    auto phase_start = staged_clock::now();
    auto lap = [&phase_start](double& acc) {
        auto now = staged_clock::now();
        acc += std::chrono::duration<double, std::milli>(now - phase_start).count();
        phase_start = now;
    };

    // ---- Request lists (who owns which of my ghosts; who wants mine) ----
    std::unordered_map<int, std::vector<int64_t>> requests_to_send;
    for (int64_t morton : lvl.ghost_id) {
        const int owner = owner_rank_for_morton(morton);
        if (owner != rank) {
            requests_to_send[owner].push_back(morton);
        }
    }
    for (auto& [peer, mortons] : requests_to_send) {
        std::sort(mortons.begin(), mortons.end());
        for (int64_t& m : mortons) {
            if (block_ghosts == nullptr || block_ghosts->count(m) != 0) m |= STAGED_REQUEST_BLOCKS_BIT;
        }
    }
    auto decode_request = [](int64_t v, bool& wants_blocks) {
        wants_blocks = (v & STAGED_REQUEST_BLOCKS_BIT) != 0;
        return v & ~STAGED_REQUEST_BLOCKS_BIT;
    };

    state.peers.clear();
    for (const auto& [peer, _] : requests_to_send) {
        state.peers.push_back(peer);
    }
    std::sort(state.peers.begin(), state.peers.end());
    const size_t num_peers = state.peers.size();

    std::vector<MPI_Request> reqs;
    std::vector<int> req_send_counts(num_peers, 0);
    std::vector<int> req_recv_counts(num_peers, 0);
    for (size_t i = 0; i < num_peers; ++i) {
        MPI_Request r;
        int ierr = MPI_Irecv(&req_recv_counts[i], 1, MPI_INT,
                             state.peers[i], 140, level_comm, &r);
        abort_on_mpi_error(ierr, "Irecv(staged-counts)", state.peers[i], sizeof(int));
        reqs.push_back(r);
    }
    for (size_t i = 0; i < num_peers; ++i) {
        req_send_counts[i] = static_cast<int>(requests_to_send[state.peers[i]].size());
        MPI_Request r;
        int ierr = MPI_Isend(&req_send_counts[i], 1, MPI_INT,
                             state.peers[i], 140, level_comm, &r);
        abort_on_mpi_error(ierr, "Isend(staged-counts)", state.peers[i], sizeof(int));
        reqs.push_back(r);
    }
    if (!reqs.empty()) {
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        reqs.clear();
    }

    std::vector<std::vector<int64_t>> requests_received(num_peers);
    for (size_t i = 0; i < num_peers; ++i) {
        if (req_recv_counts[i] > 0) {
            requests_received[i].resize(static_cast<size_t>(req_recv_counts[i]));
            MPI_Request r;
            int ierr = MPI_Irecv(requests_received[i].data(), req_recv_counts[i],
                                 MPI_INT64_T, state.peers[i], 141, level_comm, &r);
            abort_on_mpi_error(ierr, "Irecv(staged-mortons)", state.peers[i],
                               static_cast<size_t>(req_recv_counts[i]) * sizeof(int64_t));
            reqs.push_back(r);
        }
    }
    for (size_t i = 0; i < num_peers; ++i) {
        if (req_send_counts[i] > 0) {
            MPI_Request r;
            int ierr = MPI_Isend(requests_to_send[state.peers[i]].data(),
                                 req_send_counts[i], MPI_INT64_T,
                                 state.peers[i], 141, level_comm, &r);
            abort_on_mpi_error(ierr, "Isend(staged-mortons)", state.peers[i],
                               static_cast<size_t>(req_send_counts[i]) * sizeof(int64_t));
            reqs.push_back(r);
        }
    }
    if (!reqs.empty()) {
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        reqs.clear();
    }
    lap(state.t_request_ms);

    // ---- Partition by stage on both sides ----
    std::unordered_map<int64_t, int> stage_map;
    build_staged_stage_map(lvl, stage_map);

    for (int s = 0; s < STAGED_HALO_STAGES; ++s) {
        state.send_mortons[s].assign(num_peers, {});
        state.recv_mortons[s].assign(num_peers, {});
        state.send_buffers[s].assign(num_peers, {});
        state.recv_buffers[s].assign(num_peers, {});
        state.send_sizes[s].assign(num_peers, 0);
        state.recv_sizes[s].assign(num_peers, 0);
        state.recv_requests[s].clear();
        state.stage_arrived[s] = false;
    }
    for (int s = 0; s < STAGED_HALO_STAGES; ++s) {
        state.send_requests[s].clear();
    }

    // Decode the request lists: every ghost gets a shell (stage 0); only the
    // ghosts flagged by the requester get their blocks (stages 1..3).
    std::vector<std::unordered_set<int64_t>> peer_block_requests(num_peers);
    for (size_t i = 0; i < num_peers; ++i) {
        for (int64_t& v : requests_received[i]) {
            bool wants_blocks = false;
            const int64_t morton = decode_request(v, wants_blocks);
            v = morton;
            state.send_mortons[0][i].push_back(morton);
            if (wants_blocks) {
                state.send_mortons[staged_stage_of(stage_map, morton)][i].push_back(morton);
                peer_block_requests[i].insert(morton);
            }
        }
        for (int64_t& v : requests_to_send[state.peers[i]]) {
            bool wants_blocks = false;
            const int64_t morton = decode_request(v, wants_blocks);
            v = morton;
            state.recv_mortons[0][i].push_back(morton);
            if (wants_blocks) state.recv_mortons[staged_stage_of(stage_map, morton)][i].push_back(morton);
        }
    }

    // ---- Symmetric-pair pruning decisions, one per (peer, box, block) ----
    const bool prune_on = staged_halo_prune_enabled();
    std::vector<std::unordered_set<int64_t>> peer_request_sets(num_peers);
    std::unordered_map<int64_t, int> owner_rank_memo;
    auto owner_rank_cached = [&](int64_t m) {
        auto mit = owner_rank_memo.find(m);
        if (mit != owner_rank_memo.end()) {
            return mit->second;
        }
        const int owner = owner_rank_for_morton(m);
        owner_rank_memo.emplace(m, owner);
        return owner;
    };
    if (prune_on) {
        // A pair is omitted only when the receiver gets the partner's BLOCKS
        // (or holds the partner locally): shell-only ghosts cannot rebuild it.
        for (size_t i = 0; i < num_peers; ++i) {
            peer_request_sets[i] = peer_block_requests[i];
        }
    }
    std::vector<std::unordered_map<int64_t, std::vector<uint8_t>>>
        send_omit(num_peers);
    auto omit_for = [&](size_t peer_idx,
                        int64_t morton) -> const std::vector<uint8_t>* {
        if (!prune_on) {
            return nullptr;
        }
        auto oit = send_omit[peer_idx].find(morton);
        return oit == send_omit[peer_idx].end() ? nullptr : &oit->second;
    };

    // ---- Per-stage sizes, exchanged as one size_t[4] per peer ----
    std::vector<size_t> send_size_msgs(num_peers * STAGED_HALO_STAGES, 0);
    std::vector<size_t> recv_size_msgs(num_peers * STAGED_HALO_STAGES, 0);
    for (size_t i = 0; i < num_peers; ++i) {
        for (int s = 0; s < STAGED_HALO_STAGES; ++s) {
            size_t total = 0;
            for (int64_t morton : state.send_mortons[s][i]) {
                const int64_t local_idx = morton - lvl.local_morton_start;
                if (local_idx < 0 || local_idx >= lvl.num_boxes_local) {
                    throw std::runtime_error(
                        "initiate_staged_halo_gather: requested box " +
                        std::to_string(morton) + " not local");
                }
                const auto& box = lvl.local_boxes[local_idx];
                if (s == 0) {
                    total += sizeof(size_t) + get_staged_shell_size(box);
                    continue;
                }
                if (prune_on) {
                    auto& flags = send_omit[i][morton];
                    flags.assign(box.near_field_modified_interactions.size(), 0);
                    for (size_t b = 0; b < flags.size(); ++b) {
                        const auto& blk = box.near_field_modified_interactions[b];
                        flags[b] =
                            staged_should_omit_block(
                                box, blk, state.peers[i], peer_request_sets[i],
                                stage_map, s,
                                blk.neighbor_morton >= 0
                                    ? owner_rank_cached(blk.neighbor_morton)
                                    : -1)
                                ? 1
                                : 0;
                    }
                }
                total += sizeof(size_t) +
                         get_staged_blocks_size(box, omit_for(i, morton));
            }
            state.send_sizes[s][i] = total;
            send_size_msgs[i * STAGED_HALO_STAGES + s] = total;
        }
    }
    for (size_t i = 0; i < num_peers; ++i) {
        MPI_Request r;
        int ierr = MPI_Irecv(&recv_size_msgs[i * STAGED_HALO_STAGES],
                             STAGED_HALO_STAGES, MPI_UINT64_T,
                             state.peers[i], 142, level_comm, &r);
        abort_on_mpi_error(ierr, "Irecv(staged-sizes)", state.peers[i],
                           STAGED_HALO_STAGES * sizeof(size_t));
        reqs.push_back(r);
    }
    for (size_t i = 0; i < num_peers; ++i) {
        MPI_Request r;
        int ierr = MPI_Isend(&send_size_msgs[i * STAGED_HALO_STAGES],
                             STAGED_HALO_STAGES, MPI_UINT64_T,
                             state.peers[i], 142, level_comm, &r);
        abort_on_mpi_error(ierr, "Isend(staged-sizes)", state.peers[i],
                           STAGED_HALO_STAGES * sizeof(size_t));
        reqs.push_back(r);
    }
    if (!reqs.empty()) {
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        reqs.clear();
    }
    for (size_t i = 0; i < num_peers; ++i) {
        for (int s = 0; s < STAGED_HALO_STAGES; ++s) {
            state.recv_sizes[s][i] = recv_size_msgs[i * STAGED_HALO_STAGES + s];
        }
    }
    lap(state.t_sizes_ms);

    // ---- Post every stage's receives up front ----
    for (int s = 0; s < STAGED_HALO_STAGES; ++s) {
        for (size_t i = 0; i < num_peers; ++i) {
            if (state.recv_sizes[s][i] > 0) {
                state.recv_buffers[s][i].resize(state.recv_sizes[s][i]);
                int ierr = MPI_Irecv_large(state.recv_buffers[s][i].data(),
                                           state.recv_sizes[s][i], MPI_CHAR,
                                           state.peers[i], 143 + s,
                                           level_comm, state.recv_requests[s]);
                abort_on_mpi_error(ierr, "Irecv_large(staged-payload)",
                                   state.peers[i], state.recv_sizes[s][i]);
            }
        }
    }

    // ---- Pack and send, stage-priority order ----
    for (int s = 0; s < STAGED_HALO_STAGES; ++s) {
        for (size_t i = 0; i < num_peers; ++i) {
            if (state.send_sizes[s][i] == 0) {
                continue;
            }
            const auto& mortons = state.send_mortons[s][i];
            auto& buffer = state.send_buffers[s][i];
            buffer.resize(state.send_sizes[s][i]);

            std::vector<size_t> offsets(mortons.size());
            std::vector<size_t> box_sizes(mortons.size());
            size_t running = 0;
            for (size_t j = 0; j < mortons.size(); ++j) {
                const auto& box =
                    lvl.local_boxes[mortons[j] - lvl.local_morton_start];
                offsets[j] = running;
                box_sizes[j] =
                    (s == 0 ? get_staged_shell_size(box)
                            : get_staged_blocks_size(box, omit_for(i, mortons[j])));
                running += sizeof(size_t) + box_sizes[j];
            }
            if (running != state.send_sizes[s][i]) {
                throw std::runtime_error(
                    "initiate_staged_halo_gather: size bookkeeping mismatch");
            }

            std::atomic<bool> pack_error{false};
            #pragma omp parallel for schedule(dynamic)
            for (int64_t j = 0; j < static_cast<int64_t>(mortons.size()); ++j) {
                const auto& box = lvl.local_boxes[
                    mortons[static_cast<size_t>(j)] - lvl.local_morton_start];
                char* p = reinterpret_cast<char*>(buffer.data()) +
                          offsets[static_cast<size_t>(j)];
                const size_t box_size = box_sizes[static_cast<size_t>(j)];
                std::memcpy(p, &box_size, sizeof(size_t));
                p += sizeof(size_t);
                char* end =
                    (s == 0
                         ? serialize_staged_shell(box, p)
                         : serialize_staged_blocks(
                               box, p,
                               omit_for(i, mortons[static_cast<size_t>(j)])));
                if (end != p + box_size) {
                    pack_error.store(true, std::memory_order_relaxed);
                }
            }
            if (pack_error.load()) {
                throw std::runtime_error(
                    "initiate_staged_halo_gather: pack size mismatch for peer " +
                    std::to_string(state.peers[i]));
            }

            int ierr = MPI_Isend_large(buffer.data(), buffer.size(), MPI_CHAR,
                                       state.peers[i], 143 + s,
                                       level_comm, state.send_requests[s]);
            abort_on_mpi_error(ierr, "Isend_large(staged-payload)",
                               state.peers[i], buffer.size());
        }
    }

    state.active = true;
    lap(state.t_pack_ms);

    // Assisting point data (small, synchronous) while payloads are in flight.
    gather_CA_assisting_boxes_factorization(
        tree, level, level_comm, 411);
    lap(state.t_assist_ms);
}

/**
 * @brief Complete one stage: wait for its payloads and deserialize them.
 * Stage 0 materializes the ghost box shells; stages 1-3 attach the entry
 * blocks of that geometric group's ghosts.
 */
template<typename CoordType, typename DataType>
void staged_halo_wait_stage(
    TreeLevel<CoordType, DataType>& lvl,
    StagedHaloState<CoordType, DataType>& state,
    int s) {

    if (!state.active || state.stage_arrived[s]) {
        return;
    }

    // Opportunistically release earlier stages' send buffers: their receivers
    // have consumed them by the time we advance to a later stage.
    for (int t = 0; t < s; ++t) {
        if (!state.send_requests[t].empty()) {
            int done = 0;
            MPI_Testall(static_cast<int>(state.send_requests[t].size()),
                        state.send_requests[t].data(), &done, MPI_STATUSES_IGNORE);
            if (done) {
                state.send_requests[t].clear();
                for (auto& buf : state.send_buffers[t]) {
                    std::vector<StagedByte>().swap(buf);
                }
            }
        }
    }

    auto wait_start = std::chrono::high_resolution_clock::now();
    if (!state.recv_requests[s].empty()) {
        MPI_Waitall(static_cast<int>(state.recv_requests[s].size()),
                    state.recv_requests[s].data(), MPI_STATUSES_IGNORE);
        state.recv_requests[s].clear();
    }
    auto deser_start = std::chrono::high_resolution_clock::now();
    state.t_wait_ms[s] += std::chrono::duration<double, std::milli>(
        deser_start - wait_start).count();

    if (s == 0) {
        lvl.ghost_boxes.resize(lvl.ghost_id.size());
    }

    // Sequential header scan collects every box's span, then boxes
    // deserialize in parallel (distinct targets).
    {
        struct StagedUnpackItem {
            const char* ptr;
            size_t box_size;
            int64_t ghost_idx;
        };
        std::vector<StagedUnpackItem> items;

        for (size_t i = 0; i < state.peers.size(); ++i) {
            if (state.recv_sizes[s][i] == 0) {
                continue;
            }
            const char* ptr =
                reinterpret_cast<const char*>(state.recv_buffers[s][i].data());
            const char* end = ptr + state.recv_sizes[s][i];
            for (int64_t morton : state.recv_mortons[s][i]) {
                size_t box_size = 0;
                std::memcpy(&box_size, ptr, sizeof(size_t));
                ptr += sizeof(size_t);
                items.push_back(StagedUnpackItem{
                    ptr, box_size, lvl.ghost_id_to_index.at(morton)});
                ptr += box_size;
            }
            if (ptr != end) {
                throw std::runtime_error(
                    "staged_halo_wait_stage: buffer size mismatch from peer " +
                    std::to_string(state.peers[i]));
            }
        }

        std::atomic<bool> unpack_error{false};
        #pragma omp parallel for schedule(dynamic)
        for (int64_t j = 0; j < static_cast<int64_t>(items.size()); ++j) {
            const auto& item = items[static_cast<size_t>(j)];
            auto& box = lvl.ghost_boxes[item.ghost_idx];
            const char* box_end =
                (s == 0 ? deserialize_staged_shell(box, item.ptr)
                        : deserialize_staged_blocks_into(box, item.ptr));
            if (box_end != item.ptr + item.box_size) {
                unpack_error.store(true, std::memory_order_relaxed);
            }
        }
        if (unpack_error.load()) {
            throw std::runtime_error(
                "staged_halo_wait_stage: span mismatch during parallel unpack");
        }

        for (size_t i = 0; i < state.peers.size(); ++i) {
            state.recv_buffers[s][i].clear();
            state.recv_buffers[s][i].shrink_to_fit();
        }
    }
    // Without overlap every stage completes before elimination, so all shipped
    // views are still at entry state: rebuild the pruned ones by transpose
    // once the last stage is in. (Overlap mode rebuilds per-stage in
    // staged_halo_merge_stage instead, via entry transposes and mirrors.)
    if (s == STAGED_HALO_STAGES - 1 && !lvl.staged_overlap_on) {
        reconstruct_CA_symmetric_ghost_halo(lvl);
    }
    state.stage_arrived[s] = true;
    state.t_deser_ms[s] += std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - deser_start).count();

    // In overlap mode the diverted-update merge for this stage runs in
    // staged_halo_merge_stage (factorization.hpp), which needs the kernel to
    // materialize entry-fresh bases (e.g. leaf-level near blocks).
}

template<typename CoordType, typename DataType, typename KernelType>
void staged_halo_merge_stage(
    TreeLevel<CoordType, DataType>& level,
    StagedHaloState<CoordType, DataType>& state,
    int stage,
    KernelType* kernel) {
    if (!level.staged_overlap_on || stage < 1) {
        return;
    }

    std::unordered_set<int64_t> newly_arrived;
    for (size_t peer = 0; peer < state.peers.size(); ++peer) {
        for (int64_t morton : state.recv_mortons[stage][peer]) {
            newly_arrived.insert(morton);
            level.staged_unarrived.erase(morton);
        }
    }
    if (level.staged_pending == nullptr) {
        throw std::runtime_error(
            "staged_halo_merge_stage: pending state missing");
    }

    auto resolve_box = [&level](int64_t morton)
        -> BoxData<CoordType, DataType>* {
        if (auto* box = level.find_local_box(morton)) return box;
        return level.find_ghost_box(morton);
    };
    auto transpose_into_owned = [](
        const ModifiedBlock<DataType>& source,
        ModifiedBlock<DataType>& target) {
        const int64_t rows = source.a_ns_cols();
        const int64_t cols = source.a_ns_rows();
        std::vector<DataType> values(
            static_cast<size_t>(rows * cols));
        for (int64_t col = 0; col < cols; ++col) {
            for (int64_t row = 0; row < rows; ++row) {
                values[static_cast<size_t>(row + col * rows)] =
                    source.a_ns(col, row);
            }
        }
        target.symmetric_A_NS.reset();
        target.symmetric_A_NS_orientation = false;
        target.A_NS = MatrixStorage<DataType>();
        target.set_a_ns_owned(
            rows, cols, std::move(values),
            MatrixStorage<DataType>::FULL);
    };
    auto share_current_pair = [](
        ModifiedBlock<DataType>& partner,
        ModifiedBlock<DataType>& target) {
        if (partner.symmetric_A_NS &&
            partner.symmetric_A_NS == target.symmetric_A_NS) {
            return;
        }

        // A staged mirror deliberately supersedes the target's entry-state
        // view.  An earlier color may already have reduced the partner to its
        // skeleton, so the stale target and current partner transpose need
        // not have the same dimensions before this rebind.
        const int64_t expected_rows = partner.a_ns_cols();
        const int64_t expected_cols = partner.a_ns_rows();
        target.symmetric_A_NS.reset();
        target.symmetric_A_NS_orientation = false;
        target.A_NS = MatrixStorage<DataType>();
        target.share_transposed_a_ns_from(partner);
        if (target.a_ns_rows() != expected_rows ||
            target.a_ns_cols() != expected_cols) {
            throw std::runtime_error(
                "staged_halo_merge_stage: mirror rebind dimension mismatch");
        }
    };

    const int dimension = level.dimension;
    auto& pending = *level.staged_pending;
    std::lock_guard<std::mutex> lock(pending.mutex);

    // Recreate payloads omitted by symmetric staged pruning. Same-stage
    // partners remain independent until both sides' queued deltas are replayed.
    for (size_t peer = 0; peer < state.peers.size(); ++peer) {
        for (int64_t morton : state.recv_mortons[stage][peer]) {
            auto* arrived_box = level.find_ghost_box(morton);
            if (arrived_box == nullptr) {
                throw std::runtime_error(
                    "staged_halo_merge_stage: arrived ghost missing");
            }
            for (auto& block :
                 arrived_box->near_field_modified_interactions) {
                if (!block.reconstruct_A_NS_from_symmetric_pair) continue;
                block.reconstruct_A_NS_from_symmetric_pair = false;
                const int64_t partner_morton = block.neighbor_morton;
                auto* partner = resolve_box(partner_morton);
                if (partner == nullptr) {
                    throw std::runtime_error(
                        "staged_halo_merge_stage: pruned pair partner missing");
                }
                const auto partner_it =
                    partner->near_field_interaction_map.find(morton);
                if (partner_it ==
                    partner->near_field_interaction_map.end()) {
                    throw std::runtime_error(
                        "staged_halo_merge_stage: pruned reciprocal missing");
                }
                auto& partner_block =
                    partner->near_field_modified_interactions[
                        static_cast<size_t>(partner_it->second)];
                if (!partner_block.a_ns_is_allocated()) {
                    throw std::runtime_error(
                        "staged_halo_merge_stage: pruned payload missing");
                }
                if (newly_arrived.find(partner_morton) ==
                    newly_arrived.end()) {
                    pending.mirrors.insert(
                        {morton, partner_morton});
                } else {
                    transpose_into_owned(partner_block, block);
                    pending.mirrors.insert(
                        {morton, partner_morton});
                }
            }
        }
    }

    size_t kept = 0;
    for (size_t index = 0; index < pending.deltas.size(); ++index) {
        auto& delta = pending.deltas[index];
        if (level.staged_unarrived.find(delta.target_morton) !=
            level.staged_unarrived.end()) {
            if (kept != index) {
                pending.deltas[kept] = std::move(delta);
            }
            ++kept;
            continue;
        }

        if (!delta.is_schur) {
            const bool partner_unarrived =
                level.staged_unarrived.find(delta.neighbor_morton) !=
                level.staged_unarrived.end();
            const bool partner_new =
                newly_arrived.find(delta.neighbor_morton) !=
                newly_arrived.end();
            const bool partner_assisting =
                level.assisting_box_points_for_kernel_evaluation.find(
                    delta.neighbor_morton) !=
                level.assisting_box_points_for_kernel_evaluation.end();
            if (!partner_unarrived && !partner_new &&
                !partner_assisting) {
                pending.mirrors.insert(
                    {delta.target_morton, delta.neighbor_morton});
                continue;
            }
        }

        DeferredXnnAccumulatedTarget<DataType> target_state;
        target_state.target.box_morton = delta.target_morton;
        target_state.target.neighbor_morton = delta.neighbor_morton;
        target_state.target.kind = delta.is_schur
            ? DeferredXnnTargetKind::SCHUR
            : DeferredXnnTargetKind::NEAR_A_NS;
        target_state.rows = delta.rows;
        target_state.cols = delta.cols;
        target_state.data =
            materialize_deferred_xnn_target_matrix_for_accumulation(
                target_state.target, delta.rows, delta.cols,
                level, kernel, dimension);
        accumulate_deferred_xnn_matrix_in_place(
            target_state.data, delta.delta);
        flush_deferred_xnn_target_matrix_from_accumulation(
            target_state, level);
    }
    pending.deltas.resize(kept);

    for (auto it = pending.mirrors.begin();
         it != pending.mirrors.end();) {
        const int64_t target_morton = it->first;
        const int64_t partner_morton = it->second;
        if (level.staged_unarrived.find(target_morton) !=
                level.staged_unarrived.end() ||
            level.staged_unarrived.find(partner_morton) !=
                level.staged_unarrived.end()) {
            ++it;
            continue;
        }

        auto* target = resolve_box(target_morton);
        auto* partner = resolve_box(partner_morton);
        if (target == nullptr || partner == nullptr) {
            throw std::runtime_error(
                "staged_halo_merge_stage: mirror endpoint missing");
        }
        const auto partner_it =
            partner->near_field_interaction_map.find(target_morton);
        if (partner_it == partner->near_field_interaction_map.end()) {
            throw std::runtime_error(
                "staged_halo_merge_stage: mirror partner view missing");
        }
        auto& partner_block =
            partner->near_field_modified_interactions[
                static_cast<size_t>(partner_it->second)];
        if (!partner_block.a_ns_is_allocated()) {
            throw std::runtime_error(
                "staged_halo_merge_stage: mirror partner storage missing");
        }

        auto target_it =
            target->near_field_interaction_map.find(partner_morton);
        if (target_it == target->near_field_interaction_map.end()) {
            ModifiedBlock<DataType> block;
            block.neighbor_morton = partner_morton;
            const int64_t block_index = static_cast<int64_t>(
                target->near_field_modified_interactions.size());
            target->near_field_modified_interactions.push_back(
                std::move(block));
            target->near_field_interaction_map[partner_morton] =
                block_index;
            target_it =
                target->near_field_interaction_map.find(partner_morton);
        }
        auto& target_block =
            target->near_field_modified_interactions[
                static_cast<size_t>(target_it->second)];
        share_current_pair(partner_block, target_block);
        target->num_near_field_interactions = static_cast<int64_t>(
            target->near_field_modified_interactions.size());
        it = pending.mirrors.erase(it);
    }
}

/// Drain outstanding sends and release the staged buffers.
template<typename CoordType, typename DataType>
void staged_halo_finish(StagedHaloState<CoordType, DataType>& state) {
    if (!state.active) {
        return;
    }
    for (int s = 0; s < STAGED_HALO_STAGES; ++s) {
        if (!state.send_requests[s].empty()) {
            MPI_Waitall(static_cast<int>(state.send_requests[s].size()),
                        state.send_requests[s].data(), MPI_STATUSES_IGNORE);
            state.send_requests[s].clear();
        }
    }
    for (int s = 0; s < STAGED_HALO_STAGES; ++s) {
        for (auto& buf : state.send_buffers[s]) {
            std::vector<StagedByte>().swap(buf);
        }
        for (auto& buf : state.recv_buffers[s]) {
            std::vector<StagedByte>().swap(buf);
        }
    }
    state.active = false;
}

} // namespace fmm
