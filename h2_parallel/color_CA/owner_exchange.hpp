#pragma once
// ---------------------------------------------------------------------------
// owner_exchange.hpp — component-owner elimination, step 1: faces only.
//
// Design: concise-algorithm-ref/component-owner-dataflow-design.md, §10.
//
// Legacy mode 1: blue and purple stay replicated exactly as today; each
// green face patch is eliminated by ONE of its sharers (dataflow.hpp
// assignment).
// Legacy mode 2: every boundary component (blue cluster incl. orange,
// purple run, green face) is eliminated by one owner — no replicated
// elimination at all.  Today's variable-depth halo already holds every box
// of every component a rank shares, so no new fetch is needed.
// After every owned-colour wave's per-box region the owners ship the
// post-`compute_and_modify` state of the sources they eliminated to every
// rank whose halo holds those boxes; receivers install it and then run the
// unchanged candidates / owner / mirror / finalize passes.  In the production
// (deferred, symmetric) path `compute_and_modify` writes only the box's own
// fields, so "install + replay the passes" reproduces the replicated run
// bit for bit.
//
// This is the synchronous stepping stone (one exchange per green wave).  The
// counter-driven, per-component, overlapped schedule (dataflow::Runtime)
// replaces the per-wave exchange in the next step.
// ---------------------------------------------------------------------------

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <array>
#include <mpi.h>
#include <omp.h>

#include "tree_impl.hpp"
#include "serialization.hpp"
#include "dataflow.hpp"

namespace fmm {

inline int owner_component_mode() {
    return ca_owner_component_mode();
}

inline bool box_checksum_enabled() {
    return false;
}

inline void owner_mpi_check(int ierr, const char* what) {
    if (ierr != MPI_SUCCESS) {
        throw std::runtime_error(std::string("owner_exchange: ") + what + " failed with MPI error " +
                                 std::to_string(ierr));
    }
}

// ===========================================================================
// Eliminated-source state (post compute_and_modify, pre finalize)
// ===========================================================================
// Everything the per-wave passes, the transition and the solve read from a
// source box: the full BoxData serializer (indices, T, X_RR + pivots/format,
// X_RS, X_SR, Schur, X_RN, X_NR, near/far blocks + maps, use_full_set) plus
// the fields it does not cover: X_RR_full (lazy far-field) and the deferred
// temp2 rows with their per-neighbor row counts (consumed by the owner pass,
// moved into X_NR by finalize).

template<typename CoordType, typename DataType>
size_t eliminated_state_size(const BoxData<CoordType, DataType>& box) {
    size_t size = get_serialized_size(box);
    size += get_serialized_size(box.X_RR_full);
    size += sizeof(size_t) + box.deferred_xnn_temp2.size() * sizeof(DataType);
    size += sizeof(size_t) + box.deferred_xnn_neighbor_point_counts.size() * sizeof(int64_t);
    return size;
}

template<typename CoordType, typename DataType>
char* serialize_eliminated_state(const BoxData<CoordType, DataType>& box, char* buffer) {
    buffer = serialize(box, buffer);
    buffer = serialize(box.X_RR_full, buffer);
    size_t n = box.deferred_xnn_temp2.size();
    std::memcpy(buffer, &n, sizeof(size_t));
    buffer += sizeof(size_t);
    if (n > 0) {
        std::memcpy(buffer, box.deferred_xnn_temp2.data(), n * sizeof(DataType));
        buffer += n * sizeof(DataType);
    }
    n = box.deferred_xnn_neighbor_point_counts.size();
    std::memcpy(buffer, &n, sizeof(size_t));
    buffer += sizeof(size_t);
    if (n > 0) {
        std::memcpy(buffer, box.deferred_xnn_neighbor_point_counts.data(), n * sizeof(int64_t));
        buffer += n * sizeof(int64_t);
    }
    return buffer;
}

template<typename CoordType, typename DataType>
const char* deserialize_eliminated_state(BoxData<CoordType, DataType>& box, const char* buffer) {
    buffer = deserialize(box, buffer);
    buffer = deserialize(box.X_RR_full, buffer);
    size_t n = 0;
    std::memcpy(&n, buffer, sizeof(size_t));
    buffer += sizeof(size_t);
    box.deferred_xnn_temp2.resize(n);
    if (n > 0) {
        std::memcpy(box.deferred_xnn_temp2.data(), buffer, n * sizeof(DataType));
        buffer += n * sizeof(DataType);
    }
    std::memcpy(&n, buffer, sizeof(size_t));
    buffer += sizeof(size_t);
    box.deferred_xnn_neighbor_point_counts.resize(n);
    if (n > 0) {
        std::memcpy(box.deferred_xnn_neighbor_point_counts.data(), buffer, n * sizeof(int64_t));
        buffer += n * sizeof(int64_t);
    }
    return buffer;
}

// ===========================================================================
// Per-level state
// ===========================================================================

template<typename CoordType, typename DataType>
struct OwnerFacesLevelState {
    bool active = false;
    dataflow::ProcessGrid grid;
    dataflow::ComponentGraph graph;
    uint32_t my_pid = 0;
    int my_rank = -1;
    std::vector<dataflow::Color> owned_colors;     ///< colours under ownership (mode 1: green; mode 2: all)
    std::unordered_set<int64_t> eliminate_here;   ///< owned-colour boxes this rank eliminates
    int64_t present = 0;                           ///< owned-colour boxes present (local + ghost)
    int64_t owned_local = 0;                       ///< of those, eliminated here and home here
    int64_t owned_remote = 0;                      ///< eliminated here, home elsewhere
    std::vector<int> peer_ranks;                   ///< MPI ranks of active process-grid neighbours
    std::vector<uint32_t> peer_pids;
    /// ship_to[i]: boxes this rank eliminates that peer i holds (peer i asked
    /// for them at setup) — the exact receiver set, no geometric guessing.
    std::vector<std::unordered_set<int64_t>> ship_to;
    // timing / volume
    double t_pack_ms = 0, t_sizes_ms = 0, t_wire_ms = 0, t_install_ms = 0;
    size_t bytes_sent = 0, bytes_recv = 0;
    int exchanges = 0;

    /// Whether the CA colour group named `color_name` ("blue", "purple",
    /// "green"; "boundary"/"interior" are never owned) is under ownership.
    bool owns_color(const std::string& color_name) const {
        if (!active) return false;
        dataflow::Color c;
        if (color_name == "blue") c = dataflow::Color::BLUE;
        else if (color_name == "purple") c = dataflow::Color::PURPLE;
        else if (color_name == "green") c = dataflow::Color::GREEN;
        else return false;
        return std::find(owned_colors.begin(), owned_colors.end(), c) != owned_colors.end();
    }
    /// Whether this rank eliminates box `morton` of an owned colour group.
    bool eliminates(int64_t morton) const { return !active || eliminate_here.count(morton) != 0; }
};

/// Build the component graph for this level, assign face owners, and decide
/// which green boxes this rank eliminates.  Collective in effect (every rank
/// computes the identical graph); no communication.
template<typename CoordType, typename DataType>
void owner_faces_setup(TreeLevel<CoordType, DataType>& level, int32_t level_index, int dimension, int rank,
                       OwnerFacesLevelState<CoordType, DataType>& st) {
    using namespace dataflow;
    st = OwnerFacesLevelState<CoordType, DataType>{};
    if (!level.is_process_active || level.num_active_processes <= 1 || level.my_morton_id < 0) return;

    st.grid = ProcessGrid(dimension, 1u << level_index, static_cast<uint32_t>(level.num_active_processes));
    ColorOf color_of = make_ca_classifier(
        st.grid, level_index,
        [](int64_t m, int32_t lvl, int32_t dim, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_blue_box(m, lvl, dim, gs, lgs, off);
        },
        [](int64_t m, int32_t lvl, int32_t dim, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_orange_box(m, lvl, dim, gs, lgs, off);
        },
        [](int64_t m, int32_t lvl, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_purple_box(m, lvl, gs, lgs, off);
        });
    st.graph = build_component_graph(st.grid, color_of);
    assign_owners(st.graph, AssignPolicy::LPT_PER_COLOR);
    {
        const std::string e = validate_graph(st.graph);
        if (!e.empty()) throw std::runtime_error("owner_faces_setup: invalid graph: " + e);
        const std::string l = validate_against_lists(st.graph, level.blue, level.orange, level.purple, level.green);
        if (!l.empty()) throw std::runtime_error("owner_faces_setup: graph disagrees with level lists: " + l);
    }
    st.my_pid = static_cast<uint32_t>(level.my_morton_id);
    st.my_rank = rank;
    if (owner_component_mode() >= 2) st.owned_colors = boundary_colors(dimension);
    else st.owned_colors = {Color::GREEN};
    auto owned_color = [&](Color c) {
        return std::find(st.owned_colors.begin(), st.owned_colors.end(), c) != st.owned_colors.end();
    };

    // Boxes this rank eliminates: those of owned components of owned colours.
    // Every such box must be present here (local or ghost): the variable-depth
    // halo holds every box of every component the rank shares.
    for (const auto& k : st.graph.comps) {
        if (!owned_color(k.color) || k.owner != static_cast<int32_t>(st.my_pid)) continue;
        for (int64_t m : k.boxes) {
            const bool local = level.find_local_box(m) != nullptr;
            if (!local && level.find_ghost_box(m) == nullptr)
                throw std::runtime_error("owner_faces_setup: owned " + std::string(color_name(k.color)) + " box " +
                                         std::to_string(m) + " of component " + std::to_string(k.id) +
                                         " is not held by rank " + std::to_string(rank));
            st.eliminate_here.insert(m);
            if (local) ++st.owned_local; else ++st.owned_remote;
        }
    }
    // Consistency: every owned-colour box present here belongs to a component
    // whose owner will ship it to us (we are in its visible set) unless we own it.
    auto check_list = [&](const std::vector<int64_t>& list) {
        for (int64_t m : list) {
            ++st.present;
            const Component& k = st.graph.comp_of_box(m);
            if (k.owner == static_cast<int32_t>(st.my_pid)) continue;
            if (!k.is_visible(st.my_pid))
                throw std::runtime_error("owner_faces_setup: box " + std::to_string(m) + " present on rank " +
                                         std::to_string(rank) + " but its component " + std::to_string(k.id) +
                                         " does not list this rank as visible");
        }
    };
    if (owned_color(Color::BLUE)) { check_list(level.blue); check_list(level.orange); }
    if (owned_color(Color::PURPLE)) check_list(level.purple);
    if (owned_color(Color::GREEN)) check_list(level.green);

    // Peers: active ranks whose brick is within one brick of mine.
    const auto pc = st.grid.proc_coords(st.my_pid);
    const int zr = dimension == 3 ? 1 : 0;
    for (int dz = -zr; dz <= zr; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) continue;
                const int64_t px = static_cast<int64_t>(pc[0]) + dx, py = static_cast<int64_t>(pc[1]) + dy,
                              pz = static_cast<int64_t>(pc[2]) + dz;
                const int64_t P = st.grid.procs_per_dim;
                if (px < 0 || py < 0 || pz < 0 || px >= P || py >= P || pz >= P) continue;
                const uint32_t q = static_cast<uint32_t>(dimension == 2
                    ? morton::encode_2d(static_cast<uint32_t>(px), static_cast<uint32_t>(py))
                    : morton::encode_3d(static_cast<uint32_t>(px), static_cast<uint32_t>(py), static_cast<uint32_t>(pz)));
                auto it = level.morton_to_rank.find(static_cast<int>(q));
                if (it == level.morton_to_rank.end())
                    throw std::runtime_error("owner_faces_setup: no MPI rank for process id " + std::to_string(q));
                st.peer_pids.push_back(q);
                st.peer_ranks.push_back(it->second);
            }

    // Request handshake: tell each peer which of its eliminations we hold.
    // Exact by construction (the receiver enumerates what it holds), so the
    // per-wave exchange never guesses the halo's contents.
    {
        const size_t npeers = st.peer_ranks.size();
        std::vector<std::vector<int64_t>> request(npeers);
        auto collect = [&](const std::vector<int64_t>& list) {
            for (int64_t m : list) {
                const Component& k = st.graph.comp_of_box(m);
                if (k.owner == static_cast<int32_t>(st.my_pid)) continue;
                for (size_t i = 0; i < npeers; ++i)
                    if (st.peer_pids[i] == static_cast<uint32_t>(k.owner)) { request[i].push_back(m); break; }
            }
        };
        if (owned_color(Color::BLUE)) { collect(level.blue); collect(level.orange); }
        if (owned_color(Color::PURPLE)) collect(level.purple);
        if (owned_color(Color::GREEN)) collect(level.green);

        std::vector<int64_t> send_counts(npeers), recv_counts(npeers, 0);
        for (size_t i = 0; i < npeers; ++i) send_counts[i] = static_cast<int64_t>(request[i].size());
        std::vector<MPI_Request> reqs;
        reqs.reserve(4 * npeers);
        for (size_t i = 0; i < npeers; ++i) {
            MPI_Request r;
            owner_mpi_check(MPI_Irecv(&recv_counts[i], 1, MPI_INT64_T, st.peer_ranks[i], 3203, MPI_COMM_WORLD, &r),
                            "MPI_Irecv(request-counts)");
            reqs.push_back(r);
            owner_mpi_check(MPI_Isend(&send_counts[i], 1, MPI_INT64_T, st.peer_ranks[i], 3203, MPI_COMM_WORLD, &r),
                            "MPI_Isend(request-counts)");
            reqs.push_back(r);
        }
        if (!reqs.empty()) {
            owner_mpi_check(MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE),
                            "MPI_Waitall(request-counts)");
            reqs.clear();
        }
        std::vector<std::vector<int64_t>> received(npeers);
        for (size_t i = 0; i < npeers; ++i) {
            received[i].resize(static_cast<size_t>(recv_counts[i]));
            if (recv_counts[i] > 0) {
                MPI_Request r;
                owner_mpi_check(MPI_Irecv(received[i].data(), static_cast<int>(recv_counts[i]), MPI_INT64_T,
                                          st.peer_ranks[i], 3204, MPI_COMM_WORLD, &r),
                                "MPI_Irecv(request-lists)");
                reqs.push_back(r);
            }
            if (send_counts[i] > 0) {
                MPI_Request r;
                owner_mpi_check(MPI_Isend(request[i].data(), static_cast<int>(send_counts[i]), MPI_INT64_T,
                                          st.peer_ranks[i], 3204, MPI_COMM_WORLD, &r),
                                "MPI_Isend(request-lists)");
                reqs.push_back(r);
            }
        }
        if (!reqs.empty()) {
            owner_mpi_check(MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE),
                            "MPI_Waitall(request-lists)");
        }
        st.ship_to.assign(npeers, {});
        for (size_t i = 0; i < npeers; ++i) {
            for (int64_t m : received[i]) {
                if (st.eliminate_here.count(m) == 0)
                    throw std::runtime_error("owner_faces_setup: rank " + std::to_string(st.peer_ranks[i]) +
                                             " requested box " + std::to_string(m) + " which rank " +
                                             std::to_string(rank) + " does not eliminate");
                st.ship_to[i].insert(m);
            }
        }
    }
    st.active = true;
}

/// One-line summary for the level banner.
template<typename CoordType, typename DataType>
std::string owner_faces_describe(const OwnerFacesLevelState<CoordType, DataType>& st) {
    if (!st.active) return "owner-faces: inactive";
    std::string colors;
    for (dataflow::Color c : st.owned_colors) colors += std::string(colors.empty() ? "" : "+") + dataflow::color_name(c);
    return "owner-" + colors + ": " + dataflow::describe(st.graph) + " | this rank eliminates " +
           std::to_string(st.owned_local + st.owned_remote) + " of " + std::to_string(st.present) +
           " owned-colour boxes present (" + std::to_string(st.owned_remote) + " with a remote home), peers " +
           std::to_string(st.peer_ranks.size());
}

// ===========================================================================
// Per-wave exchange of eliminated-source state
// ===========================================================================
// After the per-box region of a green wave: every source of the wave that this
// rank eliminated is shipped to each peer whose halo holds it; every source of
// the wave present here but eliminated elsewhere is received from its owner
// and installed in place (local or ghost BoxData).  Sizes first (one uint64
// per peer, zero allowed), then payloads; packing is box-parallel.

template<typename CoordType, typename DataType>
void owner_faces_exchange_wave(TreeLevel<CoordType, DataType>& level,
                               const std::vector<int64_t>& wave_boxes,
                               OwnerFacesLevelState<CoordType, DataType>& st) {
    using clock = std::chrono::high_resolution_clock;
    if (!st.active) return;
    const int tag_sizes = 3201, tag_payload = 3202;
    const size_t npeers = st.peer_ranks.size();
    auto t0 = clock::now();
    auto lap = [&](double& acc) {
        auto t1 = clock::now();
        acc += std::chrono::duration<double, std::milli>(t1 - t0).count();
        t0 = t1;
    };

    // ---- what goes where -------------------------------------------------
    std::vector<std::vector<int64_t>> send_boxes(npeers);
    std::vector<std::vector<size_t>> send_box_sizes(npeers);
    std::vector<size_t> send_sizes(npeers, 0), recv_sizes(npeers, 0);
    std::vector<int64_t> expected_from(npeers, 0);
    for (int64_t m : wave_boxes) {
        const dataflow::Component& k = st.graph.comp_of_box(m);
        if (k.owner == static_cast<int32_t>(st.my_pid)) {
            BoxData<CoordType, DataType>* box = level.find_local_box(m);
            if (box == nullptr) box = level.find_ghost_box(m);
            if (box == nullptr) throw std::runtime_error("owner_faces_exchange_wave: eliminated box not found");
            const size_t sz = eliminated_state_size(*box);
            for (size_t i = 0; i < npeers; ++i) {
                if (st.ship_to[i].count(m) == 0) continue;
                send_boxes[i].push_back(m);
                send_box_sizes[i].push_back(sz);
                send_sizes[i] += sizeof(int64_t) + sizeof(size_t) + sz;
            }
        } else {
            // present here (wave boxes are), so its owner must send it
            size_t i = 0;
            for (; i < npeers; ++i) if (st.peer_pids[i] == static_cast<uint32_t>(k.owner)) break;
            if (i == npeers)
                throw std::runtime_error("owner_faces_exchange_wave: owner of box " + std::to_string(m) +
                                         " (pid " + std::to_string(k.owner) + ") is not a peer of rank " +
                                         std::to_string(st.my_rank));
            ++expected_from[i];
        }
    }
    lap(st.t_pack_ms);

    // ---- sizes -------------------------------------------------------------
    std::vector<MPI_Request> requests;
    requests.reserve(2 * npeers);
    for (size_t i = 0; i < npeers; ++i) {
        MPI_Request req;
        owner_mpi_check(MPI_Irecv(&recv_sizes[i], 1, MPI_UINT64_T, st.peer_ranks[i], tag_sizes, MPI_COMM_WORLD, &req),
                        "MPI_Irecv(sizes)");
        requests.push_back(req);
    }
    for (size_t i = 0; i < npeers; ++i) {
        MPI_Request req;
        owner_mpi_check(MPI_Isend(&send_sizes[i], 1, MPI_UINT64_T, st.peer_ranks[i], tag_sizes, MPI_COMM_WORLD, &req),
                        "MPI_Isend(sizes)");
        requests.push_back(req);
    }
    if (!requests.empty()) {
        owner_mpi_check(MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
                        "MPI_Waitall(sizes)");
        requests.clear();
    }
    lap(st.t_sizes_ms);

    // ---- payload receives -------------------------------------------------
    std::vector<std::vector<char>> recv_buffers(npeers);
    for (size_t i = 0; i < npeers; ++i) {
        if (recv_sizes[i] == 0) continue;
        recv_buffers[i].resize(recv_sizes[i]);
        owner_mpi_check(MPI_Irecv_large(recv_buffers[i].data(), recv_sizes[i], MPI_CHAR, st.peer_ranks[i],
                                        tag_payload, MPI_COMM_WORLD, requests),
                        "MPI_Irecv_large(payload)");
        st.bytes_recv += recv_sizes[i];
    }

    // ---- pack + send (double buffered) --------------------------------------
    {
        std::vector<char> send_buffers[2];
        std::vector<MPI_Request> send_requests[2];
        int parity = 0;
        for (size_t i = 0; i < npeers; ++i) {
            if (send_sizes[i] == 0) continue;
            if (!send_requests[parity].empty()) {
                owner_mpi_check(MPI_Waitall(static_cast<int>(send_requests[parity].size()),
                                            send_requests[parity].data(), MPI_STATUSES_IGNORE),
                                "MPI_Waitall(send-buffer)");
                send_requests[parity].clear();
            }
            lap(st.t_wire_ms);

            const auto& mortons = send_boxes[i];
            const auto& sizes = send_box_sizes[i];
            std::vector<size_t> offsets(mortons.size());
            size_t running = 0;
            for (size_t j = 0; j < mortons.size(); ++j) {
                offsets[j] = running;
                running += sizeof(int64_t) + sizeof(size_t) + sizes[j];
            }
            if (running != send_sizes[i]) throw std::runtime_error("owner_faces_exchange_wave: size bookkeeping");

            auto& buffer = send_buffers[parity];
            buffer.resize(send_sizes[i]);
            std::atomic<bool> pack_error{false};
            #pragma omp parallel for schedule(dynamic)
            for (int64_t j = 0; j < static_cast<int64_t>(mortons.size()); ++j) {
                const int64_t m = mortons[static_cast<size_t>(j)];
                BoxData<CoordType, DataType>* box = level.find_local_box(m);
                if (box == nullptr) box = level.find_ghost_box(m);
                char* p = buffer.data() + offsets[static_cast<size_t>(j)];
                std::memcpy(p, &m, sizeof(int64_t));
                p += sizeof(int64_t);
                std::memcpy(p, &sizes[static_cast<size_t>(j)], sizeof(size_t));
                p += sizeof(size_t);
                char* end = serialize_eliminated_state(*box, p);
                if (end != p + sizes[static_cast<size_t>(j)]) pack_error.store(true, std::memory_order_relaxed);
            }
            if (pack_error.load()) throw std::runtime_error("owner_faces_exchange_wave: pack size mismatch");
            lap(st.t_pack_ms);

            owner_mpi_check(MPI_Isend_large(buffer.data(), send_sizes[i], MPI_CHAR, st.peer_ranks[i], tag_payload,
                                            MPI_COMM_WORLD, send_requests[parity]),
                            "MPI_Isend_large(payload)");
            st.bytes_sent += send_sizes[i];
            parity ^= 1;
        }
        for (int p = 0; p < 2; ++p) {
            if (send_requests[p].empty()) continue;
            owner_mpi_check(MPI_Waitall(static_cast<int>(send_requests[p].size()), send_requests[p].data(),
                                        MPI_STATUSES_IGNORE),
                            "MPI_Waitall(send-drain)");
            send_requests[p].clear();
        }
        if (!requests.empty()) {
            owner_mpi_check(MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
                            "MPI_Waitall(payload)");
            requests.clear();
        }
        lap(st.t_wire_ms);
    }

    // ---- install -------------------------------------------------------------
    for (size_t i = 0; i < npeers; ++i) {
        int64_t installed = 0;
        const char* p = recv_buffers[i].data();
        const char* end = p + recv_buffers[i].size();
        while (p < end) {
            int64_t m = 0;
            size_t sz = 0;
            std::memcpy(&m, p, sizeof(int64_t));
            p += sizeof(int64_t);
            std::memcpy(&sz, p, sizeof(size_t));
            p += sizeof(size_t);
            BoxData<CoordType, DataType>* box = level.find_local_box(m);
            if (box == nullptr) box = level.find_ghost_box(m);
            if (box == nullptr)
                throw std::runtime_error("owner_faces_exchange_wave: received state for box " + std::to_string(m) +
                                         " which rank " + std::to_string(st.my_rank) + " does not hold");
            const char* q = deserialize_eliminated_state(*box, p);
            if (q != p + sz) throw std::runtime_error("owner_faces_exchange_wave: unpack size mismatch");
            p = q;
            ++installed;
        }
        if (installed != expected_from[i])
            throw std::runtime_error("owner_faces_exchange_wave: rank " + std::to_string(st.my_rank) + " expected " +
                                     std::to_string(expected_from[i]) + " sources from rank " +
                                     std::to_string(st.peer_ranks[i]) + ", installed " + std::to_string(installed));
    }
    ++st.exchanges;
    lap(st.t_install_ms);
}

template<typename CoordType, typename DataType>
void owner_faces_print_timing(const OwnerFacesLevelState<CoordType, DataType>& st, int level_index) {
    if (!st.active) return;
    std::printf("  [owner-faces] level %d: %d wave exchanges, sent %.1f MB, recv %.1f MB | ms: bookkeeping+pack %.1f, "
                "sizes %.1f, wire %.1f, install %.1f\n",
                level_index, st.exchanges, st.bytes_sent / 1048576.0, st.bytes_recv / 1048576.0, st.t_pack_ms,
                st.t_sizes_ms, st.t_wire_ms, st.t_install_ms);
    std::fflush(stdout);
}

// ===========================================================================
// Per-level checksum oracle (FMM_BOX_CHECKSUM=1)
// ===========================================================================
// Fixed-order sums of |entry| over the local boxes' factors; gathered to the
// print rank and summed in rank order, so two runs that agree bit for bit
// print identical lines.  Collective over comm (inactive ranks contribute 0).

template<typename CoordType, typename DataType>
void print_level_checksum(TreeLevel<CoordType, DataType>& level, int level_index, MPI_Comm comm, int print_rank,
                          int rank, bool active) {
    if (!box_checksum_enabled()) return;
    constexpr int NV = 6;
    double v[NV] = {0, 0, 0, 0, 0, 0};
    auto sum_abs = [](const std::vector<DataType>& d) {
        double s = 0.0;
        for (const auto& x : d) s += static_cast<double>(std::abs(x));
        return s;
    };
    if (active) {
        // Per-box sums in parallel (each one sequential, as before), folded
        // into v in box and block order afterwards — bitwise the same values
        // as the earlier serial loop, ~60 GB less single-threaded reading per
        // fine level (this was ~20 s of "Total level time" in every mode).
        const size_t nb = level.local_boxes.size();
        std::vector<std::array<double, 5>> per(nb);
        std::vector<std::vector<double>> blk(nb);
        #pragma omp parallel for schedule(dynamic)
        for (int64_t i = 0; i < static_cast<int64_t>(nb); ++i) {
            const auto& box = level.local_boxes[static_cast<size_t>(i)];
            auto& p = per[static_cast<size_t>(i)];
            p[0] = static_cast<double>(box.skeleton_indices.size());
            p[1] = sum_abs(box.interpolation_matrix.data);
            p[2] = sum_abs(box.X_RR.data);
            p[3] = sum_abs(box.X_NR.data);
            p[4] = sum_abs(box.schur_complement.data);
            auto& b = blk[static_cast<size_t>(i)];
            b.reserve(box.near_field_modified_interactions.size());
            for (const auto& m : box.near_field_modified_interactions) b.push_back(sum_abs(m.A_NS.data));
        }
        for (size_t i = 0; i < nb; ++i) {
            v[0] += per[i][0];
            v[1] += per[i][1];
            v[2] += per[i][2];
            v[3] += per[i][3];
            for (double x : blk[i]) v[4] += x;
            v[5] += per[i][4];
        }
    }
    int size = 1;
    MPI_Comm_size(comm, &size);
    std::vector<double> all(static_cast<size_t>(NV) * static_cast<size_t>(size), 0.0);
    MPI_Gather(v, NV, MPI_DOUBLE, all.data(), NV, MPI_DOUBLE, print_rank, comm);
    if (rank != print_rank) return;
    double tot[NV] = {0, 0, 0, 0, 0, 0};
    for (int r = 0; r < size; ++r)
        for (int j = 0; j < NV; ++j) tot[j] += all[static_cast<size_t>(r) * NV + static_cast<size_t>(j)];
    std::printf("  [checksum] level %d: skel %.0f | T %.17g | X_RR %.17g | X_NR %.17g | near %.17g | schur %.17g\n",
                level_index, tot[0], tot[1], tot[2], tot[3], tot[4], tot[5]);
    std::fflush(stdout);
}

}  // namespace fmm
