#pragma once
// A component-owner CA level (H2_CA_owner_component=3) on the device
// (h2_gpu/GPU_CA_PLAN.md, M5).
//
// The host schedule (core/owner_schedule.hpp: the component graph, the
// owner assignment, the component-aware wave partition, the payload
// recipients, the sketch and interest filters) is built as for the CPU
// path; this file turns it into the eliminator's plan (owner_plan.hpp),
// fetches the entry blocks of the owned components' remote boxes between
// the devices when a device transition kept the level's blocks, and runs
// the level: the synchronous colour schedule (per boundary colour, the
// owned components' boxes in their waves, then the colour's transport,
// which ships the owned boxes' state to its recipients and installs the
// colour's arrivals), then the interior as a replicated CA level runs it.
// The factors of this rank's boxes eliminated elsewhere land here (FULL),
// so the solves and multiplies run the level as a Color level
// (TreeLevel::ca_home_factors).
#include "factorization_driver.hpp"
#include "owner_plan.hpp"
#include "core/owner_schedule.hpp"
#include "gpu_common/bpack_env.hpp"

#ifdef H2_HAVE_GPU

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// The eliminator's view of the schedule.
template<typename CoordType, typename DataType>
OwnerPlanView build_owner_plan(ParallelTree<CoordType, DataType>* tree, int lvl,
                               const OwnerScheduleState<CoordType, DataType>& st) {
    using dataflow::Component;
    const auto& level = tree->levels[static_cast<size_t>(lvl)];
    if (!st.active) throw std::runtime_error("build_owner_plan: the owner schedule is not active");
    OwnerPlanView v;
    v.num_boundary_colors = st.num_boundary_colors;
    v.num_waves = st.num_waves;
    v.interior_seq0 = st.num_boundary_colors * st.num_waves;
    auto rank_of_pid = [&](uint32_t pid) {
        const auto it = level.morton_to_rank.find(static_cast<int>(pid));
        if (it == level.morton_to_rank.end()) throw std::runtime_error("build_owner_plan: no rank for region " + std::to_string(pid));
        return it->second;
    };
    auto pid_to_rank = std::make_shared<std::unordered_map<int, uint32_t>>();  // tree rank -> region
    for (const auto& [pid, rank] : level.morton_to_rank) (*pid_to_rank)[rank] = static_cast<uint32_t>(pid);
    const uint32_t my_pid = st.my_pid;
    for (const Component& k : st.graph.comps) {
        const int32_t ci = st.comp_seq_color[static_cast<size_t>(k.id)];
        for (int64_t m : k.boxes) {
            v.comp_of_box[m] = k.id;
            const auto w = st.box_wave.find(m);   // (the waves of the components eliminated or installed here)
            if (w != st.box_wave.end()) v.seq_of_box[m] = ci * st.num_waves + w->second;
        }
    }
    v.colours.resize(static_cast<size_t>(st.num_boundary_colors));
    for (auto& c : v.colours) c.owned_waves.resize(static_cast<size_t>(st.num_waves));
    std::vector<std::set<int>> send(static_cast<size_t>(st.num_boundary_colors)), recv(static_cast<size_t>(st.num_boundary_colors));
    for (int32_t id : st.rt.owned()) {
        const Component& k = st.graph.comps[static_cast<size_t>(id)];
        const int32_t ci = st.comp_seq_color[static_cast<size_t>(id)];
        for (int64_t m : k.boxes) {
            if (st.grid.proc_of_box(m) != my_pid) v.owned_ghosts.insert(m);
            v.colours[static_cast<size_t>(ci)].owned_waves[static_cast<size_t>(st.wave_of(m))].push_back(m);
            std::vector<std::pair<int, OwnerPayload>> out;
            for (const auto& [pid, kind] : dataflow::payload_recipients(st.graph, my_pid, m)) {
                const int rank = rank_of_pid(pid);
                out.emplace_back(rank, static_cast<OwnerPayload>(kind));
                send[static_cast<size_t>(ci)].insert(rank);
            }
            if (!out.empty()) v.recipients[m] = std::move(out);
        }
    }
    for (auto& c : v.colours) {
        for (auto& w : c.owned_waves) std::sort(w.begin(), w.end());
    }
    for (int32_t id : st.rt.shared_not_owned()) {
        const Component& k = st.graph.comps[static_cast<size_t>(id)];
        const int32_t ci = st.comp_seq_color[static_cast<size_t>(id)];
        recv[static_cast<size_t>(ci)].insert(rank_of_pid(static_cast<uint32_t>(k.owner)));
        for (int64_t m : k.boxes) {
            if (st.grid.proc_of_box(m) == my_pid) v.foreign_local.insert(m);
        }
    }
    for (int32_t id : dataflow::auxiliary_components(st.graph, my_pid)) {
        const Component& k = st.graph.comps[static_cast<size_t>(id)];
        const int32_t ci = st.comp_seq_color[static_cast<size_t>(id)];
        recv[static_cast<size_t>(ci)].insert(rank_of_pid(static_cast<uint32_t>(k.owner)));
    }
    for (int c = 0; c < st.num_boundary_colors; ++c) {
        v.colours[static_cast<size_t>(c)].send_peers.assign(send[static_cast<size_t>(c)].begin(), send[static_cast<size_t>(c)].end());
        v.colours[static_cast<size_t>(c)].recv_peers.assign(recv[static_cast<size_t>(c)].begin(), recv[static_cast<size_t>(c)].end());
    }
    const OwnerScheduleState<CoordType, DataType>* sp = &st;
    v.interest_of = [sp, pid_to_rank](int rank, int64_t m) {
        const auto it = pid_to_rank->find(rank);
        if (it == pid_to_rank->end()) throw std::runtime_error("owner plan: rank " + std::to_string(rank) + " has no region");
        return sp->interest_of(it->second, m);
    };
    v.set_sketch_context = [](int32_t comp) { owner_current_sketch_comp = comp; };
    // the interior: the local boxes off the brick's faces, by Morton parity
    const int parity_waves = 1 << tree->dimension;
    v.interior_waves.assign(static_cast<size_t>(parity_waves), {});
    for (int64_t local_idx : level.interior_id) {
        const int64_t m = local_idx + level.local_morton_start;
        v.interior_waves[static_cast<size_t>(m & (parity_waves - 1))].push_back(m);
    }
    return v;
}

// The entry blocks of the owned components' remote boxes, between the
// devices (a device transition kept the level's blocks): the lists of the
// fetch plan, both sides derived from the shared graph.
template<typename CoordType, typename DataType>
CaHaloStats exchange_ca_owner_blocks(ParallelTree<CoordType, DataType>* tree, int lvl,
                                     DeviceLevelBlocks<DataType>& blocks, MPI_Comm level_comm,
                                     const OwnerScheduleState<CoordType, DataType>& st) {
    const auto& level = tree->levels[static_cast<size_t>(lvl)];
    auto rank_of_pid = [&](uint32_t pid) { return level.morton_to_rank.at(static_cast<int>(pid)); };
    const dataflow::FetchPlan plan = dataflow::build_fetch_plan(st.graph, st.my_pid);
    std::map<int, CaHaloPeer> peer_of;
    const auto* graph = &st.graph;
    for (const auto& [color, by_owner] : plan.serve) {
        for (const auto& [owner_pid, boxes] : by_owner) {
            CaHaloPeer& p = peer_of[rank_of_pid(owner_pid)];
            p.rank = rank_of_pid(owner_pid);
            p.requested.insert(p.requested.end(), boxes.begin(), boxes.end());
            const uint32_t pid = owner_pid;
            p.fetches = [graph, pid](int64_t x) {
                const int32_t c = graph->comp_id_of_box(x);
                return c >= 0 && graph->comps[static_cast<size_t>(c)].owner == static_cast<int32_t>(pid) &&
                       graph->grid.proc_of_box(x) != pid;
            };
        }
    }
    for (const auto& [color, by_home] : plan.fetch) {
        for (const auto& [home_pid, boxes] : by_home) {
            CaHaloPeer& p = peer_of[rank_of_pid(home_pid)];
            p.rank = rank_of_pid(home_pid);
            p.mine.insert(p.mine.end(), boxes.begin(), boxes.end());
        }
    }
    for (auto& kv : peer_of) {
        std::sort(kv.second.requested.begin(), kv.second.requested.end());
        std::sort(kv.second.mine.begin(), kv.second.mine.end());
    }
    const uint32_t my_pid = st.my_pid;
    auto i_fetch = [graph, my_pid](int64_t x) {
        const int32_t c = graph->comp_id_of_box(x);
        return c >= 0 && graph->comps[static_cast<size_t>(c)].owner == static_cast<int32_t>(my_pid) &&
               graph->grid.proc_of_box(x) != my_pid;
    };
    return exchange_ca_blocks_for(tree, lvl, blocks, level_comm, peer_of, i_fetch, "CA owner fetch");
}

// Level start: the plan, the fetch of the entry blocks when they are on the
// device, the eliminator.  Returns whether the device runs the level.
template<typename CoordType, typename DataType, typename KernelType>
bool start_owner_level_gpu(ParallelTree<CoordType, DataType>* tree, int lvl,
                           ColorGpuDriver<CoordType, DataType, KernelType>& driver,
                           OwnerScheduleState<CoordType, DataType>& st, OwnerPlanView& view, MPI_Comm level_comm,
                           bool announce) {
    using clock = std::chrono::high_resolution_clock;
    auto& level = tree->levels[static_cast<size_t>(lvl)];
    const auto t0 = clock::now();
    owner_schedule_check_presence(level, st);
    view = build_owner_plan(tree, lvl, st);
    const double t_plan = std::chrono::duration<double>(clock::now() - t0).count();
    if (driver.has_device_blocks()) {
        const CaHaloStats h = exchange_ca_owner_blocks(tree, lvl, *driver.device_blocks(), level_comm, st);
        if (announce) {
            std::printf("  [gpu] level %d owner fetch: %.2f s (plan %.2f, pack %.2f, MPI %.2f, unpack %.2f), "
                        "sent %.2f GB, received %.2f GB in %lld blocks%s\n",
                        lvl, h.total, h.plan, h.pack, h.mpi, h.unpack, h.bytes_sent / 1e9, h.bytes_received / 1e9,
                        static_cast<long long>(h.blocks_received), h.staged ? " (through the host)" : "");
            std::fflush(stdout);
        }
    }
    const auto t1 = clock::now();
    const bool on = driver.start_owner_level(lvl, level_comm, &view, announce);
    if (announce) {
        std::printf("  [gpu] level %d owner start: plan %.2f s, eliminator %.2f s\n", lvl, t_plan,
                    std::chrono::duration<double>(clock::now() - t1).count());
        std::fflush(stdout);
    }
    return on;
}

// The level's elimination on the device.
template<typename CoordType, typename DataType, typename KernelType>
void factorize_owner_level_gpu(ParallelTree<CoordType, DataType>* tree, int lvl,
                               ColorGpuDriver<CoordType, DataType, KernelType>& driver,
                               OwnerScheduleState<CoordType, DataType>& st, const OwnerPlanView& view,
                               int64_t& total_skeleton, int64_t& total_redundant, int64_t& local_max_skel,
                               bool print_detail, int level_print_rank) {
    using clock = std::chrono::high_resolution_clock;
    auto& level = tree->levels[static_cast<size_t>(lvl)];
    const bool print = print_detail && tree->mpi_rank == level_print_rank;
    const bool timeline = env::trace("exchange");
    if (timeline) MPI_Barrier(st.comm);
    const auto t_level = clock::now();
    std::string marks;
    auto mark = [&](const char* what) {
        if (!timeline) return;
        char buf[64];
        std::snprintf(buf, sizeof(buf), " %s %.2f", what, std::chrono::duration<double>(clock::now() - t_level).count());
        marks += buf;
    };
    level.sketch_eliminated_filter = &st.sketch_filter;
    owner_current_sketch_comp = -1;
    auto run_wave = [&](const std::vector<int64_t>& wave, int32_t seq) {
        if (wave.empty()) return;
        driver.eliminate_wave(wave, static_cast<int>(seq));
        for (int64_t m : wave) {
            level.eliminated_boxes.insert(m);
            level.elimination_wave[m] = seq;
        }
    };
    static const char* const names[3] = {"blue/orange", "purple", "green"};
    double t_waves = 0.0, t_exchange = 0.0;
    size_t owned_boxes = 0;
    for (int c = 0; c < view.num_boundary_colors; ++c) {
        const OwnerPlanView::Colour& col = view.colours[static_cast<size_t>(c)];
        size_t boxes = 0, waves = 0;
        for (const auto& w : col.owned_waves) {
            boxes += w.size();
            waves += w.empty() ? 0 : 1;
        }
        owned_boxes += boxes;
        if (print) {
            std::cout << "  Processing owner colour " << (c < 3 ? names[c] : "?") << " (" << boxes << " owned boxes in "
                      << waves << " waves, device)..." << std::endl;
        }
        const auto t0 = clock::now();
        for (int w = 0; w < view.num_waves; ++w) {
            run_wave(col.owned_waves[static_cast<size_t>(w)], static_cast<int32_t>(c * view.num_waves + w));
        }
        const auto t1 = clock::now();
        mark(c == 0 ? "blue-waves" : c == 1 ? "purple-waves" : "green-waves");
        driver.owner_exchange(c);
        const auto t2 = clock::now();
        mark(c == 0 ? "blue-transport" : c == 1 ? "purple-transport" : "green-transport");
        t_waves += std::chrono::duration<double>(t1 - t0).count();
        t_exchange += std::chrono::duration<double>(t2 - t1).count();
        if (print) {
            std::printf("  [gpu] level %d colour %s: waves %.2f s, transport %.2f s\n", lvl, c < 3 ? names[c] : "?",
                        std::chrono::duration<double>(t1 - t0).count(), std::chrono::duration<double>(t2 - t1).count());
            std::fflush(stdout);
        }
    }
    {
        size_t boxes = 0;
        for (const auto& w : view.interior_waves) boxes += w.size();
        if (print) std::cout << "  Processing CA interior group (" << boxes << " boxes, device)..." << std::endl;
        const auto t0 = clock::now();
        for (size_t w = 0; w < view.interior_waves.size(); ++w) {
            run_wave(view.interior_waves[w], view.interior_seq0 + static_cast<int32_t>(w));
        }
        t_waves += std::chrono::duration<double>(clock::now() - t0).count();
        mark("interior");
    }
    if (timeline) {
        std::printf("[owner] rank %d level %d timeline (s since the level's start):%s\n", tree->mpi_rank, lvl, marks.c_str());
        std::fflush(stdout);
    }
    owner_current_sketch_comp = -1;
    level.sketch_eliminated_filter = nullptr;
    // the solves run the level in this order (the level's cleanup drops elimination_wave)
    level.ca_home_order.clear();
    for (const auto& box : level.local_boxes) {
        const auto it = level.elimination_wave.find(box.morton_index);
        if (it == level.elimination_wave.end()) {
            throw std::runtime_error("owner level: box " + std::to_string(box.morton_index) + " was never eliminated");
        }
        level.ca_home_order[box.morton_index] = it->second;
    }
    if (print) {
        std::printf("  [gpu] level %d owner schedule: waves %.2f s, colour transports %.2f s (owned boxes %zu in %zu "
                    "components, installed components %zu)\n",
                    lvl, t_waves, t_exchange, owned_boxes, st.rt.owned().size(), st.rt.shared_not_owned().size());
        std::fflush(stdout);
    }
    for (const auto& box : level.local_boxes) {
        total_skeleton += static_cast<int64_t>(box.skeleton_indices.size());
        total_redundant += static_cast<int64_t>(box.redundant_indices.size());
        local_max_skel = std::max<int64_t>(local_max_skel, box.skeleton_indices.size());
    }
}

}  // namespace gpu
}  // namespace fmm
#endif  // H2_HAVE_GPU
