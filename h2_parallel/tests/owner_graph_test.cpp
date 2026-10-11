// Test and timing of the component-owner schedule's graph (core/dataflow.hpp;
// h2_gpu/GPU_CA_PLAN.md, M5e): the local construction (every rank from its
// own neighbourhood plus one allgather of part descriptors) against the
// global one (the whole grid's flood fill), for 16^3 bricks on the 2^3 and
// 4^3 process grids (Laplace 160^3 on 8 ranks, 320^3 on 64: every brick's
// local graph is compared) and timed on the 8^3 grid (640^3 on 512 ranks).
//
//   make h2_owner_graph_test && ./h2_owner_graph_test     (in a build directory), or
//   g++ -std=c++17 -O2 -I../core owner_graph_test.cpp -o owner_graph_test && ./owner_graph_test
//
// (no MPI calls: the allgather is the concatenation of every brick's stream)

#include "dataflow.hpp"
#include "halo_geometry.hpp"

#include <chrono>
#include <cstdio>
#include <string>

using namespace fmm;
using namespace fmm::dataflow;

static int failures = 0;

static void fail(const std::string& what, int level, uint32_t procs, uint32_t pid) {
    ++failures;
    std::printf("FAIL level %d, %u bricks, brick %u: %s\n", level, procs, pid, what.c_str());
}

static double seconds_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

static ColorOf classifier(const ProcessGrid& grid, int32_t level) {
    return make_ca_classifier(
        grid, level,
        [](int64_t m, int32_t lvl, int32_t dim, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_blue_box(m, lvl, dim, gs, lgs, off);
        },
        [](int64_t m, int32_t lvl, int32_t dim, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_orange_box(m, lvl, dim, gs, lgs, off);
        },
        [](int64_t m, int32_t lvl, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_purple_box(m, lvl, gs, lgs, off);
        });
}

template <typename T>
static bool same(const std::vector<T>& a, const std::vector<T>& b) { return a == b; }

// the local graph of brick `pid` against the global one
static void compare(const ComponentGraph& G, const ComponentGraph& L, int level, uint32_t procs, uint32_t pid) {
    if (G.comps.size() != L.comps.size()) {
        fail("component count " + std::to_string(L.comps.size()) + " against " + std::to_string(G.comps.size()),
             level, procs, pid);
        return;
    }
    for (size_t i = 0; i < G.comps.size(); ++i) {
        const Component& g = G.comps[i];
        const Component& l = L.comps[i];
        const std::string who = "component " + std::to_string(i);
        if (l.id != g.id) fail(who + ": id", level, procs, pid);
        if (l.color != g.color) fail(who + ": colour", level, procs, pid);
        if (l.first_morton != g.first_morton) fail(who + ": first Morton index", level, procs, pid);
        if (l.cost != g.cost) fail(who + ": cost", level, procs, pid);
        if (l.owner != g.owner) fail(who + ": owner " + std::to_string(l.owner) + " against " + std::to_string(g.owner), level, procs, pid);
        if (l.upward_rank != g.upward_rank) fail(who + ": upward rank", level, procs, pid);
        if (!same(l.sharers, g.sharers)) fail(who + ": sharers", level, procs, pid);
        if (!same(l.adjacent, g.adjacent)) fail(who + ": adjacent", level, procs, pid);
        if (!same(l.visible, g.visible)) fail(who + ": visible", level, procs, pid);
        if (!same(l.preds, g.preds)) fail(who + ": preds", level, procs, pid);
        if (!same(l.succs, g.succs)) fail(who + ": succs", level, procs, pid);
        if (!same(l.preds2, g.preds2)) fail(who + ": preds2", level, procs, pid);
        if (!same(l.succs2, g.succs2)) fail(who + ": succs2", level, procs, pid);
        // the boxes: complete for the components touching the brick's
        // neighbourhood entirely, a subset otherwise
        const bool touches = g.has_sharer(pid) || g.is_adjacent(pid);
        if (touches && !l.boxes_complete) fail(who + ": touches the brick but its boxes are not complete", level, procs, pid);
        if (l.boxes_complete) {
            if (!same(l.boxes, g.boxes)) fail(who + ": boxes", level, procs, pid);
        } else {
            for (int64_t m : l.boxes)
                if (!std::binary_search(g.boxes.begin(), g.boxes.end(), m)) fail(who + ": a box not in the component", level, procs, pid);
        }
    }
    // every boundary box of the neighbourhood, with the same component
    for (uint32_t q : L.local_bricks)
        for (int64_t m : G.grid.brick_boundary_boxes(q))
            if (L.comp_id_of_box(m) != G.comp_id_of_box(m)) { fail("box " + std::to_string(m) + " of brick " + std::to_string(q), level, procs, pid); break; }
    if (L.box_to_comp.size() != [&] { size_t n = 0; for (uint32_t q : L.local_bricks) n += G.grid.brick_boundary_boxes(q).size(); return n; }())
        fail("box_to_comp size", level, procs, pid);
    if (!same(L.interior[pid].preds, G.interior[pid].preds)) fail("interior node", level, procs, pid);
    // what the schedule derives from the graph
    if (!same(auxiliary_components(L, pid), auxiliary_components(G, pid))) fail("auxiliary components", level, procs, pid);
    const FetchPlan fl = build_fetch_plan(L, pid), fg = build_fetch_plan(G, pid);
    if (fl.fetch != fg.fetch || fl.serve != fg.serve || fl.foreign_eliminated != fg.foreign_eliminated) fail("fetch plan", level, procs, pid);
    for (const auto& k : G.comps) {
        if (k.owner != static_cast<int32_t>(pid)) continue;
        for (int64_t m : k.boxes)
            if (payload_recipients(L, pid, m) != payload_recipients(G, pid, m)) { fail("payload recipients of box " + std::to_string(m), level, procs, pid); break; }
    }
    const std::string e = validate_graph(L);
    if (!e.empty()) fail("validate_graph: " + e, level, procs, pid);
    if (assignment_hash(L) != assignment_hash(G)) fail("assignment hash", level, procs, pid);
}

static void run(int32_t level, uint32_t procs, bool compare_all) {
    const ProcessGrid grid(3, 1u << level, procs);
    const ColorOf color_of = classifier(grid, level);
    auto t0 = std::chrono::steady_clock::now();
    ComponentGraph G = build_component_graph(grid, color_of);
    const double t_global = seconds_since(t0);
    // the owner assignment: the indexed polish against the reference one
    {
        ComponentGraph Gr = G;
        t0 = std::chrono::steady_clock::now();
        assign_owners(Gr, AssignPolicy::LPT_PER_COLOR, nullptr, true);
        const double t_ref = seconds_since(t0);
        ComponentGraph Gc = G;
        assign_owners(Gc, AssignPolicy::LPT_CUMULATIVE, nullptr, true);
        t0 = std::chrono::steady_clock::now();
        assign_owners(G, AssignPolicy::LPT_PER_COLOR);
        const double t_new = seconds_since(t0);
        ComponentGraph Gc2 = Gr;
        assign_owners(Gc2, AssignPolicy::LPT_CUMULATIVE);
        for (size_t i = 0; i < G.comps.size(); ++i) {
            if (G.comps[i].owner != Gr.comps[i].owner) fail("indexed polish differs from the reference (per colour) at component " + std::to_string(i), level, procs, 0);
            if (Gc2.comps[i].owner != Gc.comps[i].owner) fail("indexed polish differs from the reference (cumulative) at component " + std::to_string(i), level, procs, 0);
        }
        std::printf("  assign_owners: reference polish %.3f s, indexed %.3f s\n", t_ref, t_new);
    }
    // every brick's parts (what each rank computes), their streams concatenated
    std::vector<LocalParts> parts(procs);
    std::vector<int64_t> all;
    double t_parts_max = 0.0;
    for (uint32_t pid = 0; pid < procs; ++pid) {
        t0 = std::chrono::steady_clock::now();
        parts[pid] = local_component_parts(grid, color_of, pid);
        t_parts_max = std::max(t_parts_max, seconds_since(t0));
        all.insert(all.end(), parts[pid].stream.begin(), parts[pid].stream.end());
    }
    double t_assemble_max = 0.0;
    size_t compared = 0;
    for (uint32_t pid = 0; pid < procs; ++pid) {
        if (!compare_all && pid != 0 && pid != procs / 2 && pid != procs - 1) continue;
        t0 = std::chrono::steady_clock::now();
        ComponentGraph L = assemble_component_graph(grid, pid, parts[pid], all);
        assign_owners(L, AssignPolicy::LPT_PER_COLOR);
        t_assemble_max = std::max(t_assemble_max, seconds_since(t0));
        compare(G, L, level, procs, pid);
        {   // the cumulative policy (the CPU's fine levels: its tie-break reads the boxes at home) is the same on the local graph
            ComponentGraph Lc = assemble_component_graph(grid, pid, parts[pid], all);
            assign_owners(Lc, AssignPolicy::LPT_CUMULATIVE);
            ComponentGraph Gc = G;
            assign_owners(Gc, AssignPolicy::LPT_CUMULATIVE);
            for (size_t i = 0; i < G.comps.size(); ++i)
                if (Lc.comps[i].owner != Gc.comps[i].owner) { fail("cumulative owner differs on the local graph at component " + std::to_string(i), level, procs, pid); break; }
        }
        ++compared;
    }
    std::printf("level %d, %u bricks of %u^3: %zu components | global build %.3f s | local: parts %.3f s, assemble+assign %.3f s "
                "(max over bricks), stream %zu int64 per brick, %.1f MB gathered | %zu bricks compared, %d failures so far\n",
                level, procs, grid.brick, G.comps.size(), t_global, t_parts_max, t_assemble_max,
                parts[0].stream.size(), all.size() * 8.0 / 1e6, compared, failures);
    std::fflush(stdout);
}

int main() {
    run(5, 8, true);
    run(6, 64, true);
    run(7, 512, false);
    // 2D, and 8^3 bricks (the 160^3 case's level 4)
    {
        const ProcessGrid grid(2, 1u << 6, 16);
        const ColorOf color_of = classifier(grid, 6);
        ComponentGraph G = build_component_graph(grid, color_of);
        assign_owners(G, AssignPolicy::LPT_PER_COLOR);
        std::vector<LocalParts> parts(16);
        std::vector<int64_t> all;
        for (uint32_t pid = 0; pid < 16; ++pid) {
            parts[pid] = local_component_parts(grid, color_of, pid);
            all.insert(all.end(), parts[pid].stream.begin(), parts[pid].stream.end());
        }
        for (uint32_t pid = 0; pid < 16; ++pid) {
            ComponentGraph L = assemble_component_graph(grid, pid, parts[pid], all);
            assign_owners(L, AssignPolicy::LPT_PER_COLOR);
            compare(G, L, 6, 16, pid);
        }
        std::printf("2D level 6, 16 bricks of %u^2: %zu components, %d failures so far\n", grid.brick, G.comps.size(), failures);
    }
    run(4, 8, true);
    std::printf(failures == 0 ? "OK\n" : "%d FAILURES\n", failures);
    return failures == 0 ? 0 : 1;
}
