#pragma once
// A component-owner CA level (H2_CA_owner_component=3) as the device
// eliminator sees it (h2_gpu/GPU_CA_PLAN.md, M5).  ca_owner.hpp builds it
// from the host schedule (core/owner_schedule.hpp: the component graph, the
// owner assignment, the component-aware wave partition, the payload
// recipients); level_eliminator.hpp reads it and knows nothing of the graph.
//
// Box kinds on an owner level, besides the local boxes eliminated here as on
// any level: owned ghosts (the owned components' boxes of other ranks,
// eliminated here, their state shipped home), foreign local boxes
// (eliminated by another rank, their FULL state lands here), installed
// copies (ghosts eliminated elsewhere whose generators arrive as COMPACT
// payloads), and skeleton-only boxes (two hops away: index sets).
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fmm {
namespace gpu {

enum class OwnerPayload : int8_t { FULL = 0, COMPACT = 1, SKELETON = 2 };

struct OwnerPlanView {
    int num_boundary_colors = 3;
    int num_waves = 8;  // elimination classes per colour (the component partition's, <= 8)
    // ghost boxes eliminated here: the owned components' boxes of other ranks
    std::unordered_set<int64_t> owned_ghosts;
    // local boxes eliminated by another rank (their FULL state arrives)
    std::unordered_set<int64_t> foreign_local;
    // boundary box -> its elimination sequence number (colour * num_waves +
    // wave; the same on every rank) and its component
    std::unordered_map<int64_t, int32_t> seq_of_box;
    std::unordered_map<int64_t, int32_t> comp_of_box;
    // owned box -> the tree ranks its state goes to after its elimination
    std::unordered_map<int64_t, std::vector<std::pair<int, OwnerPayload>>> recipients;
    // whether tree rank `rank` holds a pair of interest with box m: its own
    // box, or a box of a component it owns (what a COMPACT payload's X_NR
    // rows cover)
    std::function<bool(int, int64_t)> interest_of;
    // the sketch's view of eliminated boxes while a box is planned: the
    // box's component (-1: an interior box)
    std::function<void(int32_t)> set_sketch_context;
    struct Colour {
        std::vector<int> send_peers, recv_peers;         // tree ranks, sorted
        std::vector<std::vector<int64_t>> owned_waves;  // by wave: the owned boxes (local and ghosts)
    };
    std::vector<Colour> colours;
    std::vector<std::vector<int64_t>> interior_waves;  // the local interior boxes, by Morton parity
    int32_t interior_seq0 = 0;                         // num_boundary_colors * num_waves

    bool owned_ghost(int64_t m) const { return owned_ghosts.count(m) != 0; }
    bool foreign(int64_t m) const { return foreign_local.count(m) != 0; }
    int32_t comp_of(int64_t m) const {
        const auto it = comp_of_box.find(m);
        return it == comp_of_box.end() ? -1 : it->second;
    }
};

}  // namespace gpu
}  // namespace fmm
