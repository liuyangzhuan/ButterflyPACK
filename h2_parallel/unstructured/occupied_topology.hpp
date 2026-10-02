#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace butterfly {
namespace color_unstructured {

struct OccupiedTopology {
    // Sorted Morton indices. Storage is proportional to occupied boxes rather
    // than the volume of the embedding grid.
    std::vector<std::vector<int64_t>> occupied;
    std::vector<int64_t> occupied_count;

    bool contains(int level, int64_t morton) const {
        if (level < 0 || level >= static_cast<int>(occupied.size()) ||
            morton < 0) {
            return false;
        }
        const auto& mortons = occupied[static_cast<size_t>(level)];
        return std::binary_search(mortons.begin(), mortons.end(), morton);
    }
};

template<typename CoordType, typename DataType>
OccupiedTopology discover_occupied_topology(
    const fmm::ParallelTree<CoordType, DataType>* tree) {
    if (tree == nullptr || tree->num_levels <= 0) {
        throw std::invalid_argument(
            "color_unstructured: cannot discover occupancy from an empty tree");
    }

    OccupiedTopology topology;
    topology.occupied.resize(static_cast<size_t>(tree->num_levels));
    topology.occupied_count.assign(static_cast<size_t>(tree->num_levels), 0);

    const int leaf_level = tree->num_levels - 1;
    const int64_t global_leaf_boxes =
        tree->levels[leaf_level].num_boxes_global;
    auto& leaf_mortons = topology.occupied[static_cast<size_t>(leaf_level)];

    std::vector<int64_t> local_leaf_mortons;
    int64_t local_points = 0;
    for (const auto& box : tree->levels[leaf_level].local_boxes) {
        if (box.num_points <= 0) continue;
        if (box.morton_index < 0 ||
            box.morton_index >= global_leaf_boxes) {
            throw std::runtime_error(
                "color_unstructured: leaf Morton index is out of range");
        }
        local_leaf_mortons.push_back(box.morton_index);
        local_points += box.num_points;
    }

    int64_t global_points = local_points;
    MPI_Allreduce(
        MPI_IN_PLACE, &global_points, 1, MPI_INT64_T, MPI_SUM, tree->comm);
    if (global_points != tree->num_points) {
        throw std::runtime_error(
            "color_unstructured: occupied leaves contain " +
            std::to_string(global_points) + " points, expected " +
            std::to_string(tree->num_points));
    }

    if (local_leaf_mortons.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            "color_unstructured: local occupied-leaf count exceeds MPI int range");
    }

    int comm_size = 1;
    MPI_Comm_size(tree->comm, &comm_size);
    const int local_count = static_cast<int>(local_leaf_mortons.size());
    std::vector<int> counts(static_cast<size_t>(comm_size), 0);
    MPI_Allgather(
        &local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, tree->comm);

    std::vector<int> displacements(static_cast<size_t>(comm_size), 0);
    int64_t total_count = 0;
    for (int rank = 0; rank < comm_size; ++rank) {
        if (counts[static_cast<size_t>(rank)] < 0 ||
            total_count >
                static_cast<int64_t>(std::numeric_limits<int>::max()) -
                    counts[static_cast<size_t>(rank)]) {
            throw std::overflow_error(
                "color_unstructured: global occupied-leaf count exceeds MPI int range");
        }
        displacements[static_cast<size_t>(rank)] =
            static_cast<int>(total_count);
        total_count += counts[static_cast<size_t>(rank)];
    }

    leaf_mortons.resize(static_cast<size_t>(total_count));
    MPI_Allgatherv(
        local_leaf_mortons.data(), local_count, MPI_INT64_T,
        leaf_mortons.data(), counts.data(), displacements.data(),
        MPI_INT64_T, tree->comm);
    std::sort(leaf_mortons.begin(), leaf_mortons.end());
    leaf_mortons.erase(
        std::unique(leaf_mortons.begin(), leaf_mortons.end()),
        leaf_mortons.end());

    const int children_per_parent =
        morton::children_per_box(tree->dimension);
    for (int level = leaf_level; level > 0; --level) {
        const auto& children = topology.occupied[static_cast<size_t>(level)];
        auto& parents = topology.occupied[static_cast<size_t>(level - 1)];
        parents.reserve(children.size());
        for (int64_t child : children) {
            parents.push_back(child / children_per_parent);
        }
        parents.erase(
            std::unique(parents.begin(), parents.end()), parents.end());
    }

    for (int level = 0; level < tree->num_levels; ++level) {
        topology.occupied_count[static_cast<size_t>(level)] =
            static_cast<int64_t>(
                topology.occupied[static_cast<size_t>(level)].size());
    }
    if (topology.occupied_count.front() != 1) {
        throw std::runtime_error(
            "color_unstructured: nonempty problem must have one occupied root");
    }
    return topology;
}

template<typename CoordType, typename DataType>
std::vector<int64_t> occupied_local_indices(
    const fmm::TreeLevel<CoordType, DataType>& level) {
    std::vector<int64_t> indices;
    indices.reserve(level.boundary_id.size() + level.interior_id.size());
    indices.insert(
        indices.end(), level.boundary_id.begin(), level.boundary_id.end());
    indices.insert(
        indices.end(), level.interior_id.begin(), level.interior_id.end());
    return indices;
}

template<typename CoordType, typename DataType>
std::vector<int64_t> occupied_local_indices_from_topology(
    const fmm::TreeLevel<CoordType, DataType>& level,
    int level_index,
    const OccupiedTopology& topology) {
    std::vector<int64_t> indices;
    if (!level.is_process_active || level.local_boxes.empty()) return indices;

    const int64_t first = level.local_boxes.front().morton_index;
    const int64_t last = level.local_boxes.back().morton_index;
    const auto& occupied = topology.occupied.at(
        static_cast<size_t>(level_index));
    auto begin = std::lower_bound(occupied.begin(), occupied.end(), first);
    auto end = std::upper_bound(begin, occupied.end(), last);
    indices.reserve(static_cast<size_t>(std::distance(begin, end)));
    for (auto it = begin; it != end; ++it) {
        const int64_t local_index = *it - first;
        if (local_index < 0 ||
            local_index >= static_cast<int64_t>(level.local_boxes.size()) ||
            level.local_boxes[static_cast<size_t>(local_index)].morton_index !=
                *it) {
            throw std::runtime_error(
                "color_unstructured: occupied Morton index is outside the "
                "local contiguous box slab");
        }
        indices.push_back(local_index);
    }
    return indices;
}

template<typename CoordType, typename DataType>
void compute_occupied_neighbor_lists(
    fmm::TreeLevel<CoordType, DataType>& level,
    int level_index,
    int dimension,
    const OccupiedTopology& topology) {
    const uint32_t grid_size = 1u << level_index;
    const auto occupied_indices = occupied_local_indices_from_topology(
        level, level_index, topology);

    for (int64_t local_index : occupied_indices) {
        auto& box = level.local_boxes[static_cast<size_t>(local_index)];

        box.one_hop.clear();
        for (uint64_t neighbor : morton::neighbors_nd(
                 dimension, box.morton_index, grid_size)) {
            const int64_t neighbor_morton = static_cast<int64_t>(neighbor);
            if (topology.contains(level_index, neighbor_morton)) {
                box.one_hop.push_back(neighbor_morton);
            }
        }
        box.use_full_set.assign(box.one_hop.size(), 1);

        box.two_hop.clear();
        for (uint64_t neighbor : morton::neighbors_2hop_nd(
                 dimension, box.morton_index, grid_size)) {
            const int64_t neighbor_morton = static_cast<int64_t>(neighbor);
            if (topology.contains(level_index, neighbor_morton)) {
                box.two_hop.push_back(neighbor_morton);
            }
        }
    }
}

template<typename CoordType, typename DataType>
std::vector<int64_t> occupied_box_indices_from_topology(
    const std::vector<fmm::BoxData<CoordType, DataType>>& boxes,
    int level_index,
    const OccupiedTopology& topology) {
    std::vector<int64_t> indices;
    if (boxes.empty()) return indices;

    const int64_t first = boxes.front().morton_index;
    const int64_t last = boxes.back().morton_index;
    const auto& occupied = topology.occupied.at(
        static_cast<size_t>(level_index));
    auto begin = std::lower_bound(occupied.begin(), occupied.end(), first);
    auto end = std::upper_bound(begin, occupied.end(), last);
    indices.reserve(static_cast<size_t>(std::distance(begin, end)));
    for (auto it = begin; it != end; ++it) {
        const int64_t local_index = *it - first;
        if (local_index < 0 ||
            local_index >= static_cast<int64_t>(boxes.size()) ||
            boxes[static_cast<size_t>(local_index)].morton_index != *it) {
            throw std::runtime_error(
                "color_unstructured: occupied Morton index is outside the "
                "contiguous box slab");
        }
        indices.push_back(local_index);
    }
    return indices;
}

template<typename CoordType, typename DataType>
void compute_occupied_neighbor_lists(
    std::vector<fmm::BoxData<CoordType, DataType>>& boxes,
    int level_index,
    int dimension,
    const OccupiedTopology& topology) {
    const uint32_t grid_size = 1u << level_index;
    const auto occupied_indices = occupied_box_indices_from_topology(
        boxes, level_index, topology);

    for (int64_t local_index : occupied_indices) {
        auto& box = boxes[static_cast<size_t>(local_index)];
        box.one_hop.clear();
        for (uint64_t neighbor : morton::neighbors_nd(
                 dimension, box.morton_index, grid_size)) {
            const int64_t neighbor_morton = static_cast<int64_t>(neighbor);
            if (topology.contains(level_index, neighbor_morton)) {
                box.one_hop.push_back(neighbor_morton);
            }
        }
        box.use_full_set.assign(box.one_hop.size(), 1);

        box.two_hop.clear();
        for (uint64_t neighbor : morton::neighbors_2hop_nd(
                 dimension, box.morton_index, grid_size)) {
            const int64_t neighbor_morton = static_cast<int64_t>(neighbor);
            if (topology.contains(level_index, neighbor_morton)) {
                box.two_hop.push_back(neighbor_morton);
            }
        }
    }
}

template<typename CoordType, typename DataType>
void filter_box_neighbors(fmm::BoxData<CoordType, DataType>& box,
                          const OccupiedTopology& topology) {
    if (!topology.contains(box.level, box.morton_index)) {
        box.one_hop.clear();
        box.two_hop.clear();
        box.use_full_set.clear();
        return;
    }

    std::vector<int64_t> one_hop;
    std::vector<int64_t> use_full_set;
    one_hop.reserve(box.one_hop.size());
    use_full_set.reserve(box.one_hop.size());
    for (size_t index = 0; index < box.one_hop.size(); ++index) {
        const int64_t neighbor = box.one_hop[index];
        if (!topology.contains(box.level, neighbor)) continue;
        one_hop.push_back(neighbor);
        use_full_set.push_back(
            index < box.use_full_set.size() ? box.use_full_set[index] : 1);
    }
    box.one_hop.swap(one_hop);
    box.use_full_set.swap(use_full_set);

    box.two_hop.erase(
        std::remove_if(
            box.two_hop.begin(), box.two_hop.end(),
            [&](int64_t neighbor) {
                return !topology.contains(box.level, neighbor);
            }),
        box.two_hop.end());
}

template<typename CoordType, typename DataType>
void filter_level_topology(fmm::ParallelTree<CoordType, DataType>* tree,
                           int level_index,
                           const OccupiedTopology& topology) {
    auto& level = tree->levels[static_cast<size_t>(level_index)];

    level.boundary_id.clear();
    level.interior_id.clear();
    const auto occupied_indices = occupied_local_indices_from_topology(
        level, level_index, topology);
    for (int64_t local_index : occupied_indices) {
        auto& box = level.local_boxes[static_cast<size_t>(local_index)];
        filter_box_neighbors(box, topology);
        if (box.on_boundary) {
            level.boundary_id.push_back(local_index);
        } else {
            level.interior_id.push_back(local_index);
        }
    }

    std::vector<int64_t> assisting_mortons;
    assisting_mortons.reserve(
        level.assisting_box_points_for_kernel_evaluation.size());
    for (const auto& entry :
         level.assisting_box_points_for_kernel_evaluation) {
        if (topology.contains(level_index, entry.first)) {
            assisting_mortons.push_back(entry.first);
        }
    }
    std::sort(assisting_mortons.begin(), assisting_mortons.end());
    level.assisting_box_points_for_kernel_evaluation.clear();
    for (size_t index = 0; index < assisting_mortons.size(); ++index) {
        level.assisting_box_points_for_kernel_evaluation[
            assisting_mortons[index]] = static_cast<int64_t>(index);
    }
    level.assisting_boxes.clear();
    level.assisting_boxes.resize(assisting_mortons.size());

    std::vector<std::pair<int64_t, int>> solve_mortons;
    solve_mortons.reserve(
        level.ghost_and_assisting_box_points_for_solve_map.size());
    for (const auto& entry :
         level.ghost_and_assisting_box_points_for_solve_map) {
        if (!topology.contains(level_index, entry.first)) continue;
        int is_ghost = 0;
        if (entry.second >= 0 &&
            entry.second < static_cast<int64_t>(level.is_ghost_solve.size())) {
            is_ghost = level.is_ghost_solve[static_cast<size_t>(entry.second)];
        }
        solve_mortons.emplace_back(entry.first, is_ghost);
    }
    std::sort(solve_mortons.begin(), solve_mortons.end());
    level.ghost_and_assisting_box_points_for_solve_map.clear();
    level.is_ghost_solve.clear();
    for (size_t index = 0; index < solve_mortons.size(); ++index) {
        level.ghost_and_assisting_box_points_for_solve_map[
            solve_mortons[index].first] = static_cast<int64_t>(index);
        level.is_ghost_solve.push_back(solve_mortons[index].second);
    }
    level.ghost_and_assisting_boxes_for_solve.clear();
}

template<typename CoordType, typename DataType>
void add_empty_parent_transition_assisting_slots(
    fmm::ParallelTree<CoordType, DataType>* tree,
    int child_level_index,
    const OccupiedTopology& topology) {
    auto& child_level =
        tree->levels[static_cast<size_t>(child_level_index)];
    if (!child_level.is_process_active || child_level.local_boxes.empty() ||
        child_level_index <= 0) {
        return;
    }

    const int children_per_parent =
        morton::children_per_box(tree->dimension);
    const int parent_level_index = child_level_index - 1;
    const uint32_t parent_grid_size = 1u << parent_level_index;
    const int64_t local_parent_start =
        child_level.local_boxes.front().morton_index / children_per_parent;
    const int64_t local_parent_end =
        child_level.local_boxes.back().morton_index / children_per_parent;

    std::vector<int64_t> empty_children;
    const auto& occupied_parents =
        topology.occupied[static_cast<size_t>(parent_level_index)];
    auto parent_begin = std::lower_bound(
        occupied_parents.begin(), occupied_parents.end(), local_parent_start);
    auto parent_end = std::upper_bound(
        parent_begin, occupied_parents.end(), local_parent_end);
    for (auto parent_it = parent_begin;
         parent_it != parent_end; ++parent_it) {
        const int64_t parent = *parent_it;
        const auto neighbors = morton::neighbors_nd(
            tree->dimension, parent, parent_grid_size);
        for (uint64_t neighbor_value : neighbors) {
            const int64_t neighbor_parent =
                static_cast<int64_t>(neighbor_value);
            if (neighbor_parent >= local_parent_start &&
                neighbor_parent <= local_parent_end) {
                continue;
            }
            const int64_t first_child =
                neighbor_parent * children_per_parent;
            for (int child = 0; child < children_per_parent; ++child) {
                const int64_t child_morton = first_child + child;
                if (child_level
                        .assisting_box_points_for_kernel_evaluation
                        .count(child_morton) != 0) {
                    continue;
                }
                if (topology.contains(child_level_index, child_morton)) {
                    throw std::runtime_error(
                        "color_unstructured: occupied remote child " +
                        std::to_string(child_morton) +
                        " is missing before parent assembly");
                }
                empty_children.push_back(child_morton);
            }
        }
    }

    std::sort(empty_children.begin(), empty_children.end());
    empty_children.erase(
        std::unique(empty_children.begin(), empty_children.end()),
        empty_children.end());
    for (int64_t child_morton : empty_children) {
        const int64_t slot =
            static_cast<int64_t>(child_level.assisting_boxes.size());
        fmm::PointDataRequest<CoordType> empty;
        empty.morton_index = child_morton;
        child_level.assisting_boxes.push_back(std::move(empty));
        child_level.assisting_box_points_for_kernel_evaluation[
            child_morton] = slot;
    }
}

template<typename CoordType, typename DataType>
OccupiedTopology prepare_occupied_topology(
    fmm::ParallelTree<CoordType, DataType>* tree,
    int verbosity) {
    OccupiedTopology topology = discover_occupied_topology(tree);
    for (int level = 0; level < tree->num_levels; ++level) {
        filter_level_topology(tree, level, topology);
    }

    if (verbosity >= 0 && tree->mpi_rank == 0) {
        std::cout << "color_unstructured occupied boxes:" << std::endl;
        for (int level = 0; level < tree->num_levels; ++level) {
            const int64_t full = tree->levels[level].num_boxes_global;
            const int64_t occupied =
                topology.occupied_count[static_cast<size_t>(level)];
            std::cout << "  level " << level << ": " << occupied
                      << " / " << full << std::endl;
        }
    }
    return topology;
}

}  // namespace color_unstructured
}  // namespace butterfly
