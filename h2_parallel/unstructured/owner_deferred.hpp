#pragma once

#include "core/butterfly_types.hpp"

namespace butterfly {
namespace color_unstructured {
using namespace fmm;

template<typename CoordType, typename DataType, typename KernelType>
void apply_unstructured_owner_deferred_xnn_updates_for_candidate_box(
    int64_t candidate_morton,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    const std::unordered_set<int64_t>& wave_box_set,
    DeferredXnnOwnerScratch<DataType>& scratch,
    std::vector<DeferredXnnTargetKey>& mirror_targets,
    PendingFactorUpdates<DataType>* pending,
    bool include_ghosts = false,
    const DeferredPairFilter* pair_filter = nullptr) {
    mirror_targets.clear();
    // (entries through a collective callback: the remote additions held in
    // the candidate until its requests are all fetched, so that an abandoned
    // candidate changes nothing outside itself)
    PendingFactorUpdates<DataType> candidate_adds;
    PendingFactorUpdates<DataType>* add_target =
        (kernel_entries(kernel, 0) != nullptr && pending != nullptr) ? &candidate_adds : pending;

    DeferredXnnEndpoint<CoordType, DataType> candidate_endpoint =
        resolve_deferred_xnn_endpoint(level, candidate_morton);
    BoxData<CoordType, DataType>* candidate_box = candidate_endpoint.box;
    const bool candidate_is_local = candidate_endpoint.is_local;
    const bool candidate_is_writable =
        candidate_is_local || (include_ghosts && candidate_endpoint.is_ghost);
    const bool candidate_is_assisting = candidate_endpoint.is_assisting;

    if (!candidate_is_writable && !candidate_is_assisting) {
        return;
    }
    const bool generated_near = generator_near_enabled() &&
        lazy_far_field_mode() == LazyFarFieldMode::LAZY;
    if (candidate_is_assisting && generated_near) {
        return;
    }
    if (candidate_is_writable &&
        (candidate_box == nullptr || candidate_box->one_hop.empty() ||
         deferred_xnn_should_skip_owner_candidate_box(
             candidate_box, level, include_ghosts))) {
        return;
    }

    const int dimension = level.dimension;
    scratch.accumulated_targets.clear();
    scratch.accumulated_target_indices.clear();
    scratch.preallocated_mirror_targets.clear();

    const bool candidate_unarrived =
        staged_halo_box_unarrived(level, candidate_morton);
    std::vector<typename TreeLevel<CoordType, DataType>::StagedPendingDelta>
        local_pending;
    std::vector<std::pair<int64_t, int64_t>> local_mirrors;

    auto elimination_sequence = [&](int64_t morton) -> int32_t {
        if (level.eliminated_boxes.find(morton) ==
            level.eliminated_boxes.end()) {
            return std::numeric_limits<int32_t>::min();
        }
        const auto it = level.elimination_wave.find(morton);
        return it == level.elimination_wave.end()
            ? std::numeric_limits<int32_t>::min()
            : it->second;
    };

    std::vector<int64_t> candidate_sources;
    if (candidate_is_writable) {
        for (int64_t source_morton : candidate_box->one_hop) {
            if (wave_box_set.find(source_morton) != wave_box_set.end()) {
                candidate_sources.push_back(source_morton);
            }
        }
    } else {
        for (int64_t source_morton : wave_box_set) {
            BoxData<CoordType, DataType>* source_box =
                resolve_deferred_xnn_box(level, source_morton);
            if (source_box == nullptr || source_box->deferred_xnn_temp2.empty()) {
                continue;
            }
            if (std::find(
                    source_box->one_hop.begin(),
                    source_box->one_hop.end(),
                    candidate_morton) != source_box->one_hop.end()) {
                candidate_sources.push_back(source_morton);
            }
        }
    }

    if (candidate_sources.size() > 1) {
        std::stable_sort(
            candidate_sources.begin(), candidate_sources.end(),
            [&](int64_t first, int64_t second) {
                return elimination_sequence(first) <
                    elimination_sequence(second);
            });
    }
    const int32_t candidate_sequence =
        elimination_sequence(candidate_morton);

    for (int64_t source_morton : candidate_sources) {
        // The source keeps both deferred temp2 and the original X_NR until this
        // post-wave pass consumes them.
        BoxData<CoordType, DataType>* source_box =
            resolve_deferred_xnn_box(level, source_morton);
        if (source_box == nullptr || source_box->deferred_xnn_temp2.empty()) {
            continue;
        }
        const int32_t source_sequence =
            elimination_sequence(source_morton);
        if (candidate_sequence > source_sequence) {
            continue;
        }
        if (!source_box->X_NR.is_allocated()) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: source X_NR missing");
        }

        const auto& source_neighbor_counts = source_box->deferred_xnn_neighbor_point_counts;
        if (source_neighbor_counts.size() != source_box->one_hop.size()) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: neighbor count size mismatch");
        }

        size_t source_candidate_idx = source_box->one_hop.size();
        int64_t candidate_col_offset = 0;
        for (size_t idx = 0; idx < source_box->one_hop.size(); ++idx) {
            if (source_box->one_hop[idx] == candidate_morton) {
                source_candidate_idx = idx;
                break;
            }
            // Source X_NR stores neighbor blocks contiguously in one_hop order.
            candidate_col_offset += source_neighbor_counts[idx];
        }

        if (source_candidate_idx == source_box->one_hop.size()) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: candidate not found in source one_hop");
        }

        const int64_t n_candidate = source_neighbor_counts[source_candidate_idx];
        if (n_candidate == 0) {
            continue;
        }

        // Step 7 is fully deferred in store=true mode. If the source pair target
        // is local, mirror the already-updated Step-5 block later without a GEMM.
        // If the source pair target is remote, serialize that transpose now as a
        // REPLACE payload for the target owner.
        const bool source_pair_wanted =
            pair_filter == nullptr ||
            (*pair_filter)(source_morton, candidate_morton, source_morton);
        if (candidate_is_writable && source_pair_wanted) {
            DeferredXnnTargetKey source_mirror_target;
            source_mirror_target.box_morton = candidate_morton;
            source_mirror_target.neighbor_morton = source_morton;
            source_mirror_target.kind = DeferredXnnTargetKind::NEAR_A_NS;
            auto source_mirror_insert =
                scratch.preallocated_mirror_targets.insert(source_mirror_target);
            if (source_mirror_insert.second) {
                if (candidate_unarrived) {
                    local_mirrors.emplace_back(
                        candidate_morton, source_morton);
                } else {
                    ensure_symmetric_owner_deferred_source_pair_target_matrix_storage(
                        source_mirror_target,
                        static_cast<int64_t>(source_box->skeleton_indices.size()),
                        n_candidate,
                        level);
                    mirror_targets.push_back(source_mirror_target);
                }
            }
        } else if (!candidate_is_writable && source_pair_wanted) {
            if (pending == nullptr) {
                throw std::runtime_error(
                    "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: remote source pair requires pending updates");
            }
            deferred_xnn_cache_remote_source_pair_replace_from_source_box(
                source_box, candidate_morton, n_candidate, *pending);
        }

        const int64_t total_neighbor_points = source_box->X_NR.rows;
        const int64_t r = source_box->X_NR.cols;
        if (static_cast<int64_t>(source_box->deferred_xnn_temp2.size()) != total_neighbor_points * r) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: deferred temp2 size mismatch");
        }

        scratch.owned_row_blocks.clear();
        int64_t packed_total_rows = 0;
        int64_t row_offset = 0;
        const bool candidate_is_eliminated =
            include_ghosts &&
            level.eliminated_boxes.find(candidate_morton) !=
                level.eliminated_boxes.end();
        for (size_t row_idx = 0; row_idx < source_box->one_hop.size(); ++row_idx) {
            const int64_t neighbor_morton = source_box->one_hop[row_idx];
            const int64_t n_neighbor = source_neighbor_counts[row_idx];
            const bool is_diagonal = (neighbor_morton == candidate_morton);

            if (pair_filter != nullptr && n_neighbor > 0 &&
                !(*pair_filter)(
                    source_morton, candidate_morton, neighbor_morton)) {
                row_offset += n_neighbor;
                continue;
            }

            DeferredXnnTargetKind target_kind = DeferredXnnTargetKind::SCHUR;
            bool accumulate_locally = false;
            bool emit_remote_add = false;

            if (include_ghosts) {
                if (is_diagonal) {
                    accumulate_locally = true;
                } else if (n_neighbor > 0) {
                    DeferredXnnEndpoint<CoordType, DataType> neighbor_endpoint =
                        resolve_deferred_xnn_endpoint(level, neighbor_morton);
                    if (neighbor_endpoint.box == nullptr &&
                        !neighbor_endpoint.is_assisting) {
                        throw std::runtime_error(
                            "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: "
                            "CA neighbor box missing");
                    }

                    const bool neighbor_is_eliminated =
                        level.eliminated_boxes.find(neighbor_morton) !=
                            level.eliminated_boxes.end();
                    const int32_t neighbor_sequence =
                        elimination_sequence(neighbor_morton);
                    if (neighbor_is_eliminated &&
                        neighbor_sequence > source_sequence) {
                        row_offset += n_neighbor;
                        continue;
                    }
                    if (candidate_endpoint.is_ghost &&
                        neighbor_endpoint.is_ghost &&
                        candidate_is_eliminated && neighbor_is_eliminated) {
                        row_offset += n_neighbor;
                        continue;
                    }

                    if (neighbor_endpoint.is_assisting) {
                        const auto candidate_key =
                            deferred_xnn_global_ownership_key_for_morton(
                                level, candidate_morton);
                        const auto assisting_key =
                            deferred_xnn_global_ownership_key_for_morton(
                                level, neighbor_morton);
                        if (assisting_key < candidate_key) {
                            std::ostringstream oss;
                            oss << "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: "
                                   "assisting endpoint would own CA pair"
                                << " candidate_box=" << candidate_morton
                                << " assisting_box=" << neighbor_morton
                                << " candidate_color=" << candidate_key.color_order
                                << " assisting_color=" << assisting_key.color_order;
                            throw std::runtime_error(oss.str());
                        }
                    }

                    target_kind = deferred_xnn_boxes_are_one_hop(
                        dimension, candidate_morton, neighbor_morton) ?
                            DeferredXnnTargetKind::NEAR_A_NS :
                            DeferredXnnTargetKind::FAR_A_NS;

                    // CA lazy mode regenerates two-hop fill at its readers;
                    // neither the canonical owner nor its reciprocal view
                    // should materialize a far block here.
                    if (target_kind == DeferredXnnTargetKind::FAR_A_NS &&
                        lazy_far_field_mode() == LazyFarFieldMode::LAZY) {
                        row_offset += n_neighbor;
                        continue;
                    }

                    const bool candidate_owns_pair =
                        deferred_xnn_first_box_owns_pair(
                            level, candidate_morton, neighbor_morton,
                            neighbor_endpoint.is_assisting,
                            /*use_CA_ownership=*/true);
                    if (!candidate_owns_pair &&
                        !neighbor_endpoint.is_assisting) {
                        DeferredXnnTargetKey mirror_target;
                        mirror_target.box_morton = candidate_morton;
                        mirror_target.neighbor_morton = neighbor_morton;
                        mirror_target.kind = target_kind;

                        auto insert_result =
                            scratch.preallocated_mirror_targets.insert(
                                mirror_target);
                        if (insert_result.second) {
                            if (candidate_unarrived) {
                                if (!staged_halo_box_unarrived(
                                        level, neighbor_morton)) {
                                    local_mirrors.emplace_back(
                                        candidate_morton,
                                        neighbor_morton);
                                }
                            } else {
                                ensure_symmetric_owner_deferred_xnn_target_matrix_storage(
                                    mirror_target, n_neighbor,
                                    n_candidate, level);
                                mirror_targets.push_back(mirror_target);
                            }
                        }
                    }
                    accumulate_locally = candidate_owns_pair;
                }
            } else if (is_diagonal) {
                // The color candidate Schur update is local only when the
                // candidate is local. A remote-only diagonal is transported.
                accumulate_locally = candidate_is_local;
                emit_remote_add = !candidate_is_local;
            } else if (n_neighbor > 0) {
                DeferredXnnEndpoint<CoordType, DataType> neighbor_endpoint =
                    resolve_deferred_xnn_endpoint(level, neighbor_morton);
                if (neighbor_endpoint.box == nullptr && !neighbor_endpoint.is_assisting) {
                    throw std::runtime_error(
                        "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: neighbor box missing");
                }

                target_kind = deferred_xnn_boxes_are_one_hop(
                    dimension, candidate_morton, neighbor_morton) ?
                        DeferredXnnTargetKind::NEAR_A_NS :
                        DeferredXnnTargetKind::FAR_A_NS;

                if (target_kind == DeferredXnnTargetKind::FAR_A_NS &&
                    lazy_far_field_mode() == LazyFarFieldMode::LAZY) {
                    row_offset += n_neighbor;
                    continue;
                }

                if (candidate_is_local) {
                    if (neighbor_endpoint.is_local) {
                        const bool candidate_owns_pair = deferred_xnn_first_box_owns_pair(
                            level,
                            candidate_morton,
                            neighbor_morton,
                            false);

                        if (!candidate_owns_pair) {
                            DeferredXnnTargetKey mirror_target;
                            mirror_target.box_morton = candidate_morton;
                            mirror_target.neighbor_morton = neighbor_morton;
                            mirror_target.kind = target_kind;

                            auto insert_result =
                                scratch.preallocated_mirror_targets.insert(mirror_target);
                            if (insert_result.second) {
                                ensure_symmetric_owner_deferred_xnn_target_matrix_storage(
                                    mirror_target, n_neighbor, n_candidate, level);
                                mirror_targets.push_back(mirror_target);
                            }
                        }

                        // Case 1: local-local pair. The canonical local owner keeps
                        // the block here and the mirror pass reconstructs the other
                        // local endpoint later.
                        accumulate_locally = candidate_owns_pair;
                    } else {
                        // Case 2: local-remote pair. The local candidate owns the
                        // local storage update regardless of Morton order, and the
                        // same slice is also exported as a canonical remote ADD.
                        accumulate_locally = true;
                        emit_remote_add = !generated_near;
                    }
                } else {
                    if (!neighbor_endpoint.is_local) {
                        // Case 3: remote-remote pair. There is no local storage to
                        // update, so only one assisting endpoint emits the canonical
                        // ADD. Smaller Morton is the deterministic emitter.
                        emit_remote_add = (candidate_morton < neighbor_morton);
                    }
                    // candidate remote + neighbor local: skip here because the local
                    // candidate will rebuild its own storage and emit the remote ADD.
                }
            }

            if ((accumulate_locally || emit_remote_add) && n_neighbor > 0) {
                DeferredXnnOwnedRowBlock block;
                block.local_order = static_cast<uint64_t>(row_idx);
                block.neighbor_morton = neighbor_morton;
                block.kind = is_diagonal ?
                    DeferredXnnTargetKind::SCHUR : target_kind;
                block.source_row_offset = row_offset;
                block.packed_row_offset = packed_total_rows;
                block.rows = n_neighbor;
                block.accumulate_locally = accumulate_locally;
                block.emit_remote_add = emit_remote_add;
                scratch.owned_row_blocks.push_back(block);
                packed_total_rows += n_neighbor;
            }

            row_offset += n_neighbor;
        }

        if (row_offset != total_neighbor_points) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: row offset mismatch");
        }
        if (packed_total_rows == 0) {
            continue;
        }

        for (const auto& block : scratch.owned_row_blocks) {
            std::vector<DataType> block_update(
                static_cast<size_t>(block.rows * n_candidate));
            int m = static_cast<int>(block.rows);
            int k = static_cast<int>(r);
            int lda = static_cast<int>(total_neighbor_points);
            int ldb = static_cast<int>(total_neighbor_points);
            int ldc = static_cast<int>(block.rows);
            DataType alpha = 1.0;
            DataType beta = 0.0;
            const DataType* source_rows_ptr =
                source_box->deferred_xnn_temp2.data() + block.source_row_offset;

            const DataType* candidate_rows =
                source_box->X_NR.data.data() + candidate_col_offset;
            split_columns_chunked(
                n_candidate, scratch.split_threads,
                [&](int64_t column_begin, int64_t column_end) {
                    int split_columns =
                        static_cast<int>(column_end - column_begin);
                    gemm_("N", "T", &m, &split_columns, &k,
                        &alpha,
                        source_rows_ptr, &lda,
                        candidate_rows + column_begin, &ldb,
                        &beta,
                        block_update.data() + column_begin * block.rows,
                        &ldc);
                });

            if (candidate_unarrived) {
                if (block.kind == DeferredXnnTargetKind::FAR_A_NS) {
                    throw std::runtime_error(
                        "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: "
                        "staged halo overlap requires lazy far-field mode");
                }

                const bool is_schur_delta =
                    block.kind == DeferredXnnTargetKind::SCHUR;
                if (!is_schur_delta &&
                    staged_halo_box_unarrived(
                        level, block.neighbor_morton)) {
                    typename TreeLevel<CoordType, DataType>::StagedPendingDelta
                        transpose_delta;
                    transpose_delta.target_morton = block.neighbor_morton;
                    transpose_delta.neighbor_morton = candidate_morton;
                    transpose_delta.is_schur = false;
                    transpose_delta.rows = n_candidate;
                    transpose_delta.cols = block.rows;
                    transpose_delta.delta.resize(
                        static_cast<size_t>(n_candidate * block.rows));
                    for (int64_t row = 0; row < block.rows; ++row) {
                        for (int64_t col = 0; col < n_candidate; ++col) {
                            transpose_delta.delta[static_cast<size_t>(
                                col + row * n_candidate)] =
                                block_update[static_cast<size_t>(
                                    row + col * block.rows)];
                        }
                    }
                    local_pending.push_back(
                        std::move(transpose_delta));
                }

                typename TreeLevel<CoordType, DataType>::StagedPendingDelta
                    delta;
                delta.target_morton = candidate_morton;
                delta.neighbor_morton = is_schur_delta
                    ? candidate_morton
                    : block.neighbor_morton;
                delta.is_schur = is_schur_delta;
                delta.rows = block.rows;
                delta.cols = n_candidate;
                delta.delta = std::move(block_update);
                local_pending.push_back(std::move(delta));
                continue;
            }

            if (block.accumulate_locally) {
                DeferredXnnTargetKey target;
                target.box_morton = candidate_morton;
                target.neighbor_morton = block.neighbor_morton;
                target.kind = block.kind;

                auto state_it = scratch.accumulated_target_indices.find(target);
                if (state_it == scratch.accumulated_target_indices.end()) {
                    size_t new_idx = scratch.accumulated_targets.size();
                    scratch.accumulated_target_indices.emplace(target, new_idx);

                    DeferredXnnAccumulatedTarget<DataType> target_state;
                    target_state.target = target;
                    target_state.rows = block.rows;
                    target_state.cols = n_candidate;
                    target_state.data =
                        materialize_deferred_xnn_target_matrix_for_accumulation(
                            target, block.rows, n_candidate, level, kernel, dimension);
                    scratch.accumulated_targets.push_back(std::move(target_state));
                    state_it = scratch.accumulated_target_indices.find(target);
                } else {
                    auto& existing = scratch.accumulated_targets[state_it->second];
                    if (existing.rows != block.rows || existing.cols != n_candidate) {
                        throw std::runtime_error(
                            "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: inconsistent target dimensions");
                    }
                }

                auto& target_state = scratch.accumulated_targets[state_it->second];
                accumulate_deferred_xnn_matrix_in_place(
                    target_state.data, block_update);
            }

            if (block.emit_remote_add) {
                if (pending == nullptr) {
                    throw std::runtime_error(
                        "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: remote ADD requires pending updates");
                }

                const EdgeKind edge_kind =
                    (block.kind == DeferredXnnTargetKind::SCHUR) ? EdgeKind::Diag :
                    (block.kind == DeferredXnnTargetKind::NEAR_A_NS) ? EdgeKind::Near :
                                                                       EdgeKind::Far;
                deferred_xnn_accumulate_canonical_edge_delta_from_slice(
                    *add_target,
                    candidate_morton,
                    block.neighbor_morton,
                    edge_kind,
                    block_update.data(),
                    block.rows,
                    0,
                    block.rows,
                    n_candidate);
            }
        }
    }

    // (entries through a collective callback: every request made; abandoned
    // here if one was not fetched, before anything outside the candidate
    // changes; then its remote additions, held back until now)
    entry_checkpoint();
    if (add_target != pending) {
        for (auto& [key, block] : candidate_adds.accumulated_deltas) {
            auto& destination = pending->accumulated_deltas[key];
            if (destination.data.empty()) {
                destination = std::move(block);
            } else {
                if (destination.rows != block.rows || destination.cols != block.cols ||
                    destination.data.size() != block.data.size()) {
                    throw std::runtime_error(
                        "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: "
                        "inconsistent remote addition");
                }
                for (size_t i = 0; i < block.data.size(); ++i) destination.data[i] += block.data[i];
            }
        }
    }

    for (auto& target_state : scratch.accumulated_targets) {
        flush_deferred_xnn_target_matrix_from_accumulation(target_state, level);
    }

    if (!local_pending.empty() || !local_mirrors.empty()) {
        if (level.staged_pending == nullptr) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_for_candidate_box: "
                "staged pending state missing");
        }
        auto& staged_pending = *level.staged_pending;
        std::lock_guard<std::mutex> lock(staged_pending.mutex);
        for (auto& delta : local_pending) {
            staged_pending.deltas.push_back(std::move(delta));
        }
        for (const auto& mirror : local_mirrors) {
            staged_pending.mirrors.insert(mirror);
        }
    }
}


// A nonsymmetric matrix (lazy Schur mode 2): a local candidate C takes, from
// each source B of the wave among its neighbors, its two views of B (B's
// views of C at B's skeleton: A_NS[B] = K(S_B, C) from B's A_SN[C], A_SN[B]
// = K(C, S_B) from B's A_NS[C]), the Schur update temp2_B[C] X_RN_B[:, C],
// and for each neighbor D of B one hop from C its views A_SN[D] +=
// temp2_B[C] X_RN_B[:, D] and A_NS[D] += temp2_B[D] X_RN_B[:, C] (X_RN_B
// the source's original, before finalize_deferred_xnn_source_box).  The far
// pairs are the lazy fill's; a remote box forms its own updates from the
// source's generator (form_near_updates_from_generators_nonsymmetric).  A
// candidate writes only its own views, so the candidates run in parallel.
// Of a pair of local candidates the lower one updates its views (C's A_SN[D]
// and D's A_NS[C] are the same block), the other copies them afterwards
// (mirror_targets, apply_unstructured_owner_deferred_xnn_mirror_nonsymmetric).
template<typename CoordType, typename DataType, typename KernelType>
void apply_unstructured_owner_deferred_xnn_updates_nonsymmetric(
    int64_t candidate_morton,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    const std::unordered_set<int64_t>& wave_box_set,
    DeferredXnnOwnerScratch<DataType>& scratch,
    std::vector<DeferredXnnTargetKey>& mirror_targets) {
    mirror_targets.clear();
    BoxData<CoordType, DataType>* candidate = level.find_local_box(candidate_morton);
    if (candidate == nullptr || candidate->one_hop.empty()) return;
    const int dimension = level.dimension;
    std::vector<DataType> delta;

    // the pass's sources among the candidate's neighbors, and the candidate's
    // position in each source's neighbor list (counts.size(): none)
    auto source_of = [&](int64_t source_morton, size_t& c_idx) -> BoxData<CoordType, DataType>* {
        if (wave_box_set.count(source_morton) == 0) return nullptr;
        BoxData<CoordType, DataType>* source = level.find_local_box(source_morton);
        if (source == nullptr || source->deferred_xnn_temp2.empty()) return nullptr;
        c_idx = source->one_hop.size();
        for (size_t idx = 0; idx < source->one_hop.size(); ++idx) {
            if (source->one_hop[idx] == candidate_morton) c_idx = idx;
        }
        return source;
    };
    // Entries through a collective callback: the kernel blocks the pass
    // starts the candidate's Schur block and views from (those not stored
    // yet, nor made earlier in the pass) asked for first, the same blocks at
    // the same points; abandoned here if one was not fetched, before the
    // candidate changes (the pass adds to what it finds).
    if (kernel_entries(kernel, 0) != nullptr) {
        bool schur = candidate->schur_complement.is_allocated();
        std::set<std::pair<int64_t, bool>> made;
        auto stored = [&](int64_t neighbor, bool a_sn) {
            if (made.count({neighbor, a_sn}) != 0) return true;
            const auto* mb = find_directed_block(*candidate, neighbor, EdgeKind::Near, a_sn);
            return mb != nullptr && (a_sn ? mb->A_SN.is_allocated() : mb->a_ns_is_allocated());
        };
        for (int64_t source_morton : candidate->one_hop) {
            size_t c_idx = 0;
            BoxData<CoordType, DataType>* source = source_of(source_morton, c_idx);
            if (source == nullptr) continue;
            const auto& counts = source->deferred_xnn_neighbor_point_counts;
            if (counts.size() != source->one_hop.size() || c_idx >= counts.size()) continue;  // (the pass reports it)
            const int64_t n_c = counts[c_idx];
            if (n_c == 0) continue;
            made.insert({source_morton, true});
            made.insert({source_morton, false});
            if (source->X_NR.cols == 0) continue;
            for (size_t d_idx = 0; d_idx < counts.size(); ++d_idx) {
                const int64_t n_d = counts[d_idx];
                if (n_d == 0) continue;
                const int64_t neighbor = source->one_hop[d_idx];
                if (d_idx == c_idx) {
                    if (!schur) {
                        (void)directed_view_kernel_block(level, kernel, candidate_morton, n_c, candidate_morton, n_c);
                    }
                    schur = true;
                    continue;
                }
                if (!deferred_xnn_boxes_are_one_hop(dimension, candidate_morton, neighbor)) continue;
                if (neighbor < candidate_morton && level.find_local_box(neighbor) != nullptr) continue;
                if (!stored(neighbor, true)) {
                    (void)directed_view_kernel_block(level, kernel, candidate_morton, n_c, neighbor, n_d);
                    made.insert({neighbor, true});
                }
                if (!stored(neighbor, false)) {
                    (void)directed_view_kernel_block(level, kernel, neighbor, n_d, candidate_morton, n_c);
                    made.insert({neighbor, false});
                }
            }
        }
        entry_checkpoint();
    }

    for (int64_t source_morton : candidate->one_hop) {
        if (wave_box_set.count(source_morton) == 0) continue;
        BoxData<CoordType, DataType>* source = level.find_local_box(source_morton);
        if (source == nullptr || source->deferred_xnn_temp2.empty()) continue;
        const auto& counts = source->deferred_xnn_neighbor_point_counts;
        if (counts.size() != source->one_hop.size() || !source->X_RN.is_allocated() ||
            source->X_NR.rows != source->X_RN.cols) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_nonsymmetric: incomplete source " +
                std::to_string(source_morton));
        }
        std::vector<int64_t> offsets(counts.size() + 1, 0);
        size_t c_idx = counts.size();
        for (size_t idx = 0; idx < counts.size(); ++idx) {
            offsets[idx + 1] = offsets[idx] + counts[idx];
            if (source->one_hop[idx] == candidate_morton) c_idx = idx;
        }
        if (c_idx == counts.size()) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_nonsymmetric: candidate not a neighbor of its source");
        }
        const int64_t n_c = counts[c_idx];
        if (n_c == 0) continue;

        // the candidate's views of the source, the source's of it (Step 7)
        const auto* source_ns = find_directed_block(*source, candidate_morton, EdgeKind::Near, false);
        const auto* source_sn = find_directed_block(*source, candidate_morton, EdgeKind::Near, true);
        if (source_ns == nullptr || !source_ns->a_ns_is_allocated() ||
            source_sn == nullptr || !source_sn->A_SN.is_allocated()) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_updates_nonsymmetric: source views missing");
        }
        upsert_directed_block(*candidate, source_morton, EdgeKind::Near, true)
            .A_SN.set_owned(source_ns->a_ns_rows(), source_ns->a_ns_cols(), source_ns->copy_a_ns_data(),
                            MatrixStorage<DataType>::FULL);
        upsert_directed_block(*candidate, source_morton, EdgeKind::Near, false)
            .set_a_ns_owned(source_sn->A_SN.rows, source_sn->A_SN.cols, source_sn->copy_a_sn_data(),
                            MatrixStorage<DataType>::FULL);

        const int64_t total = source->X_NR.rows;
        const int64_t r = source->X_NR.cols;
        if (r == 0) continue;
        const DataType* temp2 = source->deferred_xnn_temp2.data();  // total x r
        const DataType* x_rn = source->X_RN.data.data();            // r x total
        const int64_t x_rn_ld = source->X_RN.lda;
        // delta = temp2[rows of a] X_RN[:, cols of b]
        auto product = [&](int64_t a_off, int64_t a_rows, int64_t b_off, int64_t b_cols) {
            delta.resize(static_cast<size_t>(a_rows * b_cols));
            int m = static_cast<int>(a_rows);
            int k = static_cast<int>(r);
            int lda = static_cast<int>(total);
            int ldb = static_cast<int>(x_rn_ld);
            int ldc = static_cast<int>(a_rows);
            DataType alpha = 1.0;
            DataType beta = 0.0;
            split_columns_chunked(b_cols, scratch.split_threads, [&](int64_t c0, int64_t c1) {
                int n = static_cast<int>(c1 - c0);
                gemm_("N", "N", &m, &n, &k, &alpha, temp2 + a_off, &lda,
                      x_rn + (b_off + c0) * x_rn_ld, &ldb, &beta, delta.data() + c0 * a_rows, &ldc);
            });
        };
        for (size_t d_idx = 0; d_idx < counts.size(); ++d_idx) {
            const int64_t n_d = counts[d_idx];
            if (n_d == 0) continue;
            const int64_t neighbor = source->one_hop[d_idx];
            if (d_idx == c_idx) {
                product(offsets[c_idx], n_c, offsets[c_idx], n_c);
                add_to_schur_block(level, kernel, *candidate, n_c, delta.data(), n_c);
                continue;
            }
            if (!deferred_xnn_boxes_are_one_hop(dimension, candidate_morton, neighbor)) continue;
            if (neighbor < candidate_morton && level.find_local_box(neighbor) != nullptr) {
                DeferredXnnTargetKey mirror;
                mirror.box_morton = candidate_morton;
                mirror.neighbor_morton = neighbor;
                mirror.kind = DeferredXnnTargetKind::NEAR_A_NS;
                if (std::find(mirror_targets.begin(), mirror_targets.end(), mirror) == mirror_targets.end()) {
                    mirror_targets.push_back(mirror);
                }
                continue;
            }
            product(offsets[c_idx], n_c, offsets[d_idx], n_d);
            add_to_directed_view(level, kernel, *candidate, neighbor, EdgeKind::Near, true,
                                 n_c, n_d, delta.data(), n_c);
            product(offsets[d_idx], n_d, offsets[c_idx], n_c);
            add_to_directed_view(level, kernel, *candidate, neighbor, EdgeKind::Near, false,
                                 n_c, n_d, delta.data(), n_d);
        }
    }
}

// Before the mirrors run in parallel: every view a mirror writes made in its
// box, so that the mirrors only fill blocks.  A view made inside the parallel
// pass grows its box's block list (and map) while another mirror reads a view
// of that box: the reader then copies freed memory.
template<typename CoordType, typename DataType>
void prepare_unstructured_owner_deferred_xnn_mirror_nonsymmetric(
    TreeLevel<CoordType, DataType>& level,
    const std::vector<std::vector<DeferredXnnTargetKey>>& mirror_targets) {
    for (const auto& targets : mirror_targets) {
        for (const auto& mirror : targets) {
            BoxData<CoordType, DataType>* box = level.find_local_box(mirror.box_morton);
            if (box == nullptr) continue;  // (the mirror reports it)
            upsert_directed_block(*box, mirror.neighbor_morton, EdgeKind::Near, true);
            upsert_directed_block(*box, mirror.neighbor_morton, EdgeKind::Near, false);
        }
    }
}

// After every candidate's pass: a candidate's views of a lower local
// candidate, copied from that one's (its A_SN[D] = D's A_NS[C], its A_NS[D]
// = D's A_SN[C]).  Reads only views no mirror writes; the views it writes are
// made beforehand (prepare_unstructured_owner_deferred_xnn_mirror_nonsymmetric).
template<typename CoordType, typename DataType>
void apply_unstructured_owner_deferred_xnn_mirror_nonsymmetric(
    TreeLevel<CoordType, DataType>& level,
    const std::vector<DeferredXnnTargetKey>& mirror_targets) {
    for (const auto& mirror : mirror_targets) {
        BoxData<CoordType, DataType>* box = level.find_local_box(mirror.box_morton);
        BoxData<CoordType, DataType>* owner = level.find_local_box(mirror.neighbor_morton);
        if (box == nullptr || owner == nullptr) {
            throw std::runtime_error("apply_unstructured_owner_deferred_xnn_mirror_nonsymmetric: box missing");
        }
        const auto* owner_ns = find_directed_block(*owner, mirror.box_morton, EdgeKind::Near, false);
        const auto* owner_sn = find_directed_block(*owner, mirror.box_morton, EdgeKind::Near, true);
        if (owner_ns == nullptr || !owner_ns->a_ns_is_allocated() ||
            owner_sn == nullptr || !owner_sn->A_SN.is_allocated()) {
            throw std::runtime_error(
                "apply_unstructured_owner_deferred_xnn_mirror_nonsymmetric: owner views missing");
        }
        upsert_directed_block(*box, mirror.neighbor_morton, EdgeKind::Near, true)
            .A_SN.set_owned(owner_ns->a_ns_rows(), owner_ns->a_ns_cols(), owner_ns->copy_a_ns_data(),
                            MatrixStorage<DataType>::FULL);
        upsert_directed_block(*box, mirror.neighbor_morton, EdgeKind::Near, false)
            .set_a_ns_owned(owner_sn->A_SN.rows, owner_sn->A_SN.cols, owner_sn->copy_a_sn_data(),
                            MatrixStorage<DataType>::FULL);
    }
}

}  // namespace color_unstructured
}  // namespace butterfly
