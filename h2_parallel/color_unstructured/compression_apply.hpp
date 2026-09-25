#pragma once

namespace butterfly {
namespace color_unstructured {

template<typename CoordType, typename DataType>
void apply_compressed_interactions_serial(
    TreeLevel<CoordType, DataType>& level,
    const std::vector<SolveDataRequest<CoordType, DataType>>& source_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& target_data,
    const std::unordered_map<int64_t, std::vector<DataType>>& remote_sources) {

    for (size_t box_index = 0;
         box_index < level.local_boxes.size();
         ++box_index) {
        const auto& target = level.local_boxes[box_index];
        if (target.num_points == 0) continue;
        auto& target_vector = target_data[box_index].left_side;
        for (const auto& block : target.h2_interaction_blocks) {
            std::vector<DataType> source;
            if (level.find_local_box(block.source_morton) != nullptr) {
                source = h2_local_vector_for_morton(
                    level, source_data, block.source_morton, true);
            } else {
                const auto it = remote_sources.find(block.source_morton);
                if (it == remote_sources.end()) {
                    throw std::runtime_error(
                        "color_unstructured compressed multiply: "
                        "missing remote multipole vector");
                }
                source = it->second;
            }

            std::vector<DataType> contribution;
            const int64_t nrhs = target_data[box_index].nrhs;
            h2_matrix_vector_product(
                block.matrix, source, contribution, nrhs, 'N');
            const int64_t skeleton_count =
                static_cast<int64_t>(target.skeleton_indices.size());
            if (contribution.size() !=
                static_cast<size_t>(skeleton_count * nrhs)) {
                throw std::runtime_error(
                    "color_unstructured compressed multiply: "
                    "interaction target rank mismatch");
            }
            const int64_t target_rows =
                target_data[box_index].num_points;
            for (int64_t column = 0; column < nrhs; ++column) {
                for (int64_t i = 0; i < skeleton_count; ++i) {
                    target_vector.at(static_cast<size_t>(
                        target.skeleton_indices[static_cast<size_t>(i)] +
                        column * target_rows)) +=
                        contribution[static_cast<size_t>(
                            i + column * skeleton_count)];
                }
            }
        }
    }
}

template<typename CoordType, typename DataType>
void apply_compressed_leaf_near_serial(
    TreeLevel<CoordType, DataType>& leaf,
    const std::vector<SolveDataRequest<CoordType, DataType>>& source_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& target_data,
    const std::unordered_map<int64_t, std::vector<DataType>>& remote_sources) {

    for (size_t box_index = 0;
         box_index < leaf.local_boxes.size();
         ++box_index) {
        const auto& target = leaf.local_boxes[box_index];
        if (target.num_points == 0) continue;
        auto& output = target_data[box_index].left_side;
        for (const auto& block : target.h2_near_blocks) {
            std::vector<DataType> source;
            if (leaf.find_local_box(block.source_morton) != nullptr) {
                source = h2_local_vector_for_morton(
                    leaf, source_data, block.source_morton, false);
            } else {
                const auto it = remote_sources.find(block.source_morton);
                if (it == remote_sources.end()) {
                    throw std::runtime_error(
                        "color_unstructured compressed multiply: "
                        "missing remote near-field vector");
                }
                source = it->second;
            }

            std::vector<DataType> contribution;
            const int64_t nrhs = target_data[box_index].nrhs;
            h2_matrix_vector_product(
                block.matrix, source, contribution, nrhs, 'N');
            if (contribution.size() != output.size()) {
                throw std::runtime_error(
                    "color_unstructured compressed multiply: "
                    "near-field target size mismatch");
            }
            for (size_t i = 0; i < contribution.size(); ++i) {
                output[i] += contribution[i];
            }
        }
    }
}

template<typename CoordType, typename DataType>
void compressed_multiply(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& input,
    std::vector<DataType>& output,
    int nrhs,
    bool verbose) {

    if (nrhs <= 0 || input.size() % static_cast<size_t>(nrhs) != 0) {
        throw std::invalid_argument(
            "color_unstructured compressed multiply: "
            "invalid batched input dimensions");
    }
    const int64_t local_points =
        static_cast<int64_t>(input.size() / static_cast<size_t>(nrhs));

    const int rank = tree->mpi_rank;
    const int leaf_level = tree->num_levels - 1;
    const int first_h2_level = std::min(2, leaf_level);
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> source_data(
        static_cast<size_t>(tree->num_levels));
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> target_data(
        static_cast<size_t>(tree->num_levels));

    int64_t input_offset = 0;
    for (int level_number = first_h2_level;
         level_number <= leaf_level;
         ++level_number) {
        auto& level = tree->levels[level_number];
        if (!level.is_process_active) continue;
        source_data[level_number].resize(level.local_boxes.size());
        target_data[level_number].resize(level.local_boxes.size());
        for (size_t box_index = 0;
             box_index < level.local_boxes.size();
             ++box_index) {
            const auto& box = level.local_boxes[box_index];
            auto& source = source_data[level_number][box_index];
            auto& target = target_data[level_number][box_index];
            source.initialize(box.morton_index, rank, box.num_points, nrhs);
            target.initialize(box.morton_index, rank, box.num_points, nrhs);
            source.skeleton_indices = box.skeleton_indices;
            source.redundant_indices = box.redundant_indices;
            target.skeleton_indices = box.skeleton_indices;
            target.redundant_indices = box.redundant_indices;

            if (level_number == leaf_level) {
                for (int column = 0; column < nrhs; ++column) {
                    for (int64_t i = 0; i < box.num_points; ++i) {
                        source.left_side[static_cast<size_t>(
                            i + static_cast<int64_t>(column) *
                                box.num_points)] =
                            input[static_cast<size_t>(
                                input_offset + i +
                                static_cast<int64_t>(column) *
                                    local_points)];
                    }
                }
                input_offset += box.num_points;
                source.right_side = source.left_side;
            }
        }
    }
    if (input_offset != local_points) {
        throw std::runtime_error(
            "color_unstructured compressed multiply: "
            "local input length does not match leaf DOFs");
    }

    if (leaf_level >= 2) {
        for (int level_number = leaf_level;
             level_number >= 2;
             --level_number) {
            auto& level = tree->levels[level_number];
            if (level.is_process_active) {
                for (size_t box_index = 0;
                     box_index < level.local_boxes.size();
                     ++box_index) {
                    if (level.local_boxes[box_index].num_points == 0) continue;
                    apply_h2_upward_projection(
                        level.local_boxes[box_index],
                        source_data[level_number][box_index]);
                }
            }
            if (level_number > 2) {
                gather_skeleton_to_parent(
                    level, tree->levels[level_number - 1],
                    source_data[level_number],
                    source_data[level_number - 1],
                    tree->dimension, tree->comm);
            }
        }

        for (int level_number = 2;
             level_number <= leaf_level;
             ++level_number) {
            auto& level = tree->levels[level_number];
            std::vector<int64_t> needed;
            if (level.is_process_active) {
                for (const auto& box : level.local_boxes) {
                    if (box.num_points == 0) continue;
                    for (const auto& block : box.h2_interaction_blocks) {
                        if (level.find_local_box(
                                block.source_morton) == nullptr) {
                            needed.push_back(block.source_morton);
                        }
                    }
                }
            }
            const auto remote_sources = exchange_h2_vectors_onehop(
                tree, level_number, source_data[level_number],
                std::move(needed), true, 800 + 8 * level_number);

            if (level.is_process_active) {
                apply_compressed_interactions_serial(
                    level, source_data[level_number],
                    target_data[level_number], remote_sources);
                for (size_t box_index = 0;
                     box_index < level.local_boxes.size();
                     ++box_index) {
                    if (level.local_boxes[box_index].num_points == 0) continue;
                    apply_h2_downward_interpolation(
                        level.local_boxes[box_index],
                        target_data[level_number][box_index]);
                }
            }

            if (level_number < leaf_level) {
                scatter_solution_to_children(
                    tree->levels[level_number + 1], level,
                    target_data[level_number + 1],
                    target_data[level_number],
                    tree->dimension, tree->comm);
            }
        }
    }

    auto& leaf = tree->levels[leaf_level];
    std::vector<int64_t> near_needed;
    if (leaf.is_process_active) {
        for (const auto& box : leaf.local_boxes) {
            if (box.num_points == 0) continue;
            for (const auto& block : box.h2_near_blocks) {
                if (leaf.find_local_box(block.source_morton) == nullptr) {
                    near_needed.push_back(block.source_morton);
                }
            }
        }
    }
    const auto remote_near_sources = exchange_h2_vectors_onehop(
        tree, leaf_level, source_data[leaf_level],
        std::move(near_needed), false, 1800 + 8 * leaf_level);
    if (leaf.is_process_active) {
        apply_compressed_leaf_near_serial(
            leaf, source_data[leaf_level], target_data[leaf_level],
            remote_near_sources);
    }

    output.assign(input.size(), DataType{0});
    if (leaf.is_process_active) {
        int64_t local_row = 0;
        for (const auto& box_data : target_data[leaf_level]) {
            for (int column = 0; column < nrhs; ++column) {
                for (int64_t i = 0; i < box_data.num_points; ++i) {
                    output[static_cast<size_t>(
                        local_row + i +
                        static_cast<int64_t>(column) * local_points)] =
                        box_data.left_side[static_cast<size_t>(
                            i + static_cast<int64_t>(column) *
                                box_data.num_points)];
                }
            }
            local_row += box_data.num_points;
        }
        if (local_row != local_points) {
            throw std::runtime_error(
                "color_unstructured compressed multiply: "
                "local output row mismatch");
        }
    }

    if (verbose && rank == smallest_active_rank(leaf)) {
        std::cout << "Unstructured H2 compression-only multiply complete"
                  << std::endl;
    }
}

template<typename CoordType, typename DataType, typename KernelType>
double compression_quick_verification(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    int num_src = 20,
    unsigned seed = 12345,
    bool verbose = true) {

    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<
            std::is_same_v<DataType, std::complex<float>>, float,
            DataType>>;

    auto mag2 = [](const DataType& value) -> RealType {
        const RealType magnitude = std::abs(value);
        return magnitude * magnitude;
    };

    const auto verification = make_sparse_mvp_verification_data(
        tree, kernel, num_src, seed);
    std::vector<DataType> h2_output;
    const double start = MPI_Wtime();
    compressed_multiply(tree, verification.input, h2_output, 1, false);
    double elapsed = MPI_Wtime() - start;
    if (h2_output.size() != verification.exact_output.size()) {
        throw std::runtime_error(
            "color_unstructured compression verification: "
            "local output length mismatch");
    }

    RealType local_norms[3] = {0, 0, 0};
    for (size_t i = 0; i < h2_output.size(); ++i) {
        local_norms[0] += mag2(h2_output[i]);
        local_norms[1] += mag2(verification.exact_output[i]);
        local_norms[2] += mag2(
            h2_output[i] - verification.exact_output[i]);
    }

    RealType global_norms[3] = {0, 0, 0};
    const MPI_Datatype real_mpi =
        std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(
        local_norms, global_norms, 3, real_mpi, MPI_SUM, tree->comm);
    MPI_Allreduce(
        MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, tree->comm);

    const RealType norm_h2 = std::sqrt(global_norms[0]);
    const RealType norm_exact = std::sqrt(global_norms[1]);
    const RealType acc_mvp = global_norms[1] > RealType(0)
        ? std::sqrt(global_norms[2] / global_norms[1])
        : RealType(0);
    if (tree->mpi_rank == 0 && verbose) {
        std::cout << "H2_CheckError(compression quick): fnorm: "
                  << std::scientific << std::setprecision(7)
                  << norm_h2 << "  " << norm_exact
                  << "  acc_mvp: " << acc_mvp
                  << "  time_mvp: " << elapsed << std::endl;
    }
    return static_cast<double>(acc_mvp);
}

template<typename CoordType, typename DataType>
void compressed_bicgstab(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& rhs,
    std::vector<DataType>& solution,
    int nrhs,
    double tolerance,
    int max_iterations,
    int* completed_iterations,
    double* final_relative_residual,
    bool verbose) {

    if (tolerance <= 0.0) {
        throw std::invalid_argument(
            "color_unstructured BiCGSTAB: tolerance must be positive");
    }
    if (max_iterations <= 0) {
        throw std::invalid_argument(
            "color_unstructured BiCGSTAB: max_iterations must be positive");
    }
    if (nrhs <= 0 || rhs.size() % static_cast<size_t>(nrhs) != 0) {
        throw std::invalid_argument(
            "color_unstructured BiCGSTAB: invalid batched RHS dimensions");
    }

    const size_t value_count = rhs.size();
    const size_t local_rows = value_count / static_cast<size_t>(nrhs);
    solution.assign(value_count, DataType{0});
    std::vector<DataType> residual = rhs;
    std::vector<DataType> shadow = residual;
    std::vector<DataType> direction(value_count, DataType{0});
    std::vector<DataType> image(value_count, DataType{0});
    std::vector<DataType> intermediate(value_count, DataType{0});
    std::vector<DataType> intermediate_image(value_count, DataType{0});

    const std::vector<double> rhs_norm =
        h2_global_column_norms(rhs, nrhs, tree->comm);
    std::vector<double> relative_residual(
        static_cast<size_t>(nrhs), 0.0);
    std::vector<bool> active(static_cast<size_t>(nrhs), false);
    std::vector<int> column_iterations(static_cast<size_t>(nrhs), 0);
    for (int column = 0; column < nrhs; ++column) {
        if (rhs_norm[static_cast<size_t>(column)] > 0.0) {
            active[static_cast<size_t>(column)] = true;
            relative_residual[static_cast<size_t>(column)] = 1.0;
        }
    }

    std::vector<DataType> rho_previous(
        static_cast<size_t>(nrhs), DataType{1});
    std::vector<DataType> alpha(static_cast<size_t>(nrhs), DataType{1});
    std::vector<DataType> omega(static_cast<size_t>(nrhs), DataType{1});
    auto unusable_scalar = [](const DataType& value) {
        const double magnitude = std::abs(value);
        return magnitude == 0.0 || !std::isfinite(magnitude);
    };
    auto any_active = [&]() {
        return std::any_of(
            active.begin(), active.end(), [](bool value) {
                return value;
            });
    };
    auto clear_column = [&](std::vector<DataType>& values, int column) {
        const size_t offset = static_cast<size_t>(column) * local_rows;
        std::fill(
            values.begin() + static_cast<ptrdiff_t>(offset),
            values.begin() + static_cast<ptrdiff_t>(
                offset + local_rows),
            DataType{0});
    };

    for (int iteration = 1;
         iteration <= max_iterations && any_active();
         ++iteration) {
        const std::vector<DataType> rho =
            h2_global_column_dots(
                shadow, residual, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) {
                clear_column(direction, column);
                continue;
            }
            if (unusable_scalar(rho[static_cast<size_t>(column)])) {
                throw std::runtime_error(
                    "color_unstructured BiCGSTAB breakdown: "
                    "rho is zero for RHS " + std::to_string(column));
            }
            const size_t offset =
                static_cast<size_t>(column) * local_rows;
            if (iteration == 1) {
                std::copy_n(
                    residual.begin() + static_cast<ptrdiff_t>(offset),
                    local_rows,
                    direction.begin() + static_cast<ptrdiff_t>(offset));
            } else {
                if (unusable_scalar(
                        omega[static_cast<size_t>(column)])) {
                    throw std::runtime_error(
                        "color_unstructured BiCGSTAB breakdown: "
                        "omega is zero for RHS " +
                        std::to_string(column));
                }
                const DataType beta =
                    (rho[static_cast<size_t>(column)] /
                     rho_previous[static_cast<size_t>(column)]) *
                    (alpha[static_cast<size_t>(column)] /
                     omega[static_cast<size_t>(column)]);
                for (size_t row = 0; row < local_rows; ++row) {
                    const size_t index = offset + row;
                    direction[index] = residual[index] + beta *
                        (direction[index] -
                         omega[static_cast<size_t>(column)] *
                             image[index]);
                }
            }
        }

        compressed_multiply(tree, direction, image, nrhs, false);
        const std::vector<DataType> denominator =
            h2_global_column_dots(
                shadow, image, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) {
                clear_column(intermediate, column);
                continue;
            }
            if (unusable_scalar(
                    denominator[static_cast<size_t>(column)])) {
                throw std::runtime_error(
                    "color_unstructured BiCGSTAB breakdown: "
                    "alpha denominator is zero for RHS " +
                    std::to_string(column));
            }
            alpha[static_cast<size_t>(column)] =
                rho[static_cast<size_t>(column)] /
                denominator[static_cast<size_t>(column)];
            const size_t offset =
                static_cast<size_t>(column) * local_rows;
            for (size_t row = 0; row < local_rows; ++row) {
                const size_t index = offset + row;
                intermediate[index] = residual[index] -
                    alpha[static_cast<size_t>(column)] * image[index];
            }
        }

        const std::vector<double> intermediate_norm =
            h2_global_column_norms(
                intermediate, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) continue;
            relative_residual[static_cast<size_t>(column)] =
                intermediate_norm[static_cast<size_t>(column)] /
                rhs_norm[static_cast<size_t>(column)];
            if (relative_residual[static_cast<size_t>(column)] <=
                tolerance) {
                const size_t offset =
                    static_cast<size_t>(column) * local_rows;
                for (size_t row = 0; row < local_rows; ++row) {
                    const size_t index = offset + row;
                    solution[index] +=
                        alpha[static_cast<size_t>(column)] *
                        direction[index];
                }
                active[static_cast<size_t>(column)] = false;
                column_iterations[static_cast<size_t>(column)] =
                    iteration;
                clear_column(intermediate, column);
            }
        }
        if (!any_active()) break;

        compressed_multiply(
            tree, intermediate, intermediate_image, nrhs, false);
        const std::vector<DataType> omega_numerator =
            h2_global_column_dots(
                intermediate_image, intermediate, nrhs, tree->comm);
        const std::vector<DataType> omega_denominator =
            h2_global_column_dots(
                intermediate_image, intermediate_image, nrhs,
                tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) continue;
            if (unusable_scalar(
                    omega_denominator[
                        static_cast<size_t>(column)])) {
                throw std::runtime_error(
                    "color_unstructured BiCGSTAB breakdown: "
                    "omega denominator is zero for RHS " +
                    std::to_string(column));
            }
            omega[static_cast<size_t>(column)] =
                omega_numerator[static_cast<size_t>(column)] /
                omega_denominator[static_cast<size_t>(column)];
            const size_t offset =
                static_cast<size_t>(column) * local_rows;
            for (size_t row = 0; row < local_rows; ++row) {
                const size_t index = offset + row;
                solution[index] +=
                    alpha[static_cast<size_t>(column)] *
                        direction[index] +
                    omega[static_cast<size_t>(column)] *
                        intermediate[index];
                residual[index] = intermediate[index] -
                    omega[static_cast<size_t>(column)] *
                        intermediate_image[index];
            }
        }

        const std::vector<double> residual_norm =
            h2_global_column_norms(
                residual, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) continue;
            relative_residual[static_cast<size_t>(column)] =
                residual_norm[static_cast<size_t>(column)] /
                rhs_norm[static_cast<size_t>(column)];
            if (relative_residual[static_cast<size_t>(column)] <=
                tolerance) {
                active[static_cast<size_t>(column)] = false;
                column_iterations[static_cast<size_t>(column)] =
                    iteration;
            } else {
                rho_previous[static_cast<size_t>(column)] =
                    rho[static_cast<size_t>(column)];
            }
        }
    }

    const int iterations_done = *std::max_element(
        column_iterations.begin(), column_iterations.end());
    const double max_relative_residual = *std::max_element(
        relative_residual.begin(), relative_residual.end());
    if (completed_iterations) {
        *completed_iterations = iterations_done;
    }
    if (final_relative_residual) {
        *final_relative_residual = max_relative_residual;
    }
    if (any_active()) {
        throw std::runtime_error(
            "H2 BiCGSTAB did not converge in " +
            std::to_string(max_iterations) +
            " iterations; maximum relative residual=" +
            std::to_string(max_relative_residual));
    }

    if (verbose && tree->mpi_rank == smallest_active_rank(
            tree->levels[tree->num_levels - 1])) {
        std::cout << "H2 BiCGSTAB converged in "
                  << iterations_done
                  << " iterations, maximum relative residual="
                  << max_relative_residual << std::endl;
    }
}

}  // namespace color_unstructured
}  // namespace butterfly
