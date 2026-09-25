#pragma once

namespace butterfly {
namespace color_unstructured {

inline void validate_options(const ProgramOptions& options) {
    if (options.unstructured != 1) {
        throw std::invalid_argument(
            "color_unstructured requires H2_unstructured=1");
    }
    validate_h2_backend_selection(options);
}

template<typename CoordType, typename DataType>
void validate_tree(const H2<CoordType, DataType>* solver) {
    if (solver == nullptr || solver->tree == nullptr) {
        throw std::invalid_argument(
            "color_unstructured: solver tree has not been initialized");
    }
    validate_options(solver->options);
    for (int level = 0; level < solver->tree->num_levels; ++level) {
        if (solver->tree->level_requests_CA(level)) {
            throw std::runtime_error(
                "color_unstructured: tree contains a CA level");
        }
    }
}

}  // namespace color_unstructured
}  // namespace butterfly
