#pragma once

// =============================================================================
// Umbrella header for the H2 (format-7) ButterflyPACK integration.
//
// The implementation is split into modules (paths from h2_parallel/), each
// header including what it uses:
//
//   - core/                         the shared core of every H2 path: tree,
//                                   per-box factorization (Color, CA0, CA3),
//                                   solve, multiply, H2 types and init
//     - butterfly_types.hpp         traits, H2Kernel, H2, ProgramOptions, SparseTestVector
//     - butterfly_init.hpp          level calc, parse_program_options, h2_initiate
//   - structured/                   the structured-grid backend (H2_unstructured=0)
//     - butterfly_solve.hpp         gather_local_solution, hierarchical_solve/mul_parallel
//     - butterfly_compression.hpp   ID-only H2 construction and nested-basis matvec
//     - butterfly_verification.hpp  verify_solution_direct, h2_direct_verification, h2_quick_verification
//     - butterfly_factorization.hpp logdet, hierarchical_factorization_parallel, butterfly_factorization_parallel
//   - unstructured/                 the unstructured-grid backend (H2_unstructured=1, Color only)
//   - h2_backend_dispatch.hpp       each operation to one backend or the other
//
// The GPU backend of H2 is GPU_BACKEND/h2_gpu/.  C_BPACK_wrapper.cpp includes
// ONLY this file.
// =============================================================================

#include "core/butterfly_types.hpp"
#include "core/butterfly_init.hpp"
#include "structured/butterfly_solve.hpp"
#include "structured/butterfly_compression.hpp"
#include "structured/butterfly_verification.hpp"
#include "structured/butterfly_factorization.hpp"
#include "h2_backend_dispatch.hpp"
