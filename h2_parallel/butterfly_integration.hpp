#pragma once

// =============================================================================
// Umbrella header for the H2 (format-7) ButterflyPACK integration.
//
// The implementation is split into modules. The include order below encodes the
// dependency DAG, so every cross-module call is defined before it is used (no
// forward declarations needed):
//
//     types  ->  init  ->  solve  ->  verification  ->  factorization
//
//   - butterfly_types.hpp         traits, H2Kernel, H2, ProgramOptions, SparseTestVector
//   - butterfly_init.hpp          level calc, parse_program_options, h2_initiate
//   - butterfly_solve.hpp         gather_local_solution, hierarchical_solve/mul_parallel
//   - butterfly_verification.hpp  verify_solution_direct, h2_direct_verification, h2_verification_smvp
//   - butterfly_factorization.hpp logdet, hierarchical_factorization_parallel, butterfly_factorization_parallel
//
// C_BPACK_wrapper.cpp includes ONLY this file.
// The previous monolithic version is preserved as butterfly_integration_old.cpp.
// =============================================================================

#include "butterfly_types.hpp"
#include "butterfly_init.hpp"
#include "butterfly_solve.hpp"
#include "butterfly_verification.hpp"
#include "butterfly_factorization.hpp"
