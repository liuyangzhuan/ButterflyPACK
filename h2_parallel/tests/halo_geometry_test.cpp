// Self-test of halo_geometry.hpp: for 2D and 3D levels of several sizes and
// region counts, the ghost and assisting sets by region (what the halo
// gather's senders compute) against the reference enumeration from the
// boundary boxes (the original compute_ghost_boxes_and_colors /
// compute_assisting_boxes), for every region.
//
//   make h2_halo_geometry_test && ./h2_halo_geometry_test     (in a build directory), or
//   g++ -std=c++17 -O2 -I../core halo_geometry_test.cpp -o halo_geometry_test && ./halo_geometry_test
//
// (no MPI calls: the header needs only morton.hpp; about 4 minutes)

#include "halo_geometry.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace fmm;

static int failures = 0;

static void check(bool ok, const char* what, int dim, int level, int procs, unsigned q, unsigned r) {
    if (ok) return;
    ++failures;
    std::printf("FAIL %s: dim %d level %d procs %d region %u in region %u\n", what, dim, level, procs, q, r);
}

static void run(int dim, int level, int procs) {
    const RegionLayout L = RegionLayout::of(dim, level, procs);
    if (!L.uniform) {
        std::printf("skip dim %d level %d procs %d: not uniform\n", dim, level, procs);
        return;
    }
    const uint32_t radius = halo_region_radius(L);
    long pairs = 0;
    for (uint32_t q = 0; q < L.regions; ++q) {
        const RegionGhosts ref = ghosts_of_region_reference(L, q);
        const RegionGhosts geo = ghosts_of_region(L, q);
        check(ref.ghosts == geo.ghosts, "ghosts", dim, level, procs, q, q);
        check(ref.blue == geo.blue, "blue ghosts", dim, level, procs, q, q);
        check(ref.orange == geo.orange, "orange completions", dim, level, procs, q, q);
        const RegionAssisting ref_a = assisting_of_region_reference(L, q, ref.ghosts, ref.ghosts);
        // by region, over the candidate radius: the union must be the whole set
        std::vector<int64_t> ghosts_union, assisting_union;
        for_each_region_within(L, q, radius, [&](uint32_t r) {
            const std::vector<int64_t> g = ghosts_of_region_in(L, q, r);
            const std::vector<int64_t> a = assisting_of_region_in(L, q, r);
            // each part is exactly the reference set restricted to r
            std::vector<int64_t> g_ref, a_ref;
            for (int64_t m : ref.ghosts) if (L.region_of(m) == r) g_ref.push_back(m);
            for (int64_t m : ref_a.all) if (L.region_of(m) == r) a_ref.push_back(m);
            check(g == g_ref, "ghosts in region", dim, level, procs, q, r);
            check(a == a_ref, "assisting in region", dim, level, procs, q, r);
            ghosts_union.insert(ghosts_union.end(), g.begin(), g.end());
            assisting_union.insert(assisting_union.end(), a.begin(), a.end());
            ++pairs;
        });
        std::sort(ghosts_union.begin(), ghosts_union.end());
        std::sort(assisting_union.begin(), assisting_union.end());
        check(ghosts_union == ref.ghosts, "ghost union over the radius", dim, level, procs, q, q);
        check(assisting_union == ref_a.all, "assisting union over the radius", dim, level, procs, q, q);
        // no region beyond the radius holds any of them
        for (int64_t m : ref.ghosts) check(L.region_distance(q, L.region_of(m)) <= radius, "ghost radius", dim, level, procs, q, L.region_of(m));
        for (int64_t m : ref_a.all) check(L.region_distance(q, L.region_of(m)) <= radius, "assisting radius", dim, level, procs, q, L.region_of(m));
    }
    std::printf("dim %d level %d procs %d (width %u, radius %u): %u regions, %ld region pairs checked\n", dim, level,
                procs, L.width, radius, L.regions, pairs);
}

int main() {
    // 3D: widths 2, 4, 8, 16 at several levels; 2D: widths 2, 4, 8
    run(3, 2, 8);
    run(3, 3, 8);
    run(3, 3, 64);
    run(3, 4, 8);
    run(3, 4, 64);
    run(3, 4, 512);
    run(3, 5, 64);
    run(3, 5, 512);
    run(3, 6, 512);
    run(3, 7, 512);
    run(2, 2, 4);
    run(2, 3, 4);
    run(2, 3, 16);
    run(2, 4, 16);
    run(2, 5, 16);
    run(2, 5, 64);
    if (failures != 0) {
        std::printf("%d failures\n", failures);
        return 1;
    }
    std::printf("all checks passed\n");
    return 0;
}
