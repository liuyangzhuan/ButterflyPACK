// Column-independence of the device solve's per-box kernels (solve_kernels.cu).
// They take the right-hand sides in passes of several columns, reading a
// box's factors once a pass, and a column's arithmetic is the single-column
// order: every column of a many-column launch must be bitwise the result of
// that column launched alone.  Here random boxes - small (x_N in shared
// memory), leaf-sized (x_N in the work buffer) and large (fewer columns a
// pass) - go through each kernel with several column counts, including
// counts that leave a partial last pass, and each column of the vectors (and
// of the work buffer where it is an output: the neighbors' updates of the
// solve's forward and the multiply's backward step; elsewhere it is scratch
// for x_N when x_N does not fit in shared memory) is compared bitwise with
// its own single-column run.  Boxes beyond kSolveSharedLimit (one column a
// block, in the card's opt-in shared memory) and beyond the card's shared
// memory (global memory) are among them.  Then the global-memory kernels
// against the shared-memory ones, bitwise: each launch again with box bounds
// beyond the card's shared memory, its scratch as the solve sizes it
// (solve_global_points) and one block a box (each block taking every column
// in turn).
//
//   make h2_gpu_solve_columns
//   srun -n 1 --gpus 1 ./h2_gpu_solve_columns
//   srun -n 1 --gpus 1 ./h2_gpu_solve_columns --time [boxes]
//
// --time: each kernel's time on a wave of leaf-sized boxes (default 1024),
// one column and 16, in double; with LD_PRELOAD of another libbutterflypack
// (the launches are the same functions) it times that library's kernels.
//
// Kernels: launch_solve_forward / _backward / _diagonal and launch_mul_forward
// / _diagonal / _backward, in double and complex.

#include "../solve_kernels.hpp"
#include "../gpu_runtime.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace fmm::gpu;

namespace {

template<typename T>
T random_value(std::mt19937_64& g) {
    std::normal_distribution<double> d(0.0, 1.0);
    if constexpr (std::is_same_v<T, double>) {
        return d(g);
    } else {
        const double re = d(g), im = d(g);
        return dcomplex(re, im);
    }
}

template<typename T>
T* upload(const std::vector<T>& h, std::vector<void*>& allocs) {
    T* d = nullptr;
    check_cuda(cudaMalloc(&d, std::max<size_t>(h.size(), 1) * sizeof(T)), "cudaMalloc");
    if (!h.empty()) check_cuda(cudaMemcpy(d, h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice), "upload");
    allocs.push_back(d);
    return d;
}

template<typename T>
void download(const T* d, std::vector<T>& h) {
    if (!h.empty()) check_cuda(cudaMemcpy(h.data(), d, h.size() * sizeof(T), cudaMemcpyDeviceToHost), "download");
}

struct Shape {
    const char* name;
    int boxes;            // launched; as many again are their neighbours
    int nmin, nmax;       // points a box
    double kfrac;         // skeleton share
    int slots_min, slots_max;
};

// A wave of boxes with their neighbours (read through the slots only) and
// ghosts (another rank's boxes, in the ghost area).
template<typename T>
struct Wave {
    std::vector<SolveBox> boxes;   // the launched boxes first
    std::vector<int> n_of;         // points of every box and every ghost (for the column layout)
    std::vector<int64_t> vec_of;   // their offsets (level vector, ghost area)
    std::vector<int64_t> work_of;  // the launched boxes' work offsets
    int64_t points = 0, ghost_points = 0, work_points = 0;
    int count = 0, ghosts = 0, max_n = 0, max_r = 0, max_ntot = 0;
    SolveBox* d_boxes = nullptr;
    SolveSlot* d_slots = nullptr;
    int* d_wave = nullptr;
    std::vector<void*> allocs;
    ~Wave() {
        for (void* p : allocs) cudaFree(p);
    }
};

template<typename T>
void build(Wave<T>& w, const Shape& s, std::mt19937_64& g) {
    std::uniform_int_distribution<int> npts(s.nmin, s.nmax), nsl(s.slots_min, s.slots_max);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const int nb = 2 * s.boxes;
    const int ng = s.boxes;  // ghosts
    std::vector<std::vector<int>> skel(nb + ng), red(nb + ng);
    std::vector<const int*> d_skel(nb + ng);
    std::vector<int> k_of(nb + ng);
    w.count = s.boxes;
    w.ghosts = ng;
    for (int b = 0; b < nb + ng; ++b) {
        const int n = npts(g);
        const int k = std::max(1, std::min(n - 1, static_cast<int>(s.kfrac * n * (0.8 + 0.4 * u(g)))));
        std::vector<int> perm(static_cast<size_t>(n));
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), g);
        skel[b].assign(perm.begin(), perm.begin() + k);
        red[b].assign(perm.begin() + k, perm.end());
        k_of[b] = k;
        d_skel[b] = upload(skel[b], w.allocs);
        w.n_of.push_back(n);
        if (b < nb) {
            w.vec_of.push_back(w.points);
            w.points += n;
        } else {
            w.vec_of.push_back(w.ghost_points);
            w.ghost_points += n;
        }
    }
    std::vector<SolveSlot> slots;
    std::uniform_int_distribution<int> pick_nb(s.boxes, nb - 1), pick_g(nb, nb + ng - 1);
    for (int b = 0; b < nb; ++b) {
        SolveBox sb{};
        sb.vec = w.vec_of[b];
        sb.n = w.n_of[b];
        sb.k = k_of[b];
        sb.r = sb.n - sb.k;
        sb.skel = d_skel[b];
        sb.red = upload(red[b], w.allocs);
        sb.slot0 = static_cast<int>(slots.size());
        sb.nslots = 0;
        sb.ntot = 0;
        if (b < s.boxes) {
            const int ns = nsl(g);
            for (int a = 0; a < ns; ++a) {
                const bool ghost = u(g) < 0.3;
                const int m = ghost ? pick_g(g) : pick_nb(g);
                SolveSlot sl{};
                sl.row0 = sb.ntot;
                sl.full = u(g) < 0.5 ? 1 : 0;
                sl.count = sl.full ? w.n_of[m] : k_of[m];
                sl.ghost = ghost ? 1 : 0;
                sl.vec = w.vec_of[m];
                sl.n = w.n_of[m];
                sl.skel = d_skel[m];
                slots.push_back(sl);
                sb.ntot += sl.count;
                ++sb.nslots;
            }
            const int k = sb.k, r = sb.r, nt = sb.ntot;
            auto random_block = [&](size_t count) {
                std::vector<T> h(count);
                for (auto& x : h) x = random_value<T>(g);
                return h;
            };
            // a step left out now and then (the kernels' null branches)
            sb.T = u(g) < 0.1 ? nullptr : upload(random_block(static_cast<size_t>(k) * r), w.allocs);
            sb.xsr = u(g) < 0.1 ? nullptr : upload(random_block(static_cast<size_t>(k) * r), w.allocs);
            sb.xnr = upload(random_block(static_cast<size_t>(nt) * r), w.allocs);
            std::vector<T> lu = random_block(static_cast<size_t>(r) * r);
            for (int i = 0; i < r; ++i) lu[static_cast<size_t>(i) * r + i] = lu[static_cast<size_t>(i) * r + i] + T(2.0 * r);
            sb.lu = upload(lu, w.allocs);
            std::vector<int> ipiv(static_cast<size_t>(r));
            for (int i = 0; i < r; ++i) ipiv[i] = i + 1 + std::uniform_int_distribution<int>(0, r - 1 - i)(g);
            sb.ipiv = upload(ipiv, w.allocs);
            sb.work = w.work_points;
            w.work_of.push_back(w.work_points);
            w.work_points += nt;
            w.max_n = std::max(w.max_n, sb.n);
            w.max_r = std::max(w.max_r, sb.r);
            w.max_ntot = std::max(w.max_ntot, nt);
        }
        w.boxes.push_back(sb);
    }
    w.d_boxes = upload(w.boxes, w.allocs);
    w.d_slots = upload(slots, w.allocs);
    std::vector<int> wave(static_cast<size_t>(s.boxes));
    std::iota(wave.begin(), wave.end(), 0);
    w.d_wave = upload(wave, w.allocs);
}

// Column c of a buffer of blocks (block i: `len[i]` points at `off[i]`, its
// nrhs columns at nrhs * off[i] + c * len[i]) as a one-column buffer.
template<typename T>
std::vector<T> column(const std::vector<T>& all, const std::vector<int64_t>& off, const std::vector<int>& len,
                      int nrhs, int c, int64_t points) {
    std::vector<T> one(static_cast<size_t>(points));
    for (size_t i = 0; i < off.size(); ++i) {
        const T* src = all.data() + off[i] * nrhs + static_cast<int64_t>(c) * len[i];
        std::copy(src, src + len[i], one.data() + off[i]);
    }
    return one;
}

enum class Kernel { forward, backward, diagonal, mul_forward, mul_diagonal, mul_backward };
const char* name_of(Kernel k) {
    switch (k) {
        case Kernel::forward: return "solve forward";
        case Kernel::backward: return "solve backward";
        case Kernel::diagonal: return "solve diagonal";
        case Kernel::mul_forward: return "multiply forward";
        case Kernel::mul_diagonal: return "multiply diagonal";
        case Kernel::mul_backward: return "multiply backward";
    }
    return "?";
}

// A launch's global memory (scratch, points: none when null) and the box
// bounds it is sized for (0: the wave's own).  Bounds beyond the card's
// shared memory send boxes that fit it through the global-memory kernels.
// trans: A^T's diagonal steps (the boxes have no X_RS^T or X_RN^T for the
// transposed sweeps).
template<typename T>
struct Launch {
    T* scratch = nullptr;
    int64_t points = 0;
    int max_n = 0, max_r = 0;
    bool trans = false;
};

template<typename T>
void launch(Kernel kind, const Wave<T>& w, T* vec, T* ghost, T* work, int nrhs, cudaStream_t stream,
            const Launch<T>& how = Launch<T>{}) {
    const int max_n = how.max_n > 0 ? how.max_n : w.max_n;
    const int max_r = how.max_r > 0 ? how.max_r : w.max_r;
    switch (kind) {
        case Kernel::forward:
            launch_solve_forward(w.d_boxes, w.d_wave, w.count, vec, work, nrhs, max_n, stream, false, how.scratch,
                                 how.points);
            break;
        case Kernel::backward:
            launch_solve_backward(w.d_boxes, w.d_slots, w.d_wave, w.count, vec, static_cast<const T*>(ghost), work, nrhs,
                                  max_n, w.max_ntot, stream, false, how.scratch, how.points);
            break;
        case Kernel::diagonal:
            launch_solve_diagonal(w.d_boxes, w.d_wave, w.count, vec, nrhs, max_r, stream, how.trans, how.scratch,
                                  how.points);
            break;
        case Kernel::mul_forward:
            launch_mul_forward(w.d_boxes, w.d_slots, w.d_wave, w.count, vec, static_cast<const T*>(ghost), work, nrhs,
                               max_n, w.max_ntot, stream, false, how.scratch, how.points);
            break;
        case Kernel::mul_diagonal:
            launch_mul_diagonal(w.d_boxes, w.d_wave, w.count, vec, nrhs, max_r, stream, how.trans, how.scratch,
                                how.points);
            break;
        case Kernel::mul_backward:
            launch_mul_backward(w.d_boxes, w.d_wave, w.count, vec, work, nrhs, max_n, stream, false, how.scratch,
                                how.points);
            break;
    }
    check_cuda(cudaStreamSynchronize(stream), name_of(kind));
}

// The work buffer holds the kernel's output (the neighbors' updates), not
// scratch for x_N.
bool work_is_output(Kernel k) { return k == Kernel::forward || k == Kernel::mul_backward; }

// Returns the number of columns that differ from their single-column run.
template<typename T>
long check(Kernel kind, const Wave<T>& w, int nrhs, std::mt19937_64& g) {
    cudaStream_t stream = Context::instance().stream();
    const int nb = static_cast<int>(w.boxes.size());
    std::vector<int64_t> box_off(w.vec_of.begin(), w.vec_of.begin() + nb);
    std::vector<int> box_len(w.n_of.begin(), w.n_of.begin() + nb);
    std::vector<int64_t> ghost_off(w.vec_of.begin() + nb, w.vec_of.end());
    std::vector<int> ghost_len(w.n_of.begin() + nb, w.n_of.end());
    std::vector<int> work_len;
    for (int b = 0; b < w.count; ++b) work_len.push_back(w.boxes[static_cast<size_t>(b)].ntot);
    auto random_vector = [&](int64_t n) {
        std::vector<T> h(static_cast<size_t>(n));
        for (auto& x : h) x = random_value<T>(g);
        return h;
    };
    const std::vector<T> vec0 = random_vector(w.points * nrhs);
    const std::vector<T> ghost0 = random_vector(w.ghost_points * nrhs);
    const std::vector<T> work0 = random_vector(std::max<int64_t>(w.work_points, 1) * nrhs);
    std::vector<void*> allocs;
    // boxes beyond the card's shared memory: their global memory, as the solve sizes it
    Launch<T> how;
    how.points = solve_global_points(w.count, w.max_n, nrhs, sizeof(T));
    if (how.points > 0) how.scratch = upload(std::vector<T>(static_cast<size_t>(how.points)), allocs);
    // all the columns at once
    T* d_vec = upload(vec0, allocs);
    T* d_ghost = upload(ghost0, allocs);
    T* d_work = upload(work0, allocs);
    launch(kind, w, d_vec, d_ghost, d_work, nrhs, stream, how);
    std::vector<T> vec_all(vec0.size()), work_all(work0.size());
    download(d_vec, vec_all);
    download(d_work, work_all);
    long differ = 0;
    for (int c = 0; c < nrhs; ++c) {
        // column c alone
        std::vector<void*> one_allocs;
        T* v1 = upload(column(vec0, box_off, box_len, nrhs, c, w.points), one_allocs);
        T* g1 = upload(column(ghost0, ghost_off, ghost_len, nrhs, c, w.ghost_points), one_allocs);
        T* w1 = upload(column(work0, w.work_of, work_len, nrhs, c, std::max<int64_t>(w.work_points, 1)), one_allocs);
        launch(kind, w, v1, g1, w1, 1, stream, how);
        std::vector<T> vec_one(static_cast<size_t>(w.points)), work_one(static_cast<size_t>(std::max<int64_t>(w.work_points, 1)));
        download(v1, vec_one);
        download(w1, work_one);
        for (void* p : one_allocs) cudaFree(p);
        const std::vector<T> vec_c = column(vec_all, box_off, box_len, nrhs, c, w.points);
        const std::vector<T> work_c = column(work_all, w.work_of, work_len, nrhs, c, std::max<int64_t>(w.work_points, 1));
        const bool same_vec = std::memcmp(vec_c.data(), vec_one.data(), vec_one.size() * sizeof(T)) == 0;
        const bool same_work = !work_is_output(kind) || w.work_points == 0 ||
                               std::memcmp(work_c.data(), work_one.data(), static_cast<size_t>(w.work_points) * sizeof(T)) == 0;
        if (!same_vec || !same_work) {
            if (++differ <= 2) {
                std::printf("    %s, nrhs %d: column %d differs from its single-column run (%s%s)\n", name_of(kind), nrhs,
                            c, same_vec ? "" : "vector ", same_work ? "" : "work");
            }
        }
    }
    for (void* p : allocs) cudaFree(p);
    return differ;
}

// The launch with box bounds beyond the card's shared memory (the
// global-memory kernels) against the launch as it is, bitwise: its scratch
// as the solve sizes it, then one block a box.  Returns the launches that
// differ (of two).
template<typename T>
long check_global(Kernel kind, const Wave<T>& w, int nrhs, bool trans, std::mt19937_64& g) {
    cudaStream_t stream = Context::instance().stream();
    const int beyond = static_cast<int>(solve_shared_limit() / sizeof(T)) + 1;
    auto random_vector = [&](int64_t n) {
        std::vector<T> h(static_cast<size_t>(n));
        for (auto& x : h) x = random_value<T>(g);
        return h;
    };
    const std::vector<T> vec0 = random_vector(w.points * nrhs);
    const std::vector<T> ghost0 = random_vector(w.ghost_points * nrhs);
    const std::vector<T> work0 = random_vector(std::max<int64_t>(w.work_points, 1) * nrhs);
    auto run = [&](const Launch<T>& how, std::vector<T>& vec, std::vector<T>& work) {
        std::vector<void*> allocs;
        T* d_vec = upload(vec0, allocs);
        T* d_ghost = upload(ghost0, allocs);
        T* d_work = upload(work0, allocs);
        launch(kind, w, d_vec, d_ghost, d_work, nrhs, stream, how);
        vec.resize(vec0.size());
        work.resize(work0.size());
        download(d_vec, vec);
        download(d_work, work);
        for (void* p : allocs) cudaFree(p);
    };
    Launch<T> shared;
    shared.trans = trans;
    std::vector<T> vec_ref, work_ref;
    run(shared, vec_ref, work_ref);
    long differ = 0;
    for (int variant = 0; variant < 2; ++variant) {
        Launch<T> global = shared;
        global.max_n = global.max_r = beyond;
        global.points = variant == 0 ? solve_global_points(w.count, beyond, nrhs, sizeof(T))
                                     : static_cast<int64_t>(w.count) * beyond;
        std::vector<void*> allocs;
        global.scratch = upload(std::vector<T>(static_cast<size_t>(global.points)), allocs);
        std::vector<T> vec, work;
        run(global, vec, work);
        for (void* p : allocs) cudaFree(p);
        const bool same_vec = std::memcmp(vec.data(), vec_ref.data(), vec.size() * sizeof(T)) == 0;
        const bool same_work = !work_is_output(kind) || w.work_points == 0 ||
                               std::memcmp(work.data(), work_ref.data(), work.size() * sizeof(T)) == 0;
        if (!same_vec || !same_work) {
            ++differ;
            std::printf("    %s%s, nrhs %d: global memory (%s) differs from shared memory (%s%s)\n", name_of(kind),
                        trans ? " (transposed)" : "", nrhs, variant == 0 ? "the solve's scratch" : "one block a box",
                        same_vec ? "" : "vector ", same_work ? "" : "work");
        }
    }
    return differ;
}

template<typename T>
int test_type(const char* type, const std::vector<int>& counts, unsigned seed) {
    std::mt19937_64 g(seed);
    // (elements of shared memory: kSolveSharedLimit, the card's opt-in limit)
    const int pass = static_cast<int>(kSolveSharedLimit / sizeof(T));
    const int card = static_cast<int>(solve_shared_limit() / sizeof(T));
    const Shape shapes[] = {
        {"small boxes (x_N in shared memory)", 64, 20, 64, 0.3, 2, 6},
        {"leaf-sized boxes (x_N in the work buffer)", 48, 180, 216, 0.2, 14, 24},
        {"large boxes (fewer columns a pass)", 12, 600, 1000, 0.15, 2, 5},
        {"boxes beyond kSolveSharedLimit (one column a block)", 2, pass + 200, pass + 1200, 0.1, 1, 2},
        {"boxes beyond the card's shared memory (global memory)", 1, card + card / 20, card + card / 10, 0.02, 1, 1},
    };
    const Kernel kernels[] = {Kernel::forward, Kernel::backward, Kernel::diagonal,
                              Kernel::mul_forward, Kernel::mul_diagonal, Kernel::mul_backward};
    const std::vector<int> few = {2, 5};  // (the two shapes of large boxes)
    int failures = 0;
    for (size_t si = 0; si < std::size(shapes); ++si) {
        const Shape& s = shapes[si];
        const bool large = si >= 3;
        Wave<T> w;
        build(w, s, g);
        std::printf("%s, %s: %d boxes, max n %d, max r %d, max ntot %d\n", type, s.name, w.count, w.max_n, w.max_r,
                    w.max_ntot);
        for (Kernel k : kernels) {
            long bad = 0, total = 0;
            for (int nrhs : large ? few : counts) {
                bad += check(k, w, nrhs, g);
                total += nrhs;
            }
            std::printf("  %-18s %ld of %ld columns differ\n", name_of(k), bad, total);
            if (bad != 0) ++failures;
        }
        if (large) continue;
        // the same boxes on global memory
        for (Kernel k : kernels) {
            const bool diagonal = k == Kernel::diagonal || k == Kernel::mul_diagonal;
            long bad = 0, total = 0;
            for (int trans = 0; trans < (diagonal ? 2 : 1); ++trans) {
                for (int nrhs : counts) {
                    bad += check_global(k, w, nrhs, trans != 0, g);
                    total += 2;
                }
            }
            std::printf("  %-18s %ld of %ld global-memory launches differ from shared memory\n", name_of(k), bad,
                        total);
            if (bad != 0) ++failures;
        }
    }
    return failures;
}

// Each kernel's time (cudaEvent, the mean of five launches after a first)
template<typename T>
void time_kernels(int boxes) {
    std::mt19937_64 g(4242);
    const Shape leaf{"leaf-sized boxes", boxes, 180, 216, 0.2, 14, 24};
    Wave<T> w;
    build(w, leaf, g);
    std::printf("timing: %d leaf-sized boxes, max n %d, max r %d, max ntot %d\n", w.count, w.max_n, w.max_r,
                w.max_ntot);
    cudaStream_t stream = Context::instance().stream();
    cudaEvent_t e0, e1;
    check_cuda(cudaEventCreate(&e0), "event");
    check_cuda(cudaEventCreate(&e1), "event");
    const Kernel kernels[] = {Kernel::forward, Kernel::backward, Kernel::diagonal,
                              Kernel::mul_forward, Kernel::mul_diagonal, Kernel::mul_backward};
    for (int nrhs : {1, 16}) {
        std::vector<void*> allocs;
        auto random_vector = [&](int64_t n) {
            std::vector<T> h(static_cast<size_t>(n));
            for (auto& x : h) x = random_value<T>(g);
            return h;
        };
        T* d_vec = upload(random_vector(w.points * nrhs), allocs);
        T* d_ghost = upload(random_vector(w.ghost_points * nrhs), allocs);
        T* d_work = upload(random_vector(std::max<int64_t>(w.work_points, 1) * nrhs), allocs);
        for (Kernel k : kernels) {
            launch(k, w, d_vec, d_ghost, d_work, nrhs, stream);
            float total = 0.0f;
            for (int rep = 0; rep < 5; ++rep) {
                check_cuda(cudaEventRecord(e0, stream), "event");
                launch(k, w, d_vec, d_ghost, d_work, nrhs, stream);
                check_cuda(cudaEventRecord(e1, stream), "event");
                check_cuda(cudaEventSynchronize(e1), "event");
                float ms = 0.0f;
                check_cuda(cudaEventElapsedTime(&ms, e0, e1), "event");
                total += ms;
            }
            std::printf("  %2d column%s  %-18s %8.3f ms\n", nrhs, nrhs == 1 ? " " : "s", name_of(k), total / 5.0f);
        }
        for (void* p : allocs) cudaFree(p);
    }
    cudaEventDestroy(e0);
    cudaEventDestroy(e1);
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Context::instance().activate();
    if (argc > 1 && std::string(argv[1]) == "--time") {
        time_kernels<double>(argc > 2 ? std::atoi(argv[2]) : 1024);
        MPI_Finalize();
        return 0;
    }
    std::printf("solve columns: each column of a many-column launch against that column alone, bitwise "
                "(the card's shared memory a block: %zu bytes)\n", solve_shared_limit());
    int failures = 0;
    failures += test_type<double>("double", {2, 3, 8, 11, 17}, 12345);
    failures += test_type<dcomplex>("complex", {2, 3, 4, 5, 9}, 67890);
    std::printf("%s\n", failures == 0 ? "all columns bitwise equal" : "FAILED");
    MPI_Finalize();
    return failures == 0 ? 0 : 1;
}
