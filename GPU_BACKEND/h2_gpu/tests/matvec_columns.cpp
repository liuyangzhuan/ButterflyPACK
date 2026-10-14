// Column-independence of the device H2 matvec's kernels (h2_matvec_kernels.cu).
// With several right-hand sides they take the columns in passes, reading a
// box's blocks once a pass, and a column's arithmetic is the single-column
// order: every column of a many-column launch must be bitwise the result of
// that column launched alone.  Here random boxes, blocks and pairs go
// through each kernel with several column counts (including counts that
// leave a partial last pass), and each column of the outputs is compared
// bitwise with its own single-column run.  The accumulation covers boxes of
// every warp-kernel size (strips of 32 to 256 rows) and larger ones (the
// block kernel), blocks with K and precomputed vectors, from the local,
// ghost and partial arrays.
//
//   make h2_gpu_matvec_columns
//   srun -n 1 --gpus 1 ./h2_gpu_matvec_columns
//   srun -n 1 --gpus 1 ./h2_gpu_matvec_columns --time
//
// --time: each kernel's time on leaf-like shapes, one column and 16, in
// double; with LD_PRELOAD of another libbutterflypack (the launches are the
// same functions) it times that library's kernels.

#include "../h2_matvec_kernels.hpp"
#include "../gpu_runtime.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <functional>
#include <map>
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

template<typename T>
std::vector<T> random_vector(size_t n, std::mt19937_64& g) {
    std::vector<T> h(n);
    for (auto& x : h) x = random_value<T>(g);
    return h;
}

// Source arrays of the accumulation: segments (base in points, len points;
// column c at base * nrhs + c * len), as the compression lays them out.
struct Segments {
    std::vector<int64_t> base;
    std::vector<int> len;
    int64_t points = 0;
    int add(int len_) {
        base.push_back(points);
        len.push_back(len_);
        points += len_;
        return static_cast<int>(base.size()) - 1;
    }
};

// Column c of a buffer of segments as a one-column buffer.
template<typename T>
std::vector<T> column(const std::vector<T>& all, const Segments& sg, int nrhs, int c) {
    std::vector<T> one(static_cast<size_t>(std::max<int64_t>(sg.points, 1)));
    for (size_t i = 0; i < sg.base.size(); ++i) {
        const T* src = all.data() + sg.base[i] * nrhs + static_cast<int64_t>(c) * sg.len[i];
        std::copy(src, src + sg.len[i], one.data() + sg.base[i]);
    }
    return one;
}

// A level of boxes: their vectors (column c at c * ld + vec), interpolations,
// and coupling or near blocks over three source arrays.
template<typename T>
struct Level {
    std::vector<H2Box> boxes;
    std::vector<H2BlockRef> blocks;
    Segments src[3];                  // local, ghost, partial
    int64_t ld = 0;
    int max_k = 0, max_r = 0, max_rows = 0, max_cols = 0;
    H2Box* d_boxes = nullptr;
    H2BlockRef* d_blocks = nullptr;
    std::vector<void*> allocs;
    ~Level() {
        for (void* p : allocs) cudaFree(p);
    }
};

// `near`: blocks of n rows (else k); rows of a box drawn from [nmin, nmax];
// `pool` distinct K matrices of each shape shared by the blocks (0: one each).
template<typename T>
void build_level(Level<T>& lv, int boxes, int nmin, int nmax, double kfrac, int bmin, int bmax, bool near,
                 std::mt19937_64& g, int pool = 0) {
    std::uniform_int_distribution<int> npts(nmin, nmax), nbl(bmin, bmax);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::map<std::pair<int, int>, std::vector<const void*>> pools;
    for (int b = 0; b < boxes; ++b) {
        H2Box hb{};
        hb.n = npts(g);
        hb.k = std::max(1, std::min(hb.n - 1, static_cast<int>(kfrac * hb.n * (0.8 + 0.4 * u(g)))));
        hb.r = hb.n - hb.k;
        std::vector<int> perm(static_cast<size_t>(hb.n));
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), g);
        hb.skel = upload(std::vector<int>(perm.begin(), perm.begin() + hb.k), lv.allocs);
        hb.red = upload(std::vector<int>(perm.begin() + hb.k, perm.end()), lv.allocs);
        hb.T = u(g) < 0.1 ? nullptr : upload(random_vector<T>(static_cast<size_t>(hb.k) * hb.r, g), lv.allocs);
        hb.vec = lv.ld;
        lv.ld += hb.n;
        const int m = near ? hb.n : hb.k;
        const int first = static_cast<int>(lv.blocks.size());
        const int nb = nbl(g);
        for (int j = 0; j < nb; ++j) {
            H2BlockRef br{};
            br.ghost = static_cast<int>(u(g) * 3.0);
            const bool vector_only = u(g) < 0.1;  // a precomputed vector (m entries)
            br.cols = vector_only ? m : std::max(1, npts(g) / (near ? 1 : 3));
            const int len = (vector_only ? m : br.cols) + static_cast<int>(u(g) * 8);  // a segment, read at an offset
            const int sgi = lv.src[br.ghost].add(len);
            br.base = lv.src[br.ghost].base[static_cast<size_t>(sgi)];
            br.stride = len;
            br.off = len - (vector_only ? m : br.cols);
            if (vector_only) {
                br.K = nullptr;
            } else {
                auto& pl = pools[{m, br.cols}];
                if (pool > 0 && static_cast<int>(pl.size()) >= pool) {
                    br.K = pl[std::uniform_int_distribution<size_t>(0, pl.size() - 1)(g)];
                } else {
                    br.K = upload(random_vector<T>(static_cast<size_t>(m) * br.cols, g), lv.allocs);
                    if (pool > 0) pl.push_back(br.K);
                }
                lv.max_cols = std::max(lv.max_cols, br.cols);
            }
            lv.blocks.push_back(br);
        }
        if (near) {
            hb.near0 = first;
            hb.nnear = nb;
        } else {
            hb.block0 = first;
            hb.nblocks = nb;
        }
        lv.max_k = std::max(lv.max_k, hb.k);
        lv.max_r = std::max(lv.max_r, hb.r);
        lv.max_rows = std::max(lv.max_rows, m);
        lv.boxes.push_back(hb);
    }
    lv.d_boxes = upload(lv.boxes, lv.allocs);
    lv.d_blocks = upload(lv.blocks, lv.allocs);
}

enum class Kernel { upward, downward, accumulate, pairs };
const char* name_of(Kernel k) {
    switch (k) {
        case Kernel::upward: return "upward";
        case Kernel::downward: return "downward";
        case Kernel::accumulate: return "accumulate";
        case Kernel::pairs: return "pairs";
    }
    return "?";
}

// One launch of a level kernel on the level vector `vec` (ld per column) and
// the source arrays; returns after it is done.
template<typename T>
void launch_level(Kernel kind, const Level<T>& lv, bool near, T* vec, T* const src[3], int nrhs) {
    cudaStream_t stream = Context::instance().stream();
    const int count = static_cast<int>(lv.boxes.size());
    switch (kind) {
        case Kernel::upward:
            launch_h2_upward(lv.d_boxes, count, vec, lv.ld, nrhs, lv.max_r, stream);
            break;
        case Kernel::downward:
            launch_h2_downward(lv.d_boxes, count, vec, lv.ld, nrhs, lv.max_k, stream);
            break;
        case Kernel::accumulate:
            launch_h2_accumulate(lv.d_boxes, count, lv.d_blocks, vec, lv.ld, static_cast<const T*>(src[0]),
                                 static_cast<const T*>(src[1]), static_cast<const T*>(src[2]), nrhs, near,
                                 lv.max_rows, std::max(lv.max_cols, 1), stream);
            break;
        default:
            break;
    }
    check_cuda(cudaStreamSynchronize(stream), name_of(kind));
}

// Columns of a level kernel that differ from their single-column run.
template<typename T>
long check_level(Kernel kind, const Level<T>& lv, bool near, int nrhs, std::mt19937_64& g) {
    std::vector<void*> allocs;
    const std::vector<T> vec0 = random_vector<T>(static_cast<size_t>(lv.ld) * nrhs, g);
    std::vector<T> src0[3];
    T* d_src[3];
    for (int a = 0; a < 3; ++a) {
        src0[a] = random_vector<T>(static_cast<size_t>(std::max<int64_t>(lv.src[a].points, 1)) * nrhs, g);
        d_src[a] = upload(src0[a], allocs);
    }
    T* d_vec = upload(vec0, allocs);
    launch_level(kind, lv, near, d_vec, d_src, nrhs);
    std::vector<T> vec_all(vec0.size());
    download(d_vec, vec_all);
    long differ = 0;
    for (int c = 0; c < nrhs; ++c) {
        std::vector<void*> one_allocs;
        const std::vector<T> v1(vec0.begin() + static_cast<int64_t>(c) * lv.ld,
                                vec0.begin() + static_cast<int64_t>(c + 1) * lv.ld);
        T* d_v1 = upload(v1, one_allocs);
        T* d_s1[3];
        for (int a = 0; a < 3; ++a) d_s1[a] = upload(column(src0[a], lv.src[a], nrhs, c), one_allocs);
        launch_level(kind, lv, near, d_v1, d_s1, 1);
        std::vector<T> one(static_cast<size_t>(lv.ld));
        download(d_v1, one);
        for (void* p : one_allocs) cudaFree(p);
        if (std::memcmp(one.data(), vec_all.data() + static_cast<int64_t>(c) * lv.ld, one.size() * sizeof(T)) != 0) {
            if (++differ <= 2) {
                std::printf("    %s, nrhs %d: column %d differs from its single-column run\n", name_of(kind), nrhs, c);
            }
        }
    }
    for (void* p : allocs) cudaFree(p);
    return differ;
}

// Pairs of boxes t > s of a symmetric kernel: their vectors in x (column c
// at c * ldx + offset), their partial vectors in `partial` (c * ldp + offset).
template<typename T>
struct Pairs {
    std::vector<H2Pair> pairs;
    int64_t ldx = 0, ldp = 0;
    int max_cols = 0;
    H2Pair* d_pairs = nullptr;
    std::vector<void*> allocs;
    ~Pairs() {
        for (void* p : allocs) cudaFree(p);
    }
};

template<typename T>
void build_pairs(Pairs<T>& ps, int count, int nmin, int nmax, std::mt19937_64& g, int pool = 0) {
    std::uniform_int_distribution<int> npts(nmin, nmax);
    std::vector<std::pair<std::pair<int, int>, const void*>> made;
    for (int p = 0; p < count; ++p) {
        H2Pair pr{};
        pr.rows = npts(g);
        pr.cols = npts(g);
        if (pool > 0 && static_cast<int>(made.size()) >= pool) {  // (timing: reuse a matrix of a pool)
            const auto& e = made[std::uniform_int_distribution<size_t>(0, made.size() - 1)(g)];
            pr.rows = e.first.first;
            pr.cols = e.first.second;
            pr.K = e.second;
        } else {
            pr.K = upload(random_vector<T>(static_cast<size_t>(pr.rows) * pr.cols, g), ps.allocs);
            if (pool > 0) made.push_back({{pr.rows, pr.cols}, pr.K});
        }
        pr.xt = ps.ldx;
        ps.ldx += pr.rows;
        pr.xs = ps.ldx;
        ps.ldx += pr.cols;
        pr.pa = ps.ldp;
        ps.ldp += pr.rows;
        pr.pb = ps.ldp;
        ps.ldp += pr.cols;
        ps.max_cols = std::max(ps.max_cols, pr.cols);
        ps.pairs.push_back(pr);
    }
    ps.d_pairs = upload(ps.pairs, ps.allocs);
}

template<typename T>
void launch_pairs(const Pairs<T>& ps, const T* x, T* partial, int nrhs) {
    cudaStream_t stream = Context::instance().stream();
    launch_h2_pairs(ps.d_pairs, static_cast<int>(ps.pairs.size()), x, ps.ldx, partial, ps.ldp, nrhs, ps.max_cols,
                    stream);
    check_cuda(cudaStreamSynchronize(stream), "pairs");
}

template<typename T>
long check_pairs(const Pairs<T>& ps, int nrhs, std::mt19937_64& g) {
    std::vector<void*> allocs;
    const std::vector<T> x0 = random_vector<T>(static_cast<size_t>(ps.ldx) * nrhs, g);
    const std::vector<T> p0 = random_vector<T>(static_cast<size_t>(ps.ldp) * nrhs, g);
    T* d_x = upload(x0, allocs);
    T* d_p = upload(p0, allocs);
    launch_pairs(ps, d_x, d_p, nrhs);
    std::vector<T> p_all(p0.size());
    download(d_p, p_all);
    long differ = 0;
    for (int c = 0; c < nrhs; ++c) {
        std::vector<void*> one_allocs;
        T* d_x1 = upload(std::vector<T>(x0.begin() + c * ps.ldx, x0.begin() + (c + 1) * ps.ldx), one_allocs);
        T* d_p1 = upload(std::vector<T>(p0.begin() + c * ps.ldp, p0.begin() + (c + 1) * ps.ldp), one_allocs);
        launch_pairs(ps, d_x1, d_p1, 1);
        std::vector<T> one(static_cast<size_t>(ps.ldp));
        download(d_p1, one);
        for (void* p : one_allocs) cudaFree(p);
        if (std::memcmp(one.data(), p_all.data() + c * ps.ldp, one.size() * sizeof(T)) != 0) {
            if (++differ <= 2) std::printf("    pairs, nrhs %d: column %d differs from its single-column run\n", nrhs, c);
        }
    }
    for (void* p : allocs) cudaFree(p);
    return differ;
}

template<typename T>
int test_type(const char* type, const std::vector<int>& counts, unsigned seed) {
    std::mt19937_64 g(seed);
    int failures = 0;
    auto report = [&](const char* what, long bad, long total) {
        std::printf("  %-44s %ld of %ld columns differ\n", what, bad, total);
        if (bad != 0) ++failures;
    };
    struct LevelShape { const char* name; int boxes, nmin, nmax; double kfrac; int bmin, bmax; bool near; };
    const LevelShape shapes[] = {
        {"coupling, up to 32 rows", 64, 20, 90, 0.3, 4, 30, false},
        {"coupling, up to 64 rows", 64, 60, 200, 0.3, 4, 30, false},
        {"near, up to 128 rows", 48, 60, 128, 0.3, 4, 20, true},
        {"near, up to 256 rows (leaf-sized)", 32, 180, 256, 0.2, 10, 27, true},
        {"near, up to 600 rows (the block kernel too)", 12, 150, 600, 0.2, 3, 8, true},
    };
    for (const LevelShape& sh : shapes) {
        Level<T> lv;
        build_level(lv, sh.boxes, sh.nmin, sh.nmax, sh.kfrac, sh.bmin, sh.bmax, sh.near, g);
        std::printf("%s, %s: %d boxes, max k %d, max r %d, max rows %d, max cols %d\n", type, sh.name, sh.boxes,
                    lv.max_k, lv.max_r, lv.max_rows, lv.max_cols);
        for (Kernel k : {Kernel::upward, Kernel::downward, Kernel::accumulate}) {
            long bad = 0, total = 0;
            for (int nrhs : counts) {
                bad += check_level(k, lv, sh.near, nrhs, g);
                total += nrhs;
            }
            report(name_of(k), bad, total);
        }
    }
    for (int nmax : {64, 216, 512}) {
        Pairs<T> ps;
        build_pairs(ps, 96, 16, nmax, g);
        std::printf("%s, pairs of up to %d rows and columns: %zu pairs\n", type, nmax, ps.pairs.size());
        long bad = 0, total = 0;
        for (int nrhs : counts) {
            bad += check_pairs(ps, nrhs, g);
            total += nrhs;
        }
        report("pairs", bad, total);
    }
    return failures;
}

// Each kernel's time (cudaEvent, the mean of five launches after a first)
template<typename T>
float time_ms(const std::function<void()>& launch) {
    cudaStream_t stream = Context::instance().stream();
    cudaEvent_t e0, e1;
    check_cuda(cudaEventCreate(&e0), "event");
    check_cuda(cudaEventCreate(&e1), "event");
    launch();
    float total = 0.0f;
    for (int rep = 0; rep < 5; ++rep) {
        check_cuda(cudaEventRecord(e0, stream), "event");
        launch();
        check_cuda(cudaEventRecord(e1, stream), "event");
        check_cuda(cudaEventSynchronize(e1), "event");
        float ms = 0.0f;
        check_cuda(cudaEventElapsedTime(&ms, e0, e1), "event");
        total += ms;
    }
    cudaEventDestroy(e0);
    cudaEventDestroy(e1);
    return total / 5.0f;
}

template<typename T>
void time_kernels() {
    std::mt19937_64 g(4242);
    // leaf-like: 512 boxes of ~216 points, 20-27 near blocks from a pool of
    // 512 matrices (beyond L2); coupling blocks of the skeletons, ~64 rows;
    // 2048 pairs of 216 x 216 from a pool of 512
    Level<T> near, coupling;
    build_level(near, 512, 200, 216, 0.25, 20, 27, true, g, 512);
    build_level(coupling, 1024, 200, 256, 0.25, 60, 100, false, g, 2048);
    Pairs<T> ps;
    build_pairs(ps, 2048, 200, 216, g, 512);
    std::printf("timing: near %zu boxes (max rows %d, %zu blocks), coupling %zu boxes (max rows %d, %zu blocks), "
                "%zu pairs\n", near.boxes.size(), near.max_rows, near.blocks.size(), coupling.boxes.size(),
                coupling.max_rows, coupling.blocks.size(), ps.pairs.size());
    for (int nrhs : {1, 16}) {
        std::vector<void*> allocs;
        auto level_buffers = [&](const Level<T>& lv, T*& vec, T* src[3]) {
            vec = upload(random_vector<T>(static_cast<size_t>(lv.ld) * nrhs, g), allocs);
            for (int a = 0; a < 3; ++a) {
                src[a] = upload(random_vector<T>(static_cast<size_t>(std::max<int64_t>(lv.src[a].points, 1)) * nrhs,
                                                 g), allocs);
            }
        };
        T *vn, *vc, *sn[3], *sc[3];
        level_buffers(near, vn, sn);
        level_buffers(coupling, vc, sc);
        T* x = upload(random_vector<T>(static_cast<size_t>(ps.ldx) * nrhs, g), allocs);
        T* p = upload(random_vector<T>(static_cast<size_t>(ps.ldp) * nrhs, g), allocs);
        const char* cols = nrhs == 1 ? " 1 column " : "16 columns";
        std::printf("  %s  upward            %8.3f ms\n", cols,
                    time_ms<T>([&] { launch_level(Kernel::upward, coupling, false, vc, sc, nrhs); }));
        std::printf("  %s  downward          %8.3f ms\n", cols,
                    time_ms<T>([&] { launch_level(Kernel::downward, coupling, false, vc, sc, nrhs); }));
        std::printf("  %s  accumulate near   %8.3f ms\n", cols,
                    time_ms<T>([&] { launch_level(Kernel::accumulate, near, true, vn, sn, nrhs); }));
        std::printf("  %s  accumulate coupl. %8.3f ms\n", cols,
                    time_ms<T>([&] { launch_level(Kernel::accumulate, coupling, false, vc, sc, nrhs); }));
        std::printf("  %s  pairs             %8.3f ms\n", cols, time_ms<T>([&] { launch_pairs(ps, x, p, nrhs); }));
        for (void* q : allocs) cudaFree(q);
    }
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Context::instance().activate();
    if (argc > 1 && std::string(argv[1]) == "--time") {
        time_kernels<double>();
        MPI_Finalize();
        return 0;
    }
    std::printf("matvec columns: each column of a many-column launch against that column alone, bitwise\n");
    int failures = 0;
    failures += test_type<double>("double", {2, 3, 8, 11, 17}, 12345);
    failures += test_type<dcomplex>("complex", {2, 3, 4, 5, 9}, 67890);
    std::printf("%s\n", failures == 0 ? "all columns bitwise equal" : "FAILED");
    MPI_Finalize();
    return failures == 0 ? 0 : 1;
}
