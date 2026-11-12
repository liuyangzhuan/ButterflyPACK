// Round-trip test of the device butterfly layout: build a few levels of
// blocks with known values (including absent local blocks), upload, zero the
// host copies, download, and compare.  Nothing but the layout is exercised.
#include "bf_gpu/bf_upload.hpp"
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <random>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  using bf::gpu::BlockDesc; using bf::gpu::LevelDesc; using bf::gpu::BfDevice;
  // three levels: V (1 x 3), one interior (2 x 3), U (2 x 1); one absent block
  std::vector<LevelDesc> lv(3);
  lv[0].nr = 1; lv[0].nc = 3; lv[0].inc_r = 1; lv[0].inc_c = 1; lv[0].idx_r = 1; lv[0].idx_c = 1;
  lv[1].nr = 2; lv[1].nc = 3; lv[1].inc_r = 1; lv[1].inc_c = 1; lv[1].idx_r = 1; lv[1].idx_c = 1;
  lv[2].nr = 2; lv[2].nc = 1; lv[2].inc_r = 1; lv[2].inc_c = 1; lv[2].idx_r = 1; lv[2].idx_c = 1;
  std::vector<std::vector<std::pair<int,int>>> shp(3);
  shp[0] = {{7,4},{5,4},{6,4}};
  shp[1] = {{4,3},{4,3},{0,0},{4,3},{4,3},{4,3}};   // the third is absent
  shp[2] = {{9,3},{8,3}};
  std::vector<std::vector<std::vector<double>>> store(3);
  std::vector<std::vector<const double*>> up(3);
  std::vector<std::vector<double*>> down(3);
  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<double> uni(-1.0, 1.0);
  for (int l = 0; l < 3; ++l) {
    store[l].resize(shp[l].size());
    for (size_t b = 0; b < shp[l].size(); ++b) {
      size_t n = size_t(shp[l][b].first) * shp[l][b].second;
      store[l][b].resize(n);
      for (auto& v : store[l][b]) v = uni(rng);
      up[l].push_back(n ? store[l][b].data() : nullptr);
      down[l].push_back(n ? store[l][b].data() : nullptr);
    }
  }
  std::vector<std::vector<std::vector<double>>> gold = store;
  BfDevice dev;
  if (!bf::gpu::upload(lv, up, shp, 1, 0, 17, 12, 17, 12, dev)) {
    std::printf("FAIL: upload refused (no arena room)\n"); MPI_Finalize(); return 1;
  }
  std::printf("uploaded %lld elements over %zu blocks\n",
              (long long)dev.elements, dev.blocks.size());
  // the layout must address blocks the way the host does
  const BlockDesc& b11 = dev.block(1, 1, 1);
  if (b11.rows != 4 || b11.cols != 3) {
    std::printf("FAIL: block(1,1,1) is %dx%d, expected 4x3\n", b11.rows, b11.cols);
    MPI_Finalize(); return 1;
  }
  for (auto& l : store) for (auto& b : l) std::fill(b.begin(), b.end(), 0.0);
  bf::gpu::download(dev, down);
  double worst = 0.0;
  for (int l = 0; l < 3; ++l)
    for (size_t b = 0; b < store[l].size(); ++b)
      for (size_t i = 0; i < store[l][b].size(); ++i)
        worst = std::max(worst, std::abs(store[l][b][i] - gold[l][b][i]));
  bf::gpu::release(dev);
  std::printf("%s: worst round-trip difference %.3e; released (resident=%d)\n",
              worst == 0.0 ? "PASS" : "FAIL", worst, (int)dev.resident());
  MPI_Finalize();
  return worst == 0.0 ? 0 : 1;
}
