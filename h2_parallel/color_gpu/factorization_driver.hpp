#pragma once
// The GPU side of the Color factorization loop (H2_use_gpu), shared by the
// full-grid loop (butterfly_factorization.hpp) and the unstructured one
// (color_unstructured/factorization_impl.hpp).  Per level: the device box
// path (level_eliminator.hpp) and its transports, the level end, the level
// transition, and the root.  The loops keep their host code and their own
// wave schedules; each call here replaces the host step it names.
//
// A tree of occupied boxes (`occupancy`, the unstructured backend) keeps
// empty boxes in the local slabs; the eliminator skips them, and its
// transports follow the unstructured host transport's protocol (every rank
// in the assisting exchange, generators also to the assisting boxes'
// requesters).  The loop's host_transport callback sets the host
// transport's flags when the device exchange is off.

#ifdef H2_HAVE_GPU

#include "level_eliminator.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace fmm {
namespace gpu {

template<typename CoordType, typename DataType, typename KernelType>
class ColorGpuDriver {
public:
    using Tree = ParallelTree<CoordType, DataType>;
    using Level = TreeLevel<CoordType, DataType>;
    using Box = BoxData<CoordType, DataType>;
    using Duration = std::chrono::high_resolution_clock::duration;

    struct Options {
        double tolerance = 0.0;
        FactorizationMethod method = FactorizationMethod::LU;
        int use_sketch = 0;
        int lazy_schur = 0;
        bool is_symmetric = true;
        bool is_hermitian = false;
        bool occupancy = false;
    };

    ColorGpuDriver(Tree* tree, KernelType* kernel, const Options& options)
        : tree_(tree), kernel_(kernel), opt_(options) {}

    // ---- where the blocks of a level go next (tested where levels start)

    // The device box path runs level `lvl`.
    bool box_path_runs(int lvl) const {
        if (!color_gpu_enabled() || lvl <= 1 || tree_->level_uses_CA(lvl)) return false;
        if (!tree_->levels[static_cast<size_t>(lvl)].is_process_active) return false;
        const bool streamed =
            opt_.use_sketch == 2 && opt_.is_symmetric && !opt_.is_hermitian && tree_->id_proxy_mode != 2;
        if (!streamed || opt_.lazy_schur == 0) return false;
        return level_eliminator_would_run(tree_, lvl, kernel_, opt_.method);
    }
    // Level 1 is not eliminated: its transition builds the root, which is
    // factored on the device when level 1 runs on one rank.
    bool root_path_runs() const {
        if (!color_gpu_enabled() || tree_->num_levels < 2) return false;
        const auto& l1 = tree_->levels[1];
        const auto& l0 = tree_->levels[0];
        if (!l1.is_process_active || l1.num_active_processes != 1) return false;
        if (!l0.is_process_active || l0.num_active_processes != 1) return false;
        if (opt_.method != FactorizationMethod::LU) return false;
        return level_eliminator_would_run(tree_, 1, kernel_, opt_.method);
    }
    bool keeps_blocks_of(int lvl) const {
        return lvl == 0 ? true : (lvl == 1 ? root_path_runs() : box_path_runs(lvl));
    }

    // ---- a Color level that eliminates

    // Level start (all ranks of level_comm): the device box path when it
    // covers the level (streamed sketches, lazy far fill, a symmetric
    // kernel).  Its sketches run on the host when a box of any rank is
    // wider than the device sketch's packed lists take.  `announce`: this
    // rank prints the decision.
    void start_level(int lvl, bool use_streamed_level, MPI_Comm level_comm, bool announce) {
        if (!(color_gpu_enabled() && use_streamed_level && lazy_far_field_mode() == LazyFarFieldMode::LAZY &&
              opt_.is_symmetric && !opt_.is_hermitian)) {
            return;
        }
        int fits = 1;
        for (const Box& box : tree_->levels[static_cast<size_t>(lvl)].local_boxes) {
            if (box.num_points > kEntryDestMask) fits = 0;
        }
        MPI_Allreduce(MPI_IN_PLACE, &fits, 1, MPI_INT, MPI_MIN, level_comm);
        std::string reason;
        level_ = make_level_eliminator(tree_, lvl, kernel_, opt_.tolerance, opt_.method, &reason,
                                       std::move(blocks_), opt_.occupancy, fits != 0);
        level_index_ = lvl;
        if (announce) {
            std::cout << "  GPU box path: "
                      << (level_ ? std::string(tensor_core_gemm() ? "on (FP64 tensor-core GEMMs)" : "on")
                                 : "off (" + reason + ")")
                      << (level_ && !fits ? ", host sketches (boxes over " + std::to_string(kEntryDestMask) +
                                                " points)"
                                          : std::string())
                      << std::endl;
        }
    }
    bool on() const { return level_ != nullptr; }
    // one active rank: the level has nothing to transport
    bool local() const {
        return level_ != nullptr && tree_->levels[static_cast<size_t>(level_index_)].num_active_processes == 1;
    }

    // A transport of a device level on several ranks: the device exchange,
    // or the generators of the last wave, the host transport
    // (host_transport(&installed) installs the remote generators and returns
    // its MPI time), and the device updates from the installed generators.
    template<typename HostTransport>
    Duration transport(PendingFactorUpdates<DataType>& pending, HostTransport&& host_transport) {
        if (level_->device_exchange()) return level_->exchange();
        const auto t_emit = std::chrono::steady_clock::now();
        level_->emit_generators(pending);
        const auto t_exchange = std::chrono::steady_clock::now();
        std::vector<int64_t> installed;
        const Duration duration = host_transport(&installed);
        eliminator_stats().emit += std::chrono::duration<double>(t_exchange - t_emit).count();
        eliminator_stats().exchange +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_exchange).count();
        level_->receive_remote(installed);
        return duration;
    }

    // Box region, owner pass, mirror and finalize of a wave on the device;
    // returns its number of boundary boxes.  The loop records the wave's
    // boxes as eliminated.
    int eliminate_wave(const std::vector<int64_t>& wave, int counter) {
        return level_->eliminate_wave(wave, counter);
    }

    // Level end, after the final transport (all ranks of level_comm): the
    // factors on the host, the blocks kept for the transition, and whether
    // the solve's factors stay on the device (only if every rank of the
    // level still holds all of its own).
    void finish_level(int lvl, const PendingFactorUpdates<DataType>& pending, MPI_Comm level_comm, bool announce) {
        if (level_) {
            level_->finish();
            if (local() && (!pending.replace_blocks.empty() || !pending.accumulated_deltas.empty() ||
                            !pending.generators.empty())) {
                throw std::runtime_error("device level left updates for other ranks");
            }
        }
        if (device_solve_enabled() && device_solve_keep()) {
            int kept = level_ && level_->solve_factors_kept() ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &kept, 1, MPI_INT, MPI_MIN, level_comm);
            if (level_) level_->commit_solve_factors(kept != 0);
            if (announce) {
                std::printf("  [gpu] level %d solve factors: %s\n", lvl,
                            kept ? "kept on the device" : "to be uploaded at the first solve");
                std::fflush(stdout);
            }
        }
    }

    // Level 1 (active, not eliminated) adopts the blocks of the level-2
    // transition; its own transition builds the root on the device.
    void adopt_level(int lvl) {
        if (!blocks_) return;
        std::string reason;
        level_ = make_level_eliminator(tree_, lvl, kernel_, opt_.tolerance, opt_.method, &reason,
                                       std::move(blocks_), opt_.occupancy);
        level_index_ = lvl;
        level_->adopt_without_elimination();
    }

    // ---- the level transition

    // Device transition: this level's blocks never leave the device, and the
    // parent's go to the next level's eliminator (or to the host when that
    // level runs there).  Each rank builds the parents of its own boxes, with
    // its copies of the blocks shared with other ranks; at a process
    // reduction they then go to the host, for their new owner.
    // `build_structure()` returns the parent boxes without their blocks.
    // Returns whether the parents were built here; otherwise the level's
    // blocks are back in the host BoxData for the host transition.
    template<typename StructureBuilder>
    bool transition(const Level& level, const Level& parent_level, std::vector<Box>& parents,
                    StructureBuilder&& build_structure, bool announce) {
        bool on_device = false;
        const bool reduction_ahead = parent_level.num_active_processes != level.num_active_processes;
        if (level_ && level.is_process_active && level_->can_build_parent() &&
            (reduction_ahead || parent_level.is_process_active) && opt_.is_symmetric && !opt_.is_hermitian) {
            parents = build_structure();
            blocks_ = level_->build_parent(parents);
            level_.reset();
            if (reduction_ahead || !keeps_blocks_of(level_index_ - 1)) {
                download_level_blocks(*blocks_, parents);
                blocks_.reset();
            }
            on_device = true;
            if (announce) {
                const auto& e = eliminator_stats();
                std::printf("  [gpu] level %d device transition: %.2f s (%lld chunks: plan %.2f [restore %.2f, %.2f GB], "
                            "blocks %.2f, P %.2f, fill %.2f), %lld fill GEMMs (%.1f GFLOP), blocks %s\n",
                            level_index_, e.transition, static_cast<long long>(e.tr_chunks), e.tr_plan,
                            e.tr_restore, e.tr_restore_bytes / 1e9, e.tr_blocks, e.tr_p, e.tr_fill,
                            static_cast<long long>(e.transition_fill_gemms), e.transition_flops / 1e9,
                            blocks_ ? "kept on the device" : "copied to the host");
                std::fflush(stdout);
            }
        }
        if (level_) {
            level_->download_blocks();  // host transition
            level_.reset();
        }
        return on_device;
    }

    // ---- the root

    // The level-1 transition assembled the root block on the device.
    bool holds_root_block() const { return blocks_ != nullptr; }

    // The root's LU on the device: of the block the level-1 transition
    // assembled there, or of the host-assembled block (multi-rank runs)
    // when the heap has room.  Returns false when the host factors it.  A
    // symmetric problem factors the block's symmetric part (as the boxes'
    // X_RR, see compute_and_modify).
    bool factor_root(Box& root, bool announce) {
        const int64_t n = root.num_points;
        const bool symmetric = opt_.is_symmetric && !opt_.is_hermitian;
        if (blocks_) {
            if (announce) std::cout << "  Schur complement size: " << n << " × " << n << std::endl;
            factor_root_on_device(*blocks_, root, symmetric);
            blocks_.reset();
            if (announce) std::cout << "  ✓ Root LU factorization complete (GPU)" << std::endl;
            return true;
        }
        if (color_gpu_enabled() && opt_.method == FactorizationMethod::LU &&
            factor_host_root_on_device(root, symmetric)) {
            if (announce) {
                std::cout << "  Schur complement size: " << n << " × " << n << std::endl;
                std::cout << "  ✓ Root LU factorization complete (GPU)" << std::endl;
            }
            return true;
        }
        return false;
    }

    // Statistics of the level's device box path on this rank, if a wave ran
    // there (H2_PHASE_REPORT=1 resets them per level); returns whether it
    // printed.
    static bool report_level(int lvl) {
        const auto& e = eliminator_stats();
        if (e.boxes == 0) return false;
        std::printf(
            "  [gpu] level %d box path (rank 0): begin %.2f s, sketch %.2f s, ID %.2f s, plan %.2f s, "
            "device %.2f s, download %.2f s, store %.2f s, finish %.2f s | up %.2f GB, "
            "down %.2f GB | owner %lld gemms in %lld batches, %lld new targets | heap peak %.2f GB, %lld reclaims\n",
            lvl, e.begin, e.sketch, e.id, e.plan, e.device, e.download, e.store, e.finish,
            e.bytes_up / 1e9, e.bytes_down / 1e9, static_cast<long long>(e.owner_gemms),
            static_cast<long long>(e.owner_batches), static_cast<long long>(e.new_targets),
            static_cast<double>(e.heap_peak) / 1e9, static_cast<long long>(e.heap_reclaims));
        std::printf("  [gpu] level %d detail: sketch plan %.2f s, sketch device %.2f s "
                    "(meta build %.2f, upload %.2f, rows %.2f, stored %.2f, P %.2f, fill %.2f, "
                    "ID %.2f, ranks down %.2f), background copies: busy %.2f s, level-end wait %.2f s\n",
                    lvl, e.sketch_plan, e.sketch_gpu, e.sk_meta, e.sk_upload, e.sk_rows,
                    e.sk_stored, e.sk_p, e.sk_fill, e.sk_id, e.sk_download, e.finish_store, e.finish_sources);
        std::printf("  [gpu] level %d background copier: busy %.2f s, of which waiting for device data "
                    "%.2f s\n", lvl, e.finish_store, e.copier_wait);
        std::printf("  [gpu] level %d elimination device: fills %.2f, X_RR/X_SR %.2f, LU %.2f, X_NR %.2f, "
                    "solves %.2f, Schur+near %.2f, owner targets %.2f, owner GEMMs %.2f s\n",
                    lvl, e.el[0], e.el[1], e.el[2], e.el[3], e.el[4], e.el[5], e.el[6], e.el[7]);
        std::printf("  [gpu] level %d host: plan boxes %.2f s, owner pass %.2f s (overlapped), launch %.2f s, "
                    "sketch wait %.2f s, heap-reclaim wait %.2f s, exchange-buffer wait %.2f s | "
                    "GF/s: owner %.0f, solves %.0f\n",
                    lvl, e.plan_boxes, e.plan_owner, e.launch, e.sk_wait, e.reclaim_wait,
                    e.exchange_wait,
                    e.owner_flops / std::max(e.el[7], 1e-9) / 1e9,
                    e.solve_flops / std::max(e.el[4], 1e-9) / 1e9);
        if (e.remote_generators > 0 || e.exchange > 0.0) {
            const auto& tt = transport_timers();
            std::printf("  [gpu] level %d ranks (%s memory): emit %.2f s, exchange %.2f s [sizes %.2f, payload %.2f "
                        "(serialize %.2f), deserialize %.2f, assisting %.2f, install %.2f; sent %.2f GB, "
                        "received %.2f GB], receive %.2f s (%lld remote generators, %lld buffers "
                        "outside the exchange arena)\n",
                        lvl, e.device_exchange ? "device" : "host", e.emit, e.exchange,
                        tt.sizes, tt.payload, tt.serialize,
                        tt.deserialize, tt.assisting, tt.install, tt.bytes_sent / 1e9,
                        tt.bytes_received / 1e9, e.remote, static_cast<long long>(e.remote_generators),
                        static_cast<long long>(e.exchange_fallbacks));
        }
        std::fflush(stdout);
        return true;
    }

private:
    Tree* tree_;
    KernelType* kernel_;
    Options opt_;
    int level_index_ = -1;
    std::unique_ptr<LevelEliminatorBase<CoordType, DataType>> level_;  // this level's device box path
    std::unique_ptr<DeviceLevelBlocks<DataType>> blocks_;              // the next level's blocks
};

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
