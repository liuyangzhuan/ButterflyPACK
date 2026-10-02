#pragma once
// GPU execution of the deferred X_NN owner pass of a Color wave.
//
// The host records every candidate box's pass (apply_owner_deferred_xnn_
// updates_for_candidate_box with a DeferredXnnOwnerRecord): the GEMM tasks and
// the target blocks, loaded or synthesized as the CPU pass would.  The device
// then runs the GEMMs and the host writes the targets back and emits the
// remote ADDs.  Updates to one target are applied in its recorded (canonical)
// order: sub-batch j holds the j-th task of every target, so no two GEMMs of a
// launch write the same block and no atomics are needed.
//
// Data stays host-resident here: sources and targets travel over PCIe every
// wave.  This is the first milestone of the backend (infrastructure and the
// batching); device-resident level data comes with the box region.

#ifdef H2_HAVE_GPU

#include "gpu_runtime.hpp"

#include <omp.h>

#include <chrono>
#include <cstring>
#include <exception>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fmm {
namespace gpu {

// Per-level totals of the GPU owner pass (seconds and bytes, this rank).
struct OwnerPassStats {
    double record = 0.0;
    double pack = 0.0;
    double transfer_in = 0.0;
    double gemm = 0.0;
    double transfer_out = 0.0;
    double write_back = 0.0;
    double bytes_in = 0.0;
    double bytes_out = 0.0;
    int64_t gemms = 0;
    int64_t launches = 0;
    int64_t chunks = 0;
};

inline OwnerPassStats& owner_pass_stats() {
    static OwnerPassStats stats;
    return stats;
}

namespace owner_detail {

inline double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// Reusable host and device storage of the pass.
struct Workspace {
    DeviceBuffer device;
    PinnedBuffer pinned;
    std::vector<magma_int_t> sizes;
    std::vector<const void*> a_ptrs;
    std::vector<const void*> b_ptrs;
    std::vector<void*> c_ptrs;
};

inline Workspace& workspace() {
    static Workspace ws;
    return ws;
}

// One variable-size GEMM batch laid out in the device metadata block.
struct Launch {
    size_t first = 0;   // first entry in the pointer arrays
    size_t count = 0;
    size_t sizes = 0;   // first entry of its six size arrays (count + 1 each)
    bool accumulate = true;
};

}  // namespace owner_detail

template<typename CoordType, typename DataType, typename KernelType>
void run_owner_pass(
    const std::vector<int64_t>& candidates,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    const std::unordered_set<int64_t>& wave_box_set,
    std::vector<std::vector<DeferredXnnTargetKey>>& mirror_targets,
    std::vector<PendingFactorUpdates<DataType>>& thread_pending,
    int split_threads) {
    using clock = std::chrono::steady_clock;
    using Record = DeferredXnnOwnerRecord<DataType>;
    auto& stats = owner_pass_stats();
    const size_t num_candidates = candidates.size();
    if (num_candidates == 0) {
        return;
    }

    // ---- 1. record every candidate's pass (the CPU pass minus its GEMMs)
    auto t0 = clock::now();
    std::vector<Record> records(num_candidates);
    {
        std::exception_ptr failure;
        std::mutex failure_mutex;
        #pragma omp parallel default(shared)
        {
            const int tid = omp_get_thread_num();
            DeferredXnnOwnerScratch<DataType> scratch;
            scratch.split_threads = split_threads;
            #pragma omp for schedule(dynamic)
            for (int64_t i = 0; i < static_cast<int64_t>(num_candidates); ++i) {
                try {
                    apply_owner_deferred_xnn_updates_for_candidate_box(
                        candidates[static_cast<size_t>(i)], level, kernel, wave_box_set, scratch,
                        mirror_targets[static_cast<size_t>(i)],
                        &thread_pending[static_cast<size_t>(tid)],
                        /*include_ghosts=*/false, /*pair_filter=*/nullptr,
                        &records[static_cast<size_t>(i)]);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(failure_mutex);
                    if (!failure) failure = std::current_exception();
                }
            }
        }
        if (failure) std::rethrow_exception(failure);
    }
    stats.record += owner_detail::seconds_since(t0);

    Context& ctx = Context::instance();
    ctx.activate();
    owner_detail::Workspace& ws = owner_detail::workspace();
    const size_t elem = sizeof(DataType);

    auto target_bytes = [&](const Record& rec) {
        size_t bytes = 0;
        for (const auto& t : rec.targets) bytes += align_up(static_cast<size_t>(t.rows * t.cols) * elem);
        for (const auto& task : rec.tasks)
            if (task.emit_remote_add) bytes += align_up(static_cast<size_t>(task.rows * task.cols) * elem);
        return bytes;
    };

    // Device memory budget of one chunk (targets, deltas, sources, metadata).
    size_t free_bytes = 0, total_bytes = 0;
    check_cuda(cudaMemGetInfo(&free_bytes, &total_bytes), "cudaMemGetInfo");
    const size_t budget = static_cast<size_t>(0.85 * static_cast<double>(free_bytes + ws.device.capacity()));

    size_t begin = 0;
    while (begin < num_candidates) {
        // ---- 2. chunk of candidates whose targets, deltas and sources fit
        size_t end = begin;
        size_t chunk_targets = 0;
        size_t chunk_sources = 0;
        std::unordered_map<int64_t, size_t> source_slot;     // morton -> slot
        std::vector<const typename Record::Task*> source_task;  // one task per slot (pointers, ld, r)
        while (end < num_candidates) {
            const Record& rec = records[end];
            size_t extra_sources = 0;
            std::unordered_set<int64_t> fresh;
            for (const auto& task : rec.tasks)
                if (!source_slot.count(task.source_morton) && fresh.insert(task.source_morton).second)
                    extra_sources += 2 * align_up(static_cast<size_t>(task.ld * task.r) * elem);
            const size_t extra_targets = target_bytes(rec);
            if (end > begin && chunk_targets + chunk_sources + extra_targets + extra_sources > budget) break;
            for (const auto& task : rec.tasks)
                if (source_slot.emplace(task.source_morton, source_task.size()).second) source_task.push_back(&task);
            chunk_targets += extra_targets;
            chunk_sources += extra_sources;
            ++end;
        }
        ++stats.chunks;

        // ---- 3. layout: [sources: temp2, x_nr][targets][deltas]
        t0 = clock::now();
        std::vector<size_t> source_offset(source_task.size());
        size_t data_bytes = 0;
        for (size_t s = 0; s < source_task.size(); ++s) {
            source_offset[s] = data_bytes;
            data_bytes += 2 * align_up(static_cast<size_t>(source_task[s]->ld * source_task[s]->r) * elem);
        }
        const size_t targets_begin = data_bytes;
        std::vector<std::vector<size_t>> target_offset(end - begin);
        std::vector<std::vector<size_t>> delta_offset(end - begin);
        for (size_t c = begin; c < end; ++c) {
            const Record& rec = records[c];
            auto& toff = target_offset[c - begin];
            toff.resize(rec.targets.size());
            for (size_t t = 0; t < rec.targets.size(); ++t) {
                toff[t] = data_bytes;
                data_bytes += align_up(static_cast<size_t>(rec.targets[t].rows * rec.targets[t].cols) * elem);
            }
        }
        const size_t deltas_begin = data_bytes;
        for (size_t c = begin; c < end; ++c) {
            const Record& rec = records[c];
            auto& doff = delta_offset[c - begin];
            doff.assign(rec.tasks.size(), 0);
            for (size_t k = 0; k < rec.tasks.size(); ++k) {
                if (!rec.tasks[k].emit_remote_add) continue;
                doff[k] = data_bytes;
                data_bytes += align_up(static_cast<size_t>(rec.tasks[k].rows * rec.tasks[k].cols) * elem);
            }
        }

        // ---- 4. GEMM launches: sub-batch j = j-th task of every target, then
        //         one independent batch for the remote-ADD deltas
        char* d_data = nullptr;  // assigned after the metadata size is known
        ws.a_ptrs.clear();
        ws.b_ptrs.clear();
        ws.c_ptrs.clear();
        ws.sizes.clear();
        std::vector<owner_detail::Launch> launches;
        std::vector<std::vector<size_t>> by_target;  // task indices per (candidate, target)
        std::vector<std::pair<size_t, size_t>> target_owner;  // (candidate, target) of by_target entry
        for (size_t c = begin; c < end; ++c) {
            const Record& rec = records[c];
            const size_t first = by_target.size();
            by_target.resize(first + rec.targets.size());
            for (size_t t = 0; t < rec.targets.size(); ++t) target_owner.emplace_back(c, t);
            for (size_t k = 0; k < rec.tasks.size(); ++k)
                if (rec.tasks[k].target >= 0) by_target[first + static_cast<size_t>(rec.tasks[k].target)].push_back(k);
        }
        size_t depth = 0;
        for (const auto& list : by_target) depth = std::max(depth, list.size());

        // Pointers are recorded as byte offsets from d_data and converted below.
        struct Entry { size_t a, b, c; magma_int_t m, n, k, ld_a, ld_b, ld_c; };
        std::vector<Entry> entries;
        auto source_ptr = [&](const typename Record::Task& task, bool x_nr, int64_t row_offset) {
            const size_t s = source_slot.at(task.source_morton);
            const size_t base = source_offset[s] +
                (x_nr ? align_up(static_cast<size_t>(task.ld * task.r) * elem) : 0);
            return base + static_cast<size_t>(row_offset) * elem;
        };
        auto push_launch = [&](size_t first_entry, bool accumulate) {
            owner_detail::Launch launch;
            launch.first = first_entry;
            launch.count = entries.size() - first_entry;
            launch.accumulate = accumulate;
            if (launch.count > 0) launches.push_back(launch);
        };
        for (size_t j = 0; j < depth; ++j) {
            const size_t first_entry = entries.size();
            for (size_t e = 0; e < by_target.size(); ++e) {
                if (j >= by_target[e].size()) continue;
                const size_t c = target_owner[e].first;
                const size_t t = target_owner[e].second;
                const auto& task = records[c].tasks[by_target[e][j]];
                const auto& tgt = records[c].targets[t];
                entries.push_back(Entry{source_ptr(task, false, task.a_row_offset),
                                        source_ptr(task, true, task.b_row_offset),
                                        target_offset[c - begin][t],
                                        static_cast<magma_int_t>(task.rows), static_cast<magma_int_t>(task.cols),
                                        static_cast<magma_int_t>(task.r), static_cast<magma_int_t>(task.ld),
                                        static_cast<magma_int_t>(task.ld), static_cast<magma_int_t>(tgt.rows)});
            }
            push_launch(first_entry, true);
        }
        {
            const size_t first_entry = entries.size();
            for (size_t c = begin; c < end; ++c) {
                const Record& rec = records[c];
                for (size_t k = 0; k < rec.tasks.size(); ++k) {
                    const auto& task = rec.tasks[k];
                    if (!task.emit_remote_add) continue;
                    entries.push_back(Entry{source_ptr(task, false, task.a_row_offset),
                                            source_ptr(task, true, task.b_row_offset),
                                            delta_offset[c - begin][k],
                                            static_cast<magma_int_t>(task.rows), static_cast<magma_int_t>(task.cols),
                                            static_cast<magma_int_t>(task.r), static_cast<magma_int_t>(task.ld),
                                            static_cast<magma_int_t>(task.ld), static_cast<magma_int_t>(task.rows)});
                }
            }
            push_launch(first_entry, false);
        }

        // metadata: three pointer arrays, then six size arrays per launch
        const size_t num_entries = entries.size();
        const size_t pointers_begin = align_up(data_bytes);
        const size_t sizes_begin = align_up(pointers_begin + 3 * num_entries * sizeof(void*));
        size_t num_sizes = 0;
        for (auto& launch : launches) {
            launch.sizes = num_sizes;
            num_sizes += 6 * (launch.count + 1);
        }
        const size_t total_bytes_needed = sizes_begin + num_sizes * sizeof(magma_int_t);
        d_data = static_cast<char*>(ws.device.reserve(total_bytes_needed));
        char* h_data = static_cast<char*>(ws.pinned.reserve(total_bytes_needed));

        // pack sources and targets (and the metadata) into the pinned buffer
        #pragma omp parallel for schedule(dynamic)
        for (int64_t s = 0; s < static_cast<int64_t>(source_task.size()); ++s) {
            const auto* task = source_task[static_cast<size_t>(s)];
            const size_t n = static_cast<size_t>(task->ld * task->r) * elem;
            std::memcpy(h_data + source_offset[static_cast<size_t>(s)], task->temp2, n);
            std::memcpy(h_data + source_offset[static_cast<size_t>(s)] + align_up(n), task->x_nr, n);
        }
        #pragma omp parallel for schedule(dynamic)
        for (int64_t c = static_cast<int64_t>(begin); c < static_cast<int64_t>(end); ++c) {
            const Record& rec = records[static_cast<size_t>(c)];
            for (size_t t = 0; t < rec.targets.size(); ++t)
                std::memcpy(h_data + target_offset[static_cast<size_t>(c) - begin][t], rec.targets[t].data.data(),
                            rec.targets[t].data.size() * elem);
        }
        {
            void** h_ptrs = reinterpret_cast<void**>(h_data + pointers_begin);
            for (size_t e = 0; e < num_entries; ++e) {
                h_ptrs[e] = d_data + entries[e].a;
                h_ptrs[num_entries + e] = d_data + entries[e].b;
                h_ptrs[2 * num_entries + e] = d_data + entries[e].c;
            }
            magma_int_t* h_sizes = reinterpret_cast<magma_int_t*>(h_data + sizes_begin);
            for (const auto& launch : launches) {
                magma_int_t* m = h_sizes + launch.sizes;
                magma_int_t* n = m + (launch.count + 1);
                magma_int_t* k = n + (launch.count + 1);
                magma_int_t* la = k + (launch.count + 1);
                magma_int_t* lb = la + (launch.count + 1);
                magma_int_t* lc = lb + (launch.count + 1);
                for (size_t i = 0; i < launch.count; ++i) {
                    const Entry& en = entries[launch.first + i];
                    m[i] = en.m; n[i] = en.n; k[i] = en.k; la[i] = en.ld_a; lb[i] = en.ld_b; lc[i] = en.ld_c;
                }
                m[launch.count] = n[launch.count] = k[launch.count] = 0;
                la[launch.count] = lb[launch.count] = lc[launch.count] = 0;
            }
        }
        stats.pack += owner_detail::seconds_since(t0);

        // ---- 5. upload, GEMMs, download of targets and deltas
        cudaStream_t stream = ctx.stream();
        t0 = clock::now();
        check_cuda(cudaMemcpyAsync(d_data, h_data, deltas_begin, cudaMemcpyHostToDevice, stream), "upload data");
        check_cuda(cudaMemcpyAsync(d_data + pointers_begin, h_data + pointers_begin,
                                   total_bytes_needed - pointers_begin, cudaMemcpyHostToDevice, stream),
                   "upload metadata");
        // the remote-ADD deltas are written with beta = 0; start them from zeros
        // rather than trusting the GEMM never to read uninitialized C
        check_cuda(cudaMemsetAsync(d_data + deltas_begin, 0, data_bytes - deltas_begin, stream), "zero deltas");
        check_cuda(cudaStreamSynchronize(stream), "upload");
        stats.transfer_in += owner_detail::seconds_since(t0);
        stats.bytes_in += static_cast<double>(deltas_begin + (total_bytes_needed - pointers_begin));

        t0 = clock::now();
        {
            void** d_ptrs = reinterpret_cast<void**>(d_data + pointers_begin);
            magma_int_t* d_sizes = reinterpret_cast<magma_int_t*>(d_data + sizes_begin);
            for (const auto& launch : launches) {
                magma_int_t* m = d_sizes + launch.sizes;
                magma_int_t* n = m + (launch.count + 1);
                magma_int_t* k = n + (launch.count + 1);
                magma_int_t* la = k + (launch.count + 1);
                magma_int_t* lb = la + (launch.count + 1);
                magma_int_t* lc = lb + (launch.count + 1);
                gemm_nt_vbatched<DataType>(
                    m, n, k, DataType(1.0),
                    reinterpret_cast<const DataType* const*>(d_ptrs + launch.first), la,
                    reinterpret_cast<const DataType* const*>(d_ptrs + num_entries + launch.first), lb,
                    launch.accumulate ? DataType(1.0) : DataType(0.0),
                    reinterpret_cast<DataType**>(d_ptrs + 2 * num_entries + launch.first), lc,
                    static_cast<magma_int_t>(launch.count), ctx.queue());
                stats.gemms += static_cast<int64_t>(launch.count);
                ++stats.launches;
            }
            check_cuda(cudaStreamSynchronize(stream), "owner-pass GEMMs");
        }
        stats.gemm += owner_detail::seconds_since(t0);

        t0 = clock::now();
        check_cuda(cudaMemcpyAsync(h_data + targets_begin, d_data + targets_begin, data_bytes - targets_begin,
                                   cudaMemcpyDeviceToHost, stream),
                   "download");
        check_cuda(cudaStreamSynchronize(stream), "download");
        stats.transfer_out += owner_detail::seconds_since(t0);
        stats.bytes_out += static_cast<double>(data_bytes - targets_begin);

        // ---- 6. write the targets back and emit the remote ADDs
        t0 = clock::now();
        {
            std::exception_ptr failure;
            std::mutex failure_mutex;
            #pragma omp parallel for schedule(dynamic)
            for (int64_t c = static_cast<int64_t>(begin); c < static_cast<int64_t>(end); ++c) {
                try {
                    Record& rec = records[static_cast<size_t>(c)];
                    for (size_t t = 0; t < rec.targets.size(); ++t) {
                        auto& tgt = rec.targets[t];
                        std::memcpy(tgt.data.data(), h_data + target_offset[static_cast<size_t>(c) - begin][t],
                                    tgt.data.size() * elem);
                        flush_deferred_xnn_target_matrix_from_accumulation(tgt, level);
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lock(failure_mutex);
                    if (!failure) failure = std::current_exception();
                }
            }
            if (failure) std::rethrow_exception(failure);
        }
        for (size_t c = begin; c < end; ++c) {
            const Record& rec = records[c];
            for (size_t k = 0; k < rec.tasks.size(); ++k) {
                const auto& task = rec.tasks[k];
                if (!task.emit_remote_add) continue;
                const EdgeKind edge_kind =
                    (task.kind == DeferredXnnTargetKind::SCHUR) ? EdgeKind::Diag :
                    (task.kind == DeferredXnnTargetKind::NEAR_A_NS) ? EdgeKind::Near : EdgeKind::Far;
                deferred_xnn_accumulate_canonical_edge_delta_from_slice(
                    thread_pending[0], rec.candidate_morton, task.neighbor_morton, edge_kind,
                    reinterpret_cast<const DataType*>(h_data + delta_offset[c - begin][k]),
                    task.rows, 0, task.rows, task.cols);
            }
        }
        stats.write_back += owner_detail::seconds_since(t0);
        begin = end;
    }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
