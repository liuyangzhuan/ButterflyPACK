#pragma once
// H2_GPU_WAVE_TRACE=1: a host timeline of the waves of each device level
// (Color or CA), built from the eliminator's per-phase totals taken around
// each wave.  The level's first rank prints it at the level end: the level
// start, every wave (its host phases and its device times), the time between
// waves (the loop, and a Color level's transports), and the level finish.
//
// Per wave, the host runs: the sketch plan, the sketch launch, a wait for the
// previous wave's elimination and this wave's sketch and ID (sketch wait),
// the ID results into the boxes (ID store), the elimination plan, its
// launch, the owner pass (plan and launch), and the host store.  The device
// meanwhile runs the previous wave's elimination, then this wave's sketch and
// ID; it is idle from the end of the ID until this wave's first launch.

#ifdef H2_HAVE_GPU

#include "level_eliminator.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace fmm {
namespace gpu {

inline bool wave_trace_enabled() {
    static const bool enabled = [] {
        const char* v = std::getenv("H2_GPU_WAVE_TRACE");
        return v != nullptr && std::atoi(v) != 0;
    }();
    return enabled;
}

class WaveTrace {
public:
    using clock = std::chrono::steady_clock;

    // Level start (start_level / start_ca_level): entry and exit.
    void start_begin() {
        records_.clear();
        t_level_ = clock::now();
    }
    void start_end() {
        start_s_ = seconds(t_level_, clock::now());
        last_end_ = clock::now();
    }
    void wave_begin() {
        before_ = eliminator_stats();
        t_wave_ = clock::now();
        gap_ = seconds(last_end_, t_wave_);
    }
    void wave_end(int seq, size_t boxes) {
        Record r;
        r.seq = seq;
        r.boxes = static_cast<int>(boxes);
        r.gap = gap_;
        r.wall = seconds(t_wave_, clock::now());
        r.d = Delta::of(eliminator_stats(), before_);
        records_.push_back(r);
        last_end_ = clock::now();
    }
    void finish_begin() {
        t_finish_ = clock::now();
        tail_ = seconds(last_end_, t_finish_);
        before_ = eliminator_stats();
    }
    void finish_end() {
        finish_s_ = seconds(t_finish_, clock::now());
        finish_d_ = Delta::of(eliminator_stats(), before_);
    }
    bool empty() const { return records_.empty(); }

    // This rank's totals, for the spread over the level's ranks: waves (wall),
    // boxes, sketch wait, ID store, plan, device sketch, device elimination.
    static constexpr int kSummary = 7;
    std::array<double, kSummary> summary() const {
        Delta sum;
        double wall = 0.0, boxes = 0.0;
        for (const Record& r : records_) {
            sum.add(r.d);
            wall += r.wall;
            boxes += r.boxes;
        }
        return {wall, boxes, sum.sk_wait, sum.id, sum.plan, sum.sk_device(), sum.el_sum() + finish_d_.el_sum()};
    }

    void print(int lvl, int rank) const {
        if (records_.empty()) return;
        Delta sum;
        double wall = 0.0, gaps = 0.0;
        for (size_t i = 0; i < records_.size(); ++i) {
            sum.add(records_[i].d);
            wall += records_[i].wall;
            if (i > 0) gaps += records_[i].gap;
        }
        auto ms = [](double s) { return 1e3 * s; };
        auto other = [](const Delta& d, double w) { return w - d.sketch - d.id - d.plan - d.launch - d.store; };
        std::printf("  [gpu] level %d wave trace (rank %d), ms: start %.0f, before wave 0 %.0f | %zu waves %.0f "
                    "(sketch plan %.0f, sketch launch %.0f, sketch wait %.0f, ID store %.0f [loop %.0f], plan %.0f, "
                    "launch %.0f, store %.0f, other %.0f) | between waves %.0f | before the finish %.0f, finish %.0f\n",
                    lvl, rank, ms(start_s_), ms(records_[0].gap), records_.size(), ms(wall), ms(sum.sketch_plan),
                    ms(sum.sketch_gpu - sum.sk_wait), ms(sum.sk_wait), ms(sum.id), ms(sum.id_loop), ms(sum.plan),
                    ms(sum.launch), ms(sum.store), ms(other(sum, wall)), ms(gaps), ms(tail_), ms(finish_s_));
        std::printf("  [gpu] level %d wave trace (rank %d), device ms: sketch %.0f (upload %.0f, rows %.0f, stored %.0f, "
                    "P %.0f, fill %.0f, ID %.0f, download %.0f), elimination %.0f\n",
                    lvl, rank, ms(sum.sk_device()), ms(sum.sk_upload), ms(sum.sk_rows), ms(sum.sk_stored),
                    ms(sum.sk_p), ms(sum.sk_fill), ms(sum.sk_id), ms(sum.sk_download),
                    ms(sum.el_sum() + finish_d_.el_sum()));
        std::printf("  [gpu]   wave  seq boxes |   gap   wall | skplan sklnch skwait idstor   plan launch  store  "
                    "other | dev sketch     ID  elim   (rank %d)\n", rank);
        for (size_t i = 0; i < records_.size(); ++i) {
            const Record& r = records_[i];
            // a wave's elimination is timed when the next wave (or the finish) collects it
            const double elim = i + 1 < records_.size() ? records_[i + 1].d.el_sum() : finish_d_.el_sum();
            std::printf("  [gpu]   %4zu %4d %5d | %5.1f %6.1f | %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f | "
                        "%10.1f %6.1f %5.1f\n",
                        i, r.seq, r.boxes, ms(i > 0 ? r.gap : 0.0), ms(r.wall), ms(r.d.sketch_plan),
                        ms(r.d.sketch_gpu - r.d.sk_wait), ms(r.d.sk_wait), ms(r.d.id), ms(r.d.plan), ms(r.d.launch),
                        ms(r.d.store), ms(other(r.d, r.wall)), ms(r.d.sk_device()), ms(r.d.sk_id), ms(elim));
        }
        std::fflush(stdout);
    }

private:
    // the eliminator totals a wave adds (seconds)
    struct Delta {
        double sketch = 0, sketch_plan = 0, sketch_gpu = 0, sk_wait = 0, sk_meta = 0;
        double sk_upload = 0, sk_rows = 0, sk_stored = 0, sk_p = 0, sk_fill = 0, sk_id = 0, sk_download = 0;
        double el[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        double id = 0, id_loop = 0, plan = 0, launch = 0, store = 0;
        static Delta of(const EliminatorStats& a, const EliminatorStats& b) {
            Delta d;
            d.sketch = a.sketch - b.sketch;
            d.sketch_plan = a.sketch_plan - b.sketch_plan;
            d.sketch_gpu = a.sketch_gpu - b.sketch_gpu;
            d.sk_wait = a.sk_wait - b.sk_wait;
            d.sk_meta = a.sk_meta - b.sk_meta;
            d.sk_upload = a.sk_upload - b.sk_upload;
            d.sk_rows = a.sk_rows - b.sk_rows;
            d.sk_stored = a.sk_stored - b.sk_stored;
            d.sk_p = a.sk_p - b.sk_p;
            d.sk_fill = a.sk_fill - b.sk_fill;
            d.sk_id = a.sk_id - b.sk_id;
            d.sk_download = a.sk_download - b.sk_download;
            for (int i = 0; i < 8; ++i) d.el[i] = a.el[i] - b.el[i];
            d.id = a.id - b.id;
            d.id_loop = a.id_loop - b.id_loop;
            d.plan = a.plan - b.plan;
            d.launch = a.launch - b.launch;
            d.store = a.store - b.store;
            return d;
        }
        void add(const Delta& o) {
            sketch += o.sketch;
            sketch_plan += o.sketch_plan;
            sketch_gpu += o.sketch_gpu;
            sk_wait += o.sk_wait;
            sk_meta += o.sk_meta;
            sk_upload += o.sk_upload;
            sk_rows += o.sk_rows;
            sk_stored += o.sk_stored;
            sk_p += o.sk_p;
            sk_fill += o.sk_fill;
            sk_id += o.sk_id;
            sk_download += o.sk_download;
            for (int i = 0; i < 8; ++i) el[i] += o.el[i];
            id += o.id;
            id_loop += o.id_loop;
            plan += o.plan;
            launch += o.launch;
            store += o.store;
        }
        double sk_device() const { return sk_upload + sk_rows + sk_stored + sk_p + sk_fill + sk_id + sk_download; }
        double el_sum() const {
            double t = 0;
            for (double x : el) t += x;
            return t;
        }
    };
    struct Record {
        int seq = 0, boxes = 0;
        double gap = 0, wall = 0;
        Delta d;
    };
    static double seconds(clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    }

    std::vector<Record> records_;
    EliminatorStats before_;
    clock::time_point t_level_, t_wave_, t_finish_, last_end_;
    double start_s_ = 0, gap_ = 0, tail_ = 0, finish_s_ = 0;
    Delta finish_d_;
};

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
