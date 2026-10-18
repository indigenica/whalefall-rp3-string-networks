"""U4b: quick benchmark of the lattice kernel on this machine (~5 min).  Prints, for U(1) and RP3 at N = 256, the time of
one leapfrog step and of one measurement for several OpenMP thread counts, to choose U4_THREADS for phase 2.
Usage (from the repository):  python -B network/bench_threads.py            (run it when no other large job is running)"""
import json, math, os, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent

if len(sys.argv) > 1 and sys.argv[1] == "--child":
    sys.path.insert(0, str(HERE))
    import numpy as np
    import lattice as Lt
    out = {}
    for model in (0, 1):
        N = 256; m = math.sqrt(Lt.M2[model]); s = 1.0 / (m * N / 2)
        f, p = Lt.initial(model, N, 6.0, 61); scratch = np.empty_like(f); Lt.lib()
        Lt.lib().evolve(model, N, f, p, scratch, 20.0, 0.25, s, 2)                 # warm-up
        t = time.time(); Lt.lib().evolve(model, N, f, p, scratch, 20.5, 0.25, s, 8); step = (time.time() - t) / 8
        t = time.time(); Lt.measure(model, N, f, p, 22.5, s, spectra=True); meas = time.time() - t
        out["RP3" if model else "U1"] = dict(step_s=step, measure_s=meas)
        del f, p, scratch
    print(json.dumps(out)); sys.exit(0)

cores = os.cpu_count() or 2
res = {}
for th in sorted({max(1, cores // 4), max(1, cores // 2), max(1, (3 * cores) // 4), cores}):
    env = dict(os.environ, OMP_NUM_THREADS=str(th))
    r = subprocess.run([sys.executable, "-B", __file__, "--child"], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-2000:]); sys.exit(1)
    res[th] = json.loads(r.stdout.strip().splitlines()[-1])
    print(f"threads={th:3d}  U1: step {res[th]['U1']['step_s']:.3f}s measure {res[th]['U1']['measure_s']:.1f}s   "
          f"RP3: step {res[th]['RP3']['step_s']:.3f}s measure {res[th]['RP3']['measure_s']:.1f}s", flush=True)
best = min(res, key=lambda t: res[t]["RP3"]["step_s"])
print(f"\nRP3 step fastest with {best} threads.  Phase 2: make u4b-phase2 U4_THREADS={best}")
(HERE / "results" / "u4").mkdir(parents=True, exist_ok=True)
(HERE / "results" / "u4" / "bench_threads.json").write_text(json.dumps(dict(cores=cores, results=res, best_RP3=best), indent=1))
