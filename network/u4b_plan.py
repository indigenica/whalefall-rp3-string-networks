"""U4b (ESPEC_U4.md): choose the lattice size from the machine (RAM, cores) and print the plan with estimated cost.
Memory model measured in the cloud (peak RSS, U2/U3 driver with the slab distance transform): U(1) ~ 82 B/cell, RP3 ~ 206 B/cell,
+0.15 GB.  Time model (2 cores, OpenMP): U(1) 99 s and RP3 267 s at N = 192 (tau_end = 96), scaling ~ N^4 (cells x steps), ~1/cores^0.8.
Usage: python u4b_plan.py [--ram-frac 0.75] [--N N]      (prints JSON; produce uses the same function)"""
import json, math, os, sys
from pathlib import Path

B_CELL = {0: 82.0, 1: 206.0}; BASE_GB = 0.15
T192 = {0: 99.0, 1: 267.0}; P_TIME = 4.0   # N^3 cells x N steps (tau_end = N/2); conservative at large N
CANDIDATES = (1024, 896, 768, 640, 576, 512, 448, 384, 320, 256)


def machine():
    ram = None
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"): ram = int(line.split()[1]) / 1e6; break
    except OSError:
        pass
    return dict(ram_available_GB=ram, cores=os.cpu_count())


def mem_GB(model, N): return BASE_GB + B_CELL[model] * N ** 3 / 1e9


BENCH = Path(__file__).resolve().parent / "results" / "u4" / "bench_threads.json"


def time_h(model, N, cores):
    """with results/u4/bench_threads.json (bench_threads.py on this machine): ~2N leapfrog steps + 30 measurements, both
    scaled as N^3 from N = 256 at the fastest thread count.  Otherwise the cloud model (it underestimated the WSL machine
    ~8x for U(1) at N = 768: 2.4 h predicted for four runs, 5.2-7 h measured for each)."""
    if BENCH.exists():
        b = json.loads(BENCH.read_text()); key = "RP3" if model else "U1"
        best = min(b["results"].values(), key=lambda r: r[key]["step_s"])[key]
        return ((2 * N) * best["step_s"] + 30 * best["measure_s"]) * (N / 256) ** 3 / 3600
    return T192[model] * (N / 192) ** P_TIME * (2 / max(cores, 1)) ** 0.8 / 3600


def plan(ram_frac=0.80, N_force=None):
    m = machine(); budget = (m["ram_available_GB"] or 8.0) * ram_frac; cores = m["cores"] or 2
    pick = lambda model: next((N for N in CANDIDATES if mem_GB(model, N) <= budget and (N_force is None or N <= N_force)), None)
    N1, N2 = pick(0), pick(1)
    out = dict(machine=m, budget_GB=budget,
               phase1=dict(model="U1", N=N1, ell_end=math.log(N1 / 2) if N1 else None, runs=4,
                           mem_GB=mem_GB(0, N1) if N1 else None, hours_total=4 * time_h(0, N1, cores) if N1 else None),
               phase2=dict(model="RP3 (+U1 at the same N if different)", N=N2, ell_end=math.log(N2 / 2) if N2 else None, runs=2,
                           mem_GB=mem_GB(1, N2) if N2 else None,
                           hours_total=(2 * time_h(1, N2, cores) + (2 * time_h(0, N2, cores) if N2 != N1 else 0)) if N2 else None))
    out["note"] = ("ell_end = ln(N/2): 256 -> 4.85, 512 -> 5.55, 768 -> 5.95, 1024 -> 6.24. Gates at ell_ref = ell_end - 0.3. "
                   "If phase-1 N <= 256 the test cannot improve on U3b.")
    return out


if __name__ == "__main__":
    a = sys.argv
    rf = float(a[a.index("--ram-frac") + 1]) if "--ram-frac" in a else 0.80
    Nf = int(a[a.index("--N") + 1]) if "--N" in a else None
    print(json.dumps(plan(rf, Nf), indent=1))
