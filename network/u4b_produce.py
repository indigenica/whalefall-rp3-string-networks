"""U4b production (ESPEC_U4.md).  usage: python u4b_produce.py PHASE [--N N] [--only TAG]
PHASE in 1 (U(1) gates), 2 (RP3 + U(1)), 3 (U5d: U(1) seeds 63, 64 at the given N, ESPEC_U5), 4 (U6: RP3 seeds 63-66 + U(1) 63, 64, spectra every snapshot), 5 (U9: RP3 + U(1) with xi_i = 2, seeds 61, 62), 6 (U10a: drained networks + exchange diagnostic, ESPEC_U10a), 7 (U11: U6 drained, N = 576 fixed, ESPEC_U11), 8 (U10b: dense late starts from above, N = 576 fixed, ESPEC_U10b).  Runs already present in results/u4/runs are skipped, so the job can be restarted.
Phase 2: the RP3 runs go one after the other (each needs most of the RAM); the two U(1) runs then go at the same time, as
two processes with half the threads each, when both fit in the RAM budget (the kernel scales poorly beyond ~20 threads)."""
import json, os, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt
from u4b_plan import plan

OUT = Path(__file__).resolve().parent / "results" / "u4" / "runs"


def sets(phase, N):
    if phase == 1:   # (model, N, tau_end, seed, xi_init)
        return [(0, N, N // 2, 61, 0.25), (0, N, N // 2, 61, 2.0), (0, N, N // 2, 61, 0.5), (0, N, N // 2, 62, 0.5)]
    if phase == 3:   # U5d (ESPEC_U5): two more U(1) seeds at the phase-1 size, for the slope d xi / d ell
        return [(0, N, N // 2, 63, 0.5), (0, N, N // 2, 64, 0.5)]
    if phase == 4:   # U6 (ESPEC_U6): four more RP3 seeds + two U(1) at the phase-2 size, spectra at every snapshot
        return [(1, N, N // 2, sd, 0.5) for sd in (63, 64, 65, 66)] + [(0, N, N // 2, 63, 0.5), (0, N, N // 2, 64, 0.5)]
    if phase == 5:   # U9 (ESPEC_U9): dense initial conditions, RP3 and U(1), seeds 61, 62, spectra at every snapshot
        return [(1, N, N // 2, 61, 2.0), (1, N, N // 2, 62, 2.0), (0, N, N // 2, 61, 2.0), (0, N, N // 2, 62, 2.0)]
    return [(1, N, N // 2, 61, 0.5), (1, N, N // 2, 62, 0.5), (0, N, N // 2, 61, 0.5), (0, N, N // 2, 62, 0.5)]


# U10a (ESPEC_U10a, amendment 1): groups run in this order; a group of two U(1) runs goes as two parallel processes.
# (model, xi_init, seed, damp_ell, exch, tag suffix).  "_dr" = drained until ell = 3.8; "_x" = undrained, with the exchange diagnostic.
U10A_GROUPS = [[(1, 2.0, 61, 3.8, False, "_dr")],
               [(0, 2.0, 61, 3.8, True, "_dr"), (0, 0.5, 61, 3.8, True, "_dr")],
               [(1, 0.5, 61, 3.8, False, "_dr")],
               [(1, 2.0, 61, None, True, "_x")],
               [(0, 2.0, 61, None, True, "_x")],
               [(1, 2.0, 62, 3.8, False, "_dr")],
               [(1, 0.5, 62, 3.8, False, "_dr")],
               [(0, 2.0, 62, 3.8, True, "_dr"), (0, 0.5, 62, 3.8, True, "_dr")]]


# U11 (ESPEC_U11): U6 drained (same seeds, xi_i = 0.5), N = 576 fixed; exchange X only for U(1) (quantitative after V5).
U11_GROUPS = [[(1, 0.5, 61, 3.8, False, "_dr")],
              [(0, 0.5, 61, 3.8, True, "_dr"), (0, 0.5, 62, 3.8, True, "_dr")],
              [(1, 0.5, 62, 3.8, False, "_dr")], [(1, 0.5, 63, 3.8, False, "_dr")],
              [(0, 0.5, 63, 3.8, True, "_dr"), (0, 0.5, 64, 3.8, True, "_dr")],
              [(1, 0.5, 64, 3.8, False, "_dr")], [(1, 0.5, 65, 3.8, False, "_dr")], [(1, 0.5, 66, 3.8, False, "_dr")]]
U11_MIN_FREE_GB = 42.0
# U10b (ESPEC_U10b): dense networks started at ell_i = 3.5 (approach from above), drained until 3.8; 7th field = ell_i.
U10B_GROUPS = [[(1, 3.0, 71, 3.8, False, "_ab", 3.5)],
               [(0, 3.0, 71, 3.8, True, "_ab", 3.5), (0, 3.0, 72, 3.8, True, "_ab", 3.5)],
               [(1, 3.0, 72, 3.8, False, "_ab", 3.5)],
               [(0, 6.0, 73, 3.8, True, "_ab", 3.5)],
               [(1, 6.0, 73, 3.8, False, "_ab", 3.5)]]


def u10a_tag(model, N, sd, xi, suf): return f"{'RP3' if model else 'U1'}_N{N}_s{sd}_xi{xi:g}{suf}"


def u10a(N, only, lg, groups=None, phase=6):
    from u4b_plan import machine, mem_GB
    for grp in (groups or U10A_GROUPS):
        jobs = [(m, xi, sd, de, ex, u10a_tag(m, N, sd, xi, suf), (rest[0] if rest else 1.0)) for m, xi, sd, de, ex, suf, *rest in grp]
        if only: jobs = [j for j in jobs if j[5] == only]
        jobs = [j for j in jobs if not ((OUT / f"{j[5]}.json").exists() and (lg(f"skip {j[5]}") or True))]
        if not jobs: continue
        if phase in (7, 8):                  # U11, U10b: N is fixed; stop instead of shrinking N if the RAM is not there
            free = machine()["ram_available_GB"] or 0.0; need = sum(mem_GB(j[0], N) for j in jobs) + 2.0
            if free < max(need, U11_MIN_FREE_GB if any(j[0] == 1 for j in jobs) else need):
                lg(f"STOP: {free:.1f} GB free, {max(need, U11_MIN_FREE_GB):.1f} GB needed for {[j[5] for j in jobs]}. Free RAM and relaunch (done runs are skipped)."); sys.exit(4)
        if len(jobs) == 2 and not only:
            th = str(max(1, Lt.lib().nthreads() // 2)); lg(f"U10a pair in parallel, {th} threads each: {[j[5] for j in jobs]}")
            procs = [subprocess.Popen([sys.executable, "-B", __file__, str(phase), "--N", str(N), "--only", j[5]], env=dict(os.environ, OMP_NUM_THREADS=th))
                     for j in jobs]
            lg(f"U10a pair exit codes {[q.wait() for q in procs]}"); continue
        for m, xi, sd, de, ex, tag, li in jobs:
            r = Lt.run(m, N, N // 2, sd, sfac=1.0, xi_init=xi, dtau=0.25, nsnap=30, spectra_every=1, log=lg, damp_ell=de, exch=ex, ell_i=li)
            (OUT / f"{tag}.json").write_text(json.dumps(r, default=float)); lg(f"DONE {tag} {r['runtime_s']:.0f}s")


if __name__ == "__main__":
    phase = int(sys.argv[1]); Nf = int(sys.argv[sys.argv.index("--N") + 1]) if "--N" in sys.argv else None
    P = plan(N_force=Nf); N = (P["phase1"] if phase in (1, 3) else P["phase2"])["N"]
    if phase == 3 and Nf:                # U5d: exactly the phase-1 size (other jobs may share the RAM)
        N = Nf
    if phase in (4, 5):                  # U6, U9: the phase-2 size
        N = Nf or P["phase2"]["N"]
    spectra_every = 1 if phase in (4, 5) else 3
    if "--only" in sys.argv and Nf:      # a child of the U(1) pair: exactly the N it was given (the RAM is shared now)
        N = Nf
        spectra_every = 1 if phase in (4, 5) else 3
    OUT.mkdir(parents=True, exist_ok=True)
    log = open(OUT.parent / f"produce_phase{phase}.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    lg("PLAN " + json.dumps(P))
    if phase not in (7, 8) and (not N or N <= 256):
        lg("N <= 256: the test cannot improve on U3b (ESPEC_U4). Stop."); sys.exit(3)
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    if phase == 6:
        N = Nf or P["phase2"]["N"]; u10a(N, only, lg); sys.exit(0)
    if phase == 7:
        N = Nf or 576; u10a(N, only, lg, groups=U11_GROUPS, phase=7); sys.exit(0)
    if phase == 8:
        N = Nf or 576; u10a(N, only, lg, groups=U10B_GROUPS, phase=8); sys.exit(0)
    todo = []
    for model, n, te, sd, xi in sets(phase, N):
        tag = f"{'RP3' if model else 'U1'}_N{n}_s{sd}_xi{xi:g}"
        if only and tag != only: continue
        if (OUT / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
        todo.append((model, n, te, sd, xi, tag))
    u1 = [t for t in todo if t[0] == 0]
    from u4b_plan import mem_GB
    pair = phase in (2, 4, 5) and not only and len(u1) == 2 and 2 * mem_GB(0, N) <= P["budget_GB"]
    for model, n, te, sd, xi, tag in todo:
        if pair and model == 0: continue
        r = Lt.run(model, n, te, sd, sfac=1.0, xi_init=xi, dtau=0.25, nsnap=30, spectra_every=spectra_every, log=lg)
        r["plan"] = P
        (OUT / f"{tag}.json").write_text(json.dumps(r, default=float)); lg(f"DONE {tag} {r['runtime_s']:.0f}s")
    if pair:
        th = str(max(1, Lt.lib().nthreads() // 2)); lg(f"U(1) pair in parallel, {th} threads each: {[t[5] for t in u1]}")
        procs = [subprocess.Popen([sys.executable, "-B", __file__, str(phase), "--N", str(N), "--only", t[5]],
                                  env=dict(os.environ, OMP_NUM_THREADS=th)) for t in u1]
        rc = [q.wait() for q in procs]; lg(f"U(1) pair exit codes {rc}")
