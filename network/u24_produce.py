"""U24 (docs/preregistros/ESPEC_U24.md) and U25 (ESPEC_U25.md): drain-sensitivity and time-step networks, 384^3.
Protocol of U11 at N = 384, tau_end = 192, sfac = 1, xi_i = 0.5, ell_i = 1, 45 snapshots, spectra every 3.
U24 variants: D0 (ell_d 3.8, Gamma_d m_r), Dearly (3.4), Dlate (4.1), Dweak (m_r/2), Dstrong (2 m_r); dtau = 0.25.
U25: D0 with dtau = 0.125 and tau0 / snapshot times on the 0.25 grid (its dtau = 0.25 partner is U24 D0).
Order: per seed (61, 62, 63) the ten U24 networks, then the two U25 ones; runs already present are skipped.
    nohup python -B network/u24_produce.py > network/results/u24/produce.out 2>&1 &
    python -B network/u24_produce.py --plan          # list the queue and what is done
    python -B network/u24_produce.py --only U24      # only one of the two tests (network/u25_produce.py = --only U25)"""
import json, math, sys, time
from pathlib import Path
import lattice as Lt

HERE = Path(__file__).resolve().parent
OUT24, OUT25 = HERE / "results" / "u24", HERE / "results" / "u25"
N, TAU_END, NSNAP, SPEC_EVERY = 384, 192, 45, 3
SEEDS = (61, 62, 63)
VARIANTS = {"D0": (3.8, 1.0), "Dearly": (3.4, 1.0), "Dlate": (4.1, 1.0), "Dweak": (3.8, 0.5), "Dstrong": (3.8, 2.0)}   # ell_d, Gamma_d / m_r


def queue():
    q = []
    for sd in SEEDS:
        for v in VARIANTS:
            for model in (1, 0):
                q.append(("U24", model, sd, v, OUT24 / "runs" / f"{'RP3' if model else 'U1'}_N{N}_s{sd}_{v}.json"))
        for model in (1, 0):
            q.append(("U25", model, sd, "dt0125", OUT25 / "runs" / f"{'RP3' if model else 'U1'}_N{N}_s{sd}_dt0125.json"))
    return q


def main(argv):
    if "--plan" in argv:
        for t, m, sd, v, p in queue(): print(f"{t} {'RP3' if m else 'U1'} s{sd} {v:8s} {'done' if p.exists() else '-'}")
        return
    for d in (OUT24, OUT25): (d / "runs").mkdir(parents=True, exist_ok=True)
    log = open(OUT24 / "produce.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    lg(f"START {time.strftime('%Y-%m-%d %H:%M:%S')} threads={Lt.lib().nthreads()}")
    only = argv[argv.index("--only") + 1] if "--only" in argv else None
    for t, model, sd, v, path in queue():
        if only and t != only: continue
        if path.exists(): lg(f"skip {path.name}"); continue
        m = math.sqrt(Lt.M2[model])
        if t == "U24":
            ld, gf = VARIANTS[v]
            r = Lt.run(model, N, TAU_END, sd, sfac=1.0, xi_init=0.5, dtau=0.25, nsnap=NSNAP, spectra_every=SPEC_EVERY, log=lambda x: None,
                       damp_ell=ld, damp_rate=gf * m)
        else:
            r = Lt.run(model, N, TAU_END, sd, sfac=1.0, xi_init=0.5, dtau=0.125, tau_grid=0.25, nsnap=NSNAP, spectra_every=SPEC_EVERY,
                       log=lambda x: None, damp_ell=3.8)
        r.update(spec=f"docs/preregistros/ESPEC_{t}.md", variant=v)
        tmp = path.with_suffix(".tmp"); tmp.write_text(json.dumps(r, default=float)); tmp.rename(path)     # never a half-written result
        lg(f"DONE {t} {path.name} {r['runtime_s']:.0f}s xi_end={r['rows'][-1]['xi']:.3f}")
    lg(f"END {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main(sys.argv)
