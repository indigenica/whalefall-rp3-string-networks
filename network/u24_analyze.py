"""U24 (docs/preregistros/ESPEC_U24.md), fixed before the data: sensitivity of R_Gamma to the drain (384^3).
Window W = [4.36, 5.21], 10 equispaced points, mask rm2.5, the estimators of U11 (analyze.on_grid: local eps_lat and xi).
Per variant: R_Gamma, R_xi (ratios of seed means), eps and xi means; delta_v = X^v / X^D0 - 1 with a paired bootstrap
(seeds resampled once, the same draw for both models and every variant).  As in U11, a grid point beyond the last snapshot
with a centred difference is NaN and drops out of the means.
    python -B network/u24_analyze.py      ->  network/results/u24/U24_RESULT.json"""
import json, sys
from pathlib import Path
import numpy as np
import analyze as A2

HERE = Path(__file__).resolve().parent; D = HERE / "results" / "u24"
GRID = np.linspace(4.36, 5.21, 10)
VARIANTS = ("D0", "Dearly", "Dlate", "Dweak", "Dstrong"); SEEDS = (61, 62, 63)
NBOOT, TOL, SIG = 4000, 0.05, 0.03
U11_SAME_WINDOW = 1.562            # U11's 576^3 runs on this window (enmienda 1, informative)


def load(rundir, model, tag):
    out = {}
    for sd in SEEDS:
        p = rundir / f"{model}_N384_s{sd}_{tag}.json"
        if p.exists(): out[sd] = json.loads(p.read_text())
    return out


def tables(runs_by_tag, seeds):
    """per tag: arrays [seed, grid] of eps_lat and xi for RP3 and U1, in the order of `seeds`"""
    T = {}
    for tag, (rp, u1) in runs_by_tag.items():
        T[tag] = {k: (A2.on_grid([rp[s] for s in seeds], GRID, k), A2.on_grid([u1[s] for s in seeds], GRID, k)) for k in ("eps_lat", "xi")}
    return T


def quantities(T, idx):
    """R_Gamma, R_xi, eps_RP3, eps_U1, xi_RP3, xi_U1 with the seeds idx (ratios of seed means over the window)"""
    m = lambda a: float(np.nanmean(np.nanmean(a[idx], 0)))
    (eR, eU), (xR, xU) = T["eps_lat"], T["xi"]
    return dict(R_Gamma=m(eR) / m(eU), R_xi=m(xR) / m(xU), eps_RP3=m(eR), eps_U1=m(eU), xi_RP3=m(xR), xi_U1=m(xU))


def paired(T, ref, tags, seeds, nboot=NBOOT, seed=0):
    n = len(seeds); base = np.arange(n); rng = np.random.default_rng(seed)
    est = {t: quantities(T[t], base) for t in tags + [ref]}
    draws = [rng.integers(0, n, n) for _ in range(nboot)]
    out = {}
    for t in tags:
        d = {k: est[t][k] / est[ref][k] - 1 for k in est[t]}
        bs = {k: [] for k in d}
        for ix in draws:
            a, b = quantities(T[t], ix), quantities(T[ref], ix)
            for k in d: bs[k].append(a[k] / b[k] - 1)
        # enmienda 1 (ESPEC_U24): with < 2 seeds there is no error (every draw is the same seed); informative per-seed standard error
        per = [{k: quantities(T[t], [i])[k] / quantities(T[ref], [i])[k] - 1 for k in d} for i in range(n)]
        out[t] = {k: dict(delta=float(d[k]), sigma=(float(np.nanstd(bs[k])) if n >= 2 else None),
                          se_seeds=(float(np.std([x[k] for x in per], ddof=1) / np.sqrt(n)) if n >= 2 else None)) for k in d}
    return est, out


def main():
    runs = {v: (load(D / "runs", "RP3", v), load(D / "runs", "U1", v)) for v in VARIANTS}
    seeds = [s for s in SEEDS if all(s in rp and s in u1 for rp, u1 in runs.values())]
    if not seeds: raise SystemExit("no seed with every U24 variant in both models yet")
    T = tables(runs, seeds)
    est, delta = paired(T, "D0", [v for v in VARIANTS if v != "D0"], seeds)
    out = dict(spec="docs/preregistros/ESPEC_U24.md", window=[float(GRID[0]), float(GRID[-1])], seeds=seeds,
               provisional=len(seeds) < len(SEEDS), values=est, delta_vs_D0=delta)
    dv = [delta[v]["R_Gamma"] for v in delta]
    if len(seeds) < 2: v = f"SIN_VEREDICTO (n = {len(seeds)})"                          # enmienda 1: no sigma with one seed
    elif all(abs(x["delta"]) <= TOL and x["sigma"] <= SIG for x in dv): v = "DRENAJE_ROBUSTO"
    elif any(abs(x["delta"]) > TOL and abs(x["delta"]) > 2 * x["sigma"] for x in dv): v = "DRENAJE_SENSIBLE"
    else: v = "INDETERMINADO"
    out["verdict_U24"] = v + (" (provisional: faltan semillas)" if out["provisional"] and len(seeds) >= 2 else "")
    # enmienda 1, informative: D0 at 384^3 against U11 (576^3) evaluated on the same window [4.36, 5.21]: 1.562
    out["D0_vs_U11_same_window"] = dict(U11_576=U11_SAME_WINDOW, D0_384=est["D0"]["R_Gamma"])
    out["max_abs_delta_R_Gamma"] = float(max(abs(x["delta"]) for x in dv))
    try:
        u11 = json.loads((HERE / "results" / "u4" / "U11_RESULT.json").read_text())
        out["U11_576_R_Gamma"] = u11["R_Gamma"]["value"]                    # informative only: other N and ell_end
    except Exception: pass
    (D / "U24_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"U24 seeds {seeds}{' (PROVISIONAL)' if out['provisional'] else ''}; R_Gamma(D0) = {est['D0']['R_Gamma']:.3f}, R_xi(D0) = {est['D0']['R_xi']:.3f}")
    for t, d in delta.items():
        sg = d['R_Gamma']['sigma']
        print(f"  {t:8s} R_Gamma {est[t]['R_Gamma']:.3f}  delta {d['R_Gamma']['delta']:+.3f} +- {('%.3f' % sg) if sg is not None else 'n/a'} | "
              f"R_xi {d['R_xi']['delta']:+.3f} | eps_RP3 {d['eps_RP3']['delta']:+.3f} | eps_U1 {d['eps_U1']['delta']:+.3f}")
    print("U24:", out["verdict_U24"])
    return out


if __name__ == "__main__":
    main()
