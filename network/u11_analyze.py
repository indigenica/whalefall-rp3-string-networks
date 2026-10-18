"""U11 analysis (docs/preregistros/ESPEC_U11.md), fixed before the data: U6 drained (576^3, xi_i = 0.5, same seeds as U6).
    python -B u11_analyze.py        -> results/u4/U11_RESULT.json"""
import glob, json, math, sys
from pathlib import Path
import numpy as np
import analyze as A2
from u6_analyze import q_fit
from u10a_analyze import gamma_and_X

HERE = Path(__file__).resolve().parent
RUNS = Path(sys.argv[sys.argv.index("--runs") + 1]) if "--runs" in sys.argv else HERE / "results" / "u4" / "runs"
N = int(sys.argv[sys.argv.index("--N") + 1]) if "--N" in sys.argv else 576


def load(p): return [json.loads(Path(x).read_text()) for x in sorted(glob.glob(str(RUNS / p)))]


def corr_eps(run, grid):
    """U(1) eps_lat with the measured exchange removed: (Gamma - <X>) / norm, on grid"""
    G, X = gamma_and_X(run, grid)
    d = A2.derived(run, "rm2.5"); H = np.array([r["H"] for r in run["rows"]])
    Hg = np.interp(grid, d["ell"], H); norm = 8 * math.pi * Hg ** 3 * A2.F2 * grid * A2.KAPPA * grid
    return (G - np.where(np.isfinite(X), X, 0.0)) / norm


def boot_ratio(A, B, nboot=4000, seed=0):
    f = lambda a, b: float(np.nanmean(np.nanmean(a, 0)) / np.nanmean(np.nanmean(b, 0)))
    rng = np.random.default_rng(seed)
    bs = [f(A[rng.integers(0, len(A), len(A))], B[rng.integers(0, len(B), len(B))]) for _ in range(nboot)]
    return dict(value=f(A, B), boot_std=float(np.nanstd(bs)))


def zone(v, s):
    return "LIMPIA_ZONA_PERMITIDA" if v + s <= 1.69 else "LIMPIA_ZONA_EXCLUIDA" if v - s >= 2.09 else "LIMPIA_INTERMEDIA"


def main():
    rd, ud = load(f"RP3_N{N}_s6*_xi0.5_dr.json"), load(f"U1_N{N}_s6*_xi0.5_dr.json")
    r6, u6 = load(f"RP3_N{N}_s6*_xi0.5.json"), load(f"U1_N{N}_s6*_xi0.5.json")
    if not rd or not ud: raise SystemExit(f"not enough runs yet: RP3 {len(rd)}, U1 {len(ud)}")
    l_end = min(r["rows"][-1]["ell"] for r in rd + ud); grid = np.linspace(l_end - 1.3, l_end - 0.05, 10)
    out = dict(spec="docs/preregistros/ESPEC_U11.md", N=N, n=dict(rp3=len(rd), u1=len(ud)), window=[grid[0], grid[-1]],
               provisional=not (len(rd) >= 6 and len(ud) >= 4))
    out["R_Gamma"] = A2.window_ratio(rd, ud, grid, "eps_lat"); out["R_xi"] = A2.window_ratio(rd, ud, grid, "xi")
    A = A2.on_grid(rd, grid, "eps_lat"); B = np.array([corr_eps(r, grid) for r in ud])
    out["R_Gamma_corr"] = boot_ratio(A, B)
    out["R_Gamma_U6_same_window"] = A2.window_ratio(r6, u6, grid, "eps_lat")
    q = q_fit(rd); out["q_RP3"] = dict(q0_at_ell5=q["q0_at_ell5"], q0_std=q["q0_std"], q1_per_ell=q["q1_per_ell"], q1_std=q["q1_std"], q_end=q["q_at_ell_end"])
    pairs = {}
    for name, D, U in (("RP3", rd, r6), ("U1", ud, u6)):
        by = {r["seed"]: r for r in U}
        v = {f"s{r['seed']}": A2.window_ratio([r], [by[r["seed"]]], grid, "eps_lat")["value"] for r in D if r["seed"] in by}
        pairs[name] = dict(per_seed=v, mean=float(np.mean(list(v.values()))) if v else None, std=float(np.std(list(v.values()), ddof=1)) if len(v) > 1 else None)
    out["drained_over_undrained_by_seed"] = pairs
    le = l_end - 0.05
    mr = lambda runs: float(np.mean([np.interp(le, *(lambda d: (d["ell"], d["rhoM_over_rhoa"]))(A2.derived(r, "rm2.5"))) for r in runs])) if runs else None
    out["rhoM_over_rhoa_end"] = dict(RP3_drained=mr(rd), U1_drained=mr(ud), RP3_U6=mr(r6), U1_U6=mr(u6))
    out["f_X_U1"] = {f"s{r['seed']}": (lambda G, X: float(np.nanmean(X) / np.nanmean(G)))(*gamma_and_X(r, grid)) for r in ud}
    v, s = out["R_Gamma"]["value"], out["R_Gamma"]["boot_std"] or 0.0
    vc, sc = out["R_Gamma_corr"]["value"], out["R_Gamma_corr"]["boot_std"]
    q0, qs = out["q_RP3"]["q0_at_ell5"], out["q_RP3"]["q0_std"] or 0.0
    out["verdict_R"] = zone(v, s); out["verdict_R_corr"] = zone(vc, sc)
    out["verdict_q"] = "Q_MENOR_QUE_1" if q0 + 2 * qs < 1 else "Q_MAYOR_QUE_1" if q0 - 2 * qs > 1 else "Q_COMPATIBLE_CON_1"
    if "--runs" not in sys.argv: (HERE / "results" / "u4" / "U11_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"U11 ({'PROVISIONAL' if out['provisional'] else 'final'}; RP3 {len(rd)}, U1 {len(ud)}; window {grid[0]:.2f}-{grid[-1]:.2f})")
    print(f"  R_Gamma drained = {v:.3f} +- {s:.3f} -> {out['verdict_R']};  U(1)-X corrected = {vc:.3f} +- {sc:.3f} -> {out['verdict_R_corr']}")
    print(f"  R_xi drained = {out['R_xi']['value']:.3f};  U6 (undrained, same window) R_Gamma = {out['R_Gamma_U6_same_window']['value']:.3f}")
    print(f"  q_RP3(5) drained = {q0:.3f} +- {qs:.3f} (dq/dell = {out['q_RP3']['q1_per_ell']:.3f}) -> {out['verdict_q']}")
    print(f"  drained/undrained emission by seed: {json.dumps({k: (round(d['mean'], 3) if d['mean'] else None, d['std'] and round(d['std'], 3)) for k, d in pairs.items()})}")
    print(f"  rho_M/rho_a at window end: {json.dumps({k: (None if x is None else round(x, 3)) for k, x in out['rhoM_over_rhoa_end'].items()})}")
    print(f"  f_X (U1, drained): {json.dumps({k: round(x, 3) for k, x in out['f_X_U1'].items()})}")
    return out


if __name__ == "__main__":
    main()
