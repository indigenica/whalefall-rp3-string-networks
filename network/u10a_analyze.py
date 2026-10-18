"""U10a analysis (docs/preregistros/ESPEC_U10a.md + amendment 1), fixed before the production data.
A (verdict): drained networks (massive velocity drained until ell = 3.8).
    R_dd = R_Gamma(RP3 drained dense / U1 drained dense), R_ds = same with the sparse starts, window [ell_end - 1.3, ell_end - 0.05].
    G = |R_dd - R_ds| / mean;  dq = |q_RP3(dd) - q_RP3(ds)|;  massive excess E_M = max over the two models of
    |(rho_M/rho_a)_dense / (rho_M/rho_a)_sparse - 1| at ell_end - 0.05 (drained runs, mask rm2.5).
    RELAJACION_EXPLICA if G < 0.08 and dq < 0.15 and E_M < 0.10;  RELAJACION_NO_EXPLICA if E_M < 0.10 and G >= 0.15;  else INDETERMINADO.
B (informative): exchange X (U(1) quantitative if V5 passed; RP3 qualitative after V5 failed): f_X = <X>/Gamma in the window.
Control: the undrained RP3 run with X must reproduce U9 seed 61 in xi(ell) to 1e-5 relative.
Partial reading allowed with seed 61 only ("provisional" in the output); the verdict needs both drained seeds of each kind.
    python -B u10a_analyze.py        -> results/u4/U10A_RESULT.json"""
import glob, json, math, sys
from pathlib import Path
import numpy as np
import analyze as A2
from u6_analyze import q_fit

HERE = Path(__file__).resolve().parent
RUNS = Path(sys.argv[sys.argv.index("--runs") + 1]) if "--runs" in sys.argv else HERE / "results" / "u4" / "runs"
N = int(sys.argv[sys.argv.index("--N") + 1]) if "--N" in sys.argv else 576


def load(pat): return [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / pat)))]


def ratio_at(runs, l, key):
    v = []
    for r in runs:
        d = A2.derived(r, "rm2.5"); m = np.isfinite(d[key])
        v.append(float(np.interp(l, d["ell"][m], d[key][m])))
    return float(np.mean(v)) if v else float("nan")


def gamma_and_X(run, grid):
    """Gamma (physical, rm2.5 estimator) and the mean exchange X (samples with frac >= 0.5), both on grid"""
    d = A2.derived(run, "rm2.5"); H = np.array([r["H"] for r in run["rows"]])
    Gam = d["eps_lat"] * 8 * math.pi * H ** 3 * A2.F2 * d["ell"] * A2.KAPPA * d["ell"]
    ok = d["ok"]; G = np.interp(grid, d["ell"][ok], Gam[ok], left=np.nan, right=np.nan)
    S = [(s["ell"], s["X"]) for r in run["rows"] for s in r.get("exch", []) if s["frac"] >= 0.5]
    if len(S) < 2: return G, np.full(len(grid), np.nan)
    S = np.array(sorted(S)); lb = np.round(S[:, 0], 1)                     # average the samples of each interval
    ub = np.unique(lb); Xm = np.array([S[lb == u, 1].mean() for u in ub]); lm = np.array([S[lb == u, 0].mean() for u in ub])
    return G, np.interp(grid, lm, Xm, left=np.nan, right=np.nan)


def main():
    S = {k: load(p) for k, p in dict(RP3_dd=f"RP3_N{N}_s6*_xi2_dr.json", RP3_ds=f"RP3_N{N}_s6*_xi0.5_dr.json", U1_dd=f"U1_N{N}_s6*_xi2_dr.json",
                                     U1_ds=f"U1_N{N}_s6*_xi0.5_dr.json", RP3_ux=f"RP3_N{N}_s61_xi2_x.json", U1_ux=f"U1_N{N}_s61_xi2_x.json",
                                     RP3_u6=f"RP3_N{N}_s6*_xi0.5.json", U1_u6=f"U1_N{N}_s6*_xi0.5.json", RP3_u9=f"RP3_N{N}_s6*_xi2.json",
                                     U1_u9=f"U1_N{N}_s6*_xi2.json").items()}
    have = {k: len(v) for k, v in S.items()}
    if not (S["RP3_dd"] and S["U1_dd"]): raise SystemExit(f"not enough runs yet: {have}")
    allr = [r for k in ("RP3_dd", "RP3_ds", "U1_dd", "U1_ds") for r in S[k]]
    l_end = min(r["rows"][-1]["ell"] for r in allr); grid = np.linspace(l_end - 1.3, l_end - 0.05, 10)
    out = dict(spec="docs/preregistros/ESPEC_U10a.md", N=N, runs=have, window=[grid[0], grid[-1]],
               provisional=not (have["RP3_dd"] >= 2 and have["RP3_ds"] >= 2 and have["U1_dd"] >= 2 and have["U1_ds"] >= 2))
    R = lambda a, b, key="eps_lat": A2.window_ratio(S[a], S[b], grid, key) if S[a] and S[b] else None
    out["R_Gamma"] = dict(drained_dense=R("RP3_dd", "U1_dd"), drained_sparse=R("RP3_ds", "U1_ds"),
                          undrained_sparse_U6=R("RP3_u6", "U1_u6"), undrained_dense_U9=R("RP3_u9", "U1_u9"))
    out["R_xi"] = dict(drained_dense=R("RP3_dd", "U1_dd", "xi"), drained_sparse=R("RP3_ds", "U1_ds", "xi"))
    out["emission_drained_over_undrained"] = {f"{m}_{k}": (R(f"{m}_d{k[0]}", f"{m}_u{'9' if k == 'dense' else '6'}") or {}).get("value")
                                             for m in ("RP3", "U1") for k in ("dense", "sparse")}
    def qf(runs):
        try:
            q = q_fit(runs); return dict(q0_at_ell5=q["q0_at_ell5"], q0_std=q["q0_std"], q_end=q["q_at_ell_end"])
        except ValueError:                                   # no spectral points at ell >= 4.6 (small test boxes)
            return dict(q0_at_ell5=float("nan"), q0_std=None, q_end=float("nan"))
    out["q_RP3"] = {k: qf(S[k]) for k in ("RP3_dd", "RP3_ds") if S[k]}
    le = l_end - 0.05
    mx = {m: dict(dense=ratio_at(S[f"{m}_dd"], le, "rhoM_over_rhoa"), sparse=ratio_at(S[f"{m}_ds"], le, "rhoM_over_rhoa"),
                  undrained_dense=ratio_at(S[f"{m}_u9"], le, "rhoM_over_rhoa"), undrained_sparse=ratio_at(S[f"{m}_u6"], le, "rhoM_over_rhoa"))
          for m in ("RP3", "U1")}
    for m in mx: mx[m]["excess"] = abs(mx[m]["dense"] / mx[m]["sparse"] - 1) if mx[m]["sparse"] == mx[m]["sparse"] else float("nan")
    out["massive"] = mx
    a, b = out["R_Gamma"]["drained_dense"], out["R_Gamma"]["drained_sparse"]
    if a and b:
        G = abs(a["value"] - b["value"]) / (0.5 * (a["value"] + b["value"]))
        dq = abs(out["q_RP3"]["RP3_dd"]["q0_at_ell5"] - out["q_RP3"]["RP3_ds"]["q0_at_ell5"]) if "RP3_ds" in out["q_RP3"] else float("nan")
        EM = max(mx["RP3"]["excess"], mx["U1"]["excess"])
        verdict = ("RELAJACION_EXPLICA" if G < 0.08 and dq < 0.15 and EM < 0.10 else
                   "RELAJACION_NO_EXPLICA" if EM < 0.10 and G >= 0.15 else "INDETERMINADO")
        out.update(gap_G=G, gap_U9=0.195, dq=dq, massive_excess=EM, verdict_A=verdict,
                   R_Gamma_clean=(0.5 * (a["value"] + b["value"]) if verdict == "RELAJACION_EXPLICA" else None))
    # B: exchange
    v5p = HERE / "results" / "u4" / "U10A_V5.json"
    v5 = {r["model"]: r["passed"] for r in json.loads(v5p.read_text())["V5"]} if v5p.exists() else {}
    fx = {}
    for k in ("U1_dd", "U1_ds", "U1_ux", "RP3_ux"):
        for r in S[k]:
            G_, X_ = gamma_and_X(r, grid)
            fx[f"{k}_s{r['seed']}"] = dict(f_X=float(np.nanmean(X_) / np.nanmean(G_)) if np.isfinite(X_).any() else None,
                                           quantitative=bool(v5.get(r["model"], False)))
    out["exchange"] = dict(V5=v5, f_X=fx)
    # control: undrained RP3 with X vs U9 seed 61
    u9 = [r for r in S["RP3_u9"] if r["seed"] == 61]
    if S["RP3_ux"] and u9:
        x1 = np.array([r["xi"] for r in S["RP3_ux"][0]["rows"]]); x0 = np.array([r["xi"] for r in u9[0]["rows"]])
        out["control_reproduces_U9_s61"] = dict(max_rel=float(np.max(np.abs(x1 - x0) / np.abs(x0))), passed=bool(np.max(np.abs(x1 - x0) / np.abs(x0)) <= 1e-5))
    if "--runs" not in sys.argv: (HERE / "results" / "u4" / "U10A_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    f = lambda d: "-" if not d else f"{d['value']:.3f}" + (f" +- {d['boot_std']:.3f}" if d.get("boot_std") else "")
    print(f"U10a ({'PROVISIONAL' if out['provisional'] else 'final'}; runs {have})")
    print(f"  R_Gamma drained dense {f(a)}, drained sparse {f(b)}; undrained sparse (U6) {f(out['R_Gamma']['undrained_sparse_U6'])}, dense (U9) {f(out['R_Gamma']['undrained_dense_U9'])}")
    print(f"  q_RP3: {json.dumps({k: round(v['q0_at_ell5'], 3) for k, v in out['q_RP3'].items()})}")
    print(f"  rho_M/rho_a at window end: {json.dumps({m: {k: round(v, 3) for k, v in d.items()} for m, d in mx.items()})}")
    print(f"  emission drained/undrained: {json.dumps({k: (None if v is None else round(v, 3)) for k, v in out['emission_drained_over_undrained'].items()})}")
    if "verdict_A" in out: print(f"  G = {out['gap_G']:.3f} (U9: 0.195), dq = {out['dq']:.3f}, massive excess = {out['massive_excess']:.3f} -> {out['verdict_A']}")
    print(f"  exchange f_X: {json.dumps(fx)}")
    if "control_reproduces_U9_s61" in out: print(f"  control vs U9 s61: {out['control_reproduces_U9_s61']}")
    return out


if __name__ == "__main__":
    main()
