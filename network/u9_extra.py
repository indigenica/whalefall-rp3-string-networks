"""U9, amendment 1 (docs/preregistros/ESPEC_U9.md): informative analyses added after an external review.
They do not change the preregistered classification of u9_analyze.py.

  1. per-seed distribution of xi(ell_ref) and dxi/dell, dense vs sparse RP3
  2. did any dense network approach the attractor from above?
  3. energy histories outside the cores: a^4 rho_a, a^4 rho_M, a^3 rho_M
  4. upper bound on massive -> Goldstone conversion C and eps_corr = (Gamma - C)/norm (lower bound on string emission)
  5. end of the relaxation of the dense start (rho_M/rho_a within 10 % of the sparse mean)
  6. window scan of R_Gamma (raw and corrected) and q_RP3
  7. xi(ell) models: line a + b ell vs saturation xi_inf - c/ell

    python -B u9_extra.py [--runs DIR]        (default runs dir: results/u4/runs)
Output: results/u4/U9_EXTRA.json"""
import glob, json, math, sys
from pathlib import Path
import numpy as np
import analyze as A2
from u6_analyze import q_fit

HERE = Path(__file__).resolve().parent; N = 576
RUNS = Path(sys.argv[sys.argv.index("--runs") + 1]) if "--runs" in sys.argv else HERE / "results" / "u4" / "runs"
MASK = "rm2.5"; L_INI = (4.36, 4.66, 4.96, 5.16)


def load(pat): return [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / pat)))]


def tag(r): return f"{r['model']}_N{r['N']}_s{r['seed']}_xi{r['xi_init']:g}"


def dlog(y, lt):
    d = np.full_like(y, np.nan)
    for i in range(1, len(y) - 1):                           # same 3-snapshot local regression as analyze.derived
        if np.all(np.isfinite(y[i - 1:i + 2])): d[i] = np.polyfit(lt[i - 1:i + 2], y[i - 1:i + 2], 1)[0]
    return d


def derived_x(run):
    R = run["rows"]; d = A2.derived(run, MASK)
    tau = np.array([r["tau"] for r in R]); a = np.array([r["a"] for r in R]); H = np.array([r["H"] for r in R])
    rM = np.array([r[MASK]["rho_M"] for r in R]); lt = np.log(tau)
    C = np.clip(-H * dlog(a ** 3 * rM, lt) / a ** 3, 0, None)            # upper bound on massive -> Goldstone transfer
    C_lo = np.clip(-H * dlog(a ** 4 * rM, lt) / a ** 4, 0, None)         # amendment 1b: lower bound (only losses beyond relativistic redshift)
    norm = 8 * math.pi * H ** 3 * A2.F2 * d["ell"] * A2.KAPPA * d["ell"]
    Gam = d["eps_lat"] * norm
    d.update(eps_corr=(Gam - C) / norm, eps_corr_lo=(Gam - C_lo) / norm, conv_frac=np.where(Gam > 0, C / Gam, np.nan),
             conv_frac_lo=np.where(Gam > 0, C_lo / Gam, np.nan), a4rho_a=a ** 4 * np.array([r[MASK]["rho_a"] for r in R]),
             a4rho_M=a ** 4 * rM, a3rho_M=a ** 3 * rM)
    return d


def grid_vals(runs, grid, key):
    out = []
    for r in runs:
        d = derived_x(r); m = d["ok"] if key not in ("xi", "rhoM_over_rhoa") else np.isfinite(d[key])
        out.append(np.interp(grid, d["ell"][m], d[key][m], left=np.nan, right=np.nan) if m.sum() > 1 else np.full(len(grid), np.nan))
    return np.array(out)


def ratio(A_runs, B_runs, grid, key, nboot=2000, seed=0):
    if not A_runs or not B_runs: return None
    A, B = grid_vals(A_runs, grid, key), grid_vals(B_runs, grid, key)
    f = lambda a, b: float(np.nanmean(np.nanmean(a, 0)) / np.nanmean(np.nanmean(b, 0)))
    rng = np.random.default_rng(seed)
    bs = [f(A[rng.integers(0, len(A), len(A))], B[rng.integers(0, len(B), len(B))]) for _ in range(nboot)] if (len(A) > 1 or len(B) > 1) else []
    return dict(value=f(A, B), boot_std=float(np.nanstd(bs)) if bs else None, n=[len(A), len(B)])


def xi_at(R, l): return float(np.interp(l, [x["ell"] for x in R["rows"]], [x["xi"] for x in R["rows"]]))


def slope(R, l0, l1):
    a = np.array([(x["ell"], x["xi"]) for x in R["rows"] if l0 <= x["ell"] <= l1]); return float(np.polyfit(a[:, 0], a[:, 1], 1)[0])


def from_above(R, l_win):
    rows = [x for x in R["rows"] if x["ell"] > 2.5]; e = np.array([x["ell"] for x in rows]); x = np.array([x["xi"] for x in rows])
    w = e >= l_win; dx = np.diff(x[w])
    return dict(ell_min_after_collapse=float(e[np.argmin(x)]), xi_min=float(x.min()), n_decreasing_steps_in_window=int((dx < 0).sum()),
                n_steps_in_window=int(len(dx)), approached_from_above=bool((dx < 0).sum() > len(dx) / 2))


def relax_end(dense, sparse, l_end):
    if not dense or not sparse: return None
    g = np.arange(3.5, l_end - 0.05 + 1e-9, 0.05)
    rel = np.nanmean(grid_vals(dense, g, "rhoM_over_rhoa"), 0) / np.nanmean(grid_vals(sparse, g, "rhoM_over_rhoa"), 0) - 1
    ok = np.abs(rel) < 0.10
    for i in range(len(g)):
        if ok[i:].all(): return dict(ell_relax=float(g[i]) if g[i] <= l_end - 0.3 else None, rel_at_end=float(rel[-1]),
                                     track=[[float(a), float(b)] for a, b in zip(g[::4], rel[::4])])
    return dict(ell_relax=None, rel_at_end=float(rel[-1]), track=[[float(a), float(b)] for a, b in zip(g[::4], rel[::4])])


def xi_models(runs, rel_err, l0, l1):
    """weighted fits of the seed-mean xi: line a + b ell vs saturation xi_inf - c/ell (both linear in their parameters)"""
    if not runs: return None
    e0 = np.array([x["ell"] for x in runs[0]["rows"]]); g = e0[(e0 >= l0) & (e0 <= l1)]
    X = np.array([[xi_at(r, l) for l in g] for r in runs]); m = X.mean(0)
    err = X.std(0, ddof=1) / math.sqrt(len(runs)) if len(runs) >= 3 else rel_err(g) * m / math.sqrt(len(runs))
    def fit(u):
        Wm = np.diag(1 / err ** 2); A = np.vstack([np.ones_like(u), u]).T
        p = np.linalg.solve(A.T @ Wm @ A, A.T @ Wm @ m); return p, float(np.sum(((A @ p - m) / err) ** 2))
    pl, cl = fit(g); ps, cs = fit(-1 / g)
    dchi = cs - cl
    return dict(n_runs=len(runs), ell=g.tolist(), xi_mean=m.tolist(), err=err.tolist(), line=dict(a=float(pl[0]), b=float(pl[1]), chi2=cl),
                saturation=dict(xi_inf=float(ps[0]), c=float(ps[1]), chi2=cs), delta_chi2=float(dchi), dof=len(g) - 2,
                reading=("INDISTINGUIBLE" if abs(dchi) < 4 else "LINEAL_PREFERIDO" if dchi >= 4 else "SATURACION_PREFERIDA"))


def main():
    rp_s, rp_d = load(f"RP3_N{N}_s6*_xi0.5.json"), load(f"RP3_N{N}_s6*_xi2.json")
    u1_s, u1_d = load(f"U1_N{N}_s6*_xi0.5.json"), load(f"U1_N{N}_s6*_xi2.json")
    u1_768 = load("U1_N768_s6*_xi0.5.json")
    if not rp_d: raise SystemExit("no dense RP3 runs yet")
    l_end = min(r["rows"][-1]["ell"] for r in rp_s + rp_d + u1_s + u1_d); l_ref = l_end - 0.3; w = (l_ref - 0.7, l_ref)
    out = dict(spec="docs/preregistros/ESPEC_U9.md, enmienda 1", N=N, ell_end=l_end, ell_ref=l_ref, slope_window=list(w),
               runs=dict(rp3_sparse=[tag(r) for r in rp_s], rp3_dense=[tag(r) for r in rp_d], u1_sparse=[tag(r) for r in u1_s],
                         u1_dense=[tag(r) for r in u1_d], u1_768=[tag(r) for r in u1_768]))
    # 1 per-seed distribution
    xs, ss = [xi_at(r, l_ref) for r in rp_s], [slope(r, *w) for r in rp_s]
    out["per_seed"] = dict(sparse=dict(xi_ref=xs, slope=ss, xi_ref_mean=float(np.mean(xs)), xi_ref_std=float(np.std(xs, ddof=1)),
                                       slope_mean=float(np.mean(ss)), slope_std=float(np.std(ss, ddof=1))),
                           dense=[dict(run=tag(r), xi_ref=xi_at(r, l_ref), slope=slope(r, *w),
                                       slope_z=(slope(r, *w) - np.mean(ss)) / np.std(ss, ddof=1)) for r in rp_d])
    # 2 from above
    out["from_above"] = {tag(r): from_above(r, 4.36) for r in rp_d + u1_d}
    # 3 energy histories
    hist = {}
    for r in rp_s + rp_d + u1_s + u1_d:
        d = derived_x(r); m = d["ell"] > 2.5
        hist[tag(r)] = dict(ell=d["ell"][m].tolist(), a4rho_a=d["a4rho_a"][m].tolist(), a4rho_M=d["a4rho_M"][m].tolist(),
                            a3rho_M=d["a3rho_M"][m].tolist(), rhoM_over_rhoa=d["rhoM_over_rhoa"][m].tolist(),
                            conv_frac=[None if not np.isfinite(v) else float(v) for v in d["conv_frac"][m]])
    out["energy_histories"] = hist
    # 4 conversion bound (U4b window)
    gw = np.linspace(l_ref - 1.0, l_end - 0.05, 10)
    out["conversion"] = dict(grid=gw.tolist(), **{k: dict(conv_frac_mean=np.nanmean(grid_vals(v, gw, "conv_frac"), 0).tolist(),
                                                           conv_frac_lo_mean=np.nanmean(grid_vals(v, gw, "conv_frac_lo"), 0).tolist()) if v else None
                                                 for k, v in (("rp3_sparse", rp_s), ("rp3_dense", rp_d), ("u1_sparse", u1_s), ("u1_dense", u1_d))},
                             emission_dense_over_sparse=dict(RP3_raw=ratio(rp_d, rp_s, gw, "eps_lat"), RP3_corr=ratio(rp_d, rp_s, gw, "eps_corr"),
                                                             RP3_corr_lo=ratio(rp_d, rp_s, gw, "eps_corr_lo"),
                                                             U1_raw=ratio(u1_d, u1_s, gw, "eps_lat"), U1_corr=ratio(u1_d, u1_s, gw, "eps_corr"),
                                                             U1_corr_lo=ratio(u1_d, u1_s, gw, "eps_corr_lo")))
    # 5 relaxation end
    out["relaxation"] = dict(RP3=relax_end(rp_d, rp_s, l_end), U1=relax_end(u1_d, u1_s, l_end))
    # 6 window scan
    starts = list(L_INI); lr = (out["relaxation"]["RP3"] or {}).get("ell_relax")
    if lr is not None and all(abs(lr - s) > 0.05 for s in starts): starts.append(round(lr, 2))
    scan = []
    for l0 in sorted(starts):
        g = np.linspace(l0, l_end - 0.05, 10)
        row = dict(ell_ini=l0, R_Gamma_dense=ratio(rp_d, u1_d, g, "eps_lat"), R_Gamma_dense_corr=ratio(rp_d, u1_d, g, "eps_corr"),
                   R_Gamma_dense_corr_lo=ratio(rp_d, u1_d, g, "eps_corr_lo"),
                   R_Gamma_sparse=ratio(rp_s, u1_s, g, "eps_lat"), R_Gamma_sparse_corr=ratio(rp_s, u1_s, g, "eps_corr"),
                   R_Gamma_sparse_corr_lo=ratio(rp_s, u1_s, g, "eps_corr_lo"))
        for k, runs in (("q_RP3_dense", rp_d), ("q_RP3_sparse", rp_s)):
            try:
                q = q_fit(runs, l0=l0); npts = sum(len(s) for s in q["per_seed"])
                row[k] = dict(q0_at_ell5=q["q0_at_ell5"], q0_std=q["q0_std"], q_at_ell_end=q["q_at_ell_end"], n_points=npts,
                              reliable=npts >= 3 * len(runs))
            except Exception as ex:                          # too few spectral points for a late start
                row[k] = dict(error=str(ex))
        scan.append(row)
    out["window_scan"] = scan
    # 7 xi models
    e0 = np.array([x["ell"] for x in rp_s[0]["rows"]]); X = np.array([[xi_at(r, l) for l in e0] for r in rp_s])
    relsd = X.std(0, ddof=1) / X.mean(0)
    rel_err = lambda g: np.interp(g, e0, relsd)
    l768 = min(r["rows"][-1]["ell"] for r in u1_768) if u1_768 else None
    out["xi_models"] = dict(RP3_sparse=xi_models(rp_s, rel_err, 4.36, l_end), RP3_dense=xi_models(rp_d, rel_err, 4.36, l_end),
                            U1_sparse_576=xi_models(u1_s, rel_err, 4.36, l_end), U1_dense_576=xi_models(u1_d, rel_err, 4.36, l_end),
                            U1_sparse_768=xi_models(u1_768, rel_err, 4.36, l768) if u1_768 else None)
    (HERE / "results" / "u4" / "U9_EXTRA.json").write_text(json.dumps(out, indent=1, default=float))

    ps = out["per_seed"]
    print(f"[1] sparse RP3: xi(ref) {ps['sparse']['xi_ref_mean']:.3f} +- {ps['sparse']['xi_ref_std']:.3f}, slope {ps['sparse']['slope_mean']:.3f} +- {ps['sparse']['slope_std']:.3f}")
    for d in ps["dense"]: print(f"    dense {d['run']}: xi(ref) {d['xi_ref']:.3f}, slope {d['slope']:.3f} (z = {d['slope_z']:+.2f})")
    for k, v in out["from_above"].items():
        print(f"[2] {k}: min xi {v['xi_min']:.3f} at ell {v['ell_min_after_collapse']:.2f}; decreasing steps in window {v['n_decreasing_steps_in_window']}/{v['n_steps_in_window']}")
    c = out["conversion"]
    for k in ("rp3_sparse", "rp3_dense", "u1_sparse", "u1_dense"):
        if c[k]: print(f"[4] conversion / Gamma, {k}: upper bound {np.nanmean(c[k]['conv_frac_mean']):.2f}, lower bound {np.nanmean(c[k]['conv_frac_lo_mean']):.2f}")
    for k, v in c["emission_dense_over_sparse"].items():
        if v: print(f"    emission dense/sparse {k}: {v['value']:.2f}")
    print(f"[5] relaxation: " + "; ".join(f"{k}: ell_relax {(v or {}).get('ell_relax')}, rhoM/rhoa dense/sparse - 1 at end {(v or {}).get('rel_at_end')}" for k, v in out['relaxation'].items()))
    for r in scan:
        f = lambda x: "-" if not x else f"{x['value']:.3f}"
        qf = lambda x: "n/a" if "q0_at_ell5" not in x else f"{x['q0_at_ell5']:.2f}{'' if x['reliable'] else ' (few pts)'}"
        print(f"[6] ell_ini {r['ell_ini']:.2f}: R_G dense {f(r['R_Gamma_dense'])} (corr {f(r['R_Gamma_dense_corr_lo'])}..{f(r['R_Gamma_dense_corr'])}), "
              f"sparse {f(r['R_Gamma_sparse'])} (corr {f(r['R_Gamma_sparse_corr_lo'])}..{f(r['R_Gamma_sparse_corr'])}); "
              f"q dense {qf(r['q_RP3_dense'])}, q sparse {qf(r['q_RP3_sparse'])}")
    for k, v in out["xi_models"].items():
        if v: print(f"[7] {k}: line b = {v['line']['b']:.3f}; saturation xi_inf = {v['saturation']['xi_inf']:.2f}; dchi2 = {v['delta_chi2']:+.1f} -> {v['reading']}")
    return out


if __name__ == "__main__":
    main()
