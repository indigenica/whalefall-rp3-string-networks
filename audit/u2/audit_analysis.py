"""Auditoria independiente U2 - recomputo del analisis (check 4) desde results/runs/*.json (solo lectura).
Reimplementado sin importar analyze.py.  Salida: audit_u2/audit_analysis.json"""
import glob, json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; RUNS = HERE.parent / "results" / "runs"
F2, KAP = 0.5, 0.24


def load(pat): return [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / pat)))]


def series(run, mask="rm2.5"):
    R = run["rows"]; g = lambda k: np.array([r[k] for r in R], float)
    tau, a, H, ell, xi = g("tau"), g("a"), g("H"), g("ell"), g("xi")
    rho = np.array([r[mask]["rho_a"] for r in R], float); fr = np.array([r[mask]["frac"] for r in R], float)
    rM = np.array([r[mask]["rho_M"] for r in R], float)
    y = a ** 4 * rho; lt = np.log(tau)
    # derivative: 3-point least squares slope in ln tau (independent formula: sum (x-xb)(y-yb)/sum (x-xb)^2)
    dy = np.full(len(y), np.nan)
    for i in range(1, len(y) - 1):
        xs, ys = lt[i - 1:i + 2], y[i - 1:i + 2]
        if np.all(np.isfinite(ys)):
            dy[i] = np.sum((xs - xs.mean()) * (ys - ys.mean())) / np.sum((xs - xs.mean()) ** 2)
    Gam = H * dy / a ** 4                      # Gamma = a^-4 d(a^4 rho)/dt, d/dt = H d/dln tau in RD
    base = 8 * math.pi * H ** 3 * F2 * ell
    ok = (fr >= 0.5) & np.isfinite(Gam)
    return dict(tau=tau, a=a, H=H, ell=ell, xi=xi, y=y, yM=a ** 4 * rM, fr=fr, eps=Gam / (base * KAP * ell), eta=Gam / (base * xi), ok=ok,
                s=run["s"], m=math.sqrt(run["m2"]))


def grid_vals(runs, grid, key, mask="rm2.5"):
    out = []
    for r in runs:
        d = series(r, mask); m = d["ok"] if key not in ("xi",) else np.ones(len(d["ell"]), bool)
        out.append(np.interp(grid, d["ell"][m], d[key][m], left=np.nan, right=np.nan))
    return np.array(out)


def ratio_window(rp, u1, grid, key, mask="rm2.5"):
    A, B = grid_vals(rp, grid, key, mask), grid_vals(u1, grid, key, mask)
    return float(np.nanmean(np.nanmean(A, 0)) / np.nanmean(np.nanmean(B, 0))), A, B


def eps_int(run, l0, l1, mask="rm2.5"):
    d = series(run, mask)
    if np.interp(l0, d["ell"], d["fr"]) < 0.5: return float("nan")
    s, m = d["s"], d["m"]
    # analytic: int 8 pi H^3 f2 kap ell^2 a^4 dt, a = s tau, H = 1/(s tau^2), dt = a dtau ; use ell as variable: ell = ln(m s tau^2)
    # => tau = sqrt(e^ell/(m s)), dtau = tau/2 dell ; integrand in ell: 8 pi f2 kap ell^2 * s^2 tau^2 /(s^3 tau^6) * tau/2 * ... do numerically
    ls = np.linspace(l0, l1, 20001); tau = np.sqrt(np.exp(ls) / (m * s)); a = s * tau; H = 1 / (s * tau ** 2)
    integ = 8 * math.pi * H ** 3 * F2 * KAP * ls ** 2 * a ** 4 * a * (tau / 2)
    I = np.trapezoid(integ, ls)
    ok = np.isfinite(d["y"])
    return float((np.interp(l1, d["ell"][ok], d["y"][ok]) - np.interp(l0, d["ell"][ok], d["y"][ok])) / I)


def boot(fun, A, B, n=4000, seed=0):
    rng = np.random.default_rng(seed)
    return float(np.nanstd([fun(A[rng.integers(0, len(A), len(A))], B[rng.integers(0, len(B), len(B))]) for _ in range(n)]))


def main():
    S = {"A": ("RP3_N256_s[0-9]_f1_dt0.25_xi0.5.json", "U1_N256_s[0-9]_f1_dt0.25_xi0.5.json"),
         "B": ("RP3_N192_*xi0.5.json", "U1_N192_*xi0.5.json"), "C": ("RP3_N256_s21_f0.6*.json", "U1_N256_s21_f0.6*.json"),
         "D25": ("RP3_N128_*dt0.25_*.json", "U1_N128_*dt0.25_*.json"), "D125": ("RP3_N128_*dt0.125*.json", "U1_N128_*dt0.125*.json"),
         "E25": ("RP3_N256_s31_*xi0.25.json", "U1_N256_s31_*xi0.25.json"), "E100": ("RP3_N256_s31_*xi1.json", "U1_N256_s31_*xi1.json")}
    sets = {k: (load(a), load(b)) for k, (a, b) in S.items()}
    out = dict(available={k: [len(v[0]), len(v[1])] for k, v in sets.items()})
    rp, u1 = sets["A"]
    lfin = min(r["rows"][-1]["ell"] for r in rp + u1); grid = np.linspace(3.0, lfin - 0.05, 12)
    Rxi, Ax, Bx = ratio_window(rp, u1, grid, "xi"); RG, AG, BG = ratio_window(rp, u1, grid, "eps")
    f = lambda a, b: np.nanmean(np.nanmean(a, 0)) / np.nanmean(np.nanmean(b, 0))
    out["A"] = dict(R_xi=Rxi, R_xi_boot=boot(f, Ax, Bx), R_Gamma_local=RG, R_Gamma_local_boot=boot(f, AG, BG),
                    per_ell_R_Gamma=(np.nanmean(AG, 0) / np.nanmean(BG, 0)).tolist(), grid=grid.tolist(),
                    eps_lat_rp3=float(np.nanmean(AG)), eps_lat_u1=float(np.nanmean(BG)),
                    eta_u1=float(np.nanmean(grid_vals(u1, grid, "eta"))), eta_rp3=float(np.nanmean(grid_vals(rp, grid, "eta"))),
                    last_grid_point_valid=bool(np.isfinite(np.nanmean(AG, 0)[-1])))
    # seed-level classical error (alternative to the bootstrap with only 2 RP3 seeds)
    ma, mb = np.nanmean(AG, 1), np.nanmean(BG, 1)
    rel = math.sqrt(np.var(ma, ddof=1) / len(ma) / ma.mean() ** 2 + np.var(mb, ddof=1) / len(mb) / mb.mean() ** 2)
    out["A"]["per_seed_eps_rp3"] = ma.tolist(); out["A"]["per_seed_eps_u1"] = mb.tolist()
    out["A"]["R_Gamma_local_classical_se"] = float(RG * rel)
    # integral estimator
    l1 = lfin - 0.05
    ei_rp = np.array([eps_int(r, 3.0, l1) for r in rp]); ei_u1 = np.array([eps_int(r, 3.0, l1) for r in u1])
    fi = lambda a, b: np.nanmean(a) / np.nanmean(b)
    out["A"]["eps_int_rp3"] = ei_rp.tolist(); out["A"]["eps_int_u1"] = ei_u1.tolist()
    out["A"]["R_Gamma_int"] = float(fi(ei_rp, ei_u1)); out["A"]["R_Gamma_int_boot"] = boot(fi, ei_rp, ei_u1, seed=2)
    # sub-windows of the integral estimator (trend check robust to local noise)
    sub = {}
    for lo, hi in ((3.0, 3.6), (3.6, 4.2), (4.2, l1), (3.0, 3.9), (3.9, l1)):
        a_ = np.array([eps_int(r, lo, hi) for r in rp]); b_ = np.array([eps_int(r, lo, hi) for r in u1])
        sub[f"{lo:.2f}-{hi:.2f}"] = dict(eps_rp3=float(a_.mean()), eps_u1=float(b_.mean()), R=float(a_.mean() / b_.mean()), R_boot=boot(fi, a_, b_, n=2000))
    out["A"]["R_Gamma_int_subwindows"] = sub
    # trend: slope of per-ell ratio, finite points only
    pe = np.array(out["A"]["per_ell_R_Gamma"]); mk = np.isfinite(pe)
    sl, ic = np.polyfit(grid[mk], pe[mk], 1)
    res = pe[mk] - (sl * grid[mk] + ic); se_ols = math.sqrt(np.sum(res ** 2) / (mk.sum() - 2) / np.sum((grid[mk] - grid[mk].mean()) ** 2))
    rng = np.random.default_rng(1); sls = []
    for _ in range(2000):
        a_ = np.nanmean(AG[rng.integers(0, len(AG), len(AG))], 0); b_ = np.nanmean(BG[rng.integers(0, len(BG), len(BG))], 0)
        rr = a_ / b_; m = np.isfinite(rr)
        if m.sum() > 3: sls.append(np.polyfit(grid[m], rr[m], 1)[0])
    c1, c0 = np.polyfit(1 / grid[mk], pe[mk], 1)
    out["trend"] = dict(slope=float(sl), boot_std_seeds=float(np.std(sls)), z_boot=float(abs(sl) / np.std(sls)),
                        ols_se_from_scatter=se_ols, z_ols=float(abs(sl) / se_ols),
                        significant_by_code_bootstrap=bool(abs(sl) > 2 * np.std(sls)), significant_by_ols=bool(abs(sl) > 2 * se_ols),
                        B3_fit_c0_c1=[float(c0), float(c1)], B3_R_at={str(L): float(c0 + c1 / L) for L in (5, 25, 50, 70, 87, 100)},
                        note_code_bug="analyze.py polyfit sobre per_ell con NaN en el ultimo punto -> slope NaN -> significant=False siempre")
    # systematics as preregistered (+ amendments)
    def wr(rp_, u1_, g): return ratio_window(rp_, u1_, g, "eps")[0]
    syst = {}
    for name, key in (("box_B", "B"), ("core_C", "C"), ("IC_E25", "E25"), ("IC_E100", "E100")):
        a_, b_ = sets[key]
        if a_ and b_:
            l2 = min(r["rows"][-1]["ell"] for r in a_ + b_); g = np.linspace(3.0, l2 - 0.05, 10)
            v, vA = wr(a_, b_, g), wr(rp, u1, g[g <= lfin]); syst[name] = dict(value=v, A_same=vA, delta=abs(v - vA))
            ia = np.array([eps_int(r, 3.0, l2 - 0.05) for r in a_]); ib = np.array([eps_int(r, 3.0, l2 - 0.05) for r in b_])
            iA = np.array([eps_int(r, 3.0, min(l2, lfin) - 0.05) for r in rp]); iB = np.array([eps_int(r, 3.0, min(l2, lfin) - 0.05) for r in u1])
            syst[name + "_int"] = dict(value=float(ia.mean() / ib.mean()), A_same=float(iA.mean() / iB.mean()), delta=abs(ia.mean() / ib.mean() - iA.mean() / iB.mean()))
        elif b_ or a_:
            syst[name] = dict(note=f"incompleto: RP3={len(a_)} U1={len(b_)}")
    d25, d125 = sets["D25"], sets["D125"]
    l2 = min(r["rows"][-1]["ell"] for r in d25[0] + d25[1] + d125[0] + d125[1]); g = np.linspace(3.0, l2 - 0.05, 8)
    v1, v2 = wr(*d25, g), wr(*d125, g); syst["timestep_D"] = dict(dt025=v1, dt0125=v2, delta=abs(v1 - v2))
    for mk_ in ("rm1.5", "rm4"):
        v = ratio_window(rp, u1, grid, "eps", mask=mk_)[0]; syst["mask_" + mk_] = dict(value=v, delta=abs(v - RG))
    out["systematics"] = syst
    prim = [k for k in syst if not k.endswith("_int") and "delta" in syst[k]]
    stat = out["A"]["R_Gamma_local_boot"]
    tot = math.sqrt(stat ** 2 + sum(syst[k]["delta"] ** 2 for k in prim))
    dest = abs(out["A"]["R_Gamma_int"] - RG)
    tot_final = math.sqrt(tot ** 2 + dest ** 2) if dest > tot else tot
    # alternative: classical seed error instead of 2-seed bootstrap
    tot_cl = math.sqrt(out["A"]["R_Gamma_local_classical_se"] ** 2 + sum(syst[k]["delta"] ** 2 for k in prim))
    out["R_Gamma_combined"] = dict(value=RG, stat_boot=stat, total=tot_final, included=prim, estimator_diff=dest,
                                   total_with_classical_seed_se=tot_cl)
    # U(1) gate
    xg = grid_vals(u1, np.array([4.5]), "xi")[:, 0]; xr = grid_vals(rp, np.array([4.5]), "xi")[:, 0]
    out["U1_gate"] = dict(xi_U1_4p5_per_seed=xg.tolist(), xi_U1_4p5_mean=float(xg.mean()), xi_RP3_4p5_mean=float(xr.mean()),
                          eta_U1=out["A"]["eta_u1"], pass_xi=bool(0.5 <= xg.mean() <= 2.5), pass_eta=bool(0.3 <= out["A"]["eta_u1"] <= 3))
    for key in ("B", "E25", "E100", "C"):
        a_, b_ = sets[key]
        if b_: out["U1_gate"][f"xi_U1_4p5_{key}"] = float(np.nanmean(grid_vals(b_, np.array([4.5]), "xi")))
        if b_ and key == "C": out["U1_gate"]["xi_U1_4p2_C"] = float(np.nanmean(grid_vals(b_, np.array([4.2]), "xi")))
    # attractor diagnostics: xi(ell) per set and slope dxi/dell in window
    att = {}
    pts = np.array([3.0, 3.5, 4.0, 4.5, 4.8])
    for key in ("A", "B", "C", "E25", "E100"):
        for mi, lab in ((0, "RP3"), (1, "U1")):
            runs = sets[key][mi]
            if not runs: continue
            v = np.nanmean(grid_vals(runs, pts, "xi"), 0)
            gg = np.linspace(3.2, min(r["rows"][-1]["ell"] for r in runs), 10); vv = np.nanmean(grid_vals(runs, gg, "xi"), 0)
            att[f"{lab}_{key}"] = dict(xi_at=dict(zip([str(p) for p in pts], [None if not np.isfinite(x) else round(float(x), 4) for x in v])),
                                       dxi_dell_window=float(np.polyfit(gg, vv, 1)[0]))
    # early peak and minimum (A)
    for mi, lab in ((0, "RP3"), (1, "U1")):
        for r in sets["A"][mi][:1]:
            d = series(r); m = d["ell"] < 3.5
            ip = int(np.argmax(np.where(d["ell"] < 2.6, d["xi"], -1))); im = int(np.argmin(np.where(m & (d["ell"] > 2.5), d["xi"], 9)))
            att[f"{lab}_A_transient"] = dict(peak_xi=float(d["xi"][ip]), peak_ell=float(d["ell"][ip]), min_xi=float(d["xi"][im]), min_ell=float(d["ell"][im]),
                                             yM_over_y_at_3=float(np.interp(3.0, d["ell"], d["yM"]) / np.interp(3.0, d["ell"], d["y"])))
    out["attractor"] = att
    # massive bath vs Goldstone increment over the window (A means)
    mb = {}
    for mi, lab in ((0, "RP3"), (1, "U1")):
        y0 = np.mean([np.interp(3.0, series(r)["ell"], series(r)["yM"]) for r in sets["A"][mi]])
        y1 = np.mean([np.interp(l1, series(r)["ell"], series(r)["yM"]) for r in sets["A"][mi]])
        g0 = np.mean([np.interp(3.0, series(r)["ell"], series(r)["y"]) for r in sets["A"][mi]])
        g1 = np.mean([np.interp(l1, series(r)["ell"], series(r)["y"]) for r in sets["A"][mi]])
        mb[lab] = dict(a4rhoM_3=float(y0), a4rhoM_end=float(y1), a4rhoa_3=float(g0), a4rhoa_end=float(g1), dGold=float(g1 - g0),
                       massive_over_dGold=float(y0 / (g1 - g0)))
    out["massive_bath"] = mb
    (HERE / "audit_analysis.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
