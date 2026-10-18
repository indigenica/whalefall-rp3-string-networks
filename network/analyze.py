"""U2 analysis (ESPEC_U2.md §3-§6 + amendments): xi, Goldstone emission, R_xi, R_Gamma with seeds bootstrap and systematics.
Output: results/u2/U2_RESULT.json, figures/u2_red.png"""
import glob, json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; RUNS = HERE / "results" / "u2" / "runs"
F2 = 0.5; KAPPA = 0.24
MASKS = ("rm1.5", "rm2.5", "rm4")


def load(pattern):
    return [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / pattern)))]


def derived(run, mask="rm2.5"):
    R = run["rows"]
    tau = np.array([r["tau"] for r in R]); a = np.array([r["a"] for r in R]); H = np.array([r["H"] for r in R])
    ell = np.array([r["ell"] for r in R]); xi = np.array([r["xi"] for r in R])
    rho = np.array([r[mask]["rho_a"] for r in R]); fr = np.array([r[mask]["frac"] for r in R]); rM = np.array([r[mask]["rho_M"] for r in R])
    y = a ** 4 * rho; lt = np.log(tau)
    dy = np.full_like(y, np.nan)
    for i in range(1, len(y) - 1):                         # local linear regression over 3 snapshots (centred)
        if np.all(np.isfinite(y[i - 1:i + 2])):
            dy[i] = np.polyfit(lt[i - 1:i + 2], y[i - 1:i + 2], 1)[0]
    Gam = H * dy / a ** 4
    norm = 8 * math.pi * H ** 3 * F2 * ell
    ok = (fr >= 0.5) & np.isfinite(Gam)
    return dict(ell=ell, xi=xi, eps_lat=Gam / (norm * KAPPA * ell), eta=Gam / (norm * xi), Gam_over_H3=Gam / H ** 3,
                rhoM_over_rhoa=rM / rho, frac=fr, ok=ok)


def on_grid(runs, grid, key, mask="rm2.5"):
    out = []
    for r in runs:
        d = derived(r, mask); m = d["ok"] if key != "xi" else (d["ell"] > 0)
        x, v = d["ell"][m], d[key][m]
        out.append(np.interp(grid, x, v, left=np.nan, right=np.nan) if len(x) > 1 else np.full(len(grid), np.nan))
    return np.array(out)


def window_ratio(rp, u1, grid, key, mask="rm2.5", nboot=4000, seed=0):
    A = on_grid(rp, grid, key, mask); B = on_grid(u1, grid, key, mask)
    ratio = lambda a, b: np.nanmean(np.nanmean(a, 0)) / np.nanmean(np.nanmean(b, 0))
    est = ratio(A, B)
    rng = np.random.default_rng(seed); bs = []
    if len(A) > 1 or len(B) > 1:
        for _ in range(nboot):
            bs.append(ratio(A[rng.integers(0, len(A), len(A))], B[rng.integers(0, len(B), len(B))]))
    per_ell = np.nanmean(A, 0) / np.nanmean(B, 0)
    return dict(value=float(est), boot_std=float(np.nanstd(bs)) if bs else None, per_ell=per_ell.tolist(), grid=grid.tolist(),
                rp3_mean=np.nanmean(A, 0).tolist(), u1_mean=np.nanmean(B, 0).tolist(), n_rp3=len(A), n_u1=len(B))


def eps_integral(run, l0, l1, mask="rm2.5"):
    """amendment 5: eps_int = Delta(a^4 rho_a) / int 8 pi H^3 f^2 ell 0.24 ell a^4 dt over [l0, l1] (dt = a dtau)"""
    R = run["rows"]; ell = np.array([r["ell"] for r in R]); fr = np.array([r[mask]["frac"] for r in R])
    y = np.array([r["a"] ** 4 * r[mask]["rho_a"] for r in R])
    if np.interp(l0, ell, fr) < 0.5: return float("nan")
    s = run["s"]; tau0, tau1 = math.sqrt(math.exp(l0) / (math.sqrt(run["m2"]) * s)), math.sqrt(math.exp(l1) / (math.sqrt(run["m2"]) * s))
    tt = np.linspace(tau0, tau1, 4001); a = s * tt; H = 1 / (s * tt * tt); L = np.log(math.sqrt(run["m2"]) / H)
    integrand = 8 * math.pi * H ** 3 * F2 * L * KAPPA * L * a ** 4 * a
    I = float(np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(tt)))
    return float((np.interp(l1, ell, y) - np.interp(l0, ell, y)) / I)


def int_ratio(rp, u1, l0, l1, mask="rm2.5", nboot=4000, seed=2):
    A = np.array([eps_integral(r, l0, l1, mask) for r in rp]); B = np.array([eps_integral(r, l0, l1, mask) for r in u1])
    est = float(np.nanmean(A) / np.nanmean(B)); rng = np.random.default_rng(seed)
    bs = [np.nanmean(A[rng.integers(0, len(A), len(A))]) / np.nanmean(B[rng.integers(0, len(B), len(B))]) for _ in range(nboot)] if len(A) > 1 else []
    return dict(value=est, boot_std=float(np.std(bs)) if bs else None, eps_rp3=A.tolist(), eps_u1=B.tolist())


def spectrum_q(run, x_lo=10.0, uv_frac=0.25):
    """slope q of the instantaneous emission F(k) ~ k^-q from the last two spectral snapshots, k_phys in [x_lo H, uv_frac m]"""
    R = [r for r in run["rows"] if "spectrum" in r]
    if len(R) < 2: return None
    r1, r2 = R[-2], R[-1]
    k = 0.5 * (np.array(r2["spectrum"]["k_edges"][1:]) + np.array(r2["spectrum"]["k_edges"][:-1]))
    E1 = np.array(r1["spectrum"]["dE"]) * r1["a"] ** 4; E2 = np.array(r2["spectrum"]["dE"]) * r2["a"] ** 4
    Fk = E2 - E1; a = r2["a"]; H = r2["H"]; m = math.sqrt(run["m2"])
    kp = k / a; sel = (kp > x_lo * H) & (kp < uv_frac * m) & (Fk > 0)
    if sel.sum() < 3: return dict(q=None, n=int(sel.sum()))
    sl = np.polyfit(np.log(k[sel]), np.log(Fk[sel] / (k[1] - k[0])), 1)[0]
    return dict(q=float(-sl), n=int(sel.sum()), x_range=[float(x_lo), float(uv_frac * m / H)], ell=r2["ell"])


def main():
    sets = {
        "A": (load("RP3_N256_s*_f1_dt0.25_xi0.5.json"), load("U1_N256_s*_f1_dt0.25_xi0.5.json")),
        "B": (load("RP3_N192_*xi0.5.json"), load("U1_N192_*xi0.5.json")),
        "C": (load("RP3_N256_s21_f0.6*.json"), load("U1_N256_s21_f0.6*.json")),
        "D25": (load("RP3_N128_*dt0.25*.json"), load("U1_N128_*dt0.25*.json")),
        "D125": (load("RP3_N128_*dt0.125*.json"), load("U1_N128_*dt0.125*.json")),
        "E25": (load("RP3_N256_s31_*xi0.25.json"), load("U1_N256_s31_*xi0.25.json")),
        "E100": (load("RP3_N256_s31_*xi1.json"), load("U1_N256_s31_*xi1.json")),
    }
    avail = {k: (len(v[0]), len(v[1])) for k, v in sets.items()}
    rpA, u1A = sets["A"]
    if not rpA or not u1A:
        print("set A incomplete", avail); return
    lfin = min(r["rows"][-1]["ell"] for r in rpA + u1A)
    grid = np.linspace(3.0, lfin - 0.05, 12)
    res = dict(available=avail, window=[3.0, float(lfin)])
    res["R_xi_A"] = window_ratio(rpA, u1A, grid, "xi")
    res["R_Gamma_int_A"] = int_ratio(rpA, u1A, 3.0, lfin - 0.05)
    res["R_Gamma_A"] = window_ratio(rpA, u1A, grid, "eps_lat")
    res["eps_lat_A"] = dict(rp3=float(np.nanmean(on_grid(rpA, grid, "eps_lat"))), u1=float(np.nanmean(on_grid(u1A, grid, "eps_lat"))))
    res["eta_A"] = dict(rp3=float(np.nanmean(on_grid(rpA, grid, "eta"))), u1=float(np.nanmean(on_grid(u1A, grid, "eta"))))
    res["xi_at_4p5"] = dict(rp3=float(np.nanmean(on_grid(rpA, np.array([4.5]), "xi"))), u1=float(np.nanmean(on_grid(u1A, np.array([4.5]), "xi"))))
    res["rhoM_over_rhoa_late"] = dict(rp3=float(np.nanmean(on_grid(rpA, grid[-3:], "rhoM_over_rhoa"))), u1=float(np.nanmean(on_grid(u1A, grid[-3:], "rhoM_over_rhoa"))))
    # trend of R_Gamma in ell (bootstrap of the per-ell ratio slope)
    A = on_grid(rpA, grid, "eps_lat"); B = on_grid(u1A, grid, "eps_lat"); rng = np.random.default_rng(1); slopes = []
    for _ in range(2000):
        a_ = np.nanmean(A[rng.integers(0, len(A), len(A))], 0); b_ = np.nanmean(B[rng.integers(0, len(B), len(B))], 0)
        rr = a_ / b_; m = np.isfinite(rr)
        if m.sum() > 3: slopes.append(np.polyfit(grid[m], rr[m], 1)[0])
    rr0 = np.nanmean(A, 0) / np.nanmean(B, 0); m0 = np.isfinite(rr0)   # audit fix: polyfit over NaN never reported a trend
    sl = np.polyfit(grid[m0], rr0[m0], 1)[0]
    res["R_Gamma_trend"] = dict(slope_per_ell=float(sl), boot_std=float(np.std(slopes)) if slopes else None,
                                significant=bool(slopes and abs(sl) > 2 * np.std(slopes)))
    # systematics
    syst = {}
    def wr(rp, u1, g): return window_ratio(rp, u1, g, "eps_lat")["value"] if rp and u1 else None
    for name, key, gr in (("box_B", "B", None), ("core_resolution_C", "C", None), ("IC_xi0.25", "E25", None), ("IC_xi1.0", "E100", None)):
        rp, u1 = sets[key]
        if rp and u1:
            l2 = min(r["rows"][-1]["ell"] for r in rp + u1); g = np.linspace(3.0, l2 - 0.05, 10)
            syst[name] = dict(value=wr(rp, u1, g), A_same_window=wr(rpA, u1A, g[g <= lfin]))
            syst[name]["delta"] = abs(syst[name]["value"] - syst[name]["A_same_window"])
    if sets["D25"][0] and sets["D125"][0] and sets["D25"][1] and sets["D125"][1]:
        l2 = min(r["rows"][-1]["ell"] for r in sum(sets["D25"] + sets["D125"], [])); g = np.linspace(3.0, l2 - 0.05, 8)
        v1, v2 = wr(*sets["D25"], g), wr(*sets["D125"], g); syst["timestep_D"] = dict(dt025=v1, dt0125=v2, delta=abs(v1 - v2))
    for mk in ("rm1.5", "rm4"):
        v = window_ratio(rpA, u1A, grid, "eps_lat", mask=mk)["value"]; syst[f"mask_{mk}"] = dict(value=v, delta=abs(v - res["R_Gamma_A"]["value"]))
    for name, key in (("box_B_int", "B"), ("core_resolution_C_int", "C"), ("IC_xi0.25_int", "E25"), ("IC_xi1.0_int", "E100")):
        rp, u1 = sets[key]
        if rp and u1:
            l2 = min(r["rows"][-1]["ell"] for r in rp + u1) - 0.05
            v = int_ratio(rp, u1, 3.0, l2)["value"]; vA = int_ratio(rpA, u1A, 3.0, min(l2, lfin - 0.05))["value"]
            syst[name] = dict(value=v, A_same_window=vA, delta=abs(v - vA))
    if sets["D25"][0] and sets["D125"][0] and sets["D25"][1] and sets["D125"][1]:
        l2 = min(r["rows"][-1]["ell"] for r in sum(sets["D25"] + sets["D125"], [])) - 0.05
        v1, v2 = int_ratio(*sets["D25"], 3.0, l2)["value"], int_ratio(*sets["D125"], 3.0, l2)["value"]
        syst["timestep_D_int"] = dict(dt025=v1, dt0125=v2, delta=abs(v1 - v2))
    res["systematics"] = syst
    stat = res["R_Gamma_A"]["boot_std"] or 0.0
    prim = [k for k in syst if not k.endswith("_int")]
    tot = math.sqrt(stat ** 2 + sum(syst[k]["delta"] ** 2 for k in prim if syst[k].get("delta") is not None))
    d_est = abs(res["R_Gamma_int_A"]["value"] - res["R_Gamma_A"]["value"])
    if d_est > tot: tot = math.sqrt(tot ** 2 + d_est ** 2)
    res["R_Gamma"] = dict(value=res["R_Gamma_A"]["value"], stat=stat, total=tot, estimator_difference=d_est,
                          integral_estimator=res["R_Gamma_int_A"]["value"], integral_estimator_boot_std=res["R_Gamma_int_A"]["boot_std"])
    res["spectra_q"] = {f"{r['model']}_s{r['seed']}": spectrum_q(r) for r in rpA + u1A}
    res["U1_control_gate"] = dict(eta_in_0p3_3=bool(0.3 <= res["eta_A"]["u1"] <= 3), xi_4p5_in_0p5_2p5=bool(0.5 <= res["xi_at_4p5"]["u1"] <= 2.5))
    res["U1_control_gate"]["passed"] = all(res["U1_control_gate"].values())
    (HERE / "results" / "u2" / "U2_RESULT.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: res[k] for k in ("available", "window", "R_Gamma", "eps_lat_A", "eta_A", "xi_at_4p5", "rhoM_over_rhoa_late", "R_Gamma_trend", "U1_control_gate")}, indent=1, default=float))
    print("R_xi", res["R_xi_A"]["value"], res["R_xi_A"]["boot_std"])
    print("systematics", json.dumps(syst, default=float))
    print("q", res["spectra_q"])
    figure(sets, grid, res)


def figure(sets, grid, res):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for key, ls in (("A", "-"), ("B", "--"), ("C", ":"), ("E25", "-."), ("E100", "-.")):
        rp, u1 = sets[key]
        for runs, col, lab in ((rp, "C3", "RP³"), (u1, "C0", "U(1)")):
            for i, r in enumerate(runs):
                d = derived(r); m = d["ell"] > 1.0
                ax[0].plot(d["ell"][m], d["xi"][m], color=col, ls=ls, lw=0.8, alpha=0.7, label=f"{lab} {key}" if i == 0 else None)
                ok = d["ok"]; ax[1].plot(d["ell"][ok], d["eps_lat"][ok], color=col, ls=ls, lw=0.8, alpha=0.7)
    ax[0].axhline(1.20, color="gray", lw=0.8, ls="--"); ax[0].plot(grid, 0.24 * grid, color="k", lw=0.8)
    ax[0].text(3.05, 1.25, "Hindmarsh+ ζ* = 1.20", fontsize=7, color="gray"); ax[0].text(4.0, 0.24 * 4.0 + 0.05, "0.24 ℓ", fontsize=7)
    ax[0].set_xlabel("ℓ = ln(m/H)"); ax[0].set_ylabel("ξ"); ax[0].legend(fontsize=6, ncol=2); ax[0].set_xlim(1, None)
    ax[1].axhline(1, color="k", lw=0.8); ax[1].set_xlabel("ℓ"); ax[1].set_ylabel("ε_red = Γ/(8πH³f²ℓ·0.24ℓ)")
    ax[1].set_yscale("log")
    pe = np.array(res["R_Gamma_A"]["per_ell"]); px = np.array(res["R_xi_A"]["per_ell"])
    ax[2].plot(grid, pe, "o-", color="C2", label="R_Γ (emisión)"); ax[2].plot(grid, px, "s--", color="C4", label="R_ξ (densidad)")
    v, t = res["R_Gamma"]["value"], res["R_Gamma"]["total"]
    ax[2].axhspan(v - t, v + t, color="C2", alpha=0.15); ax[2].set_xlabel("ℓ"); ax[2].set_ylabel("RP³ / U(1)")
    ax[2].set_title(f"R_Γ = {v:.2f} ± {t:.2f}", fontsize=9); ax[2].legend(fontsize=8)
    for a_ in ax: a_.grid(alpha=0.3, lw=0.3)
    fig.suptitle("U2 · red RP³ frente a U(1), mismo código, misma f y misma masa de núcleo", fontsize=10)
    fig.tight_layout(); (HERE / "figures").mkdir(exist_ok=True); fig.savefig(HERE / "figures" / "u2_red.png", dpi=140); plt.close(fig)


if __name__ == "__main__":
    main()
