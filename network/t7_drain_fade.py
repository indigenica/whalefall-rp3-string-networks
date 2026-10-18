"""Task 7 (docs/PLAN_NUBE_2026-10-11.md): does the initial massive bath fade within the window?  Informative, no verdict.

For the U11 seeds, compares each drained network with its undrained U6 twin along the window [ell_fin - 1.3, ell_fin - 0.05]:
  - emission ratio drained/undrained per ell (estimator of U6/U11: Gamma = H a^-4 d(a^4 rho_a)/d ln tau, centred differences, mask rm2.5),
    and its slope in ell per seed;
  - rho_M/rho_a outside the cores, drained and undrained.
Also prints the tree-level decay epoch ell_dec = ln(m/Gamma) of the massive Q modes from results/s3_sectors.json.

    python3 network/t7_drain_fade.py  ->  network/results/u4/T7_DRAIN_FADE.json"""
import json, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]; RUNS = ROOT / "network" / "results" / "u4" / "runs"


def series(p):
    R = json.loads(p.read_text())["rows"]
    g = lambda k: np.array([r[k] for r in R])
    return g("tau"), g("a"), g("H"), g("ell"), np.array([r["rm2.5"]["rho_a"] for r in R]), np.array([r["rm2.5"]["rho_M"] for r in R])


def emission(p):
    tau, a, H, ell, rho, rhoM = series(p); Y = a ** 4 * rho; x = np.log(tau)
    d = np.full(len(x), np.nan); d[1:-1] = (Y[2:] - Y[:-2]) / (x[2:] - x[:-2])
    return ell, H * d / a ** 4


def main():
    lf = math.log(576 / 2); grid = np.linspace(lf - 1.3, lf - 0.05, 10); out = dict(window=[grid[0], grid[-1]], models={})
    for mod, seeds in (("RP3", range(61, 67)), ("U1", range(61, 65))):
        rat, slopes, rm = [], [], {"drained": [], "undrained": []}
        for s in seeds:
            pd, pu = RUNS / f"{mod}_N576_s{s}_xi0.5_dr.json", RUNS / f"{mod}_N576_s{s}_xi0.5.json"
            (ld, gd), (lu, gu) = emission(pd), emission(pu); okd, oku = np.isfinite(gd), np.isfinite(gu)
            r = np.interp(grid, ld[okd], gd[okd], right=np.nan) / np.interp(grid, lu[oku], gu[oku], right=np.nan)
            ok = np.isfinite(r); rat.append(r); slopes.append(float(np.polyfit(grid[ok], r[ok], 1)[0]))
            for lab, p in (("drained", pd), ("undrained", pu)):
                _, _, _, ell, rho, rhoM = series(p); rm[lab].append(np.interp(grid, ell, rhoM / rho))
        sl = np.array(slopes)
        out["models"][mod] = dict(ratio_by_ell=[float(v) if np.isfinite(v) else None for v in np.nanmean(rat, 0)], slope_per_seed=slopes, slope_mean=float(sl.mean()),
                                  slope_sem=float(sl.std(ddof=1) / math.sqrt(len(sl))),
                                  rhoM_over_rhoa={k: np.mean(v, 0).tolist() for k, v in rm.items()})
        print(f"{mod}: drained/undrained from {np.nanmean(rat, 0)[0]:.2f} to {np.nanmean(rat, 0)[-2]:.2f}; "
              f"slope {sl.mean():+.3f} ± {out['models'][mod]['slope_sem']:.3f} per unit ell; rho_M/rho_a undrained "
              f"{np.mean(rm['undrained'], 0)[0]:.2f} -> {np.mean(rm['undrained'], 0)[-1]:.2f}, drained {np.mean(rm['drained'], 0)[0]:.2f} -> {np.mean(rm['drained'], 0)[-1]:.2f}")
    s3 = json.loads((ROOT / "results" / "s3_sectors.json").read_text())
    out["tree_level_decay"] = {sp["name"]: dict(m=sp["mass_code"], Gamma=sp["width_code"], ell_dec=math.log(sp["mass_code"] / sp["width_code"]))
                               for sp in s3["species"] if sp["width_code"] > 0}
    for k, v in out["tree_level_decay"].items(): print(f"{k}: Gamma/m = {v['Gamma'] / v['m']:.2e}, ell_dec = ln(m/Gamma) = {v['ell_dec']:.1f}")
    (ROOT / "network" / "results" / "u4" / "T7_DRAIN_FADE.json").write_text(json.dumps(out, indent=1, allow_nan=False))


if __name__ == "__main__":
    main()
