"""U25 (docs/preregistros/ESPEC_U25.md), fixed before the data: time step.  dtau = 0.125 runs against their dtau = 0.25
partners, the U24 D0 networks (same seeds, same initial conditions and snapshot times).  Same window and estimators as U24;
Delta_X = X(0.125) / X(0.25) - 1 with a paired bootstrap.  Enmienda 1: the drain is first order in dtau (kern.c), so the bias of
dtau = 0.25 is b = 2|Delta| (first-order Richardson); (4/3)|Delta| (second order) is reported for reference.
    python -B network/u25_analyze.py      ->  network/results/u25/U25_RESULT.json"""
import json
from pathlib import Path
import numpy as np
import u24_analyze as U

HERE = Path(__file__).resolve().parent; D = HERE / "results" / "u25"


def main():
    half = (U.load(D / "runs", "RP3", "dt0125"), U.load(D / "runs", "U1", "dt0125"))
    full = (U.load(U.D / "runs", "RP3", "D0"), U.load(U.D / "runs", "U1", "D0"))
    seeds = [s for s in U.SEEDS if all(s in x for x in half + full)]
    if not seeds: raise SystemExit("no seed with both time steps in both models yet")
    T = U.tables({"dt0125": half, "dt025": full}, seeds)
    est, delta = U.paired(T, "dt025", ["dt0125"], seeds, seed=1); d = delta["dt0125"]
    out = dict(spec="docs/preregistros/ESPEC_U25.md", window=[float(U.GRID[0]), float(U.GRID[-1])], seeds=seeds,
               provisional=len(seeds) < len(U.SEEDS), values=est, delta=d,
               bias_dt025_first_order={k: dict(b=2 * abs(x["delta"]), sigma_b=(2 * x["sigma"] if x["sigma"] is not None else None)) for k, x in d.items()},
               bias_dt025_second_order={k: 4 / 3 * abs(x["delta"]) for k, x in d.items()})
    g, x = d["R_Gamma"], d["R_xi"]
    if len(seeds) < 2:                                   # ESPEC_U24 enmienda 1, also for U25: no sigma, no verdict with one seed
        out["verdict_U25"] = f"SIN_VEREDICTO (n = {len(seeds)})"; out["systematic_for_budget"] = None
        (D / "U25_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
        print(f"U25 seeds {seeds}: Delta R_Gamma {g['delta']:+.4f}, Delta R_xi {x['delta']:+.4f}; {out['verdict_U25']}"); return out
    bg, bx = 2 * abs(g["delta"]) + 4 * g["sigma"], 2 * abs(x["delta"]) + 4 * x["sigma"]                     # b + 2 sigma_b (enmienda 1)
    if bg <= 0.02 and bx <= 0.01: v = "PASO_ROBUSTO"
    elif (abs(g["delta"]) > 0.015 and abs(g["delta"]) > 2 * g["sigma"]) or (abs(x["delta"]) > 0.0075 and abs(x["delta"]) > 2 * x["sigma"]): v = "PASO_SENSIBLE"
    else: v = "INDETERMINADO"
    out["systematic_for_budget"] = (None if v == "PASO_ROBUSTO" else
                                   dict(R_Gamma=out["bias_dt025_first_order"]["R_Gamma"], R_xi=out["bias_dt025_first_order"]["R_xi"]))
    out["verdict_U25"] = v + (" (provisional: faltan semillas)" if out["provisional"] else "")
    # informative: where the change sits in the window (per grid point, seed means)
    (eRh, eUh), (eRf, eUf) = T["dt0125"]["eps_lat"], T["dt025"]["eps_lat"]
    out["per_ell"] = dict(grid=U.GRID.tolist(), eps_RP3=(np.nanmean(eRh, 0) / np.nanmean(eRf, 0) - 1).tolist(),
                          eps_U1=(np.nanmean(eUh, 0) / np.nanmean(eUf, 0) - 1).tolist())
    (D / "U25_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"U25 seeds {seeds}{' (PROVISIONAL)' if out['provisional'] else ''}")
    for k in ("R_Gamma", "R_xi", "eps_RP3", "eps_U1", "xi_RP3", "xi_U1"):
        print(f"  {k:8s} Delta {d[k]['delta']:+.4f} +- {d[k]['sigma']:.4f}   (bias of dtau = 0.25: first order {out['bias_dt025_first_order'][k]['b']:.4f}, second {out['bias_dt025_second_order'][k]:.4f})")
    print("U25:", out["verdict_U25"])
    return out


if __name__ == "__main__":
    main()
