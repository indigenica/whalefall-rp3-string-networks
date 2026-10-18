"""U10a, post-hoc informative numbers (ESPEC_U10a, amendment 2; written AFTER the production data, not part of the verdict).
Production ran at 512^3 (RAM planner) instead of 576^3, so the undrained references (U6, U9) are at 576^3 and the comparisons
with them use the same ABSOLUTE ell window as the 512^3 runs ([4.25, 5.50]).
    python -B u10a_extra.py   -> results/u4/U10A_EXTRA.json"""
import glob, json
from pathlib import Path
import numpy as np
import analyze as A2
from u6_analyze import q_fit
from u10a_analyze import gamma_and_X

HERE = Path(__file__).resolve().parent; RUNS = HERE / "results" / "u4" / "runs"


def load(p): return [json.loads(Path(x).read_text()) for x in sorted(glob.glob(str(RUNS / p)))]


def main():
    S = {k: load(p) for k, p in dict(Rdd="RP3_N512_s6*_xi2_dr.json", Rds="RP3_N512_s6*_xi0.5_dr.json", Udd="U1_N512_s6*_xi2_dr.json",
                                     Uds="U1_N512_s6*_xi0.5_dr.json", Rux="RP3_N512_s61_xi2_x.json", Uux="U1_N512_s61_xi2_x.json",
                                     R6="RP3_N576_s6*_xi0.5.json", U6="U1_N576_s6*_xi0.5.json", R9="RP3_N576_s6*_xi2.json", U9="U1_N576_s6*_xi2.json").items()}
    le = min(r["rows"][-1]["ell"] for k in ("Rdd", "Rds", "Udd", "Uds", "Rux", "Uux") for r in S[k]); g = np.linspace(le - 1.3, le - 0.05, 10)
    W = lambda a, b, key="eps_lat", grid=g: A2.window_ratio(S[a], S[b], grid, key)
    out = dict(note="post-hoc, informative; ESPEC_U10a amendment 2", N=512, window=[g[0], g[-1]])
    out["R_Gamma"] = {f"{a}/{b}": dict(value=W(a, b)["value"], boot_std=W(a, b)["boot_std"], R_xi=W(a, b, "xi")["value"])
                      for a, b in (("Rdd", "Udd"), ("Rds", "Uds"), ("Rux", "Uux"), ("R6", "U6"), ("R9", "U9"))}
    out["R_Gamma_per_seed"] = {f"{k}_s{sd}": A2.window_ratio([r for r in S["R" + k] if r["seed"] == sd], [r for r in S["U" + k] if r["seed"] == sd], g, "eps_lat")["value"]
                               for k in ("dd", "ds") for sd in (61, 62)}
    out["emission_ratio"] = {f"{a}/{b}": W(a, b)["value"] for a, b in (("Rdd", "Rux"), ("Udd", "Uux"), ("Rds", "R6"), ("Uds", "U6"), ("Rdd", "R9"), ("Udd", "U9"))}
    q = {}
    for k in ("Rdd", "Rds", "Rux", "R6", "R9"):
        f = q_fit(S[k]); q[k] = dict(q0_at_ell5=f["q0_at_ell5"], q0_std=f["q0_std"], q_end=f["q_at_ell_end"], dq_dell=f["q1_per_ell"])
    out["q_RP3"] = q
    m = {}
    for k in S:
        v = [[float(np.interp(l, *(lambda d: (d["ell"], d["rhoM_over_rhoa"]))(A2.derived(r, "rm2.5")))) for l in (4.3, 4.8, 5.3)] for r in S[k]]
        m[k] = np.mean(v, 0).tolist()
    out["rhoM_over_rhoa_at_ell_4.3_4.8_5.3"] = m
    fx = {}
    for k in ("Udd", "Uds", "Uux", "Rux"):
        for r in S[k]:
            G, X = gamma_and_X(r, g); fx[f"{k}_s{r['seed']}"] = float(np.nanmean(X) / np.nanmean(G))
    out["f_X"] = fx
    out["U1_dense_two_methods"] = dict(retained_after_drain=out["emission_ratio"]["Udd/Uux"], one_minus_fX_undrained=1 - fx["Uux_s61"])
    (HERE / "results" / "u4" / "U10A_EXTRA.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=lambda x: round(float(x), 3)))
    return out


if __name__ == "__main__":
    main()
