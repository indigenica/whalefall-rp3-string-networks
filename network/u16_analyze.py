"""U16 analysis (docs/preregistros/ESPEC_U16.md), with the estimators of U6 (analyze.on_grid / window_ratio / int_ratio).
Window [ell_end - 1.3, ell_end - 0.05], ell_end = ln 96, 10-point grid.  Resolution arm: per seed, matched RP3 and U(1) of the same
seed; Delta_X = X(r = 0.6) - X(r = 1) per seed, sigma = standard error of the mean of the three differences.  Volume arm: r = 1,
L = 3.33 tau_end against L = 2 tau_end, unpaired.  Output: network/results/u16/U16_RESULT.json"""
import glob, json, math
from pathlib import Path
import numpy as np
import analyze as A2

HERE = Path(__file__).resolve().parent; RUNS = HERE / "results" / "u16" / "runs"
L_END = math.log(96.0); L0, L1 = L_END - 1.3, L_END - 0.05; GRID = np.linspace(L0, L1, 10)
RES = [(192, 96, 1.0), (256, 128, 0.75), (320, 160, 0.6)]; VOL = (320, 96, 1.0)


def load(model, N, te, r, arm, seed=None):
    pat = f"{model}_N{N}_t{te}_r{r:g}_s{seed if seed else '*'}_{arm}.json"
    return [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / pat)))]


def wmean(run, key):
    return float(np.nanmean(A2.on_grid([run], GRID, key, "rm2.5")))


def stats_of(rp, u1):
    return dict(R_xi=A2.window_ratio(rp, u1, GRID, "xi")["value"], R_Gamma=A2.window_ratio(rp, u1, GRID, "eps_lat")["value"],
                R_Gamma_int=A2.int_ratio(rp, u1, L0, L1)["value"],
                xi_RP3=float(np.mean([wmean(r, "xi") for r in rp])), xi_U1=float(np.mean([wmean(r, "xi") for r in u1])),
                eps_RP3=float(np.mean([wmean(r, "eps_lat") for r in rp])), eps_U1=float(np.mean([wmean(r, "eps_lat") for r in u1])),
                n_rp3=len(rp), n_u1=len(u1))


def classify(delta, sigma, ref, names):
    rel, srel = abs(delta) / abs(ref), sigma / abs(ref)
    if rel < 0.05 and srel < 0.05: return names[0]
    if rel >= 0.05 and abs(delta) >= 2 * sigma: return names[1]
    return names[2]


def main():
    out = dict(spec="docs/preregistros/ESPEC_U16.md", window=[L0, L1], resolution={}, per_seed={}, volume={})
    for N, te, r in RES:
        rp, u1 = load("RP3", N, te, r, "res"), load("U1", N, te, r, "res")
        if rp and u1: out["resolution"][f"r{r:g}"] = dict(N=N, tau_end=te, **stats_of(rp, u1))
    seeds = (61, 62, 63); keys = ("R_xi", "R_Gamma", "R_Gamma_int", "xi_RP3", "xi_U1", "eps_RP3", "eps_U1")
    for sd in seeds:
        row = {}
        for N, te, r in RES:
            rp, u1 = load("RP3", N, te, r, "res", sd), load("U1", N, te, r, "res", sd)
            if rp and u1: row[f"r{r:g}"] = stats_of(rp, u1)
        out["per_seed"][str(sd)] = row
    complete = [sd for sd in seeds if all(f"r{r:g}" in out["per_seed"][str(sd)] for _, _, r in RES)]
    out["seeds_complete"] = complete
    if len(complete) >= 2:
        D = {}
        for k in keys:
            d = np.array([out["per_seed"][str(sd)]["r0.6"][k] - out["per_seed"][str(sd)]["r1"][k] for sd in complete])
            ref = out["resolution"]["r1"][k]
            D[k] = dict(delta=float(d.mean()), sigma=float(d.std(ddof=1) / math.sqrt(len(d))), per_seed=d.tolist(), ref_r1=ref,
                        rel=float(d.mean() / ref))
            rr = np.array([1.0, 0.75, 0.6]); y = np.array([out["resolution"][f"r{x:g}"][k] for x in rr])
            D[k]["X0_r2"] = float(np.polyfit(rr ** 2, y, 1)[1]); D[k]["X0_r1"] = float(np.polyfit(rr, y, 1)[1])
        out["delta_resolution"] = D
        names = ("RESOLUCION_ROBUSTA", "SESGO_DE_RESOLUCION", "INDETERMINADO")
        out["verdict_resolution"] = {k: classify(D[k]["delta"], D[k]["sigma"], D[k]["ref_r1"], names) for k in ("R_xi", "R_Gamma")}
    # amendment 1 (informative, no new verdict): fourth resolution r = 0.5
    rp5, u15 = load("RP3", 384, 192, 0.5, "res"), load("U1", 384, 192, 0.5, "res")
    if rp5 and u15 and len(complete) >= 2:
        out["resolution"]["r0.5"] = dict(N=384, tau_end=192, **stats_of(rp5, u15))
        c5 = [sd for sd in complete if load("RP3", 384, 192, 0.5, "res", sd) and load("U1", 384, 192, 0.5, "res", sd)]
        A5 = {}
        for sd in c5: out["per_seed"][str(sd)]["r0.5"] = stats_of(load("RP3", 384, 192, 0.5, "res", sd), load("U1", 384, 192, 0.5, "res", sd))
        for k in keys:
            d = np.array([out["per_seed"][str(sd)]["r0.5"][k] - out["per_seed"][str(sd)]["r1"][k] for sd in c5])
            rr = np.array([1.0, 0.75, 0.6, 0.5]); y = np.array([out["resolution"][f"r{x:g}"][k] for x in rr])
            A5[k] = dict(delta=float(d.mean()), sigma=float(d.std(ddof=1) / math.sqrt(len(d))) if len(d) > 1 else None, per_seed=d.tolist(),
                         rel=float(d.mean() / out["resolution"]["r1"][k]), X0_r2_4pt=float(np.polyfit(rr ** 2, y, 1)[1]),
                         X0_r1_4pt=float(np.polyfit(rr, y, 1)[1]), values_by_r=y.tolist())
        out["amendment1_r0.5"] = dict(seeds=c5, delta=A5, classification_informative={k: classify(A5[k]["delta"], A5[k]["sigma"] or 0.0, out["resolution"]["r1"][k],
                                      ("RESOLUCION_ROBUSTA", "SESGO_DE_RESOLUCION", "INDETERMINADO")) for k in ("R_xi", "R_Gamma")},
                                      R_Gamma_monotonic_0p75_to_0p5=bool(out["resolution"]["r0.5"]["R_Gamma"] <= out["resolution"]["r0.6"]["R_Gamma"] <= out["resolution"]["r0.75"]["R_Gamma"]
                                                                          or out["resolution"]["r0.5"]["R_Gamma"] >= out["resolution"]["r0.6"]["R_Gamma"] >= out["resolution"]["r0.75"]["R_Gamma"]))
    rpv, u1v = load("RP3", *VOL, "vol"), load("U1", *VOL, "vol")
    if rpv and u1v and "r1" in out["resolution"]:
        rp1, u11 = load("RP3", 192, 96, 1.0, "res"), load("U1", 192, 96, 1.0, "res")
        V = {}
        for k in ("R_xi", "R_Gamma", "R_Gamma_int"):
            big, small = stats_of(rpv, u1v)[k], stats_of(rp1, u11)[k]
            key = {"R_xi": "xi", "R_Gamma": "eps_lat"}.get(k)
            if key:   # sigma from the seed scatter of each arm (ratio of seed-mean ratios, propagated)
                def sd_ratio(rp, u1):
                    a = np.array([wmean(r, key) for r in rp]); b = np.array([wmean(r, key) for r in u1])
                    return (a.mean() / b.mean()) * math.sqrt((a.std(ddof=1) / a.mean()) ** 2 / len(a) + (b.std(ddof=1) / b.mean()) ** 2 / len(b))
                sig = math.hypot(sd_ratio(rpv, u1v), sd_ratio(rp1, u11))
            else:
                sig = float("nan")
            V[k] = dict(L333=big, L2=small, delta=big - small, sigma=sig)
        out["volume"] = V
        names = ("VOLUMEN_ROBUSTO", "SESGO_DE_VOLUMEN", "INDETERMINADO")
        out["verdict_volume"] = {k: classify(V[k]["delta"], V[k]["sigma"], V[k]["L2"], names) for k in ("R_xi", "R_Gamma")}
    def finite(x):                                            # strict JSON: non-finite values (no estimate) are null
        if isinstance(x, dict): return {k: finite(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)): return [finite(v) for v in x]
        if isinstance(x, (float, np.floating)): return float(x) if math.isfinite(x) else None
        if isinstance(x, np.integer): return int(x)
        return x
    out = finite(out)
    (HERE / "results" / "u16" / "U16_RESULT.json").write_text(json.dumps(out, indent=1, allow_nan=False))
    print(json.dumps({k: out.get(k) for k in ("resolution", "delta_resolution", "verdict_resolution", "volume", "verdict_volume", "amendment1_r0.5")},
                     indent=1, default=float))


if __name__ == "__main__":
    main()
