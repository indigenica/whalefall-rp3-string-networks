"""U6 analysis (docs/preregistros/ESPEC_U6.md): RP3 at 576^3 with seeds 61-66 and U(1) with seeds 61-64.
R_Gamma and R_xi in the U4b window, and the instantaneous emission slope q(ell) of RP3 and U(1) from consecutive spectral
snapshots (same estimator as analyze.spectrum_q).  Output: results/u4/U6_RESULT.json"""
import glob, json, math
from pathlib import Path
import numpy as np
import analyze as A2

HERE = Path(__file__).resolve().parent; RUNS = HERE / "results" / "u4" / "runs"; N = 576


def qseries(run, x_lo=10.0, uv_frac=0.25):
    R = [r for r in run["rows"] if "spectrum" in r]; m = math.sqrt(run["m2"]); out = []
    for r1, r2 in zip(R[:-1], R[1:]):
        k = 0.5 * (np.array(r2["spectrum"]["k_edges"][1:]) + np.array(r2["spectrum"]["k_edges"][:-1]))
        Fk = np.array(r2["spectrum"]["dE"]) * r2["a"] ** 4 - np.array(r1["spectrum"]["dE"]) * r1["a"] ** 4
        sel = (k / r2["a"] > x_lo * r2["H"]) & (k / r2["a"] < uv_frac * m) & (Fk > 0)
        if sel.sum() >= 8:                                     # at least 8 k-bins in the fit range
            out.append((0.5 * (r1["ell"] + r2["ell"]), float(-np.polyfit(np.log(k[sel]), np.log(Fk[sel]), 1)[0])))
    return out


def q_fit(runs, l0=4.6, nboot=4000, seed=5):
    """q(ell) = q0 + q1 (ell - 5) over ell >= l0, seeds pooled; bootstrap over seeds"""
    S = [np.array(qseries(r)) for r in runs]; S = [s[s[:, 0] >= l0] for s in S if len(s)]
    def fit(ss):
        X = np.concatenate(ss); return np.polyfit(X[:, 0] - 5.0, X[:, 1], 1)[::-1]
    q0, q1 = fit(S); rng = np.random.default_rng(seed)
    bs = np.array([fit([S[i] for i in rng.integers(0, len(S), len(S))]) for _ in range(nboot)]) if len(S) > 1 else None
    lend = min(r["rows"][-1]["ell"] for r in runs)
    return dict(q0_at_ell5=float(q0), q1_per_ell=float(q1), q0_std=float(bs[:, 0].std()) if bs is not None else None,
                q1_std=float(bs[:, 1].std()) if bs is not None else None, q_at_ell_end=float(q0 + q1 * (lend - 5)), ell_end=lend,
                per_seed=[s.tolist() for s in S], n_seeds=len(S))


def main():
    rp = [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / f"RP3_N{N}_s6*_xi0.5.json")))]
    u1 = [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / f"U1_N{N}_s6*_xi0.5.json")))]
    l_end = min(r["rows"][-1]["ell"] for r in rp + u1); l_ref = l_end - 0.3; grid = np.linspace(l_ref - 1.0, l_end - 0.05, 10)
    out = dict(spec="docs/preregistros/ESPEC_U6.md", N=N, n_rp3=len(rp), n_u1=len(u1), window=[l_ref - 1.0, l_end - 0.05],
               R_Gamma=A2.window_ratio(rp, u1, grid, "eps_lat"), R_Gamma_int=A2.int_ratio(rp, u1, l_ref - 1.0, l_end - 0.05),
               R_xi=A2.window_ratio(rp, u1, grid, "xi"), q_RP3=q_fit(rp), q_U1=q_fit(u1))
    xi_slope = lambda runs: [float(np.polyfit(*np.array([(x["ell"], x["xi"]) for x in r["rows"] if l_ref - 0.7 <= x["ell"] <= l_ref]).T, 1)[0]) for r in runs]
    out["xi_slope"] = dict(RP3=xi_slope(rp), U1=xi_slope(u1))
    q1, s1 = out["q_RP3"]["q1_per_ell"], out["q_RP3"]["q1_std"]
    out["q_RP3_trend"] = ("CRECE" if s1 and q1 > 2 * s1 else "DECRECE" if s1 and q1 < -2 * s1 else "SIN_TENDENCIA_SIGNIFICATIVA")
    (HERE / "results" / "u4" / "U6_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"RP3 {len(rp)} / U1 {len(u1)} runs; R_Gamma = {out['R_Gamma']['value']:.3f} +- {out['R_Gamma']['boot_std']:.3f}; "
          f"R_xi = {out['R_xi']['value']:.3f}; q_RP3(5) = {out['q_RP3']['q0_at_ell5']:.3f} +- {out['q_RP3']['q0_std']}, "
          f"dq/dell = {q1:.3f} +- {s1} ({out['q_RP3_trend']}); q_RP3(ell_end) = {out['q_RP3']['q_at_ell_end']:.3f}; "
          f"q_U1(5) = {out['q_U1']['q0_at_ell5']:.3f}")
    return out


if __name__ == "__main__":
    main()
