"""U9 analysis (docs/preregistros/ESPEC_U9.md): dense (xi_i = 2) vs sparse (xi_i = 0.5) RP3 networks at 576^3.
Output: results/u4/U9_RESULT.json"""
import glob, json
from pathlib import Path
import numpy as np
import analyze as A2
from u6_analyze import q_fit

HERE = Path(__file__).resolve().parent; RUNS = HERE / "results" / "u4" / "runs"; N = 576


def load(pat): return [json.loads(Path(p).read_text()) for p in sorted(glob.glob(str(RUNS / pat)))]


def xi_at(R, l): return float(np.interp(l, [x["ell"] for x in R["rows"]], [x["xi"] for x in R["rows"]]))


def slope(R, l0, l1):
    a = np.array([(x["ell"], x["xi"]) for x in R["rows"] if l0 <= x["ell"] <= l1]); return float(np.polyfit(a[:, 0], a[:, 1], 1)[0])


def main():
    rp_s, rp_d = load(f"RP3_N{N}_s6*_xi0.5.json"), load(f"RP3_N{N}_s6*_xi2.json")
    u1_s, u1_d = load(f"U1_N{N}_s6*_xi0.5.json"), load(f"U1_N{N}_s6*_xi2.json")
    if not rp_d: raise SystemExit("no dense RP3 runs yet")
    l_end = min(r["rows"][-1]["ell"] for r in rp_s + rp_d + u1_s + u1_d); l_ref = l_end - 0.3; w = (l_ref - 0.7, l_ref)
    xs, xd = np.mean([xi_at(r, l_ref) for r in rp_s]), np.mean([xi_at(r, l_ref) for r in rp_d])
    M = abs(xd - xs) / (0.5 * (xd + xs))
    ss, sd = [slope(r, *w) for r in rp_s], [slope(r, *w) for r in rp_d]
    s_s, s_d = float(np.mean(ss)), float(np.mean(sd))
    cls = ("TRANSITORIO_DESDE_ABAJO" if M < 0.2 and s_d <= 0.5 * s_s else
           "CRECIMIENTO_GENUINO" if M < 0.2 and s_d >= 0.7 * s_s else "INDETERMINADO")
    out = dict(spec="docs/preregistros/ESPEC_U9.md", N=N, ell_ref=l_ref, slope_window=list(w), n=dict(rp3_sparse=len(rp_s), rp3_dense=len(rp_d), u1_sparse=len(u1_s), u1_dense=len(u1_d)),
               xi_RP3_ref=dict(sparse=float(xs), dense=float(xd), memory=float(M)), slope_RP3=dict(sparse=ss, dense=sd, mean_sparse=s_s, mean_dense=s_d),
               classification=cls,
               xi_RP3_tracks=dict(dense=[[(x["ell"], x["xi"]) for x in r["rows"]] for r in rp_d]))
    if u1_d:
        xu_s, xu_d = np.mean([xi_at(r, l_ref) for r in u1_s]), np.mean([xi_at(r, l_ref) for r in u1_d])
        out["U1_memory"] = float(abs(xu_d - xu_s) / (0.5 * (xu_d + xu_s)))
        grid = np.linspace(l_ref - 1.0, l_end - 0.05, 10)
        out["R_Gamma_dense"] = A2.window_ratio(rp_d, u1_d, grid, "eps_lat"); out["R_xi_dense"] = A2.window_ratio(rp_d, u1_d, grid, "xi")
    out["q_RP3_dense"] = q_fit(rp_d)
    (HERE / "results" / "u4" / "U9_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"U9: memory RP3 = {M:.3f} (sparse {xs:.3f}, dense {xd:.3f}); slope sparse {s_s:.3f}, dense {s_d:.3f} -> {cls}")
    if u1_d:
        print(f"    U1 memory = {out['U1_memory']:.3f}; R_Gamma(dense) = {out['R_Gamma_dense']['value']:.3f} +- {out['R_Gamma_dense']['boot_std']}; "
              f"R_xi(dense) = {out['R_xi_dense']['value']:.3f}; q_RP3(dense, ell=5) = {out['q_RP3_dense']['q0_at_ell5']:.3f}")
    return out


if __name__ == "__main__":
    main()
