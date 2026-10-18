"""U10b analysis (docs/preregistros/ESPEC_U10b.md), fixed before the data: dense networks started late (ell_i = 3.5, approach
from above, '_ab') against the drained sparse networks of U11 (approach from below), both at 576^3 and drained until ell = 3.8.
    python -B u10b_analyze.py        -> results/u4/U10B_RESULT.json"""
import glob, json, sys
from pathlib import Path
import numpy as np
import analyze as A2
from u6_analyze import q_fit

HERE = Path(__file__).resolve().parent
RUNS = Path(sys.argv[sys.argv.index("--runs") + 1]) if "--runs" in sys.argv else HERE / "results" / "u4" / "runs"
N = int(sys.argv[sys.argv.index("--N") + 1]) if "--N" in sys.argv else 576
W0 = float(sys.argv[sys.argv.index("--W0") + 1]) if "--W0" in sys.argv else 4.36   # --W0 only for small test boxes


def load(p): return [json.loads(Path(x).read_text()) for x in sorted(glob.glob(str(RUNS / p)))]


def mean_xi(runs, grid):
    return np.mean([np.interp(grid, [x["ell"] for x in r["rows"]], [x["xi"] for x in r["rows"]]) for r in runs], axis=0)


def model_block(ab, low, l_end):
    e = np.array([x["ell"] for x in ab[0]["rows"]]); g = e[(e >= W0) & (e <= l_end - 0.05 + 1e-9)]
    xa, xl = mean_xi(ab, g), mean_xi(low, g)
    sa, sl = float(np.polyfit(g, xa, 1)[0]), float(np.polyfit(g, xl, 1)[0])
    frac_down = float(np.mean(np.diff(xa) < 0))
    le = l_end - 0.05; xa_end, xl_end = float(mean_xi(ab, [le])[0]), float(mean_xi(low, [le])[0])
    D = xa_end / xl_end; GA = bool(sa < 0 and frac_down > 0.5)
    verdict = ("CONVERGENCIA" if D < 1.15 else "ATRACTOR_CRECIENTE" if sa >= 0 else
               "ATRACTOR_ACOTADO" if GA and sl > 0 else "INDETERMINADO")
    lw = g[g >= l_end - 0.75]                                 # informative: crossing of the two straight lines (last 0.7 in ell)
    pa, pl = np.polyfit(lw, mean_xi(ab, lw), 1), np.polyfit(lw, mean_xi(low, lw), 1)
    cross = None if abs(pa[0] - pl[0]) < 1e-9 else float((pl[1] - pa[1]) / (pa[0] - pl[0]))
    per_xi = {}
    for xi in sorted({r["xi_init"] for r in ab}):
        rr = [r for r in ab if r["xi_init"] == xi]
        per_xi[f"xi_i={xi:g}"] = dict(n=len(rr), xi_end=float(mean_xi(rr, [le])[0]), slope=float(np.polyfit(g, mean_xi(rr, g), 1)[0]))
    return dict(n_ab=len(ab), n_low=len(low), ell=g.tolist(), xi_ab=xa.tolist(), xi_low=xl.tolist(), slope_ab=sa, slope_low=sl,
                frac_down_ab=frac_down, gate_GA=GA, xi_ab_end=xa_end, xi_low_end=xl_end, Delta_end=D, verdict=verdict,
                lines_cross_ell=cross, lines_cross_xi=(float(np.polyval(pa, cross)) if cross is not None else None), per_xi_init=per_xi)


def main():
    S = dict(R_ab=load(f"RP3_N{N}_s7*_xi*_ab.json"), U_ab=load(f"U1_N{N}_s7*_xi*_ab.json"),
             R_low=load(f"RP3_N{N}_s6*_xi0.5_dr.json"), U_low=load(f"U1_N{N}_s6*_xi0.5_dr.json"))
    have = {k: len(v) for k, v in S.items()}
    if not all(have.values()): raise SystemExit(f"not enough runs yet: {have}")
    l_end = min(r["rows"][-1]["ell"] for v in S.values() for r in v)
    out = dict(spec="docs/preregistros/ESPEC_U10b.md", N=N, runs=have, window=[W0, l_end - 0.05],
               provisional=not (have["R_ab"] >= 3 and have["U_ab"] >= 3))
    out["RP3"] = model_block(S["R_ab"], S["R_low"], l_end); out["U1"] = model_block(S["U_ab"], S["U_low"], l_end)
    out["R_xi_bracket_end"] = dict(low=out["RP3"]["xi_low_end"] / out["U1"]["xi_low_end"], high=out["RP3"]["xi_ab_end"] / out["U1"]["xi_ab_end"])
    grid = np.linspace(l_end - 1.3, l_end - 0.05, 10)
    for tag, a, b in (("ab", "R_ab", "U_ab"), ("low", "R_low", "U_low")):
        rg = A2.window_ratio(S[a], S[b], grid, "eps_lat"); rx = A2.window_ratio(S[a], S[b], grid, "xi")
        try:
            q = q_fit(S[a]); qq = dict(q0_at_ell5=q["q0_at_ell5"], q0_std=q["q0_std"], q_end=q["q_at_ell_end"])
        except ValueError:
            qq = None
        out[f"emission_{tag}"] = dict(R_Gamma=rg["value"], R_Gamma_std=rg["boot_std"], R_xi=rx["value"], R_Gamma_over_R_xi=rg["value"] / rx["value"], q_RP3=qq)
    (HERE / "results" / "u4" / "U10B_RESULT.json").write_text(json.dumps(out, indent=1, default=float)) if "--runs" not in sys.argv else None
    print(f"U10b ({'PROVISIONAL' if out['provisional'] else 'final'}; runs {have}; window {W0}-{l_end - 0.05:.2f})")
    for m in ("RP3", "U1"):
        b = out[m]
        print(f"  {m}: xi from above {b['xi_ab_end']:.3f} (slope {b['slope_ab']:+.3f}, {b['frac_down_ab']:.0%} steps down), from below {b['xi_low_end']:.3f} "
              f"(slope {b['slope_low']:+.3f}); Delta_end = {b['Delta_end']:.2f}; G_A {'pass' if b['gate_GA'] else 'fail'} -> {b['verdict']}")
        print(f"      lines cross at ell = {b['lines_cross_ell']}, xi = {b['lines_cross_xi']}; per xi_i: {json.dumps({k: round(v['xi_end'], 3) for k, v in b['per_xi_init'].items()})}")
    print(f"  R_xi bracket at ell_end: [{out['R_xi_bracket_end']['low']:.2f}, {out['R_xi_bracket_end']['high']:.2f}]")
    for tag in ("ab", "low"):
        e = out[f"emission_{tag}"]
        print(f"  emission ({tag}): R_Gamma = {e['R_Gamma']:.3f}, R_xi = {e['R_xi']:.3f}, R_Gamma/R_xi = {e['R_Gamma_over_R_xi']:.3f}, q_RP3 = {None if not e['q_RP3'] else round(e['q_RP3']['q0_at_ell5'], 3)}")
    return out


if __name__ == "__main__":
    main()
