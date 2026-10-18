"""U3b analysis (ESPEC_U3.md): gates G1 (U(1) control) and G2 (attractor: IC memory < 20 % at ell = 4.5), and, if both pass,
R_Gamma from set F1 with the F1-F2 difference and the U2 systematics (box, core, dt, mask) added in quadrature.
Output: results/U3B_RESULT.json"""
import json, math, glob
from pathlib import Path
import numpy as np
import analyze as A2                       # same estimators as U2 (derived, on_grid, window_ratio, int_ratio)

HERE = Path(__file__).resolve().parent
RUNS = HERE / "results" / "u3" / "runs"
U2 = HERE / "results" / "u2" / "U2_RESULT.json"


def load(p): return [json.loads(Path(x).read_text()) for x in sorted(glob.glob(str(RUNS / p)))]


def main(u2_path=U2):
    F1 = (load("RP3_N256_s4[12]_w3_xi0.5.json"), load("U1_N256_s4[12]_w3_xi0.5.json"))
    F2 = (load("RP3_N256_s43_w3_xi2.json"), load("U1_N256_s43_w3_xi2.json"))
    avail = dict(F1=[len(F1[0]), len(F1[1])], F2=[len(F2[0]), len(F2[1])])
    out = dict(spec="ESPEC_U3.md (U3b)", available=avail)
    if not F1[1]:
        out["status"] = "INCOMPLETE"; (HERE / "results" / "u3" / "U3B_RESULT.json").write_text(json.dumps(out, indent=1)); print(out); return out
    g45 = np.array([4.5])
    xi = lambda runs: float(np.nanmean(A2.on_grid(runs, g45, "xi"))) if runs else float("nan")
    eta = lambda runs, grid: float(np.nanmean(A2.on_grid(runs, grid, "eta")))
    lfin = min(r["rows"][-1]["ell"] for r in F1[0] + F1[1]); grid = np.linspace(3.0, lfin - 0.05, 12)
    rp_missing = not F1[0]
    out["xi_4p5"] = dict(F1_RP3=xi(F1[0]), F1_U1=xi(F1[1]), F2_RP3=xi(F2[0]), F2_U1=xi(F2[1]))
    out["eta_U1_F1"] = eta(F1[1], grid)
    G1 = dict(xi_U1=out["xi_4p5"]["F1_U1"], eta_U1=out["eta_U1_F1"])
    G1["passed"] = bool(0.5 <= G1["xi_U1"] <= 2.5 and 0.3 <= G1["eta_U1"] <= 3)
    mem = {}
    for m, a, b in (("RP3", "F1_RP3", "F2_RP3"), ("U1", "F1_U1", "F2_U1")):
        x1, x2 = out["xi_4p5"][a], out["xi_4p5"][b]
        mem[m] = abs(x2 - x1) / (0.5 * (x1 + x2)) if np.isfinite(x2) else float("nan")
    G2 = dict(ic_memory=mem, passed=bool(all(np.isfinite(v) and v < 0.20 for v in mem.values())))
    if rp_missing:   # amendment 3: RP3 runs of F not produced; both gates are decided by the U(1) control alone
        G2 = dict(ic_memory=mem, passed=False, note="RP3 runs of F not produced (amendment 3); U(1) memory alone exceeds 20 %" if mem["U1"] >= 0.2 else "RP3 missing")
    out["G1"], out["G2"] = G1, G2
    if rp_missing:
        out["status"] = "GATES_FAILED" if not G1["passed"] or mem["U1"] >= 0.2 else "INCOMPLETE"
        out["note"] = "R_Gamma not computed: RP3 runs of set F stopped (ESPEC_U3 amendment 3). S6 keeps the U2 value."
        (HERE / "results" / "u3" / "U3B_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
        print(json.dumps({k: out[k] for k in ("available", "xi_4p5", "G1", "G2", "status")}, indent=1, default=float)); return out
    wr = A2.window_ratio(F1[0], F1[1], grid, "eps_lat"); out["R_Gamma_F1"] = wr
    out["R_Gamma_int_F1"] = A2.int_ratio(F1[0], F1[1], 3.0, lfin - 0.05)
    out["R_xi_F1"] = A2.window_ratio(F1[0], F1[1], grid, "xi")
    if F2[0] and F2[1]:
        out["R_Gamma_F2"] = A2.window_ratio(F2[0], F2[1], grid, "eps_lat")["value"]
    u2 = json.loads(Path(u2_path).read_text())
    syst_u2 = {k: v["delta"] for k, v in u2["systematics"].items() if not k.endswith("_int") and not k.startswith("IC_") and v.get("delta") is not None}
    d_ic = abs(out.get("R_Gamma_F2", wr["value"]) - wr["value"])
    tot = math.sqrt((wr["boot_std"] or 0) ** 2 + d_ic ** 2 + sum(v ** 2 for v in syst_u2.values()))
    d_est = abs(out["R_Gamma_int_F1"]["value"] - wr["value"])
    if d_est > tot: tot = math.sqrt(tot ** 2 + d_est ** 2)
    out["R_Gamma"] = dict(value=wr["value"], stat=wr["boot_std"], ic_F1_F2=d_ic, u2_systematics=syst_u2, total=tot,
                          integral_estimator=out["R_Gamma_int_F1"]["value"])
    out["status"] = "GATES_PASSED" if (G1["passed"] and G2["passed"]) else "GATES_FAILED"
    out["note"] = ("If GATES_PASSED, S6 is recomputed with this R_Gamma (chain/s6_measured.py --u3b); otherwise the network verdict stays "
                   "INCONCLUSO and S6 keeps the U2 value.")
    (HERE / "results" / "u3" / "U3B_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: out[k] for k in ("available", "xi_4p5", "G1", "G2", "R_Gamma", "status")}, indent=1, default=float))
    return out


if __name__ == "__main__":
    main()
