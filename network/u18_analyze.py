"""U18 analysis (docs/preregistros/ESPEC_U18.md): local vs integral R_Gamma with 30 (= U16) and with dense snapshots.
Output: network/results/u18/U18_RESULT.json"""
import json, math
from pathlib import Path
import numpy as np
import analyze as A2

HERE = Path(__file__).resolve().parent; RUNS = HERE / "results" / "u18" / "runs"
L_END = math.log(96.0); L0, L1 = L_END - 1.3, L_END - 0.05; GRID = np.linspace(L0, L1, 10)
CASES = {"r1": (192, 96, 1.0), "r0.6": (320, 160, 0.6)}


def load(model, N, te, r):
    p = RUNS / f"{model}_N{N}_t{te}_r{r:g}_s61_dense.json"
    return json.loads(p.read_text()) if p.exists() else None


def subset(run, dense):
    return dict(run, rows=[x for x in run["rows"] if dense or x.get("u16_snapshot")])


def eps_smooth(run, half=0.15, mask="rm2.5"):
    """eps on GRID with the derivative of a^4 rho_a from a local linear regression in ln tau over |ell - ell_g| <= half"""
    R = run["rows"]; ell = np.array([x["ell"] for x in R]); lt = np.log([x["tau"] for x in R])
    y = np.array([x["a"] ** 4 * x[mask]["rho_a"] for x in R]); fr = np.array([x[mask]["frac"] for x in R])
    out = []
    for lg in GRID:
        sel = (np.abs(ell - lg) <= half) & np.isfinite(y) & (fr >= 0.5)
        if sel.sum() < 4: out.append(np.nan); continue
        dy = np.polyfit(lt[sel], y[sel], 1)[0]
        H = 1.0 / (run["s"] * math.exp(np.interp(lg, ell, lt)) ** 2); a = run["s"] * math.exp(np.interp(lg, ell, lt))
        Gam = H * dy / a ** 4; norm = 8 * math.pi * H ** 3 * A2.F2 * lg
        out.append(Gam / (norm * A2.KAPPA * lg))
    return np.array(out)


def oscillation(run, mask="rm2.5"):
    R = run["rows"]; ell = np.array([x["ell"] for x in R]); y = np.array([x["a"] ** 4 * x[mask]["rho_a"] for x in R])
    sel = (ell >= L0) & (ell <= L1) & np.isfinite(y); e, v = ell[sel], y[sel] / np.mean(y[sel])
    res = v - np.polyval(np.polyfit(e, v, 3), e)
    noise = float(np.sqrt(np.mean((res[1:-1] - 0.5 * (res[:-2] + res[2:])) ** 2) / 1.5))
    eu = np.linspace(e[0], e[-1], 512); ru = np.interp(eu, e, res); P = np.abs(np.fft.rfft(ru - ru.mean())) ** 2
    fr = np.fft.rfftfreq(len(eu), eu[1] - eu[0]); j = 1 + int(np.argmax(P[1:]))
    return dict(rms_residual=float(np.std(res)), noise=noise, resolvable=bool(np.std(res) > 3 * noise), period_ell=float(1 / fr[j]), n_points=int(sel.sum()))


def finite(x):
    if isinstance(x, dict): return {k: finite(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [finite(v) for v in x]
    if isinstance(x, (float, np.floating)): return float(x) if math.isfinite(x) else None
    if isinstance(x, (np.bool_,)): return bool(x)
    return x


def main():
    out = dict(spec="docs/preregistros/ESPEC_U18.md", window=[L0, L1], cases={})
    u16 = json.loads((HERE / "results" / "u16" / "U16_RESULT.json").read_text())["per_seed"]["61"]
    for key, (N, te, r) in CASES.items():
        rp, u1 = load("RP3", N, te, r), load("U1", N, te, r)
        if not (rp and u1): out["cases"][key] = dict(error="missing"); continue
        c = {}
        for lab, dense in (("30", False), ("dense", True)):
            a, b = subset(rp, dense), subset(u1, dense)
            c[f"R_xi_{lab}"] = A2.window_ratio([a], [b], GRID, "xi")["value"]
            c[f"R_Gamma_local_{lab}"] = A2.window_ratio([a], [b], GRID, "eps_lat")["value"]
            c[f"R_Gamma_int_{lab}"] = A2.int_ratio([a], [b], L0, L1)["value"]
        c["R_Gamma_local_smooth"] = float(np.nanmean(eps_smooth(rp)) / np.nanmean(eps_smooth(u1)))
        c["oscillation"] = dict(RP3=oscillation(rp), U1=oscillation(u1))
        ref = u16[key]
        c["V1"] = {k: dict(u18=c[{"R_xi": "R_xi_30", "R_Gamma": "R_Gamma_local_30", "R_Gamma_int": "R_Gamma_int_30"}[k]], u16=ref[k]) for k in ("R_xi", "R_Gamma", "R_Gamma_int")}
        c["V1_pass"] = bool(all(abs(v["u18"] / v["u16"] - 1) <= 1e-6 for v in c["V1"].values()))
        out["cases"][key] = c
    if all("V1_pass" in out["cases"][k] for k in CASES):
        A, B = out["cases"]["r1"], out["cases"]["r0.6"]
        out["V1_pass"] = bool(A["V1_pass"] and B["V1_pass"])
        d = {k: B[k] - A[k] for k in ("R_Gamma_local_30", "R_Gamma_local_dense", "R_Gamma_local_smooth", "R_Gamma_int_30", "R_Gamma_int_dense", "R_xi_dense")}
        out["delta_r0.6_minus_r1"] = d
        gap30 = [abs(X["R_Gamma_local_30"] - X["R_Gamma_int_30"]) for X in (A, B)]
        gapd = [abs(X["R_Gamma_local_dense"] - X["R_Gamma_int_dense"]) for X in (A, B)]
        out["gap_local_minus_int"] = dict(n30=gap30, dense=gapd)
        if not out["V1_pass"]: v = "INVALIDO"
        elif abs(d["R_Gamma_local_dense"]) <= 0.010 and all(g2 < 0.5 * g1 for g1, g2 in zip(gap30, gapd)): v = "MUESTREO"
        elif d["R_Gamma_local_dense"] >= 0.015: v = "FISICO"
        else: v = "INDETERMINADO"
        out["verdict_U18"] = v
    (HERE / "results" / "u18" / "U18_RESULT.json").write_text(json.dumps(finite(out), indent=1, allow_nan=False))
    print(json.dumps({k: out.get(k) for k in ("delta_r0.6_minus_r1", "gap_local_minus_int", "V1_pass", "verdict_U18")}, indent=1, default=float))


if __name__ == "__main__":
    main()
