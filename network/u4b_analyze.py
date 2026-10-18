"""U4b analysis (ESPEC_U4.md): gates G1' (xi_U1 at ell_ref in [0.5, 2.5]), G2' (IC memory < 20 % at ell_ref, xi_i 0.25 vs 2.0),
G3' (informative: dxi/dell of U(1) in [ell_ref - 0.7, ell_ref]); if RP3 runs exist and G1'+G2' pass, R_Gamma in
[ell_ref - 1, ell_end - 0.05] with seed bootstrap + U2 mask/dt systematics, and the RP3 spectral slope q.
Output: results/u4/U4B_RESULT.json (one block per lattice size found)."""
import glob, json, math
from pathlib import Path
import numpy as np
import analyze as A2

HERE = Path(__file__).resolve().parent; RUNS = HERE / "results" / "u4" / "runs"
U2 = HERE / "results" / "u2" / "U2_RESULT.json"


def load(p): return [json.loads(Path(x).read_text()) for x in sorted(glob.glob(str(RUNS / p)))]


def xi_at(runs, ell):
    v = [float(np.interp(ell, [r["ell"] for r in R["rows"]], [r["xi"] for r in R["rows"]])) for R in runs]
    return float(np.mean(v)) if v else float("nan"), v


def slope(runs, l0, l1):
    out = []
    for R in runs:
        l = np.array([r["ell"] for r in R["rows"]]); x = np.array([r["xi"] for r in R["rows"]]); m = (l >= l0) & (l <= l1)
        if m.sum() >= 3: out.append(float(np.polyfit(l[m], x[m], 1)[0]))
    return out


def block(N):
    U = {xi: load(f"U1_N{N}_s61_xi{xi:g}.json") for xi in (0.25, 2.0)}
    U05 = load(f"U1_N{N}_s6[12]_xi0.5.json"); R05 = load(f"RP3_N{N}_s6[12]_xi0.5.json")
    allr = U[0.25] + U[2.0] + U05 + R05
    if not allr: return None
    l_end = min(r["rows"][-1]["ell"] for r in allr); l_ref = l_end - 0.3
    b = dict(N=N, ell_end=l_end, ell_ref=l_ref, available=dict(U1_xi025=len(U[0.25]), U1_xi2=len(U[2.0]), U1_xi05=len(U05), RP3_xi05=len(R05)))
    x05, x05s = xi_at(U05, l_ref); b["xi_U1_ref"] = dict(mean=x05, seeds=x05s)
    xa, _ = xi_at(U[0.25], l_ref); xb, _ = xi_at(U[2.0], l_ref)
    b["G1p"] = dict(xi_U1=x05, passed=bool(np.isfinite(x05) and 0.5 <= x05 <= 2.5))
    mem = abs(xb - xa) / (0.5 * (xa + xb)) if np.isfinite(xa) and np.isfinite(xb) else float("nan")
    b["G2p"] = dict(xi_025=xa, xi_2=xb, ic_memory=mem, passed=bool(np.isfinite(mem) and mem < 0.20))
    b["G3p_slope_U1"] = dict(slopes=slope(U05 + U[0.25] + U[2.0], l_ref - 0.7, l_ref), GHV_Buschmann=[0.24, 0.254], Hindmarsh=0.0)
    if R05 and U05:
        grid = np.linspace(l_ref - 1.0, l_end - 0.05, 10)
        wr = A2.window_ratio(R05, U05, grid, "eps_lat"); b["R_Gamma"] = wr
        b["R_Gamma_int"] = A2.int_ratio(R05, U05, l_ref - 1.0, l_end - 0.05)
        b["R_xi"] = A2.window_ratio(R05, U05, grid, "xi")
        u2 = json.loads(U2.read_text())
        syst = {k: u2["systematics"][k]["delta"] for k in ("mask_rm1.5", "mask_rm4", "timestep_D") if k in u2["systematics"]}
        tot = math.sqrt((wr["boot_std"] or 0) ** 2 + sum(v ** 2 for v in syst.values()))
        d_est = abs(b["R_Gamma_int"]["value"] - wr["value"])
        if d_est > tot: tot = math.sqrt(tot ** 2 + d_est ** 2)
        b["R_Gamma_summary"] = dict(value=wr["value"], total=tot, u2_systematics=syst, estimator_difference=d_est)
        b["spectra_q_RP3"] = {f"s{r['seed']}": A2.spectrum_q(r) for r in R05}
    b["status"] = "GATES_PASSED" if (b["G1p"]["passed"] and b["G2p"]["passed"]) else "GATES_FAILED"
    return b


def inherit_gates(blocks):
    """ESPEC_U4: phase 2 (RP3 + U1 at N2 < N1) runs only if phase 1 passed its gates at N1; the gates are those of the
    phase-1 block (xi at a given ell does not depend on the box, U2/U3b/U4b).  A block with RP3 runs but without its own
    xi_i = 0.25 / 2.0 runs therefore inherits the status of the largest-N phase-1 block.  For information, the IC memory of
    the phase-1 runs is also evaluated at the start, reference point and end of the phase-2 R_Gamma window."""
    p1 = [b for b in blocks if b["available"]["U1_xi025"] and b["available"]["U1_xi2"]]
    if not p1: return blocks
    ref = max(p1, key=lambda b: b["N"])
    U = {xi: load(f"U1_N{ref['N']}_s61_xi{xi:g}.json") for xi in (0.25, 2.0)}
    for b in blocks:
        if b is ref or not b["available"]["RP3_xi05"] or b["available"]["U1_xi025"]: continue
        mem = {}
        for l in (b["ell_ref"] - 1.0, b["ell_ref"], b["ell_end"] - 0.05):
            xa, _ = xi_at(U[0.25], l); xb, _ = xi_at(U[2.0], l); mem[f"{l:.2f}"] = abs(xb - xa) / (0.5 * (xa + xb))
        b["G2p"] = dict(inherited_from_N=ref["N"], passed=ref["G2p"]["passed"], ic_memory_ref=ref["G2p"]["ic_memory"],
                        ic_memory_in_this_window=mem)
        b["status_own"] = b["status"]
        b["status"] = "GATES_PASSED" if (ref["status"] == "GATES_PASSED" and b["G1p"]["passed"]) else "GATES_FAILED"
        b["gates_note"] = f"phase-2 block: gates of the phase-1 block N={ref['N']} (ESPEC_U4); G1' also at this N"
    return blocks


def main():
    Ns = sorted({int(Path(p).name.split("_N")[1].split("_")[0]) for p in glob.glob(str(RUNS / "*.json"))})
    out = dict(spec="docs/preregistros/ESPEC_U4.md (U4b)", blocks=inherit_gates([b for b in (block(N) for N in Ns) if b]))
    # IC-memory trend with N (U3b at N = 256 as the first point)
    u3 = HERE / "results" / "u3" / "U3B_RESULT.json"
    trend = [dict(N=256, ell_ref=4.5, ic_memory=json.loads(u3.read_text())["G2"]["ic_memory"]["U1"])] if u3.exists() else []
    trend += [dict(N=b["N"], ell_ref=b["ell_ref"], ic_memory=b["G2p"]["ic_memory"]) for b in out["blocks"] if "ic_memory" in b["G2p"]]
    out["ic_memory_trend"] = trend
    (HERE / "results" / "u4").mkdir(parents=True, exist_ok=True)
    (HERE / "results" / "u4" / "U4B_RESULT.json").write_text(json.dumps(out, indent=1, default=float))
    for b in out["blocks"]:
        print(f"N={b['N']} ell_ref={b['ell_ref']:.2f} G1'={b['G1p']} G2'={b['G2p'].get('ic_memory', b['G2p'].get('ic_memory_in_this_window'))} {b['G2p']['passed']} "
              f"slope={b['G3p_slope_U1']['slopes']} status={b['status']}" + (f" R_Gamma={b['R_Gamma_summary']}" if "R_Gamma_summary" in b else ""))
    print("IC memory trend:", trend)
    return out


if __name__ == "__main__":
    main()
