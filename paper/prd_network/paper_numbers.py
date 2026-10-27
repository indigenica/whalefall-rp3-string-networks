"""Every number quoted in main.tex, recomputed from the versioned result files and looked up in the manuscript itself.

    python3 -B paper/prd_network/paper_numbers.py          # exit 1 if anything does not match

Each row computes the values from the stored results, rounds them as the text does, fills them into the LaTeX
fragment where the paper quotes them, and checks that this exact fragment is in main.tex (whitespace-insensitive).
So the check fails if the results change, if the text changes, or if the two disagree. A second check lists every
occurrence of a number in the text (abstract to conclusions) that lies outside every matched fragment and is not in
SETUP, so a new quoted number cannot slip in unchecked, and a value checked in one place does not cover its other
appearances. SETUP lists only parameters, lattice sizes, windows, literature values and single digits; numbers that
are also results (0.6, 1.6, 30) are not in it, and their set-up or literature uses are matched as fragments in CONTEXT.
The export carries this script (`make numbers` in the paper repository)."""
import json, math, re, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
R = ROOT / "network" / "results"


def J(rel):
    return json.loads((R / rel).read_text())


U4B, U6, U9, U9X = J("u4/U4B_RESULT.json"), J("u4/U6_RESULT.json"), J("u4/U9_RESULT.json"), J("u4/U9_EXTRA.json")
U10A, U10AX, V5 = J("u4/U10A_RESULT.json"), J("u4/U10A_EXTRA.json"), J("u4/U10A_V5.json")
U10B, U11, T7 = J("u4/U10B_RESULT.json"), J("u4/U11_RESULT.json"), J("u4/T7_DRAIN_FADE.json")
U16, U18, VAL = J("u16/U16_RESULT.json"), J("u18/U18_RESULT.json"), J("u2/validation.json")
AUD = json.loads((ROOT / "audit" / "u4_u11" / "a1_network.json").read_text())


def pm(v, e, d=2):
    return f"{v:.{d}f}\\pm{e:.{d}f}"


def pct(x, d=0, sign=False):
    return f"{x * 100:{'+' if sign else ''}.{d}f}%"


def memory_384(l_ref=4.96):
    """U(1) initial-condition memory at 384^3 (cloud test, ell_ref = ell_end - 0.3), as in INFORME_U4"""
    rows = []
    for xi in ("0.25", "2"):
        run = J(f"u4/cloud_N384/U1_N384_s61_xi{xi}.json")
        e = np.array([r["ell"] for r in run["rows"] if "xi" in r]); x = np.array([r["xi"] for r in run["rows"] if "xi" in r])
        rows.append(np.interp(l_ref, e, x))
    return abs(rows[1] - rows[0]) / np.mean(rows)


qv = {k: v["variants"] for k, v in AUD["q_variants"].items()}
pub_conv = "a,H del instante final del par (como el publicado)"
two_seed = [q["q0_at_ell5"] for k in ("U9 q_RP3 densas", "U10a q_RP3 densas", "U10a q_RP3 poco densas") for q in qv[k].values()]
two_seed_shift = max(max(abs(q["q0_at_ell5"] - qv[k][pub_conv]["q0_at_ell5"]) for q in qv[k].values())
                     for k in ("U9 q_RP3 densas", "U10a q_RP3 densas"))
r1, r6 = U18["cases"]["r1"], U18["cases"]["r0.6"]
SAMP = math.ceil(max(abs(c["R_Gamma_local_30"] / c["R_Gamma_local_dense"] - 1) for c in (r1, r6)) * 100) / 100   # largest departure seen, rounded up
dres = U16["delta_resolution"]; vol = U16["volume"]; am1 = U16["amendment1_r0.5"]["delta"]
u1_576_slope = float(np.mean(U6["xi_slope"]["U1"]))


def r2(x):
    return f"{x:.2f}"


def P(x, d=0):
    return f"{x * 100:.{d}f}"


def SP(x, d=0):
    return f"{x * 100:+.{d}f}"


v5 = V5["V5"][1]["rel_err"]
mem384 = memory_384()
e_RP3 = U9X["conversion"]["emission_dense_over_sparse"]["RP3_raw"]["value"]
e_U1 = U9X["conversion"]["emission_dense_over_sparse"]["U1_raw"]["value"]
bath_RP3 = 1 - U11["drained_over_undrained_by_seed"]["RP3"]["mean"]
bath_U1 = 1 - U11["drained_over_undrained_by_seed"]["U1"]["mean"]
per_len = max(abs(1 - U10B[k]["R_Gamma_over_R_xi"]) for k in ("emission_low", "emission_ab"))
RG, RGs = U11["R_Gamma"]["value"], U11["R_Gamma"]["boot_std"]
q_shift = max(q["q0_at_ell5"] for k, q in qv["U11 q_RP3"].items() if k != "ajuste ponderado por F") - qv["U11 q_RP3"][pub_conv]["q0_at_ell5"]   # window conventions only
coll = [v["ell_min_after_collapse"] for v in U9X["from_above"].values()]
qm, qs = qv["U11 q_RP3"]["a,H medios geométricos (base del auditor)"], qv["U11 q_RP3"]["a,H del instante inicial del par"]
qF = qv["U11 q_RP3"]["ajuste ponderado por F"]
_like = [c for c in AUD["checks"] if not any(s in c["name"] for s in ("base del auditor", "derivada", "R_Γ^corr", "f_X")) and c.get("published")]
aud_like = max(abs(c["audit"] / c["published"] - 1) for c in _like)          # like-for-like audit checks
aud_RG = abs(AUD["checks"][0]["audit"] / AUD["checks"][0]["published"] - 1)
U24R = json.loads((R / "u24" / "U24_RESULT.json").read_text())
_u24g = max(U24R["delta_vs_D0"].values(), key=lambda d: abs(d["R_Gamma"]["delta"]))["R_Gamma"]          # largest |delta R_Gamma| (U24)
_u24x = U24R["delta_vs_D0"]["Dstrong"]["R_xi"]
U25R = json.loads((R / "u25" / "U25_RESULT.json").read_text()); _u25 = U25R["delta"]
_b25 = {k: 2 * abs(_u25[k]["delta"]) + 4 * _u25[k]["sigma"] for k in ("R_Gamma", "R_xi")}        # b + 2 sigma_b (ESPEC_U25 enmienda 1)
_b25s = {k: 2 * abs(_u25[k]["delta"]) + 4 * _u25[k]["se_seeds"] for k in ("R_Gamma", "R_xi")}     # same with the between-seed error
_v24 = U24R["values"]; _perlen = (_v24["Dstrong"]["R_Gamma"] / _v24["Dstrong"]["R_xi"]) / (_v24["D0"]["R_Gamma"] / _v24["D0"]["R_xi"]) - 1


U19B = json.loads((R / "u19b" / "U19B_RESULT.json").read_text())
_rxi_all = [U6["R_xi"]["value"], U9["R_xi_dense"]["value"], U10B["R_xi_bracket_end"]["high"], U10B["R_xi_bracket_end"]["low"], U11["R_xi"]["value"]]
_xi16 = "1.6" if all(1.55 <= x < 1.70 for x in _rxi_all) else "?"          # every protocol within 1.6 to one digit (1.68 included)


def _u1row(xi, sd):
    run = J(f"u4/runs/U1_N768_s{sd}_xi{xi}.json")
    e = np.array([r["ell"] for r in run["rows"] if "xi" in r]); x = np.array([r["xi"] for r in run["rows"] if "xi" in r])
    return tuple(f"{np.interp(l, e, x):.3f}" for l in (3.0, 4.0, 4.5, 5.0, 5.65, 5.95))


def sci1(x):
    """one significant digit, rounded up, as (mantissa, exponent) strings: 3.5e-3 -> ("4", "-3")"""
    e = math.floor(math.log10(x)); m = math.ceil(x / 10 ** e)
    if m == 10: m, e = 1, e + 1
    return str(m), str(e)

# (where, what, values from the results, LaTeX fragment of main.tex with {} for each value, source)
ROWS = [
    # abstract
    ("abstract", "R_xi from above / from below", (r2(U10B["R_xi_bracket_end"]["high"]), r2(U10B["R_xi_bracket_end"]["low"])),
     r"give nearly the same ratio, $\Rxi={}$ (from above) and ${}$ (from below)", "U10B"),
    ("abstract", "bath share, sparse starts", (P(bath_RP3), P(bath_U1)), r"It is ${}\%$ ($\RP$) and ${}\%$ (U(1))", "U11"),
    ("abstract", "R_Gamma stat / sampling", (r2(RG), r2(RGs), r2(SAMP * RG)), r"$\RG={}\pm{}\,({{\rm stat}})\pm{}\,({{\rm sampling}})$, with drained sparse and dense starts", "U11, U18"),
    ("abstract", "protocol gap, drained", (P(U10A["gap_G"]),), r"dense starts ${}\%$ apart at $512^3$ (two seeds each)", "U10A"),
    ("abstract", "per unit length", (P(per_len),), r"per unit string length the emission is within ${}\%$", "U10B"),
    ("abstract", "q and its convention shift", (r2(U11["q_RP3"]["q0_at_ell5"]), r2(U11["q_RP3"]["q0_std"]), r2(q_shift)), r"is $q={}\pm{}$ with our estimator, but an independent re-implementation shows that the convention used to map the physical fitting window onto lattice momenta shifts it by up to ${}$", "U11, audit"),
    # §2
    ("§2", "V1: maximum energy offset", (P(max(VAL["V1"][m]["max_rel_drift"] for m in ("RP3", "U1")), 1),), r"departs from its initial value by up to ${}\%$", "u2/validation.json V1"),
    ("§2", "V1: secular drift below 1e-3", ("3" if max(VAL["V1"][m]["secular_drift_after_first_block"] for m in ("RP3", "U1")) < 1e-3 else "?",), r"the subsequent secular drift is below $10^{{-{}}}$", "u2/validation.json V1"),
    ("§2", "halving dtau changes R_Gamma by", (f"{J('u2/U2_RESULT.json')['systematics']['timestep_D']['delta']:.3f}",), r"halving it in a $128^3$ test without the drain changes $\RG$ by ${}$", "U2_RESULT systematics"),
    ("§2", "V2: linear spectrum", (f"{max(v['rel_err'] for v in VAL['V2'].values()):.0e}".split("e-0")[0],), r"reproduced to ${}\times10^{{-5}}$", "u2/validation.json V2"),
    ("§2", "trace-mode mass^2 (V2: frequency reproduced)", (f"{VAL['V2']['RP3:heavy_trace']['m2']:.2f}" if VAL["V2"]["RP3:heavy_trace"]["passed"] else "?",),
     r"the trace mode $m^2={}$", "u2/validation.json V2"),
    ("§2", "V4: string-free Goldstone energy", (P(VAL["V4"]["RP3"]["comoving_energy_spread"], 1),), r"the comoving Goldstone energy is conserved to ${}\%$", "u2/validation.json V4"),
    ("§2", "U16: R_xi, 1 -> 0.6", (f"{dres['R_xi']['rel'] * 100:.2f}", f"{dres['R_xi']['sigma'] / dres['R_xi']['ref_r1'] * 100:.2f}"), r"changes $\Rxi$ by ${}\pm{}\%$ between", "U16"),
    ("§2", "U16: R_xi at 0.5", (P(am1["R_xi"]["rel"], 1),), r"and by ${}\%$ at $m_r a\,dx=0.5$", "U16 amendment 1"),
    ("§2", "U16: local R_Gamma at 0.6 / 0.5", (SP(dres["R_Gamma"]["rel"], 1), SP(am1["R_Gamma"]["rel"])), r"(${}\%$ at $0.6$, ${}\%$ at $0.5$)", "U16"),
    ("§2", "U18: dense sampling, local vs integral", ("1" if max(abs(c["R_Gamma_local_dense"] / c["R_Gamma_int_dense"] - 1) for c in (r1, r6)) < 0.01 else "?",), r"local and integral estimators agree to $<{}\%$", "U18"),
    ("§2", "U18: dense sampling, change 1 -> 0.6", (P(r6["R_Gamma_local_dense"] / r1["R_Gamma_local_dense"] - 1), P(r6["R_Gamma_int_dense"] / r1["R_Gamma_int_dense"] - 1, 1)), r"change by ${}\%$ and ${}\%$ between", "U18"),
    ("§2", "U18: 30 snapshots, local", (SP(r1["R_Gamma_local_30"] / r1["R_Gamma_local_dense"] - 1, 1), SP(r6["R_Gamma_local_30"] / r6["R_Gamma_local_dense"] - 1, 1)), r"departs from its densely sampled value by ${}\%$ and ${}\%$", "U18"),
    ("§2", "U18: 30 snapshots, integral", (P(min(c["R_Gamma_int_30"] / c["R_Gamma_int_dense"] - 1 for c in (r1, r6))),), r"the integral one by up to ${}\%$", "U18"),
    ("§2", "U18 -> sampling systematic (largest departure, rounded up)", (P(SAMP),), r"adopt, after the data, a $\pm{}\%$ sampling systematic", "U18"),
    ("§2", "U16: volume", (f"{(3.33 / 2) ** 3:.1f}", f"{vol['R_xi']['delta'] / vol['R_xi']['L2'] * 100:.1f}", f"{vol['R_xi']['sigma'] / vol['R_xi']['L2'] * 100:.1f}",
                           f"{vol['R_Gamma']['delta'] / vol['R_Gamma']['L2'] * 100:+.1f}", f"{vol['R_Gamma']['sigma'] / vol['R_Gamma']['L2'] * 100:.1f}"),
     r"A box with ${}$ times more Hubble volumes changes $\Rxi$ and $\RG$ by ${}\pm{}\%$ and ${}\pm{}\%$", "U16 volume (L = 3.33 vs 2)"),
    # §3.1
    ("§3.1", "U(1) memory at 256^3", (P(U4B["ic_memory_trend"][0]["ic_memory"]),), r"the U(1) control kept a ${}\%$ memory", "U4B"),
    ("§3.1", "U(1) xi and memory at 768^3", (r2(U4B["blocks"][1]["xi_U1_ref"]["mean"]), P(U4B["ic_memory_trend"][1]["ic_memory"], 1)), r"$\xi_{{\rm U(1)}}(5.65)={}\in[0.5,2.5]$ and a memory $|\xi(\xi_i{{=}}2)-\xi(\xi_i{{=}}0.25)|/\bar\xi={}\%<20\%$", "U4B"),
    ("§3.1", "memory trend", (P(U4B["ic_memory_trend"][0]["ic_memory"]), P(mem384), P(U4B["ic_memory_trend"][1]["ic_memory"], 1)), r"(${}\%$, ${}\%$, ${}\%$ at $256^3$, $384^3$, $768^3$)", "U4B; cloud_N384"),
    # §3.2
    ("§3.2", "R_xi sparse (U6)", (f"{U6['R_xi']['value']:.3f}", f"{U6['R_xi']['boot_std']:.3f}"), r"$\Rxi={}\pm{}$", "U6"),
    ("§3.2", "dxi/dell RP3 vs U(1)", (r2(U9["slope_RP3"]["mean_sparse"]), r2(u1_576_slope)), r"($d\xi/d\ell={}$ against ${}$)", "U9; U6"),
    ("§3.2", "dense starts", ("5" if max(U9["xi_RP3_ref"]["memory"], abs(U9["slope_RP3"]["mean_dense"] / U9["slope_RP3"]["mean_sparse"] - 1)) < 0.05 else "?",
                             r2(U9["R_xi_dense"]["value"]), r2(U9["R_xi_dense"]["boot_std"]), f"{min(coll):.1f}", f"{max(coll):.1f}"),
     r"give the same density and slope within ${}\%$ ($\Rxi={}\pm{}$), but all collapse below their final density at $\ell\approx{}$--${}$", "U9, U9_EXTRA"),
    ("§3.2", "late dense / sparse slopes", (r2(U10B["RP3"]["slope_ab"]), r2(U10B["U1"]["slope_ab"]), f"{U10B['RP3']['slope_low']:+.2f}", f"{U10B['U1']['slope_low']:+.2f}"),
     r"($d\xi/d\ell={}$ for $\RP$, ${}$ for U(1)), while the drained sparse starts increase (${}$, ${}$)", "U10B"),
    ("§3.2", "bracket at ell = 5.61", (r2(U10B["RP3"]["xi_low_end"]), r2(U10B["RP3"]["xi_ab_end"]), r2(U10B["U1"]["xi_low_end"]), r2(U10B["U1"]["xi_ab_end"]),
                                       r2(U10B["R_xi_bracket_end"]["high"]), r2(U10B["R_xi_bracket_end"]["low"])),
     r"bracketed to $[{},{}]$ for $\RP$ and $[{},{}]$ for U(1), and the ratio observed in the two families is nearly the same, $\Rxi={}$ from above and ${}$ from below", "U10B"),
    # §3.3
    ("§3.3", "dense/sparse emission and rho_M/rho_a", (r2(e_RP3), r2(e_U1)),
     r"($\times{}$ for $\RP$, $\times{}$ for U(1)), with a harder spectrum", "U9_EXTRA"),
    ("§3.3", "dense vs sparse q (undrained), relaxation", (pm(U9["q_RP3_dense"]["q0_at_ell5"], U9["q_RP3_dense"]["q0_std"]), pm(U6["q_RP3"]["q0_at_ell5"], U6["q_RP3"]["q0_std"]),
                                                     P(U9X["relaxation"]["RP3"]["rel_at_end"]), P(U9X["relaxation"]["U1"]["rel_at_end"])),
     r"($q_{{\RP}}(\ell=5)={}$ against ${}$ with our default convention, though with two dense seeds this difference is convention dependent, \S), and end the window with a ${}$--${}\%$ larger ratio", "U9, U6, U9_EXTRA"),
    ("§3.3", "exchange closes a string-free U(1) budget to", (P(v5),), r"closes the Goldstone energy budget of a string-free box to ${}\%$", "U10A_V5"),
    ("§3.3", "dense U(1): drain / exchange", (P(U10AX["U1_dense_two_methods"]["retained_after_drain"]), P(U10AX["U1_dense_two_methods"]["one_minus_fX_undrained"])),
     r"the drain leaves ${}\%$ of the emission, and the undrained exchange attributes ${}\%$ to the strings", "U10A_EXTRA"),
    ("§3.3", "dense RP3, informative (512 drained / 576 undrained)", (P(U10AX["emission_ratio"]["Rdd/R9"]),), r"it gives the same picture: ${}\%$ of the emission remains", "U10A_EXTRA Rdd/R9"),
    ("§3.3", "drained/undrained per seed", (r2(1 - bath_RP3), r2(1 - bath_U1), r2(U11["drained_over_undrained_by_seed"]["RP3"]["std"]), r2(U11["drained_over_undrained_by_seed"]["U1"]["std"]), P(bath_RP3), P(bath_U1)),
     r"the drain leaves ${}$ of the $\RP$ emission and ${}$ of the U(1) one (seed-to-seed scatter ${}$ and ${}$): the bath was ${}\%$ and ${}\%$", "U11"),
    ("§3.3", "undrained R_Gamma", (r2(U11["R_Gamma_U6_same_window"]["value"]), r2(U11["R_Gamma_U6_same_window"]["boot_std"])), r"the undrained $\RG={}\pm{}$ was biased low", "U11"),
    ("§3.3", "drained U(1) exchange below", (f"{math.ceil(max(abs(v) for v in U11['f_X_U1'].values()) * 200) / 2:.1f}",), r"the measured exchange is below ${}\%$", "U11 f_X_U1"),
    ("§3.3", "fading of the bath", (r2(T7["models"]["RP3"]["slope_mean"]), r2(T7["models"]["RP3"]["slope_sem"]), r2(T7["models"]["U1"]["slope_mean"]), r2(T7["models"]["U1"]["slope_sem"])),
     r"rises by ${}\pm{}$ per unit $\ell$ (U(1): ${}\pm{}$)", "T7"),
    # §3.4
    ("§3.4", "R_Gamma and R_xi drained", (r2(RG), r2(RGs), r2(SAMP * RG), r2(U11["R_xi"]["value"]), r2(U11["R_xi"]["boot_std"])),
     r"\RG&={}\pm{}\,({{\rm stat}})\pm{}\,({{\rm sampling}}),\\ \Rxi&={}\pm{},", "U11, U18"),
    ("§3.4", "R_Gamma total and exchange-corrected", (r2(RG), r2(math.hypot(RGs, SAMP * RG)), r2(U11["R_Gamma_corr"]["value"]), r2(U11["R_Gamma_corr"]["boot_std"])),
     r"that is $\RG={}\pm{}$ with the statistical and snapshot-sampling errors combined, and $\RG={}\pm{}\,({{\rm stat}})$", "U11"),
    ("§3.4", "512^3 drained, sparse / dense, gaps", (r2(U10A["R_Gamma"]["drained_sparse"]["value"]), r2(U10A["R_Gamma"]["drained_sparse"]["boot_std"]),
                                                    r2(U10A["R_Gamma"]["drained_dense"]["value"]), r2(U10A["R_Gamma"]["drained_dense"]["boot_std"]), P(U10A["gap_G"]), P(U10A["gap_U9"])),
     r"the drained ratio is ${}\pm{}$ for sparse and ${}\pm{}$ for dense starts, a protocol gap of ${}\%$ against ${}\%$", "U10A"),
    ("§3.4", "per unit length", (P(per_len), r2(U10B["emission_low"]["R_Gamma_over_R_xi"]), r2(U10B["emission_ab"]["R_Gamma_over_R_xi"])),
     r"emits within ${}\%$ of U(1) from both sides of the attractor ($\RG/\Rxi={}$ from below, ${}$ from above)", "U10B"),
    ("§3.4", "q and its trend", (r2(U11["q_RP3"]["q0_at_ell5"]), r2(U11["q_RP3"]["q0_std"]), r2((U11["q_RP3"]["q_end"] - U11["q_RP3"]["q0_at_ell5"]) / (U6["q_RP3"]["ell_end"] - 5)), r2(U11["q_RP3"]["q1_std"])),
     r"$q_{{\RP}}(\ell=5)={}\pm{}$, with a weak decrease with $\ell$ ($dq/d\ell={}\pm{}$)", "U11"),
    ("§3.4", "q by convention", (r2(qm["q0_at_ell5"]), r2(qm["q0_std"]), r2(qs["q0_at_ell5"]), r2(qs["q0_std"]), r2(qF["q0_at_ell5"]), r2(qF["q0_std"])),
     r"give ${}\pm{}$ and ${}\pm{}$, and a fit weighted by $F$ gives ${}\pm{}$", "audit q_variants"),
    ("§3.4", "two seeds", (f"{two_seed_shift:.1f}", f"{min(two_seed):.1f}", f"{max(two_seed):.1f}"), r"(by up to ${}$ for the dense starts), so earlier two-seed values should be read as $q\approx{}$--${}$", "audit (U9, U10a)"),
    # §4, §5, conclusions
    ("§4", "bath share in the discussion", (P(bath_U1),), r"A ${}\%$ contamination", "U11"),
    ("§5", "verdict threshold quoted", ("1.69",), r"$\RG+\sigma\le{}$", "ESPEC_U11 (threshold, not a result)"),
    ("§5", "audit: like-for-like max and R_Gamma", (*sci1(aud_like), *sci1(aud_RG)), r"to $\le{}\times10^{{{}}}$ relative ($\RG$ itself to ${}\times10^{{{}}}$", "audit checks (like-for-like)"),
    ("§5", "exchange-corrected R_Gamma: ours and audit's", (r2(U11["R_Gamma_corr"]["value"]), r2(AUD["U11_Rcorr"][0])), r"moves the exchange-corrected value from ${}$ to ${}$", "U11, audit"),
    ("§3.4", "exchange-corrected with the audit's definition", (r2(AUD["U11_Rcorr"][0]),), r"(${}$ with the audit's definition of the per-seed exchange fraction)", "audit"),
    ("abstract", "from above R_Gamma", (r2(U10B["emission_ab"]["R_Gamma"]), r2(U10B["emission_ab"]["R_Gamma_std"])), r"give $\RG={}\pm{}$, but per unit string length", "U10B emission_ab"),
    ("§3.4", "from above R_Gamma and q", (r2(U10B["emission_ab"]["R_Gamma"]), r2(U10B["emission_ab"]["R_Gamma_std"]), r2(U10B["emission_ab"]["q_RP3"]["q0_at_ell5"]), r2(U10B["emission_ab"]["q_RP3"]["q0_std"])),
     r"give $\RG={}\pm{}$ and $q_{{\RP}}(\ell=5)={}\pm{}$ in the same window", "U10B emission_ab"),
    ("§3.4", "U(1) bath fade", (r2(T7["models"]["U1"]["slope_mean"]), r2(T7["models"]["U1"]["slope_sem"])), r"no significant fading (${}\pm{}$ per unit $\ell$", "T7_DRAIN_FADE"),
    ("§3.4", "R_xi drained vs undrained; U24 drain-rate shift of R_xi", (f"{abs(U11['R_xi']['value'] / U6['R_xi']['value'] - 1) * 100:.1f}", f"{abs(_u24x['delta']) * 100:.1f}", f"{_u24x['sigma'] * 100:.1f}"),
     r"the drained and undrained $\Rxi$ agree to ${}\%$, but in the preregistered drain test doubling its rate lowers $\Rxi$ by ${}\pm{}\%$", "U11, U6, U24"),
    ("§3.4", "U24: largest drain effect on R_Gamma", (f"{abs(_u24g['delta']) * 100:.1f}", f"{_u24g['sigma'] * 100:.1f}"), r"changes $\RG$ by at most ${}\pm{}\%$, which bounds", "U24"),
    ("§2", "U25: time step", (f"{_u25['R_Gamma']['delta'] * 100:+.2f}".replace("+", ""), f"{_u25['R_Gamma']['sigma'] * 100:.2f}", f"{_u25['R_xi']['delta'] * 100:+.2f}".replace("+", ""), f"{_u25['R_xi']['sigma'] * 100:.2f}"),
     r"changes of ${}\pm{}\%$ in $\RG$ and ${}\pm{}\%$ in $\Rxi$", "U25"),
    ("§2", "U25: bias + 2 sigma (bootstrap and between-seed)", (f"{_b25['R_Gamma'] * 100:.1f}", f"{_b25['R_xi'] * 100:.1f}", f"{_b25s['R_Gamma'] * 100:.1f}", f"{_b25s['R_xi'] * 100:.1f}"),
     r"with bias plus $2\sigma$ of ${}\%$ and ${}\%$ (${}\%$ and ${}\%$ with the between-seed error)", "U25 (enmienda 1; ESPEC_U24 enmienda 1)"),
    ("§3.4", "U24 Dstrong: R_Gamma shift and per-length shift", (f"{U24R['delta_vs_D0']['Dstrong']['R_Gamma']['delta'] * 100:.1f}", f"{U24R['delta_vs_D0']['Dstrong']['R_Gamma']['sigma'] * 100:.1f}", P(_perlen)),
     r"the emission ratio moves by ${}\pm{}\%$ because both models respond alike; the emission per unit string length, $\RG/\Rxi$, therefore moves by about ${}\%$", "U24"),
    ("concl.", "per unit length, R_Gamma total", (P(per_len), r2(RG), r2(math.hypot(RGs, SAMP * RG))), r"within ${}\%$ as efficiently, giving $\RG={}\pm{}$", "U10B, U11"),
    ("concl.", "bath share", (P(bath_RP3), P(bath_U1)), r"(${}$--${}\%$ for sparse starts", "U11"),
    # section E of PLAN_REVISION_PRD (19 Oct 2026)
    ("§2", "trace mode m a dx at the end, and at m_r a dx = 0.6", (f"{math.sqrt(VAL['V2']['RP3:heavy_trace']['m2'] / 0.24):.2f}", f"{0.6 * math.sqrt(VAL['V2']['RP3:heavy_trace']['m2'] / 0.24):.2f}"),
     r"the heavier trace mode reaches $m\,a\,dx={}$ at the end; the resolution test below, which lowers it to ${}$, covers it empirically", "u2/validation.json V2"),
    ("§2", "q fitting range in k at ell = 5", (f"{math.exp(5) / 40:.1f}",), r"$10H<k/a<m_r/4$ (a factor ${}$ in $k$ at $\ell=5$)", "definition (m_r/4 over 10H = e^ell/40)"),
    ("§3.1", "U(1) 768^3 table: seeds 63 and 64 exist (U5d)", tuple(str(sd) for sd in (63, 64) if (R / "u4" / "runs" / f"U1_N768_s{sd}_xi0.5.json").exists()),
     r"Seeds {} and {} were added later to measure the growth rate", "u4/runs"),
    *[("§3.1", f"U(1) 768^3 table, xi_i = {xi} seed {sd}", (xi, str(sd), *_u1row(xi, sd)), r"${}$ ({}) & {} & {} & {} & {} & {} & {}\\", f"u4/runs/U1_N768_s{sd}_xi{xi}.json")
      for xi, sd in (("0.25", 61), ("2", 61), ("0.5", 61), ("0.5", 62), ("0.5", 63), ("0.5", 64))],
    ("abstract", "R_xi ~ 1.6 for every protocol", (_xi16,), r"The string-density ratio is $\Rxi\simeq{}$ for every protocol tested", "U6, U9, U10B, U11"),
    ("§4", "R_xi ~ 1.6 (discussion)", (_xi16,), r"Our same-code ratio to U(1), $\Rxi\simeq{}$, is consistent", "U6, U9, U10B, U11"),
    ("concl.", "R_xi ~ 1.6 (conclusions)", (_xi16,), r"is $\approx{}$ times denser than a U(1) network", "U6, U9, U10B, U11"),
    ("§4", "U(1) control density at ell ~ 5.6", (f"{U4B['blocks'][1]['xi_U1_ref']['mean']:.1f}",), r"the U(1) control density, $\xi\simeq{}$ at $\ell\simeq5.6$", "U4B"),
    ("§4", "Yamada comparison: our RP3 bracket", (r2(U10B["RP3"]["xi_low_end"]), r2(U10B["RP3"]["xi_ab_end"])), r"our $\RP$ bracket at $\ell=5.61$, $[{},{}]$, contains his fixed-comoving-core value", "U10B"),
    ("§5", "U11 verdict margin", (r2(RG + RGs), r2(RG + math.hypot(RGs, SAMP * RG)), r2(1.69 - (RG + math.hypot(RGs, SAMP * RG)))),
     r"$\RG+\sigma={}$ with the statistical error and ${}$ with the combined one, a margin of only ${}$", "U11, U18 (threshold of ESPEC_U11)"),
    ("§3.4", "U16 resolution of R_xi, repeated in the budget", (f"{dres['R_xi']['rel'] * 100:.2f}",), r"resolution moves $\Rxi$ by ${}\%$ and leaves $\RG$ within the sampling systematic", "U16"),
    ("abstract", "q weighted by F", (r2(qF["q0_at_ell5"]), r2(qF["q0_std"])), r"a fit weighted by the spectrum gives ${}\pm{}$; we therefore quote", "audit q_variants"),
    ("§4", "limitation (ii): the per-length agreement", (P(per_len),), r"and the ${}\%$ agreement with U(1) quoted above", "U10B"),
    ("concl.", "U19b: RP3/U(1) soft yield per string energy (perspective)", (r2(U19B["blocks"]["RP3_3.2"]["Y"] / U19B["blocks"]["U1_3.2"]["Y"]), str(U19B["blocks"]["RP3_3.2"]["n"])),
     r"within a ratio of ${}$, so the denser", "U19B_RESULT"),
    ("concl.", "U19b: networks per model", (str(min(U19B["blocks"]["RP3_3.2"]["n"], U19B["blocks"]["U1_3.2"]["n"])),), r"(four networks per model at $512^3$, not independently audited)" if min(U19B["blocks"]["RP3_3.2"]["n"], U19B["blocks"]["U1_3.2"]["n"]) == 4 else "?", "U19B_RESULT"),
]

# Set-up or literature uses of numbers that are also results (so they are not in SETUP); matched as fragments.
CONTEXT = [
    r"so that $m_r a\,dx=1$, $0.75$ and $0.6$ at the end of the run",
    r"between $m_r a\,dx=1$ and $0.6$",
    r"($+1.6\%$ at $0.6$, $+4\%$ at $0.5$)",
    r"with the $\approx30$ snapshots per run used for all production networks",
    r"about $30\%$ in large U(1) simulations",
    r"His $\zeta$ is our $\xi$ (string length from pierced plaquettes with the same factor $2/3$, times $t^2/V$)",
    r"(Claude Opus 5.5, Anthropic)",                                              # model version required by the APS AI policy
]

# Numbers in the text that are set-up, parameters, definitions or literature values, not results of this work.
SETUP = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "20", "70", "100",
         "128", "192", "256", "320", "384", "512", "576", "768",                     # lattice sizes
         "400",                                                                      # model constants lambda = 320, lambda_P = 192, Z_Q = 400 (ESPEC_U2)
         "0.24", "0.25", "0.5", "0.6", "0.75", "1.3", "0.05", "2.5", "3.5", "3.8", "4.6", "4.56", "3.4", "4.1",   # parameters, windows (3.4, 4.1: drain ends of U24)
         "5.61", "5.65", "5.66", "5.95", "5.6", "1.5", "2.0", "8", "9", "3.0", "4.0", "4.5", "5.0",   # ell values, literature (Yamada), table columns
         "1.69", "1e-3", "-3", "-5", "-4", "-1", "-0.05", "6e-05"}


def manuscript():
    t = (HERE / "main.tex").read_text(encoding="utf-8")
    t = t[t.index(r"\begin{abstract}"):t.index(r"\bibliographystyle")]
    t = re.sub(r"(?<!\\)%.*", "", t)
    t = re.sub(r"\\(cite|ref|label|includegraphics)(\[[^]]*\])?\{[^}]*\}", "", t)   # citation keys carry years
    return re.sub(r"\s+", " ", t)


def main():
    tex = manuscript()
    bad, spans = 0, []
    w = max(len(r[1]) for r in ROWS)

    def mark(f):
        hits = [m.span() for m in re.finditer(re.escape(f), tex)]
        spans.extend(hits)
        return bool(hits)
    for where, what, vals, frag, src in ROWS:
        f = re.sub(r"\s+", " ", frag.format(*vals))
        ok = mark(f)
        bad += not ok
        print(f"{'ok ' if ok else 'XX '} {where:8} {what:<{w}}  {' / '.join(vals):<34} [{src}]")
        if not ok:
            print(f"      not in main.tex: {f}")
    for frag in CONTEXT:
        if not mark(re.sub(r"\s+", " ", frag)):
            bad += 1; print(f"XX  context fragment not in main.tex: {frag}")
    setup = SETUP | {n.lstrip("0") or "0" for n in SETUP}
    loose = sorted({m.group() for m in re.finditer(r"(?<![A-Za-z\d.])\d+(?:\.\d+)?", tex)                      # not digits inside a command (\tfrac12)
                    if m.group() not in setup and not any(a <= m.start() and m.end() <= b for a, b in spans)}, key=float)
    print(f"\n{len(ROWS)} quoted fragments checked against the results and the manuscript: {bad} do not match")
    if loose:
        print("Numbers in main.tex outside every checked fragment and not in SETUP (add a row, a CONTEXT fragment or, if it is a parameter, list it):")
        print("  " + ", ".join(loose))
    sys.exit(1 if bad or loose else 0)


if __name__ == "__main__":
    main()
