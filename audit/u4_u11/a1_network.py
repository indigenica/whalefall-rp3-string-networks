"""Audit U4-U11, part 1: network estimators (audit/u4_u11/PLAN.md).  Re-computes R_Gamma, R_xi, q, drained/undrained,
X and the slopes with a_lib.py (independent of network/*_analyze.py) and compares with the published JSON.

    python3 audit/u4_u11/a1_network.py   ->  audit/u4_u11/a1_network.json"""
import json, math
from pathlib import Path
import numpy as np
import a_lib as L

HERE = Path(__file__).resolve().parent
PUB = {k: json.loads((L.ROOT / "network" / "results" / "u4" / f"{k}_RESULT.json").read_text()) for k in ("U6", "U9", "U10A", "U10B", "U11")}
OUT = {"plan": "audit/u4_u11/PLAN.md", "checks": []}


def check(name, audit, pub, sigma=None, rel=None):
    if pub is None or audit is None or (isinstance(pub, float) and math.isnan(pub)):
        st = "SIN_REFERENCIA"
    elif sigma:
        z = abs(audit - pub) / sigma; st = "PASA" if z <= 0.5 else ("MATIZADO" if z <= 1 else "DISCREPANCIA")
    else:
        d = abs(audit - pub) / max(abs(pub), 1e-12); st = "PASA" if d <= 0.03 else ("MATIZADO" if d <= 0.10 else "DISCREPANCIA")
    OUT["checks"].append(dict(name=name, audit=audit, published=pub, sigma=sigma, status=st))
    print(f"{st:<15} {name:<44} audit={audit!s:<22.22} pub={pub!s:<22.22} " + (f"sigma={sigma:.4f}" if sigma else ""), flush=True)


def runs(model, N, seeds, xi, suf=""):
    return [L.load(f"{model}_N{N}_s{s}_xi{xi}{suf}.json") for s in seeds]


def rg(r3, u1, **kw):
    (R1, R2, per, s1, s2), g, A, B = L.R_gamma(r3, u1, **kw)
    return R2, s2


def rx(r3, u1):
    R1, R2, per, s1, s2 = L.R_xi(r3, u1); return R2, s2


Q_VARIANTS = {"a,H del instante final del par (como el publicado)": dict(amode="r2"), "a,H medios geométricos (base del auditor)": dict(amode="geo"),
              "a,H del instante inicial del par": dict(amode="r1"), "k en el borde superior del bin": dict(amode="geo", kmode="upper"),
              "ajuste ponderado por F": dict(amode="geo", weight=True)}


def q_points_v(run, amode="geo", kmode="center", weight=False):
    R = [r for r in run["rows"] if r.get("spectrum")]; m = math.sqrt(run["m2"]); out = []
    for r1, r2 in zip(R[:-1], R[1:]):
        k = np.array(r1["spectrum"]["k_edges"]); kc = 0.5 * (k[1:] + k[:-1]) if kmode == "center" else k[1:]
        F = r2["a"] ** 4 * np.array(r2["spectrum"]["dE"]) - r1["a"] ** 4 * np.array(r1["spectrum"]["dE"])
        am = {"geo": math.sqrt(r1["a"] * r2["a"]), "r1": r1["a"], "r2": r2["a"]}[amode]
        Hm = {"geo": math.sqrt(r1["H"] * r2["H"]), "r1": r1["H"], "r2": r2["H"]}[amode]
        kp = kc / am; sel = (kp >= 10 * Hm) & (kp <= 0.25 * m) & (F > 0)
        if sel.sum() < 8: continue
        w = np.sqrt(F[sel] / F[sel].max()) if weight else None
        out.append((0.5 * (r1["ell"] + r2["ell"]), -np.polyfit(np.log(kc[sel]), np.log(F[sel]), 1, w=w)[0], int(sel.sum()), int(sel.sum())))
    return out


def q_all(rs):
    res = {}
    orig = L.q_points
    for name, kw in Q_VARIANTS.items():
        L.q_points = (lambda run, kw=kw, **_: q_points_v(run, **kw))
        q = L.q_fit(rs); q.pop("points"); res[name] = q
    L.q_points = orig
    return res


def q_summary(tag, rs, pub, verdict_rule=None):
    qa = q_all(rs); ref = qa["a,H del instante final del par (como el publicado)"]; base = qa["a,H medios geométricos (base del auditor)"]
    check(f"{tag} q0(ℓ=5), mismo convenio", ref["q0_at_ell5"], pub["q0_at_ell5"], pub["q0_std"])
    check(f"{tag} q0(ℓ=5), base del auditor", base["q0_at_ell5"], pub["q0_at_ell5"], pub["q0_std"])
    if "q1_per_ell" in pub: check(f"{tag} q1, mismo convenio", ref["q1_per_ell"], pub["q1_per_ell"], pub.get("q1_std"))
    q0s = [v["q0_at_ell5"] for v in qa.values()]
    OUT.setdefault("q_variants", {})[tag] = dict(variants=qa, spread=float(max(q0s) - min(q0s)), q0_plus_2sigma={k: v["q0_at_ell5"] + 2 * v["q0_std"] for k, v in qa.items()})
    print(f"   {tag}: q0 entre variantes {min(q0s):.3f}–{max(q0s):.3f} (dispersión {max(q0s) - min(q0s):.3f}, σ publicada {pub['q0_std']:.3f})", flush=True)
    return qa


# ---------------- U11 ----------------
r3 = runs("RP3", 576, range(61, 67), "0.5", "_dr"); u1 = runs("U1", 576, range(61, 65), "0.5", "_dr")
r3u6 = runs("RP3", 576, range(61, 67), "0.5"); u1u6 = runs("U1", 576, range(61, 65), "0.5")
P = PUB["U11"]
R, s = rg(r3, u1); check("U11 R_Γ drenada", R, P["R_Gamma"]["value"], P["R_Gamma"]["boot_std"]); OUT["U11_R"] = (R, s)
for sch in ("central_smooth", "fit5"):
    check(f"U11 R_Γ drenada, derivada {sch}", rg(r3, u1, scheme=sch)[0], P["R_Gamma"]["value"], P["R_Gamma"]["boot_std"])
check("U11 R_ξ drenada", rx(r3, u1)[0], P["R_xi"]["value"], P["R_xi"]["boot_std"])
Rc, sc = rg(r3, u1, x_corr=1.0); check("U11 R_Γ^corr (U(1) − ⟨X⟩)", Rc, P["R_Gamma_corr"]["value"], P["R_Gamma_corr"]["boot_std"]); OUT["U11_Rcorr"] = (Rc, sc)
check("U11 R_Γ de U6 en la misma ventana", rg(r3u6, u1u6)[0], P["R_Gamma_U6_same_window"]["value"], P["R_Gamma_U6_same_window"]["boot_std"])
for mod, A, B in (("RP3", r3, r3u6), ("U1", u1, u1u6)):
    per = {}
    for a, b in zip(A, B):
        g = L.window_grid(576)
        ga = L.on_grid(*L.emission(a)[::2], g); gb = L.on_grid(*L.emission(b)[::2], g); ok = np.isfinite(ga) & np.isfinite(gb)
        per[f"s{a['seed']}"] = float(ga[ok].mean() / gb[ok].mean())
        check(f"U11 drenada/U6 {mod} s{a['seed']}", per[f"s{a['seed']}"], P["drained_over_undrained_by_seed"][mod]["per_seed"][f"s{a['seed']}"])
    check(f"U11 drenada/U6 {mod} media", float(np.mean(list(per.values()))), P["drained_over_undrained_by_seed"][mod]["mean"])
for lab, rs, key in (("RP3 drenada", r3, "RP3_drained"), ("U1 drenada", u1, "U1_drained"), ("RP3 U6", r3u6, "RP3_U6"), ("U1 U6", u1u6, "U1_U6")):
    check(f"U11 ρ_M/ρ_a al final, {lab}", L.rhoM_over_rhoa_end(rs)[0], P["rhoM_over_rhoa_end"][key])
for r in u1:
    check(f"U11 f_X U(1) s{r['seed']}", L.f_X(r), P["f_X_U1"][f"s{r['seed']}"])
q11 = q_summary("U11 q_RP3", r3, P["q_RP3"])
OUT["U11_verdicts"] = dict(
    R=("LIMPIA_ZONA_PERMITIDA" if R + s <= 1.69 else "LIMPIA_ZONA_EXCLUIDA" if R - s >= 2.09 else "LIMPIA_INTERMEDIA"),
    R_corr=("LIMPIA_ZONA_PERMITIDA" if Rc + sc <= 1.69 else "LIMPIA_ZONA_EXCLUIDA" if Rc - sc >= 2.09 else "LIMPIA_INTERMEDIA"),
    q={k: ("Q_MENOR_QUE_1" if v["q0_at_ell5"] + 2 * v["q0_std"] < 1 else "Q_MAYOR_QUE_1" if v["q0_at_ell5"] - 2 * v["q0_std"] > 1 else "Q_COMPATIBLE_CON_1")
       for k, v in q11.items()})
print("U11 veredictos auditados:", json.dumps(OUT["U11_verdicts"], ensure_ascii=False, indent=1), flush=True)

# ---------------- U6 ----------------
P = PUB["U6"]
check("U6 R_Γ", rg(r3u6, u1u6)[0], P["R_Gamma"]["value"], P["R_Gamma"]["boot_std"])
check("U6 R_ξ", rx(r3u6, u1u6)[0], P["R_xi"]["value"], P["R_xi"]["boot_std"])
Ri = np.mean([L.integral_emission(r) for r in r3u6]) / np.mean([L.integral_emission(r) for r in u1u6])
check("U6 R_Γ integral", float(Ri), P["R_Gamma_int"]["value"], P["R_Gamma_int"]["boot_std"])
q_summary("U6 q_RP3", r3u6, P["q_RP3"])
q_summary("U6 q_U1", u1u6, P["q_U1"])
lr = 5.663 - 0.3
for r, v in zip(u1u6, P["xi_slope"]["U1"]): check(f"U6 dξ/dℓ U(1) s{r['seed']}", L.xi_slope(r, lr - 0.7, lr), v)

# ---------------- U9 ----------------
P = PUB["U9"]; r3d = runs("RP3", 576, (61, 62), "2"); u1d = runs("U1", 576, (61, 62), "2")
check("U9 R_Γ densas", rg(r3d, u1d)[0], P["R_Gamma_dense"]["value"], P["R_Gamma_dense"]["boot_std"])
check("U9 R_ξ densas", rx(r3d, u1d)[0], P["R_xi_dense"]["value"], P["R_xi_dense"]["boot_std"])
for r, v in zip(r3d, P["slope_RP3"]["dense"]): check(f"U9 dξ/dℓ RP3 densa s{r['seed']}", L.xi_slope(r, *P["slope_window"]), v)
check("U9 dξ/dℓ RP3 poco densa (media)", float(np.mean([L.xi_slope(r, *P["slope_window"]) for r in r3u6])), P["slope_RP3"]["mean_sparse"])
q_summary("U9 q_RP3 densas", r3d, P["q_RP3_dense"])

# ---------------- U10a ----------------
P = PUB["U10A"]
dd3 = runs("RP3", 512, (61, 62), "2", "_dr"); ds3 = runs("RP3", 512, (61, 62), "0.5", "_dr")
dd1 = runs("U1", 512, (61, 62), "2", "_dr"); ds1 = runs("U1", 512, (61, 62), "0.5", "_dr")
Rdd, sdd = rg(dd3, dd1); Rds, sds = rg(ds3, ds1)
check("U10a R_Γ drenadas densas", Rdd, P["R_Gamma"]["drained_dense"]["value"], P["R_Gamma"]["drained_dense"]["boot_std"])
check("U10a R_Γ drenadas poco densas", Rds, P["R_Gamma"]["drained_sparse"]["value"], P["R_Gamma"]["drained_sparse"]["boot_std"])
check("U10a R_ξ drenadas densas", rx(dd3, dd1)[0], P["R_xi"]["drained_dense"]["value"], P["R_xi"]["drained_dense"]["boot_std"])
check("U10a R_ξ drenadas poco densas", rx(ds3, ds1)[0], P["R_xi"]["drained_sparse"]["value"], P["R_xi"]["drained_sparse"]["boot_std"])
G = abs(Rdd - Rds) / ((Rdd + Rds) / 2); check("U10a brecha G", G, P["gap_G"])
for lab, rs, key in (("RP3 densas", dd3, ("RP3", "dense")), ("RP3 poco densas", ds3, ("RP3", "sparse")), ("U1 densas", dd1, ("U1", "dense")), ("U1 poco densas", ds1, ("U1", "sparse"))):
    check(f"U10a ρ_M/ρ_a al final, {lab}", L.rhoM_over_rhoa_end(rs)[0], P["massive"][key[0]][key[1]])
qdd = q_summary("U10a q_RP3 densas", dd3, P["q_RP3"]["RP3_dd"]); qds = q_summary("U10a q_RP3 poco densas", ds3, P["q_RP3"]["RP3_ds"])

# ---------------- U10b ----------------
P = PUB["U10B"]
ab3 = [L.load(f) for f in ("RP3_N576_s71_xi3_ab.json", "RP3_N576_s72_xi3_ab.json", "RP3_N576_s73_xi6_ab.json")]
ab1 = [L.load(f) for f in ("U1_N576_s71_xi3_ab.json", "U1_N576_s72_xi3_ab.json", "U1_N576_s73_xi6_ab.json")]
lo, hi = P["window"]
for mod, ab, low in (("RP3", ab3, r3), ("U1", ab1, u1)):
    sab = float(np.mean([L.xi_slope(r, lo, hi) for r in ab])); slow = float(np.mean([L.xi_slope(r, lo, hi) for r in low]))
    check(f"U10b dξ/dℓ {mod} lanzadas tarde", sab, P[mod]["slope_ab"]); check(f"U10b dξ/dℓ {mod} U11", slow, P[mod]["slope_low"])
    xe = float(np.mean([np.interp(hi, L.series(r)["ell"], L.series(r)["xi"]) for r in ab])); check(f"U10b ξ final {mod} lanzadas tarde", xe, P[mod]["xi_ab_end"])
Rab, sab_ = rg(ab3, ab1); check("U10b R_Γ lanzadas tarde", Rab, P["emission_ab"]["R_Gamma"], P["emission_ab"]["R_Gamma_std"])
q_summary("U10b q_RP3 lanzadas tarde", ab3, P["emission_ab"]["q_RP3"])

OUT["summary"] = {st: sum(c["status"] == st for c in OUT["checks"]) for st in ("PASA", "MATIZADO", "DISCREPANCIA", "SIN_REFERENCIA")}
print(OUT["summary"])
(HERE / "a1_network.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float))
