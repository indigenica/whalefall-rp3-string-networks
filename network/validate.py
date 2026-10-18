"""U2 validations V1-V4 (ESPEC_U2.md §2).  Output: results/validation.json"""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt

OUT = Path(__file__).resolve().parent / "results" / "u2"; OUT.mkdir(parents=True, exist_ok=True)


def total_energy(model, N, f, p, a=1.0):
    E = np.zeros(3); Lt.lib().energy(model, N, f, p, a, E); return E.sum()


def v1():
    res = {}
    for model in (1, 0):
        N = 48; f, p = Lt.initial(model, N, 2.0, 11); sc = np.empty_like(f)
        E0 = total_energy(model, N, f, p); Es = [E0]
        for k in range(8):
            Lt.lib().evolve(model, N, f, p, sc, 0.0, 0.25, -1.0, 50); Es.append(total_energy(model, N, f, p))
        drift = max(abs(e / E0 - 1) for e in Es)
        secular = max(abs(e / Es[1] - 1) for e in Es[1:])
        res["RP3" if model else "U1"] = dict(E0=E0, max_rel_drift=drift, secular_drift_after_first_block=secular,
                                             passed=bool(drift < 1e-3), note="preregistered gate: max_rel_drift < 1e-3")
    return res


def v2():
    """single plane wave in each sector around the vacuum, Minkowski; frequency vs lattice+leapfrog dispersion"""
    N = 32; dt = 0.25; nk = 2; kc = 2 * math.pi * nk / N; klat2 = 4 * math.sin(kc / 2) ** 2
    x = np.arange(N); wave = np.cos(kc * x)[None, None, :] * np.ones((N, N, 1))
    eps = 1e-4
    sectors = {"RP3": [("Goldstone_01", 1, {4: 1.0}, 0.0), ("quintet_12", 1, {7: 1.0}, 0.24),
                       ("singlet_normal", 1, {0: math.sqrt(3) / 2, 1: -0.5 / math.sqrt(3), 2: -0.5 / math.sqrt(3), 3: -0.5 / math.sqrt(3)}, 0.24),
                       ("heavy_trace", 1, {0: 0.5, 1: 0.5, 2: 0.5, 3: 0.5}, 1.84)],
               "U1": [("Goldstone", 0, {1: 1.0}, 0.0), ("radial", 0, {0: 1.0}, 0.24)]}
    out = {}
    for mname, lst in sectors.items():
        for name, model, shape, m2 in lst:
            C = Lt.NCOMP[model]
            f = np.zeros((C, N, N, N), np.float32)
            if model == 1: f[0] = 1.0
            else: f[0] = math.sqrt(0.5)
            for c, w in shape.items(): f[c] += eps * w * wave
            f = np.ascontiguousarray(f.reshape(C, -1)); p = np.zeros_like(f); sc = np.empty_like(f)
            vec = np.zeros((C, N ** 3)); [vec.__setitem__(c, w * wave.ravel()) for c, w in shape.items()]
            vec /= np.sqrt((vec ** 2).sum())
            base = np.zeros((C, N ** 3)); base[0] = 1.0 if model == 1 else math.sqrt(0.5)
            ts, ys = [], []
            w_exp = math.sqrt(klat2 + m2); nsteps = int(12 * 2 * math.pi / w_exp / dt)
            for k in range(nsteps):
                ts.append(k * dt); ys.append(float(((f - base) * vec).sum()))
                Lt.lib().evolve(model, N, f, p, sc, 0.0, dt, -1.0, 1)
            ts, ys = np.array(ts), np.array(ys)
            w_disc = 2 / dt * math.asin(min(1.0, dt * w_exp / 2))        # leapfrog dispersion
            popt, _ = curve_fit(lambda t, A, w, ph: A * np.cos(w * t + ph), ts, ys, p0=[ys[0], w_disc, 0.0])
            rel = abs(abs(popt[1]) / w_disc - 1)
            out[f"{mname}:{name}"] = dict(omega_fit=abs(popt[1]), omega_expected=w_disc, m2=m2, rel_err=rel, passed=bool(rel < 1e-3))
    return out


def v3():
    """string / antistring pair at x = N/2 -+ N/8 (periodic-safe: far-field angle jump across the box edge < pi/2)"""
    N = 32; y, x = np.meshgrid(np.arange(N) - N / 2 + 0.5, np.arange(N) - N / 2 + 0.5, indexing="ij")
    th = np.arctan2(y, x + N / 8) - np.arctan2(y, x - N / 8)
    res = {}
    n = np.zeros((4, N, N, N), np.float32); n[0] = np.cos(th / 2)[None]; n[1] = np.sin(th / 2)[None]
    c, _, _ = Lt.strings(1, N, np.ascontiguousarray(Lt.q_from_n(n).reshape(10, -1)))
    res["RP3_pair"] = dict(count=int(c), expected=2 * N, passed=bool(c == 2 * N))
    phi = np.zeros((2, N, N, N), np.float32); phi[0] = math.sqrt(0.5) * np.cos(th)[None]; phi[1] = math.sqrt(0.5) * np.sin(th)[None]
    c, _, _ = Lt.strings(0, N, np.ascontiguousarray(phi.reshape(2, -1)))
    res["U1_pair"] = dict(count=int(c), expected=2 * N, passed=bool(c == 2 * N))
    n2 = np.zeros((4, N, N, N), np.float32); n2[0] = 1
    c, _, _ = Lt.strings(1, N, np.ascontiguousarray(Lt.q_from_n(n2).reshape(10, -1)))
    res["RP3_uniform"] = dict(count=int(c), expected=0, passed=bool(c == 0))
    # n rotating by 2 pi around each defect (a U(1)-type winding of n) is topologically trivial in RP3 -> 0
    n3 = np.zeros((4, N, N, N), np.float32); n3[0] = np.cos(th)[None]; n3[1] = np.sin(th)[None]
    c, seed, _ = Lt.strings(1, N, np.ascontiguousarray(Lt.q_from_n(n3).reshape(10, -1)))
    # the plaquette that contains the singular point has 90-degree steps (ambiguous at lattice scale, see ESPEC_U2 amendment 2);
    # the test is that no Z2 flux is detected away from the two cores (distance > 2 cells)
    r1 = np.hypot(x + N / 8, y); r2 = np.hypot(x - N / 8, y); far = (np.minimum(r1, r2) > 2.0)[None].repeat(N, 0).ravel()
    cfar = int(seed[far].sum())
    res["RP3_full_winding_is_trivial"] = dict(count_raw=int(c), marked_sites_far_from_cores=cfar, expected_far=0, passed=bool(cfar == 0))
    return res


def v4():
    res = {}
    for model in (1, 0):
        N = 64; dt = 0.25; tau0, tau1 = 8.0, 64.0; s = 1.0 / (math.sqrt(0.24) * tau1)
        x = np.arange(N); k = 2 * math.pi * 8 / N
        pk = (np.exp(-0.5 * ((x - N / 2) / 6.0) ** 2) * np.cos(k * x))[None, None, :] * np.ones((N, N, 1))
        C = Lt.NCOMP[model]; f = np.zeros((C, N, N, N), np.float32)
        if model == 1: f[0] = 1.0; f[4] = 1e-3 * pk
        else: f[0] = math.sqrt(0.5); f[1] = 1e-3 * pk
        f = np.ascontiguousarray(f.reshape(C, -1)); p = np.zeros_like(f); sc = np.empty_like(f)
        tau, vals = tau0, []
        for ts in np.geomspace(16, tau1, 7):
            nst = int(round((ts - tau) / dt)); Lt.lib().evolve(model, N, f, p, sc, tau, dt, s, nst); tau += nst * dt
            r = Lt.measure(model, N, f, p, tau, s)
            Etot = (r["E_kin"] + r["E_grad"]) * r["a"] ** 4
            vals.append(dict(tau=tau, a4_rho_a=r["rm2.5"]["rho_a"] * r["a"] ** 4, a4_Ekin_x2=2 * r["E_kin"] * r["a"] ** 4, a4_Etot=Etot, n_plaq=r["n_plaq"]))
        # rho_a oscillates between kinetic and gradient for a standing packet: compare the cycle-averaged comoving energy
        e = np.array([v["a4_Etot"] for v in vals]); ra = np.array([v["a4_rho_a"] for v in vals])
        spread_E = float(e.max() / e.min() - 1)
        mean_ratio = float(ra.mean() / e.mean())
        res["RP3" if model else "U1"] = dict(rows=vals, comoving_energy_spread=spread_E, rho_a_over_E_mean=mean_ratio,
                                             no_strings=bool(all(v["n_plaq"] == 0 for v in vals)),
                                             passed=bool(spread_E < 0.03 and abs(mean_ratio - 1) < 0.15 and all(v["n_plaq"] == 0 for v in vals)))
    return res


def quick_checks(N=32):
    """fast subset for tests/run_tests.py: V2 (linear spectrum) and V3 (string detector)"""
    a, b = v2(), v3()
    return dict(V2=a, V3=b, passed=all(v["passed"] for v in list(a.values()) + list(b.values())))


if __name__ == "__main__":
    t0 = time.time(); R = {}
    for name, fn in (("V1", v1), ("V2", v2), ("V3", v3), ("V4", v4)):
        R[name] = fn(); print(name, json.dumps({k: (v.get("passed") if isinstance(v, dict) else v) for k, v in R[name].items()}), flush=True)
    R["all_passed"] = all(v["passed"] for blk in ("V1", "V2", "V3", "V4") for v in R[blk].values())
    R["runtime_s"] = time.time() - t0
    (OUT / "validation.json").write_text(json.dumps(R, indent=1, default=float))
    print("ALL PASSED" if R["all_passed"] else "SOME FAILED", f"{R['runtime_s']:.0f}s")
