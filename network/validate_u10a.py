"""U10a, V5 (docs/preregistros/ESPEC_U10a.md): closure of the massive -> Goldstone exchange X in a string-free Minkowski box.
Vacuum + small Gaussian perturbations in every component (amplitude 0.15, lc = 4 cells), N = 64, dtau = 0.1, 600 steps,
samples every 2 steps.  Criterion per model: max_t |dE_G(t) - int_0^t X| <= 0.15 max_t |dE_G(t)|.
Also checks the drain (damp_massive): with factor 0 the massive velocity vanishes and the Goldstone velocity is unchanged.
    python -B validate_u10a.py           -> results/u4/U10A_V5.json"""
import json, math, sys
from pathlib import Path
import numpy as np
import lattice as Lt

HERE = Path(__file__).resolve().parent
N, DT, STEPS, EVERY, AMP, LC = 64, 0.1, 600, 2, 0.15, 4.0


def start(model, seed=7):
    rng = np.random.default_rng(seed)
    if model == 1:
        n0 = np.zeros((4, N, N, N), np.float32); n0[0] = 1.0
        g = Lt.gaussian_field(N, 10, LC, rng); f = Lt.q_from_n(n0) + AMP * g / g.std()
    else:
        f = np.zeros((2, N, N, N), np.float32); f[0] = math.sqrt(0.5)
        g = Lt.gaussian_field(N, 2, LC, rng); f = f + AMP * g / g.std()
    f = np.ascontiguousarray(f.reshape(Lt.NCOMP[model], -1).astype(np.float32))
    return f, np.zeros_like(f)


def v5(model):
    f, p = start(model); sc = np.empty_like(f); L = Lt.lib()
    cnt0 = Lt.strings(model, N, f)[0]
    t, E, X = [], [], []
    for k in range(STEPS // EVERY + 1):
        if k: L.evolve(model, N, f, p, sc, 0.0, DT, -1.0, EVERY)
        smp = Lt.exchange_sample(model, N, f, p, 1.0, -1.0, None)
        t.append(k * EVERY * DT); E.append(smp["e_kin"] + smp["e_grad"]); X.append(smp["X"])
    t, E, X = map(np.array, (t, E, X))
    dE = E - E[0]; IX = np.concatenate([[0.0], np.cumsum(0.5 * (X[1:] + X[:-1]) * np.diff(t))])
    err = float(np.max(np.abs(dE - IX)) / max(np.max(np.abs(dE)), 1e-30))
    cnt1 = Lt.strings(model, N, f)[0]
    return dict(model="RP3" if model else "U1", strings_start=int(cnt0), strings_end=int(cnt1), max_abs_dE=float(np.max(np.abs(dE))),
                E_G0=float(E[0]), rel_err=err, passed=bool(err <= 0.15 and cnt0 == 0 and cnt1 == 0),
                t=t[::10].tolist(), dE=dE[::10].tolist(), intX=IX[::10].tolist())


def drain_check(model):
    f, p = start(model); sc = np.empty_like(f); L = Lt.lib()
    L.evolve(model, N, f, p, sc, 0.0, DT, -1.0, 20)               # give it some velocity
    n = Lt.directions(model, N, f)
    vt0 = np.empty((10 if model else 1, N ** 3), np.float32); L.tangent_velocity(model, N, f, p, n, 1.0, vt0)
    q = p.copy(); L.damp_massive(model, N, f, q, 0.0)
    vt1 = np.empty_like(vt0); L.tangent_velocity(model, N, f, q, n, 1.0, vt1)
    vm = np.empty_like(vt0); L.normal_velocity(model, N, f, q, n, 1.0, vm)
    vm0 = np.empty_like(vt0); L.normal_velocity(model, N, f, p, n, 1.0, vm0)
    return dict(model="RP3" if model else "U1", massive_rms_before=float(np.sqrt((vm0 ** 2).mean())), massive_rms_after=float(np.sqrt((vm ** 2).mean())),
                tangent_change=float(np.sqrt(((vt1 - vt0) ** 2).mean()) / max(np.sqrt((vt0 ** 2).mean()), 1e-30)))


if __name__ == "__main__":
    out = dict(spec="docs/preregistros/ESPEC_U10a.md (V5)", V5=[v5(1), v5(0)], drain=[drain_check(1), drain_check(0)])
    for r in out["V5"]: print(f"V5 {r['model']}: strings {r['strings_start']}->{r['strings_end']}, max|dE_G| = {r['max_abs_dE']:.3e} "
                              f"(E_G0 = {r['E_G0']:.3e}), rel. error of int X = {r['rel_err']:.3f} -> {'PASS' if r['passed'] else 'FAIL'}")
    for r in out["drain"]: print(f"drain {r['model']}: massive rms {r['massive_rms_before']:.3e} -> {r['massive_rms_after']:.3e}, tangent change {r['tangent_change']:.2e}")
    (HERE / "results" / "u4" / "U10A_V5.json").write_text(json.dumps(out, indent=1))
