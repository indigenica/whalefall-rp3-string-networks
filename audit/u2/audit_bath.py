"""Auditoria U2 - prueba exploratoria propia: cuanto Goldstone produce un BANO MASIVO SIN CUERDAS en la ventana de analisis.
Motivo: en el conjunto A, a^4 rho_M (2<K_M> fuera de nucleos) en l = 3 es ~1.1x el incremento Delta(a^4 rho_a) de toda la ventana en
ambos modelos; si ese bano decae / se convierte en Goldstones, Gamma_lat incluye emision que no viene de las cuerdas.
Bano gaussiano (suavizado sigma = 2 celdas) en las direcciones masivas m^2 = 0.24 (U1: radial; RP3: bloque 3x3 de Q), p = 0,
amplitud ajustada para que <2K_M> promediado en las primeras oscilaciones ~ a^4 rho_M medido en A en l = 3.  N = 48, RD igual que A.
Salida: audit_u2/audit_bath.json"""
import ctypes, json, math
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter

HERE = Path(__file__).resolve().parent
L = ctypes.CDLL(str(HERE / "_build" / "kern_audit.so"))
fp = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS"); up = np.ctypeslib.ndpointer(np.uint8, flags="C_CONTIGUOUS")
I, D = ctypes.c_int, ctypes.c_double
L.evolve.argtypes = [I, I, fp, fp, fp, D, D, D, I]; L.principal.argtypes = [I, fp, fp, ctypes.c_void_p]
L.plaquettes.argtypes = [I, I, fp, up]; L.plaquettes.restype = ctypes.c_long; L.kinetic_split.argtypes = [I, I, fp, fp, fp, D, fp, fp]
AN = json.loads((HERE / "audit_analysis.json").read_text())["massive_bath"]
m = math.sqrt(0.24); s = 1 / (m * 128); dt = 0.25; N = 48
tau0 = math.sqrt(math.exp(3.0) / (m * s))


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def run(model, amp, gamp, seed=5):
    """amp: massive bath amplitude; gamp: Goldstone seed (rms angle, 2pi-periodic units) - present in both runs"""
    rng = np.random.default_rng(seed); rg = np.random.default_rng(seed + 100)
    C = 10 if model else 2; f = np.zeros((C, N, N, N))
    G = [gaussian_filter(rg.standard_normal((N, N, N)), 2.0, mode="wrap") for _ in range(3)]; G = [g / g.std() for g in G]
    if model:   # n = e0 cos(|b|) + bhat sin(|b|), b = (gamp/2) G  (theta_2pi = 2 x angle of n)
        b = np.stack(G) * gamp / 2; bn = np.sqrt((b ** 2).sum(0)) + 1e-30
        n = np.concatenate([np.cos(bn)[None], b / bn * np.sin(bn)])
        for i in range(4): f[i] = n[i] * n[i]
        for k, (i, j) in enumerate(PAIRS): f[4 + k] = math.sqrt(2) * n[i] * n[j]
        comps = [1, 2, 3, 7, 8, 9]       # massive bath in the lower 3x3 block (rotated frame ignored: n ~ e0)
    else:
        th = gamp * G[0]; f[0] = math.sqrt(0.5) * np.cos(th); f[1] = math.sqrt(0.5) * np.sin(th); comps = None
    if amp > 0:
        if model:
            for c in comps:
                g = gaussian_filter(rng.standard_normal((N, N, N)), 2.0, mode="wrap"); f[c] += amp * g / g.std()
        else:
            g = gaussian_filter(rng.standard_normal((N, N, N)), 2.0, mode="wrap"); rr = 1 + amp / math.sqrt(0.5) * g / g.std()
            f[0] *= rr; f[1] *= rr
    ff = np.ascontiguousarray(f.astype(np.float32).reshape(C, -1)); pp = np.zeros_like(ff); sc = np.empty_like(ff)
    tau = round(tau0 / dt) * dt; rows = []
    KG = np.empty(N ** 3, np.float32); KM = np.empty(N ** 3, np.float32); seedb = np.zeros(N ** 3, np.uint8)
    for k in range(int((128 - tau) / dt / 4)):
        L.evolve(model, N, ff, pp, sc, tau, dt, s, 4); tau += 4 * dt; a = s * tau
        nv = ff if model == 0 else np.empty((4, N ** 3), np.float32)
        if model: L.principal(N, ff, nv, None)
        L.kinetic_split(model, N, ff, pp, nv, a, KG, KM)
        npl = L.plaquettes(model, N, nv, seedb)
        rows.append((tau, math.log(m * s * tau * tau), float(2 * KG.mean() * a ** 4), float(2 * KM.mean() * a ** 4), int(npl)))
    return np.array(rows)


out = {}
for model, lab in ((0, "U1"), (1, "RP3")):
    target = AN[lab]["a4rhoM_3"]; gt = AN[lab]["a4rhoa_3"]; amp = 0.1; gamp = 0.2
    for it in range(3):                 # calibrate amplitudes on early (l < 3.2) means of 2<K_M> a^4 and 2<K_G> a^4
        r = run(model, amp, gamp); e = r[r[:, 1] < 3.2]
        amp *= math.sqrt(target / e[:, 3].mean()); gamp *= math.sqrt(gt / e[:, 2].mean())
    r = run(model, amp, gamp); r0 = run(model, 0.0, gamp)
    W = lambda x, lo, hi: x[(x[:, 1] > lo) & (x[:, 1] < hi)]
    dG = W(r, 4.6, 5)[:, 2].mean() - W(r, 0, 3.2)[:, 2].mean(); dG0 = W(r0, 4.6, 5)[:, 2].mean() - W(r0, 0, 3.2)[:, 2].mean()
    out[lab] = dict(amp=amp, gamp=gamp, a4rhoM_early=float(W(r, 0, 3.2)[:, 3].mean()), a4rhoM_late=float(W(r, 4.6, 5)[:, 3].mean()),
                    a4rhoa_early=float(W(r, 0, 3.2)[:, 2].mean()), a4rhoa_late=float(W(r, 4.6, 5)[:, 2].mean()),
                    a4rhoa_early_nobath=float(W(r0, 0, 3.2)[:, 2].mean()), a4rhoa_late_nobath=float(W(r0, 4.6, 5)[:, 2].mean()),
                    dGold_with_bath=float(dG), dGold_no_bath=float(dG0), dGold_bath=float(dG - dG0),
                    dGold_network_A=AN[lab]["dGold"], fraction_of_network_dGold=float((dG - dG0) / AN[lab]["dGold"]), max_plaq=int(max(r[:, 4].max(), r0[:, 4].max())))
    print(lab, json.dumps(out[lab]), flush=True)
fr = out["RP3"]["fraction_of_network_dGold"], out["U1"]["fraction_of_network_dGold"]
out["implied_bias_on_R_Gamma_if_subtracted"] = float((1 - fr[0]) / (1 - fr[1]))
(HERE / "audit_bath.json").write_text(json.dumps(out, indent=1))
print(out["implied_bias_on_R_Gamma_if_subtracted"])
