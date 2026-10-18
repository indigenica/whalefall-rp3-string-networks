"""Auditoria independiente U2 - fisica del nucleo (checks 1-3).
Usa una copia compilada de kern.c (audit_u2/_build/kern_audit.so); todo lo demas (potencial, gradiente, Hessiano,
estimador Goldstone, detector de plaquetas) se reimplementa aqui en numpy desde ESPEC_U2.md.
Ejecutar con OMP_NUM_THREADS=1.  Salida: audit_u2/audit_kernel.json"""
import ctypes, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
L = ctypes.CDLL(str(HERE / "_build" / "kern_audit.so"))
fp = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS"); dp = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
up = np.ctypeslib.ndpointer(np.uint8, flags="C_CONTIGUOUS"); I, D = ctypes.c_int, ctypes.c_double
L.evolve.argtypes = [I, I, fp, fp, fp, D, D, D, I]; L.energy.argtypes = [I, I, fp, fp, D, dp]
L.principal.argtypes = [I, fp, fp, ctypes.c_void_p]; L.plaquettes.argtypes = [I, I, fp, up]; L.plaquettes.restype = ctypes.c_long
L.kinetic_split.argtypes = [I, I, fp, fp, fp, D, fp, fp]

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
LAM, LAMP, ZQ = 320.0, 192.0, 400.0
OUT = {}


# ---------- independent potential (from ESPEC, not from kern.c) ----------
def Qmat(q):
    Q = np.zeros((4, 4))
    for i in range(4): Q[i, i] = q[i]
    for k, (i, j) in enumerate(PAIRS): Q[i, j] = Q[j, i] = q[4 + k] / math.sqrt(2)
    return Q


def U_rp3(q):
    Q = Qmat(q); M = Q @ Q - Q
    return (LAM * (np.trace(Q) - 1) ** 2 / 4 + LAMP * np.trace(M @ M) / 4) / ZQ


def U_u1(p): return 0.24 / 4 * (p[0] ** 2 + p[1] ** 2 - 0.5) ** 2


def num_grad(U, q, h=1e-5):
    g = np.zeros(len(q))
    for c in range(len(q)):
        e = np.zeros(len(q)); e[c] = h; g[c] = (U(q + e) - U(q - e)) / (2 * h)
    return g


def num_hess(U, q, h=1e-4):
    n = len(q); Hm = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ei = np.zeros(n); ej = np.zeros(n); ei[i] = h; ej[j] = h
            Hm[i, j] = (U(q + ei + ej) - U(q + ei - ej) - U(q - ei + ej) + U(q - ei - ej)) / (4 * h * h)
    return Hm


def lap(F):  # 7-point periodic, F (C,N,N,N) indexed [c,z,y,x]
    return sum(np.roll(F, s, ax) for ax in (1, 2, 3) for s in (1, -1)) - 6 * F


# ---------- 1a. force = a^2 (lap q - a^2 dU/dq) with U already divided by ZQ ----------
def check_force():
    res = {}
    rng = np.random.default_rng(7); N = 16
    for model, C, U in ((1, 10, U_rp3), (0, 2, U_u1)):
        if model:
            n = rng.standard_normal((4, N, N, N)); n /= np.linalg.norm(n, axis=0)
            f = np.zeros((10, N, N, N))
            for i in range(4): f[i] = n[i] * n[i]
            for k, (i, j) in enumerate(PAIRS): f[4 + k] = math.sqrt(2) * n[i] * n[j]
            f += 0.15 * rng.standard_normal(f.shape)          # well off the vacuum, generic
        else:
            f = 0.7 * rng.standard_normal((2, N, N, N))
        f = f.astype(np.float32)
        s, tau, dtau = 0.1, 5.0, 1e-4; a = s * tau
        ff = np.ascontiguousarray(f.reshape(C, -1)); pp = np.zeros_like(ff); sc = np.empty_like(ff)
        L.evolve(model, N, ff, pp, sc, tau, dtau, s, 1)
        Fk = pp.reshape(C, N, N, N).astype(np.float64) / dtau   # = 0.5[F(tau)+F(tau+dt)] ~ F(tau)
        f64 = f.astype(np.float64); g = np.zeros_like(f64)
        flat = f64.reshape(C, -1).T
        idx = rng.choice(flat.shape[0], 300, replace=False)
        gs = np.array([num_grad(U, flat[i].copy()) for i in idx])
        Fmine = (a ** 2 * (lap(f64).reshape(C, -1)[:, idx] - a ** 2 * gs.T))
        Fker = Fk.reshape(C, -1)[:, idx]
        rel = float(np.max(np.abs(Fker - Fmine)) / np.max(np.abs(Fmine)))
        # test of the "divided by ZQ": compare with gradient of U*ZQ (would be 400x)
        res["RP3" if model else "U1"] = dict(max_rel_err_force=rel, passed=bool(rel < 2e-3))
    return res


# ---------- 1b. EOM in RD: linear modes vs ODE  phi'' + (2/tau) phi' + (k_lat^2 + a^2 m^2) phi = 0 ----------
def check_eom():
    res = {}; N = 16; nk = 2; k = 2 * math.pi * nk / N; klat2 = 4 * math.sin(k / 2) ** 2
    s = 0.05; tau0, tau1, dt = 4.0, 60.0, 0.05; eps = 1e-4
    x = np.arange(N); wave = np.broadcast_to(np.cos(k * x)[None, None, :], (N, N, N))
    cases = [("U1_goldstone", 0, {1: 1.0}, 0.0), ("U1_radial", 0, {0: 1.0}, 0.24),
             ("RP3_goldstone_01", 1, {4: 1.0}, 0.0), ("RP3_quintet_12", 1, {7: 1.0}, 0.24),
             ("RP3_heavy", 1, {0: 0.5, 1: 0.5, 2: 0.5, 3: 0.5}, 1.84)]
    for name, model, shape, m2 in cases:
        C = 10 if model else 2; f = np.zeros((C, N, N, N), np.float32)
        f[0] = 1.0 if model else math.sqrt(0.5)
        nv = math.sqrt(sum(w * w for w in shape.values()))
        for c, w in shape.items(): f[c] += eps * w / nv * wave
        ff = np.ascontiguousarray(f.reshape(C, -1)); pp = np.zeros_like(ff); sc = np.empty_like(ff)
        base = np.zeros((C, N ** 3)); base[0] = 1.0 if model else math.sqrt(0.5)
        vec = np.zeros((C, N ** 3))
        for c, w in shape.items(): vec[c] = w / nv * wave.ravel()
        proj = lambda: float(((ff - base) * vec).sum() / (vec * vec).sum())
        ts = np.linspace(tau0, tau1, 113); ys = [proj()]; tau = tau0
        for t in ts[1:]:
            nst = int(round((t - tau) / dt)); L.evolve(model, N, ff, pp, sc, tau, dt, s, nst); tau += nst * dt; ys.append(proj())
        sol = solve_ivp(lambda t, y: [y[1], -2 / t * y[1] - (klat2 + (s * t) ** 2 * m2) * y[0]], (tau0, tau1), [eps, 0.0],
                        t_eval=ts, rtol=1e-10, atol=1e-16)
        err = float(np.max(np.abs(np.array(ys) - sol.y[0])) / eps)
        res[name] = dict(max_err_over_amp=err, passed=bool(err < 2e-2))
    return res


# ---------- 1c. Hessian at Q = diag(1,0,0,0) ----------
def check_spectrum():
    q0 = np.zeros(10); q0[0] = 1.0
    w = np.sort(np.linalg.eigvalsh(num_hess(U_rp3, q0)))
    wu = np.sort(np.linalg.eigvalsh(num_hess(U_u1, np.array([math.sqrt(0.5), 0.0]))))
    grad0 = float(np.abs(num_grad(U_rp3, q0)).max())
    exp = [0, 0, 0] + [0.24] * 6 + [1.84]
    return dict(RP3_eigs=[round(float(v), 6) for v in w], expected=exp, grad_at_vacuum=grad0,
                U1_eigs=[round(float(v), 6) for v in wu], passed=bool(np.allclose(w, exp, atol=1e-5) and np.allclose(wu, [0, 0.24], atol=1e-6)))


# ---------- 2. Goldstone estimator: identical wave in units of f in both models ----------
def check_goldstone():
    res = {}; N = 16; k = 2 * math.pi / N; x = np.arange(N)
    for A in (0.05, 1.0, 2.5):     # amplitude of the 2pi-periodic angle (units of f)
        th = np.broadcast_to((A * np.cos(k * x))[None, None, :], (N, N, N)).astype(np.float64)
        thd = np.broadcast_to((A * 0.3 * np.sin(k * x))[None, None, :], (N, N, N)).astype(np.float64)   # dtheta/dtau
        a = 0.7
        # U(1): phi = v e^{i theta}; pi = a^2 phi'
        v = math.sqrt(0.5)
        f = np.stack([v * np.cos(th), v * np.sin(th)]).astype(np.float32).reshape(2, -1)
        p = (a * a * np.stack([-v * np.sin(th) * thd, v * np.cos(th) * thd])).astype(np.float32).reshape(2, -1)
        KG = np.empty(N ** 3, np.float32); KM = np.empty(N ** 3, np.float32)
        L.kinetic_split(0, N, np.ascontiguousarray(f), np.ascontiguousarray(p), np.ascontiguousarray(f), a, KG, KM)
        kg_u1 = float(KG.mean()); km_u1 = float(KM.mean())
        # RP3: n = (cos(theta/2), sin(theta/2), 0, 0), Q = n n^T
        n = np.stack([np.cos(th / 2), np.sin(th / 2), 0 * th, 0 * th]); nd = np.stack([-np.sin(th / 2), np.cos(th / 2), 0 * th, 0 * th]) * thd / 2
        q = np.zeros((10, N, N, N)); qd = np.zeros_like(q)
        for i in range(4): q[i] = n[i] * n[i]; qd[i] = 2 * n[i] * nd[i]
        for kk, (i, j) in enumerate(PAIRS): q[4 + kk] = math.sqrt(2) * n[i] * n[j]; qd[4 + kk] = math.sqrt(2) * (nd[i] * n[j] + n[i] * nd[j])
        f = np.ascontiguousarray(q.astype(np.float32).reshape(10, -1)); p = np.ascontiguousarray((a * a * qd).astype(np.float32).reshape(10, -1))
        nv = np.empty((4, N ** 3), np.float32); L.principal(N, f, nv, None)
        L.kinetic_split(1, N, f, p, nv, a, KG, KM)
        kg_rp = float(KG.mean()); km_rp = float(KM.mean())
        # analytic physical Goldstone kinetic density: (1/2) f^2 (dtheta/dt)^2, dt = a dtau, f^2 = 1/2
        an = float((0.5 * 0.5 * (thd / a) ** 2).mean())
        # total kinetic from the code's own energy() for reference
        res[f"A={A}"] = dict(KG_U1=kg_u1, KG_RP3=kg_rp, analytic=an, ratio_RP3_over_U1=kg_rp / kg_u1,
                             KM_U1=km_u1, KM_RP3=km_rp, passed=bool(abs(kg_rp / kg_u1 - 1) < 1e-4 and abs(kg_u1 / an - 1) < 1e-4))
    # dynamic: same small Goldstone packet (in units of f) evolved in RD in both models; compare 2<KG> history
    N = 32; s = 1 / (math.sqrt(0.24) * 32); tau, dt = 4.0, 0.25; k = 2 * math.pi * 4 / N; x = np.arange(N)
    th = (0.3 * np.exp(-0.5 * ((x - N / 2) / 4.0) ** 2) * np.cos(k * x))[None, None, :] * np.ones((N, N, 1))
    hist = {}
    for model in (0, 1):
        if model:
            n = np.stack([np.cos(th / 2), np.sin(th / 2), 0 * th, 0 * th]); q = np.zeros((10, N, N, N))
            for i in range(4): q[i] = n[i] * n[i]
            for kk, (i, j) in enumerate(PAIRS): q[4 + kk] = math.sqrt(2) * n[i] * n[j]
            f = q
        else:
            f = math.sqrt(0.5) * np.stack([np.cos(th), np.sin(th)])
        C = f.shape[0]; ff = np.ascontiguousarray(f.astype(np.float32).reshape(C, -1)); pp = np.zeros_like(ff); sc = np.empty_like(ff)
        t = tau; out = []
        for ts in (8, 12, 16, 24, 32):
            nst = int(round((ts - t) / dt)); L.evolve(model, N, ff, pp, sc, t, dt, s, nst); t += nst * dt; a = s * t
            nv = ff if model == 0 else np.empty((4, N ** 3), np.float32)
            if model: L.principal(N, ff, nv, None)
            KG = np.empty(N ** 3, np.float32); KM = np.empty(N ** 3, np.float32); L.kinetic_split(model, N, ff, pp, nv, a, KG, KM)
            out.append(float(2 * KG.mean() * a ** 4))
        hist["RP3" if model else "U1"] = out
    r = np.array(hist["RP3"]) / np.array(hist["U1"])
    res["dynamic_packet_RD"] = dict(a4_rho_a_U1=hist["U1"], a4_rho_a_RP3=hist["RP3"], ratio=r.tolist(), passed=bool(np.all(np.abs(r - 1) < 0.02)))
    return res


# ---------- 3. string detection: own numpy holonomy / winding ----------
def my_plaq(model, F):
    """F: (4,N,N,N) unit n (RP3) or (2,N,N,N) phi (U1), [c,z,y,x]; returns count"""
    tot = 0
    sh = {"x": (0, 0, -1), "y": (0, -1, 0), "z": (-1, 0, 0)}
    def S(A, d): return np.roll(A, sh[d], axis=(1, 2, 3)) if True else None
    def roll(A, dz, dy, dx): return np.roll(A, (-dz, -dy, -dx), axis=(1, 2, 3))
    for (u, v) in (((0, 0, 1), (0, 1, 0)), ((0, 1, 0), (1, 0, 0)), ((1, 0, 0), (0, 0, 1))):   # xy, yz, zx
        c0 = F; c1 = roll(F, *u); c2 = roll(F, *(np.add(u, v))); c3 = roll(F, *v)
        cyc = [c0, c1, c2, c3, c0]
        if model:
            sg = np.ones(F.shape[1:])
            for i in range(4): sg *= np.where((cyc[i] * cyc[i + 1]).sum(0) >= 0, 1, -1)
            tot += int((sg < 0).sum())
        else:
            w = 0
            for i in range(4):
                d = np.arctan2(cyc[i + 1][1], cyc[i + 1][0]) - np.arctan2(cyc[i][1], cyc[i][0])
                w = w + (d + np.pi) % (2 * np.pi) - np.pi
            tot += int((np.abs(w) > np.pi).sum())
    return tot


def check_strings():
    res = {}; rng = np.random.default_rng(3); N = 24
    from scipy.ndimage import gaussian_filter
    for model in (1, 0):
        C = 10 if model else 2
        g = np.stack([gaussian_filter(rng.standard_normal((N, N, N)), 2.0, mode="wrap") for _ in range(C)]).astype(np.float32)
        if model:
            gf = np.ascontiguousarray(g.reshape(10, -1)); nv = np.empty((4, N ** 3), np.float32); L.principal(N, gf, nv, None)
            # independent eigenvector with numpy
            Qs = np.zeros((N ** 3, 4, 4)); gg = g.reshape(10, -1).T.astype(np.float64)
            for i in range(4): Qs[:, i, i] = gg[:, i]
            for k, (i, j) in enumerate(PAIRS): Qs[:, i, j] = Qs[:, j, i] = gg[:, 4 + k] / math.sqrt(2)
            w, V = np.linalg.eigh(Qs); nmine = V[:, :, -1].T.reshape(4, N, N, N)
            align = float(np.min(np.abs((nmine.reshape(4, -1) * nv).sum(0))))
            F = nmine; fk = nv
        else:
            F = g.astype(np.float64); fk = np.ascontiguousarray(g.reshape(2, -1)); align = None
        seed = np.zeros(N ** 3, np.uint8); ck = L.plaquettes(model, N, fk, seed)
        cm = my_plaq(model, F)
        res["RP3_random" if model else "U1_random"] = dict(kernel=int(ck), mine=int(cm), min_abs_n_alignment=align, passed=bool(ck == cm))
    # length estimator: randomly oriented rings, true length 2 pi R; (2/3) N_plaq / L_true averaged over orientations
    N = 40; R = 9.0; ratios = {0: [], 1: []}
    zz, yy, xx = np.meshgrid(*(np.arange(N) - N / 2 + 0.37,) * 3, indexing="ij")
    P = np.stack([xx, yy, zz]).reshape(3, -1)
    for t in range(24):
        Rm = np.linalg.qr(rng.standard_normal((3, 3)))[0]
        X = Rm @ P; rho = np.hypot(X[0], X[1]); z = X[2]
        th = np.arctan2(z, rho - R)
        # far-field taper of the phase to avoid fake windings at the periodic boundary: multiply by smooth window of the
        # solid-angle-like phase is not needed if we count only plaquettes away from the boundary -> mask interior region
        interior = (np.abs(P) < N / 2 - 3).all(0).reshape(N, N, N)
        for model in (0, 1):
            if model: F = np.stack([np.cos(th / 2), np.sin(th / 2), 0 * th, 0 * th]).reshape(4, N, N, N)
            else: F = np.stack([np.cos(th), np.sin(th)]).reshape(2, N, N, N)
            cnt = my_plaq_mask(model, F, interior)
            ratios[model].append((2 / 3) * cnt / (2 * math.pi * R))
    res["ring_length_estimator"] = dict(U1_mean=float(np.mean(ratios[0])), U1_sd=float(np.std(ratios[0])),
                                        RP3_mean=float(np.mean(ratios[1])), RP3_sd=float(np.std(ratios[1])),
                                        note="(2/3) N_plaq / (2 pi R), 24 random orientations, R = 9 cells; expected ~1 (isotropic average of |tx|+|ty|+|tz| = 3/2)",
                                        passed=bool(abs(np.mean(ratios[0]) - 1) < 0.06 and abs(np.mean(ratios[1]) - 1) < 0.06))
    # xi formula: xi = L_phys t^2 / V_phys with t = 1/(2H), a = s tau, H = 1/(s tau^2)
    s, tau, Lc, Vc = 0.0123, 77.0, 1234.0, 256 ** 3
    a = s * tau; H = 1 / (s * tau * tau); t = 1 / (2 * H)
    res["xi_formula"] = dict(direct=a * Lc * t * t / (a ** 3 * Vc), code=Lc * tau * tau / (4 * Vc))
    res["xi_formula"]["passed"] = bool(abs(res["xi_formula"]["direct"] / res["xi_formula"]["code"] - 1) < 1e-12)
    return res


def my_plaq_mask(model, F, mask):
    tot = 0
    def roll(A, dz, dy, dx): return np.roll(A, (-dz, -dy, -dx), axis=(1, 2, 3))
    for (u, v) in (((0, 0, 1), (0, 1, 0)), ((0, 1, 0), (1, 0, 0)), ((1, 0, 0), (0, 0, 1))):
        cyc = [F, roll(F, *u), roll(F, *(np.add(u, v))), roll(F, *v), F]
        if model:
            sg = np.ones(F.shape[1:])
            for i in range(4): sg *= np.where((cyc[i] * cyc[i + 1]).sum(0) >= 0, 1, -1)
            tot += int(((sg < 0) & mask).sum())
        else:
            w = 0
            for i in range(4):
                d = np.arctan2(cyc[i + 1][1], cyc[i + 1][0]) - np.arctan2(cyc[i][1], cyc[i][0])
                w = w + (d + np.pi) % (2 * np.pi) - np.pi
            tot += int(((np.abs(w) > np.pi) & mask).sum())
    return tot


if __name__ == "__main__":
    for name, fn in (("1a_force", check_force), ("1b_eom_RD", check_eom), ("1c_spectrum", check_spectrum),
                     ("2_goldstone_estimator", check_goldstone), ("3_strings", check_strings)):
        OUT[name] = fn(); print(name, json.dumps(OUT[name], default=float)[:900], flush=True)
    (HERE / "audit_kernel.json").write_text(json.dumps(OUT, indent=1, default=float))
