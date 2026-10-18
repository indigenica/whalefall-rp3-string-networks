"""U3b driver (U2 driver + U3b core-growth preparation; with tau_c = 0 it is exactly the U2 driver). U2 driver: ctypes wrapper of kern.c, initial conditions, measurements (ESPEC_U2.md)."""
import ctypes, math, os, subprocess, json, time
from pathlib import Path
import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
F2 = 0.5                                   # Goldstone constant f^2 (code units), both models
M2 = {1: 0.24, 0: 0.24}                    # core mass^2: quintet (RP3) / radial (U1)
NCOMP = {1: 10, 0: 2}
_L = None


def lib():
    global _L
    if _L is None:
        so = HERE / "_build" / "kern.so"; so.parent.mkdir(exist_ok=True)
        if not so.exists() or so.stat().st_mtime < (HERE / "kern.c").stat().st_mtime:
            subprocess.check_call(["gcc", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC", str(HERE / "kern.c"), "-o", str(so), "-lm"])
        L = ctypes.CDLL(str(so))
        fp = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS"); dp = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
        up = np.ctypeslib.ndpointer(np.uint8, flags="C_CONTIGUOUS"); I, D = ctypes.c_int, ctypes.c_double
        L.evolve.argtypes = [I, I, fp, fp, fp, D, D, D, I]; L.evolve.restype = None
        L.energy.argtypes = [I, I, fp, fp, D, dp]; L.energy.restype = None
        L.principal.argtypes = [I, fp, fp, ctypes.c_void_p]; L.principal.restype = None
        L.plaquettes.argtypes = [I, I, fp, up]; L.plaquettes.restype = ctypes.c_long
        L.kinetic_split.argtypes = [I, I, fp, fp, fp, D, fp, fp]; L.kinetic_split.restype = None
        L.tangent_velocity.argtypes = [I, I, fp, fp, fp, D, fp]; L.tangent_velocity.restype = None
        L.nthreads.restype = ctypes.c_int
        L.set_core.argtypes = [D]; L.set_core.restype = None
        L.set_core_late.argtypes = [D]; L.set_core_late.restype = None                 # U19b
        L.edt_periodic.argtypes = [I, up, fp]; L.edt_periodic.restype = None
        L.set_damp.argtypes = [D, D]; L.set_damp.restype = None                         # U10a
        L.damp_massive.argtypes = [I, I, fp, fp, D]; L.damp_massive.restype = None
        L.exchange.argtypes = [I, I, fp, fp, D, ctypes.c_void_p, dp]; L.exchange.restype = None
        L.normal_velocity.argtypes = [I, I, fp, fp, fp, D, fp]; L.normal_velocity.restype = None
        _L = L
    return _L


# ---------- initial conditions ----------
def gaussian_field(N, ncomp, lc, rng):
    k = np.fft.fftfreq(N) * 2 * np.pi; kz = np.fft.rfftfreq(N) * 2 * np.pi
    K2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + kz[None, None, :] ** 2
    filt = np.exp(-0.25 * K2 * lc * lc)            # amplitude filter -> P(k) ~ exp(-k^2 lc^2 / 2)
    out = np.empty((ncomp, N, N, N), np.float32)
    for c in range(ncomp):
        w = rng.standard_normal((N, N, N)).astype(np.float32)
        out[c] = np.fft.irfftn(np.fft.rfftn(w) * filt, s=(N, N, N)).astype(np.float32)
    return out


def q_from_n(n):
    """n: (4, ...) -> Q = n n^T in the 10-comp Frobenius basis (diag 4, offdiag 01,02,03,12,13,23 times sqrt2)"""
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    out = np.empty((10,) + n.shape[1:], np.float32)
    for i in range(4): out[i] = n[i] * n[i]
    for k, (i, j) in enumerate(pairs): out[4 + k] = math.sqrt(2) * n[i] * n[j]
    return out


def initial(model, N, lc, seed):
    rng = np.random.default_rng(seed)
    if model == 1:
        g = gaussian_field(N, 10, lc, rng)                    # random symmetric matrix field
        n = np.empty((4, N, N, N), np.float32)
        lib().principal(N, np.ascontiguousarray(g.reshape(10, -1)), n.reshape(4, -1), None)
        f = q_from_n(n)
    else:
        g = gaussian_field(N, 2, lc, rng); r = np.sqrt(g[0] ** 2 + g[1] ** 2) + 1e-30
        f = (math.sqrt(0.5) * g / r).astype(np.float32)
    return np.ascontiguousarray(f.reshape(NCOMP[model], -1)), np.zeros((NCOMP[model], N ** 3), np.float32)


# ---------- measurements ----------
def directions(model, N, f):
    if model == 1:
        n = np.empty((4, N ** 3), np.float32); lib().principal(N, f, n, None); return n
    return f


def strings(model, N, f, n=None):
    n = directions(model, N, f) if n is None else n
    seed = np.zeros(N ** 3, np.uint8)
    cnt = lib().plaquettes(model, N, n, seed)
    return cnt, seed, n


def distance_to_strings(N, seed, pad, slab=64):
    """Euclidean distance to the nearest pierced plaquette, periodic: exact everywhere (C, OpenMP; U4b speed-up).  The
    scipy slab version below gives the same values wherever the distance is <= pad, which is all that the masks use."""
    s = np.ascontiguousarray(seed, dtype=np.uint8)
    if not s.any(): return np.full((N, N, N), 1e9, np.float32)
    out = np.empty(N ** 3, np.float32); lib().edt_periodic(N, s, out)
    return out.reshape(N, N, N)


def distance_to_strings_scipy(N, seed, pad, slab=64):
    """Euclidean distance to the nearest pierced plaquette, periodic.  Exact wherever the distance is <= pad (all that the
    masks use, since every mask radius is <= pad - 2); larger distances are only guaranteed to be > pad.  Computed in
    z-slabs of thickness `slab` (U4b memory fix; identical masks to the single full-volume transform of U2)."""
    s = seed.reshape(N, N, N)
    if not s.any(): return np.full((N, N, N), 1e9, np.float32)
    out = np.empty((N, N, N), np.float32)
    sxy = np.pad(s, ((0, 0), (pad, pad), (pad, pad)), mode="wrap") if pad > 0 else s
    for z0 in range(0, N, slab):
        z1 = min(z0 + slab, N)
        zi = np.arange(z0 - pad, z1 + pad) % N
        d = ndimage.distance_transform_edt(sxy[zi] == 0)
        out[z0:z1] = d[pad:pad + (z1 - z0), pad:pad + N, pad:pad + N]
    return out


def measure(model, N, f, p, tau, s, rmask=(1.5, 2.5, 4.0), spectra=False, tau_c=0.0, massive_spectrum=False, return_dist=False, tau_f=0.0):
    a = s * tau; m = math.sqrt(M2[model]); H = 1 / (s * tau * tau); ell = math.log(m / H)
    a_p = s * tau_c if (tau_c > 0 and tau < tau_c) else a                  # U3b: potential scale factor (core width 1/(m a_p))
    if tau_f > 0 and tau > tau_f: a_p = s * tau_f                          # U19b: late frozen core
    ell_eff = math.log(m * a_p * tau)
    E = np.zeros(3); lib().energy(model, N, f, p, a, E)
    cnt, seed, n = strings(model, N, f)
    KG = np.empty(N ** 3, np.float32); KM = np.empty(N ** 3, np.float32)
    lib().kinetic_split(model, N, f, p, n, a, KG, KM)
    rc = [r / (m * a_p) for r in rmask]                       # mask radii in comoving cells (core width uses a_p)
    RCAP = min(N / 4, 64.0)                                   # radii beyond min(N/4, 64) cells: invalid (frac -> 0); = N/4 for N <= 256 (U2)
    pad = int(math.ceil(min(max(rc), RCAP))) + 2
    dist = distance_to_strings(N, seed, pad).ravel()
    Lc = (2.0 / 3.0) * cnt                                   # comoving length (dx = 1)
    xi = Lc * tau * tau / (4.0 * N ** 3)
    out = dict(tau=tau, a=a, a_p=a_p, H=H, ell=ell, ell_eff=ell_eff, fat=bool(a_p != a), n_plaq=int(cnt), xi=xi, E_kin=E[0], E_grad=E[1], E_pot=E[2],
               KG_all=float(KG.mean(dtype=np.float64)), KM_all=float(KM.mean(dtype=np.float64)))
    for r, rcells in zip(rmask, rc):
        if rcells > RCAP:
            out[f"rm{r:g}"] = dict(frac=0.0, rho_a=float("nan"), rho_M=float("nan")); continue
        keep = dist > rcells
        fr = float(keep.mean())
        out[f"rm{r:g}"] = dict(frac=fr, rho_a=float(2 * KG[keep].mean(dtype=np.float64)) if fr > 0 else float("nan"),
                               rho_M=float(2 * KM[keep].mean(dtype=np.float64)) if fr > 0 else float("nan"))
    if spectra:
        out["spectrum"] = tangent_spectrum(model, N, f, p, n, a, dist > min(rc[1], RCAP))
        if massive_spectrum:                                  # U10a: spectrum of the massive (normal) velocity, same mask
            out["spectrum_M"] = tangent_spectrum(model, N, f, p, n, a, dist > min(rc[1], RCAP), normal=True)
    return (out, dist) if return_dist else out


_KB = {}


def _kbins(N):
    if N not in _KB:
        kk = (np.fft.fftfreq(N) * 2 * np.pi).astype(np.float32); kz = (np.fft.rfftfreq(N) * 2 * np.pi).astype(np.float32)
        K = np.sqrt(kk[:, None, None] ** 2 + kk[None, :, None] ** 2 + kz[None, None, :] ** 2)
        nb = N // 2; edges = np.linspace(0, math.pi, nb + 1)
        idx = np.clip(np.digitize(K.ravel(), edges) - 1, 0, nb - 1).astype(np.int32)
        wz = np.ones(K.shape, np.float32); wz[:, :, 1:(-1 if N % 2 == 0 else None)] = 2      # rfft double counting
        _KB[N] = (edges, idx, wz.ravel(), nb)
    return _KB[N]


def tangent_spectrum(model, N, f, p, n, a, keep, normal=False):
    import scipy.fft as sfft
    nc = 10 if model == 1 else 1
    v = np.empty((nc, N ** 3), np.float32)
    (lib().normal_velocity if normal else lib().tangent_velocity)(model, N, f, p, n, a, v)
    edges, idx, wz, nb = _kbins(N)
    P = np.zeros(nb)
    mk = keep.reshape(N, N, N).astype(np.float32)
    for c in range(nc):
        F = sfft.rfftn(v[c].reshape(N, N, N) * mk, workers=int(os.environ.get("OMP_NUM_THREADS", os.cpu_count() or 2)))
        P += np.bincount(idx, weights=((F.real ** 2 + F.imag ** 2).ravel() * wz), minlength=nb)[:nb]
        del F
    del v
    P *= 0.5 / N ** 6 / max(float(mk.mean()), 1e-12)        # energy per comoving k-bin per unit volume (sum = 1/2 <v^2>)
    return dict(k_edges=edges.tolist(), dE=P.tolist())


# ---------- run ----------
KXI = {1: 50.0, 0: 27.2}     # xi_init * lc^2 at tau_i = 18.65 (calibration of amendment 4, N=128, 5 lc values)


def exchange_sample(model, N, f, p, tau, s, keep):
    """U10a: mean X, e_kin, e_grad (physical) over keep (uint8 mask or None) and the kept fraction"""
    out = np.zeros(4); a = s * tau if s > 0 else 1.0
    lib().exchange(model, N, f, p, a, None if keep is None else keep.ctypes.data, out)
    c = max(out[3], 1.0)
    return dict(tau=tau, X=out[0] / c, e_kin=out[1] / c, e_grad=out[2] / c, frac=out[3] / N ** 3)


def run(model, N, tau_end, seed, sfac=1.0, xi_init=0.5, dtau=0.25, nsnap=24, spectra_every=3, out_dir=None, log=print,
        damp_ell=None, exch=False, exch_samples=8, exch_stride=4, ell_i=1.0):
    """amendment 4: start at ell_i = 1 (tau0 = sqrt(e tau_end / sfac)) with lc chosen so that the initial xi ~ xi_init.
    U10a (ESPEC_U10a): damp_ell -> drain of the massive velocity at physical rate m until ell = damp_ell;
    exch -> after every snapshot, exch_samples evaluations of the exchange X every exch_stride steps, on the cells farther
    than 4/(m a) + 0.5 * exch_samples * exch_stride * dtau cells from the snapshot's strings (strings move < 0.5 cell/tau),
    plus the massive-velocity spectrum at every spectral snapshot.  With damp_ell=None and exch=False this is the old driver;
    with exch only, the evolution is bit-identical (the extra calls only read the fields)."""
    m = math.sqrt(M2[model]); s = sfac / (m * tau_end)
    tau0 = round(math.sqrt(math.exp(ell_i) * tau_end / sfac) / dtau) * dtau      # U10b: ell_i > 1 starts later (default 1: unchanged)
    lc = (tau0 / 18.65) * math.sqrt(KXI[model] / xi_init)
    f, p = initial(model, N, lc, seed); scratch = np.empty_like(f)
    lib()
    tau_d = math.sqrt(math.exp(damp_ell) / (m * s)) if damp_ell else 0.0
    lib().set_damp(tau_d, m if damp_ell else 0.0)
    taus = np.unique(np.round(np.geomspace(tau0 * 1.08, tau_end, nsnap) / dtau) * dtau)
    rows, tau, t0, keepX = [], tau0, time.time(), None
    for j, ts in enumerate(taus):
        nst = int(round((ts - tau) / dtau))
        samples = []
        if exch and keepX is not None and nst > 0:
            for _ in range(min(exch_samples, nst // exch_stride)):
                lib().evolve(model, N, f, p, scratch, tau, dtau, s, exch_stride); tau += exch_stride * dtau; nst -= exch_stride
                smp = exchange_sample(model, N, f, p, tau, s, keepX); smp["ell"] = math.log(m * s * tau * tau); samples.append(smp)
            rows[-1]["exch"] = samples
        if nst > 0:
            lib().evolve(model, N, f, p, scratch, tau, dtau, s, nst); tau += nst * dtau
        spec = (j % spectra_every == spectra_every - 1 or j == len(taus) - 1)
        r, dist = measure(model, N, f, p, tau, s, spectra=spec, massive_spectrum=exch, return_dist=True)
        if exch:
            rX = 4.0 / (m * s * tau) + 0.5 * exch_samples * exch_stride * dtau
            keepX = np.ascontiguousarray((dist.ravel() > rX).astype(np.uint8)); r["exch_mask_cells"] = rX
        del dist
        rows.append(r)
        log(f"[{'RP3' if model else 'U1'} N={N} seed={seed} sfac={sfac}] tau={tau:.1f} ell={r['ell']:.2f} xi={r['xi']:.3f} "
            f"rho_a={r['rm2.5']['rho_a']:.3e} frac={r['rm2.5']['frac']:.2f} ({time.time()-t0:.0f}s)")
    lib().set_damp(0.0, 0.0)
    res = dict(model="RP3" if model else "U1", N=N, tau_end=tau_end, seed=seed, sfac=sfac, s=s, tau0=tau0, dtau=dtau, lc=lc, xi_init=xi_init,
               f2=F2, m2=M2[model], threads=lib().nthreads(), runtime_s=time.time() - t0, rows=rows,
               damp=dict(ell=damp_ell, tau=tau_d, rate=(m if damp_ell else 0.0)), exch=exch, ell_i=ell_i)
    if out_dir:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        (Path(out_dir) / f"{res['model']}_N{N}_s{seed}_f{sfac:g}.json").write_text(json.dumps(res, default=float))
    return res


def run_u3(model, N, tau_end, seed, w=3.0, xi_init=0.5, dtau=0.25, nsnap=30, spectra_every=3, log=print):
    """U3b (ESPEC_U3.md): fat (constant comoving core width w cells) until tau_c = tau_end / w, then physical strings."""
    m = math.sqrt(M2[model]); s = 1.0 / (m * tau_end); tau_c = tau_end / w
    tau0 = round(math.e * w / dtau) * dtau                     # ell_eff(tau0) = ln(m s tau_c tau0) = ln(tau0 / w) = 1
    lc = (tau0 / 18.65) * math.sqrt(KXI[model] / xi_init)
    f, p = initial(model, N, lc, seed); scratch = np.empty_like(f)
    lib().set_core(tau_c)
    taus = np.unique(np.round(np.geomspace(tau0 * 1.08, tau_end, nsnap) / dtau) * dtau)
    rows, tau, t0 = [], tau0, time.time()
    for j, ts in enumerate(taus):
        nst = int(round((ts - tau) / dtau))
        if nst > 0:
            lib().evolve(model, N, f, p, scratch, tau, dtau, s, nst); tau += nst * dtau
        r = measure(model, N, f, p, tau, s, spectra=(j % spectra_every == spectra_every - 1 or j == len(taus) - 1), tau_c=tau_c)
        rows.append(r)
        log(f"[U3b {'RP3' if model else 'U1'} N={N} seed={seed} w={w:g} xi_i={xi_init:g}] tau={tau:.1f} ell={r['ell']:.2f} "
            f"ell_eff={r['ell_eff']:.2f} fat={r['fat']} xi={r['xi']:.3f} rho_a={r['rm2.5']['rho_a']:.3e} frac={r['rm2.5']['frac']:.2f} ({time.time()-t0:.0f}s)")
    lib().set_core(0.0)
    return dict(model="RP3" if model else "U1", N=N, tau_end=tau_end, seed=seed, sfac=1.0, s=s, tau0=tau0, tau_c=tau_c, w=w, dtau=dtau, lc=lc,
                xi_init=xi_init, f2=F2, m2=M2[model], threads=lib().nthreads(), runtime_s=time.time() - t0, rows=rows)
