"""Independent re-implementation of the U4-U11 network estimators (audit/u4_u11/PLAN.md).

Written from the definitions in ESPEC_U2 (sec. 3), ESPEC_U6, ESPEC_U10a and ESPEC_U11 only; it does not read or import any
network/*_analyze.py.  Input: the run summaries in network/results/u4/runs/*.json.

Conventions read off the data (checked in a0_units.py): H = 1/(a tau), ell = ln(m_r/H) with m_r = sqrt(m2), rho_a and rho_M
are physical energy densities outside the cores, and sum(spectrum.dE) = rho_a/2 (kinetic part of the masked Goldstone velocity)."""
import json, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "network" / "results" / "u4" / "runs"
RNG = np.random.default_rng(20261011)


def load(name):
    return json.loads((RUNS / name).read_text())


def ell_fin(run): return math.log(run["N"] / 2)


def window_grid(N, n=10):
    lf = math.log(N / 2)
    return np.linspace(lf - 1.3, lf - 0.05, n)


def series(run, mask="rm2.5"):
    R = run["rows"]
    tau = np.array([r["tau"] for r in R]); a = np.array([r["a"] for r in R]); H = np.array([r["H"] for r in R])
    ell = np.array([r["ell"] for r in R]); xi = np.array([r["xi"] for r in R])
    rho = np.array([r[mask]["rho_a"] for r in R]); rhoM = np.array([r[mask]["rho_M"] for r in R])
    frac = np.array([r[mask]["frac"] for r in R])
    return dict(tau=tau, a=a, H=H, ell=ell, xi=xi, rho=rho, rhoM=rhoM, frac=frac, f2=run["f2"], m=math.sqrt(run["m2"]))


def emission(run, scheme="central", mask="rm2.5"):
    """Gamma = H a^-4 d(a^4 rho_a)/d ln tau at each instant, normalised to eps_lat = Gamma / (8 pi H^3 f^2 ell * 0.24 ell).
    scheme: 'central' (neighbours i-1, i+1), 'central_smooth' (then a 3-point running mean), 'fit5' (LSQ slope on 5 points)."""
    s = series(run, mask)
    Y = s["a"] ** 4 * s["rho"]; x = np.log(s["tau"]); n = len(x)
    d = np.full(n, np.nan)
    if scheme in ("central", "central_smooth"):
        d[1:-1] = (Y[2:] - Y[:-2]) / (x[2:] - x[:-2])
        if scheme == "central_smooth":
            dd = d.copy(); d = np.full(n, np.nan); d[2:-2] = (dd[1:-3] + dd[2:-2] + dd[3:-1]) / 3
    elif scheme == "fit5":
        for i in range(2, n - 2):
            d[i] = np.polyfit(x[i - 2:i + 3], Y[i - 2:i + 3], 1)[0]
    G = s["H"] * d / s["a"] ** 4
    eps = G / (8 * math.pi * s["H"] ** 3 * s["f2"] * s["ell"] * 0.24 * s["ell"])
    return s["ell"], G, eps


def on_grid(ell, y, grid):
    ok = np.isfinite(y)
    return np.interp(grid, ell[ok], y[ok], left=np.nan, right=np.nan)


def ratio_boot(A, B, nboot=4000):
    """A, B: arrays (seeds x grid) for RP3 and U1.  Returns (R as mean over the grid of the ratio of seed means,
    R as ratio of window means, per-ell ratio, bootstrap std of the first)."""
    A = np.asarray(A); B = np.asarray(B)
    per = A.mean(0) / B.mean(0)                                   # grid points with no estimate (NaN) are dropped
    R1 = float(np.nanmean(per)); ok = np.isfinite(per); R2 = float(A[:, ok].mean() / B[:, ok].mean())
    bs, bs2 = [], []
    for _ in range(nboot):
        ia = RNG.integers(0, len(A), len(A)); ib = RNG.integers(0, len(B), len(B))
        bs.append(np.nanmean(A[ia].mean(0) / B[ib].mean(0))); bs2.append(A[ia][:, ok].mean() / B[ib][:, ok].mean())
    return R1, R2, per, float(np.std(bs, ddof=1)), float(np.std(bs2, ddof=1))


def R_gamma(rp3_runs, u1_runs, scheme="central", x_corr=None):
    N = rp3_runs[0]["N"]; g = window_grid(N)
    A = [on_grid(*emission(r, scheme)[::2], g) for r in rp3_runs]
    B = []
    for r in u1_runs:
        ell, G, eps = emission(r, scheme)
        if x_corr is not None:                                     # Gamma_U1 - <X>, X interpolated at the same ell
            le, X = exch_snapshots(r); X = X * x_corr
            ok = np.isfinite(G)
            Xi = np.interp(ell[ok], le, X, left=np.nan, right=np.nan)
            s = series(r); fac = 8 * math.pi * s["H"][ok] ** 3 * s["f2"] * s["ell"][ok] * 0.24 * s["ell"][ok]
            eps = eps.copy(); eps[ok] = (G[ok] - Xi) / fac
        B.append(on_grid(ell, eps, g))
    return ratio_boot(A, B), g, np.array(A), np.array(B)


def R_xi(rp3_runs, u1_runs):
    g = window_grid(rp3_runs[0]["N"])
    A = [on_grid(series(r)["ell"], series(r)["xi"], g) for r in rp3_runs]
    B = [on_grid(series(r)["ell"], series(r)["xi"], g) for r in u1_runs]
    return ratio_boot(A, B)


def integral_emission(run):
    """Delta(a^4 rho_a) over the window (ESPEC_U2 6.5), interpolated in ell."""
    s = series(run); g = window_grid(run["N"]); Y = s["a"] ** 4 * s["rho"]
    return float(np.interp(g[-1], s["ell"], Y) - np.interp(g[0], s["ell"], Y))


def q_points(run, lo=10.0, hi=0.25, min_bins=8):
    """q between consecutive snapshots with spectra: F(k) = a2^4 dE2 - a1^4 dE1 (comoving k), fit F ~ k^-q on
    k_phys = k/a_mid in [lo H_mid, hi m]."""
    R = [r for r in run["rows"] if r.get("spectrum")]; m = math.sqrt(run["m2"]); out = []
    for r1, r2 in zip(R[:-1], R[1:]):
        k = np.array(r1["spectrum"]["k_edges"]); kc = 0.5 * (k[1:] + k[:-1])
        assert np.allclose(k, r2["spectrum"]["k_edges"])
        F = r2["a"] ** 4 * np.array(r2["spectrum"]["dE"]) - r1["a"] ** 4 * np.array(r1["spectrum"]["dE"])
        am = math.sqrt(r1["a"] * r2["a"]); Hm = math.sqrt(r1["H"] * r2["H"]); lm = 0.5 * (r1["ell"] + r2["ell"])
        kp = kc / am; sel = (kp >= lo * Hm) & (kp <= hi * m) & (F > 0)
        nsel = int(np.sum((kp >= lo * Hm) & (kp <= hi * m)))
        if np.sum(sel) < min_bins: continue
        sl = np.polyfit(np.log(kc[sel]), np.log(F[sel]), 1)[0]
        out.append((lm, -sl, int(np.sum(sel)), nsel))
    return out


def q_fit(runs, ell_min=4.6, nboot=4000, **kw):
    pts = [q_points(r, **kw) for r in runs]; lf = ell_fin(runs[0])
    def fit(sets):
        L = np.array([p[0] for s in sets for p in s if p[0] >= ell_min]); Q = np.array([p[1] for s in sets for p in s if p[0] >= ell_min])
        q1, q0 = np.polyfit(L - 5, Q, 1); return q0, q1
    q0, q1 = fit(pts); bs = []
    for _ in range(nboot):
        idx = RNG.integers(0, len(pts), len(pts)); bs.append(fit([pts[i] for i in idx]))
    bs = np.array(bs)
    return dict(q0_at_ell5=float(q0), q1_per_ell=float(q1), q0_std=float(np.std(bs[:, 0], ddof=1)), q1_std=float(np.std(bs[:, 1], ddof=1)),
                q_end=float(q0 + q1 * (lf - 5)), n_points=int(sum(1 for s in pts for p in s if p[0] >= ell_min)), points=pts)


def xi_slope(run, lo, hi):
    s = series(run); m = (s["ell"] >= lo) & (s["ell"] <= hi)
    return float(np.polyfit(s["ell"][m], s["xi"][m], 1)[0])


def rhoM_over_rhoa_end(runs):
    v = []
    for r in runs:
        s = series(r); le = window_grid(r["N"])[-1]
        v.append(np.interp(le, s["ell"], s["rhoM"]) / np.interp(le, s["ell"], s["rho"]))
    return float(np.mean(v)), [float(x) for x in v]


def exch_samples(run):
    """X samples (ESPEC_U10a B): every row carries a list of {tau, ell, X, e_kin, e_grad, frac}; valid if frac >= 0.5."""
    return [e for r in run["rows"] for e in (r.get("exch") or []) if isinstance(e, dict) and e.get("frac", 0) >= 0.5]


def exch_snapshots(run):
    """<X> per snapshot: mean of the valid samples taken after it, at their mean ell."""
    pts = []
    for row in run["rows"]:
        ex = [e for e in (row.get("exch") or []) if isinstance(e, dict) and e.get("frac", 0) >= 0.5]
        if ex: pts.append((np.mean([e["ell"] for e in ex]), np.mean([e["X"] for e in ex])))
    le, X = np.array(pts).T
    return le, X


def f_X(run):
    """<X>/Gamma over the window (ratio of window means on the grid)."""
    ell, G, _ = emission(run); ok = np.isfinite(G); g = window_grid(run["N"]); le, X = exch_snapshots(run)
    Gg = np.interp(g, ell[ok], G[ok], right=np.nan); Xg = np.interp(g, le, X, right=np.nan); m = np.isfinite(Gg) & np.isfinite(Xg)
    return float(Xg[m].mean() / Gg[m].mean())
