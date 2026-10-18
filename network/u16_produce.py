"""U16 production (docs/preregistros/ESPEC_U16.md): resolution and volume arms of the RP3 / U(1) networks.
  resolution arm   (N, tau_end, r) = (192, 96, 1), (256, 128, 0.75), (320, 160, 0.6): ell_end = ln 96 and L H a = 2 at the end
                   in all three; initial conditions MATCHED: one physical field per seed, white noise on 320^3 filtered with the
                   320^3 correlation length, truncated in Fourier space to 256^3 and 192^3 (the filter is < 1e-10 at the cut)
  volume arm       N = 320, tau_end = 96, r = 1 (L = 3.33 tau_end), plain lattice.initial
U6 protocol otherwise: no drain, xi_i = 0.5, ell_i = 1, 30 snapshots, spectra every 3rd.  Runs already present are skipped.
    python -B network/u16_produce.py [res|vol|all|res4] [--only TAG]   (res4: amendment 1, N = 384, r = 0.5)        ->  network/results/u16/runs/<TAG>.json"""
import json, math, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt

OUT = Path(__file__).resolve().parent / "results" / "u16" / "runs"
NF, TAU_F, R_F = 320, 160, 0.6                 # finest grid of the resolution arm: the matched field is generated here
RES = [(192, 96, 1.0), (256, 128, 0.75), (320, 160, 0.6)]
VOL = [(320, 96, 1.0)]
RES4 = [(384, 192, 0.5)]                      # amendment 1: fourth resolution (matched field zero-padded from 320^3)
SEEDS_RES, SEEDS_VOL = (61, 62, 63), (61, 62)
XI_I, DTAU = 0.5, 0.25


def lc_of(model, tau_end, r):
    """the correlation length lattice.run would use (ell_i = 1), in cells of that grid"""
    tau0 = round(math.sqrt(math.e * tau_end / r) / DTAU) * DTAU
    return (tau0 / 18.65) * math.sqrt(Lt.KXI[model] / XI_I)


def _keep(Nf, N):
    """indices of the full-FFT axis of an Nf grid that survive truncation to N (Nyquist of N dropped)"""
    h = N // 2
    return np.r_[0:h, Nf - h + 1:Nf]


def matched_gaussian(N, ncomp, seed):
    """band-limited Gaussian components: white noise on NF^3, P(k) filter of the NF grid, truncated to N^3"""
    rng = np.random.default_rng(seed); lcf = MATCH_LC
    k = np.fft.fftfreq(NF) * 2 * np.pi; kz = np.fft.rfftfreq(NF) * 2 * np.pi
    ia, iz = _keep(NF, N), np.arange(N // 2)
    out = np.empty((ncomp, N, N, N), np.float32)
    for c in range(ncomp):
        w = rng.standard_normal((NF, NF, NF))
        if N <= NF:                            # truncation (192, 256, 320)
            W = np.fft.rfftn(w)[np.ix_(ia, ia, iz)]
            K2 = k[ia][:, None, None] ** 2 + k[ia][None, :, None] ** 2 + kz[iz][None, None, :] ** 2
            W *= np.exp(-0.25 * K2 * lcf * lcf)
            h = N // 2; sel = np.r_[0:h, N - h + 1:N]; zsel = np.arange(h)
        else:                                  # amendment 1: zero-padding of the NF modes (exact: band-limited field)
            ja, jz = _keep(NF, NF), np.arange(NF // 2)
            W = np.fft.rfftn(w)[np.ix_(ja, ja, jz)]
            K2 = k[ja][:, None, None] ** 2 + k[ja][None, :, None] ** 2 + kz[jz][None, None, :] ** 2
            W *= np.exp(-0.25 * K2 * lcf * lcf)
            h = NF // 2; sel = np.r_[0:h, N - h + 1:N]; zsel = np.arange(h)
        Wn = np.zeros((N, N, N // 2 + 1), np.complex128)
        Wn[np.ix_(sel, sel, zsel)] = W
        out[c] = np.fft.irfftn(Wn, s=(N, N, N)).astype(np.float32)
        del w, W, Wn
    return out


def matched_initial(model, N, lc_ignored, seed):
    if model == 1:
        g = matched_gaussian(N, 10, seed)
        n = np.empty((4, N, N, N), np.float32)
        Lt.lib().principal(N, np.ascontiguousarray(g.reshape(10, -1)), n.reshape(4, -1), None)
        f = Lt.q_from_n(n)
    else:
        g = matched_gaussian(N, 2, seed); rr = np.sqrt(g[0] ** 2 + g[1] ** 2) + 1e-30
        f = (math.sqrt(0.5) * g / rr).astype(np.float32)
    return np.ascontiguousarray(f.reshape(Lt.NCOMP[model], -1)), np.zeros((Lt.NCOMP[model], N ** 3), np.float32)


MATCH_LC = None
_plain_initial = Lt.initial


def tag(model, N, tau_end, r, seed, arm):
    return f"{'RP3' if model else 'U1'}_N{N}_t{tau_end}_r{r:g}_s{seed}_{arm}"


def jobs(which):
    out = []
    if which in ("res", "all"):
        for sd in SEEDS_RES:
            for model in (1, 0):
                for N, te, r in RES: out.append((model, N, te, r, sd, "res"))
    if which == "res4":
        for sd in SEEDS_RES:
            for model in (1, 0):
                for N, te, r in RES4: out.append((model, N, te, r, sd, "res"))
    if which in ("vol", "all"):
        for sd in SEEDS_VOL:
            for model in (1, 0):
                for N, te, r in VOL: out.append((model, N, te, r, sd, "vol"))
    return out


def main(argv):
    global MATCH_LC
    which = argv[1] if len(argv) > 1 and not argv[1].startswith("--") else "all"
    only = argv[argv.index("--only") + 1] if "--only" in argv else None
    OUT.mkdir(parents=True, exist_ok=True)
    log = open(OUT.parent / "produce.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    for model, N, te, r, sd, arm in jobs(which):
        t = tag(model, N, te, r, sd, arm)
        if only and t != only: continue
        if (OUT / f"{t}.json").exists(): lg(f"skip {t}"); continue
        if arm == "res":
            MATCH_LC = lc_of(model, TAU_F, R_F); Lt.initial = matched_initial
        else:
            Lt.initial = _plain_initial
        t0 = time.time()
        res = Lt.run(model, N, te, sd, sfac=r, xi_init=XI_I, dtau=DTAU, nsnap=30, spectra_every=3, log=lg)
        res["u16"] = dict(arm=arm, r=r, matched_ic=(arm == "res"), match_grid=NF if arm == "res" else None,
                          match_lc_cells_on_NF=MATCH_LC if arm == "res" else None, L_over_tau_end=N / te,
                          spec="docs/preregistros/ESPEC_U16.md")
        (OUT / f"{t}.json").write_text(json.dumps(res, default=float)); lg(f"DONE {t} {time.time() - t0:.0f}s")
        Lt.initial = _plain_initial


if __name__ == "__main__":
    main(sys.argv)
