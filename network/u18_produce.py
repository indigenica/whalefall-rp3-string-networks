"""U18 production (docs/preregistros/ESPEC_U18.md): seed 61, RP3 and U(1), matched initial conditions of U16, at
(N, tau_end, r) = (192, 96, 1) and (320, 160, 0.6), with the 30 snapshot times of U16 plus 240 dense ones (same interval).
Spectra only at the U16 snapshots that had them.  Same evolution as lattice.run (the measurements only read the fields).
    python -B network/u18_produce.py            ->  network/results/u18/runs/<TAG>.json"""
import json, math, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt
import u16_produce as P16

OUT = Path(__file__).resolve().parent / "results" / "u18" / "runs"
CASES = [(192, 96, 1.0), (320, 160, 0.6)]; SEED = 61; DTAU = 0.25; NSNAP, NDENSE = 30, 240


def run(model, N, tau_end, r, log):
    m = math.sqrt(Lt.M2[model]); s = r / (m * tau_end)
    tau0 = round(math.sqrt(math.e * tau_end / r) / DTAU) * DTAU
    P16.MATCH_LC = P16.lc_of(model, P16.TAU_F, P16.R_F)
    f, p = P16.matched_initial(model, N, None, SEED); scratch = np.empty_like(f); Lt.lib(); Lt.lib().set_damp(0.0, 0.0)
    orig = np.unique(np.round(np.geomspace(tau0 * 1.08, tau_end, NSNAP) / DTAU) * DTAU)
    spec_at = {float(t) for j, t in enumerate(orig) if (j % 3 == 2 or j == len(orig) - 1)}
    dense = np.unique(np.round(np.geomspace(tau0 * 1.08, tau_end, NDENSE) / DTAU) * DTAU)
    taus = np.unique(np.concatenate([orig, dense]))
    rows, tau, t0 = [], tau0, time.time()
    for ts in taus:
        nst = int(round((ts - tau) / DTAU))
        if nst > 0:
            Lt.lib().evolve(model, N, f, p, scratch, tau, DTAU, s, nst); tau += nst * DTAU
        rr = Lt.measure(model, N, f, p, tau, s, spectra=(float(ts) in spec_at))
        rr["u16_snapshot"] = bool(np.any(np.isclose(orig, ts))); rows.append(rr)
        if rr["u16_snapshot"]:
            log(f"[U18 {'RP3' if model else 'U1'} N={N} r={r:g}] tau={tau:.2f} ell={rr['ell']:.3f} xi={rr['xi']:.3f} ({time.time() - t0:.0f}s)")
    return dict(model="RP3" if model else "U1", N=N, tau_end=tau_end, seed=SEED, sfac=r, s=s, tau0=tau0, dtau=DTAU, xi_init=P16.XI_I,
                f2=Lt.F2, m2=Lt.M2[model], threads=Lt.lib().nthreads(), runtime_s=time.time() - t0, rows=rows,
                u18=dict(n_u16_snapshots=int(len(orig)), n_total=int(len(taus)), matched_ic=True, spec="docs/preregistros/ESPEC_U18.md"))


def main():
    OUT.mkdir(parents=True, exist_ok=True); log = open(OUT.parent / "produce.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    for N, te, r in CASES:
        for model in (1, 0):
            tag = f"{'RP3' if model else 'U1'}_N{N}_t{te}_r{r:g}_s{SEED}_dense"
            if (OUT / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
            res = run(model, N, te, r, lg)
            (OUT / f"{tag}.json").write_text(json.dumps(res, default=float)); lg(f"DONE {tag} {res['runtime_s']:.0f}s")


if __name__ == "__main__":
    main()
