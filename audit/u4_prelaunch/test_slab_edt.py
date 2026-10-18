"""Independent check of network/lattice.py::distance_to_strings (slab EDT) against the exact periodic Euclidean distance
(cKDTree with boxsize = N).  Claim under test: exact wherever the true distance <= pad; elsewhere the value is > pad.
Also checks the mask decision keep = dist > rc for every rc <= pad - 2 (what measure() uses).
Run:  python audit/u4_prelaunch/test_slab_edt.py"""
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "network"))
from lattice import distance_to_strings  # noqa: E402


def truth(N, seed):
    pts = np.argwhere(seed.reshape(N, N, N)).astype(float)
    g = np.indices((N, N, N)).reshape(3, -1).T.astype(float)
    d, _ = cKDTree(pts, boxsize=N).query(g, workers=2)
    return d.reshape(N, N, N)


def strings_seed(N, rng, nlines=6, npts=40):
    s = np.zeros((N, N, N), np.uint8)
    for _ in range(nlines):                      # random closed-ish curves (periodic random walks), like string seeds
        p = rng.integers(0, N, 3)
        for _ in range(N * 3):
            p = (p + rng.integers(-1, 2, 3)) % N; s[tuple(p)] = 1
    idx = rng.integers(0, N, (npts, 3)); s[idx[:, 0], idx[:, 1], idx[:, 2]] = 1   # isolated points
    return s.ravel()


fails = 0
rng = np.random.default_rng(7)
for N, pads, slabs, nl in ((48, (3, 8, 14), (5, 16, 64), 2), (64, (4, 18, 34), (7, 16, 64), 1), (40, (12, 22, 45), (3, 64), 0)):
    for trial in range(2):
        seed = strings_seed(N, rng, nlines=nl, npts=6 if nl == 0 else 30)
        T = truth(N, seed)
        for pad in pads:
            for slab in slabs:
                D = distance_to_strings(N, seed, pad, slab=slab)
                m = T <= pad
                err = float(np.max(np.abs(D[m] - T[m]))) if m.any() else 0.0
                lower_ok = bool(np.all(D[~m] > pad)) if (~m).any() else True
                masks_ok = all(np.array_equal(D > rc, T > rc) for rc in np.linspace(0.5, max(pad - 2, 0.5), 9))
                ok = err < 1e-5 and lower_ok and masks_ok
                fails += not ok
                print(f"N={N} trial={trial} pad={pad:>2} slab={slab:>2}  frac(T<=pad)={m.mean():.3f}  max|D-T| (T<=pad)={err:.2e}  "
                      f"D>pad where T>pad: {lower_ok}  masks identical: {masks_ok}  {'OK' if ok else 'FAIL'}")
print("ALL OK" if fails == 0 else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
