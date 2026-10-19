"""Pruebas rápidas del kernel de red (~1 min; necesitan gcc con OpenMP).  Uso:  python -B network/test_network.py

Las importa también tests/run_tests.py (las tres pruebas de red de `make test` viven aquí) y se exportan tal cual al
repositorio del artículo de la red (paper/prd_network/export/).
Quick tests of the lattice kernel; also shipped unchanged with the network-paper repository."""
import json, os, shutil, sys, time, unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


class Skip(unittest.SkipTest):
    """a test that could not run here (e.g. no gcc): reported as SKIP, never as PASS"""


def need_gcc():
    if shutil.which("gcc") is None: raise Skip("no gcc")


def test_network_kernel_small():
    """lattice kernel (N = 32): linear Goldstone and quintet frequencies; Z2 string pair detection (needs gcc + OpenMP)"""
    need_gcc()
    import validate as V
    if hasattr(V, "quick_checks"):
        res = V.quick_checks(N=32)
        assert res["passed"], res


def test_network_reproduces_u2_run():
    """the lattice kernel shipped here (U2 kernel + optional U3b core growth, off by default) reproduces a stored U2 run:
    U(1), N = 128, tau_end = 64, seed 101: identical pierced-plaquette counts at every snapshot and rho_a to 1e-4 relative
    (~30 s; needs gcc + OpenMP)"""
    need_gcc()
    import lattice as Lt
    ref = json.loads((HERE / "results" / "u2" / "runs" / "U1_N128_s101_f1_dt0.25_xi0.5.json").read_text())
    r = Lt.run(0, 128, 64, 101, sfac=1.0, xi_init=0.5, dtau=0.25, log=lambda m: None)
    for a, b in zip(r["rows"], ref["rows"]):
        assert a["n_plaq"] == b["n_plaq"]
        assert abs(a["rm2.5"]["rho_a"] / b["rm2.5"]["rho_a"] - 1) < 1e-4 if b["rm2.5"]["frac"] > 0 else True


def test_distance_slab_exact():
    """U4b memory fix: slab-wise periodic distance transform == full transform wherever distance <= pad (needs gcc: the C kernel)"""
    need_gcc()
    import lattice as Lt
    from scipy import ndimage
    rng = np.random.default_rng(7)
    for N, pad, dens in ((48, 14, 0.002), (64, 30, 0.0002)):
        seed = (rng.random(N ** 3) < dens).astype(np.uint8); s = seed.reshape(N, N, N)
        d0 = ndimage.distance_transform_edt(np.pad(s, pad, mode="wrap") == 0)[pad:-pad, pad:-pad, pad:-pad]
        d1 = Lt.distance_to_strings(N, seed, pad, slab=16); m = d0 <= pad
        assert float(np.max(np.abs(d0[m] - d1[m]))) < 1e-5 and bool(np.all(d1[~m] > pad - 1e-6))


def _same_fields(a, b):
    """the evolved fields and momenta are equal bit for bit (the measured sums use OpenMP reductions, whose order is not fixed)"""
    return all(np.array_equal(x, y) for x, y in zip(a["_fields"], b["_fields"]))


def test_drain_rate_and_tau_grid_identity():
    """U24/U25 (ESPEC_U24, ESPEC_U25): with damp_rate = m given explicitly and with tau_grid = dtau the drained run is
    bit-identical to the defaults (RP3, N = 48, drain to ell = 2.5: final fields and momenta equal bit for bit)"""
    need_gcc()
    import lattice as Lt
    m = Lt.math.sqrt(Lt.M2[1])
    kw = dict(sfac=1.0, xi_init=0.5, dtau=0.25, nsnap=6, spectra_every=2, log=lambda x: None, damp_ell=2.5, return_fields=True)
    r0 = Lt.run(1, 48, 24, 61, **kw)
    assert _same_fields(r0, Lt.run(1, 48, 24, 61, damp_rate=m, **kw))
    assert _same_fields(r0, Lt.run(1, 48, 24, 61, tau_grid=0.25, **kw))
    assert not _same_fields(r0, Lt.run(1, 48, 24, 61, damp_rate=2 * m, **kw))                    # the rate does reach the kernel
    kh = dict(kw, dtau=0.125); r2 = Lt.run(1, 48, 24, 61, tau_grid=0.25, **kh)                      # U25: same snapshot times
    assert [x["tau"] for x in r2["rows"]] == [x["tau"] for x in r0["rows"]] and r2["tau0"] == r0["tau0"]


TESTS = [test_network_kernel_small, test_network_reproduces_u2_run, test_distance_slab_exact, test_drain_rate_and_tau_grid_identity]

def run_all(tests):
    """PASS / SKIP / FAIL per test; 'ALL PASSED' only if nothing failed and nothing was skipped.  With REQUIRE_ALL=1 (CI,
    where gcc exists) a skipped test is a failure."""
    fails, skips = 0, 0
    for t in tests:
        t0 = time.time()
        try:
            t(); print(f"PASS {t.__name__} ({time.time() - t0:.1f}s)")
        except Skip as e:
            skips += 1; print(f"SKIP {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            fails += 1; print(f"FAIL {t.__name__}: {type(e).__name__}: {e}")
    strict = os.environ.get("REQUIRE_ALL") == "1"
    if fails or (skips and strict):
        print(f"{fails} FAILED, {skips} SKIPPED"); return 1
    print("ALL PASSED" if not skips else f"PASSED, {skips} SKIPPED (not run here)"); return 0


if __name__ == "__main__":
    sys.exit(run_all(TESTS))
