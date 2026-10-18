"""Pruebas rápidas del kernel de red (~1 min; necesitan gcc con OpenMP).  Uso:  python -B network/test_network.py

Las importa también tests/run_tests.py (las tres pruebas de red de `make test` viven aquí) y se exportan tal cual al
repositorio del artículo de la red (paper/prd_network/export/).
Quick tests of the lattice kernel; also shipped unchanged with the network-paper repository."""
import json, shutil, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def test_network_kernel_small():
    """lattice kernel (N = 32): linear Goldstone and quintet frequencies; Z2 string pair detection (needs gcc + OpenMP)"""
    if shutil.which("gcc") is None:
        print("   (skipped: no gcc)"); return
    import validate as V
    if hasattr(V, "quick_checks"):
        res = V.quick_checks(N=32)
        assert res["passed"], res


def test_network_reproduces_u2_run():
    """the lattice kernel shipped here (U2 kernel + optional U3b core growth, off by default) reproduces a stored U2 run:
    U(1), N = 128, tau_end = 64, seed 101 (~30 s; needs gcc + OpenMP)"""
    if shutil.which("gcc") is None:
        print("   (skipped: no gcc)"); return
    import lattice as Lt
    ref = json.loads((HERE / "results" / "u2" / "runs" / "U1_N128_s101_f1_dt0.25_xi0.5.json").read_text())
    r = Lt.run(0, 128, 64, 101, sfac=1.0, xi_init=0.5, dtau=0.25, log=lambda m: None)
    for a, b in zip(r["rows"], ref["rows"]):
        assert a["n_plaq"] == b["n_plaq"]
        assert abs(a["rm2.5"]["rho_a"] / b["rm2.5"]["rho_a"] - 1) < 1e-4 if b["rm2.5"]["frac"] > 0 else True


def test_distance_slab_exact():
    """U4b memory fix: slab-wise periodic distance transform == full transform wherever distance <= pad (needs gcc: the C kernel)"""
    if shutil.which("gcc") is None:
        print("   (skipped: no gcc)"); return
    import lattice as Lt
    from scipy import ndimage
    rng = np.random.default_rng(7)
    for N, pad, dens in ((48, 14, 0.002), (64, 30, 0.0002)):
        seed = (rng.random(N ** 3) < dens).astype(np.uint8); s = seed.reshape(N, N, N)
        d0 = ndimage.distance_transform_edt(np.pad(s, pad, mode="wrap") == 0)[pad:-pad, pad:-pad, pad:-pad]
        d1 = Lt.distance_to_strings(N, seed, pad, slab=16); m = d0 <= pad
        assert float(np.max(np.abs(d0[m] - d1[m]))) < 1e-5 and bool(np.all(d1[~m] > pad - 1e-6))


TESTS = [test_network_kernel_small, test_network_reproduces_u2_run, test_distance_slab_exact]

if __name__ == "__main__":
    fails = 0
    for t in TESTS:
        t0 = time.time()
        try:
            t(); print(f"PASS {t.__name__} ({time.time() - t0:.1f}s)")
        except Exception as e:  # noqa: BLE001
            fails += 1; print(f"FAIL {t.__name__}: {type(e).__name__}: {e}")
    print("ALL PASSED" if not fails else f"{fails} FAILED"); sys.exit(1 if fails else 0)
