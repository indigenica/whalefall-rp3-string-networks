"""U2 production (ESPEC_U2.md §4 + amendment 1).  usage: python produce.py SET   (SET in D, A_U1, A_RP3, B, C)"""
import sys, json, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt
OUT = Path(__file__).resolve().parent / "results" / "u2" / "runs"
# (model, N, tau_end, seed, sfac, dtau, xi_init)
SETS = {
    "D": [(m, 128, 64, 101, 1.0, dt, 0.5) for m in (1, 0) for dt in (0.25, 0.125)],
    "A_U1": [(0, 256, 128, sd, 1.0, 0.25, 0.5) for sd in (1, 2, 3, 4)],
    "A_RP3a": [(1, 256, 128, sd, 1.0, 0.25, 0.5) for sd in (1, 2)],
    "A_RP3b": [(1, 256, 128, sd, 1.0, 0.25, 0.5) for sd in (3, 4)],
    "E": [(m, 256, 128, 31, 1.0, 0.25, xi) for m in (0, 1) for xi in (0.25, 1.0)],
    "B": [(m, 192, 96, sd, 1.0, 0.25, 0.5) for m in (0, 1) for sd in (11, 12)],
    "C": [(m, 256, 128, 21, 0.6, 0.25, 0.5) for m in (0, 1)],
}
if __name__ == "__main__":
    log = open(OUT.parent / f"produce_{sys.argv[1]}.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    for model, N, te, sd, sf, dt, xi in SETS[sys.argv[1]]:
        tag = f"{'RP3' if model else 'U1'}_N{N}_s{sd}_f{sf:g}_dt{dt:g}_xi{xi:g}"
        if (OUT / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
        r = Lt.run(model, N, te, sd, sfac=sf, xi_init=xi, dtau=dt, log=lg)
        OUT.mkdir(parents=True, exist_ok=True); (OUT / f"{tag}.json").write_text(json.dumps(r, default=float))
        lg(f"DONE {tag} {r['runtime_s']:.0f}s")
