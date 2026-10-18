"""U3b production (ESPEC_U3.md).  usage: python produce_u3.py SET   (SET in PILOT, F1_U1, F1_RP3, F2)"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lattice as Lt
OUT = Path(__file__).resolve().parent / "results" / "u3" / "runs"
# (model, N, tau_end, seed, w, xi_init)
SETS = {
    "PILOT": [(m, 128, 64, 7, 3.0, 0.5) for m in (0, 1)],
    "F1_U1": [(0, 256, 128, sd, 3.0, 0.5) for sd in (41, 42)],
    "F1_RP3": [(1, 256, 128, sd, 3.0, 0.5) for sd in (41, 42)],
    "F2": [(m, 256, 128, 43, 3.0, 2.0) for m in (0, 1)],
}
if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    log = open(OUT.parent / f"produce_{sys.argv[1]}.log", "a")
    def lg(msg): print(msg, flush=True); log.write(msg + "\n"); log.flush()
    for model, N, te, sd, w, xi in SETS[sys.argv[1]]:
        tag = f"{'RP3' if model else 'U1'}_N{N}_s{sd}_w{w:g}_xi{xi:g}"
        if (OUT / f"{tag}.json").exists(): lg(f"skip {tag}"); continue
        r = Lt.run_u3(model, N, te, sd, w=w, xi_init=xi, log=lg)
        (OUT / f"{tag}.json").write_text(json.dumps(r, default=float)); lg(f"DONE {tag} {r['runtime_s']:.0f}s")
