"""Figure of U3b: xi(ell) of the U(1) control with and without core-growth preparation (U2 sets A/E vs U3b set F)."""
import json, glob
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def curves(pattern):
    out = []
    for p in sorted(glob.glob(str(HERE / pattern))):
        r = json.loads(Path(p).read_text()); R = r["rows"]
        out.append((Path(p).stem, np.array([x["ell"] for x in R]), np.array([x["xi"] for x in R]), np.array([x.get("fat", False) for x in R])))
    return out


def main():
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for name, l, x, fat in curves("results/u2/runs/U1_N256_*.json"):
        ax.plot(l, x, color="C0", lw=0.8, alpha=0.6, label="U(1), U2 (cuerdas físicas desde ℓ = 1)" if "s1_" in name else None)
    for name, l, x, fat in curves("results/u3/runs/U1_N256_*.json"):
        ax.plot(l[~fat], x[~fat], color="C3", lw=1.2, label=f"U(1), U3b preparación w = 3 ({name.split('_')[2]}, ξ_i={name.split('xi')[-1]})")
        ax.plot(l[fat], x[fat], color="C3", lw=0.6, ls=":")
    ax.axhspan(0.5, 2.5, color="gray", alpha=0.12); ax.text(3.05, 0.52, "puerta de control G1: ξ_U1(4.5) ∈ [0.5, 2.5]", fontsize=7)
    ax.axvline(4.5, color="k", lw=0.6, ls="--"); ax.set_xlim(1, 5); ax.set_ylim(0, 3)
    ax.set_xlabel("ℓ = ln(m/H)"); ax.set_ylabel("ξ"); ax.legend(fontsize=7); ax.grid(alpha=0.3, lw=0.3)
    ax.set_title("U3b · la preparación de núcleo no lleva el control U(1) al atractor con N = 256", fontsize=9)
    fig.tight_layout(); fig.savefig(HERE / "figures" / "u3b_control.png", dpi=140); plt.close(fig); print("figure network/figures/u3b_control.png")


if __name__ == "__main__":
    main()
