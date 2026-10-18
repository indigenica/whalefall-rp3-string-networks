"""Audit U4-U11, part 0: conventions read off the run summaries (H = 1/(a tau), ell = ln(m_r/H), sum dE = rho_a/2)."""
import json, math
import numpy as np
import a_lib as L

worst = dict(H=0.0, ell=0.0, spec=0.0); n = 0
for p in sorted(L.RUNS.glob("*.json")):
    run = json.loads(p.read_text()); m = math.sqrt(run["m2"])
    for r in run["rows"]:
        n += 1
        worst["H"] = max(worst["H"], abs(r["H"] * r["a"] * r["tau"] - 1))
        worst["ell"] = max(worst["ell"], abs(r["ell"] - math.log(m / r["H"])))
        if r.get("spectrum") and r["ell"] >= 4.0:
            worst["spec"] = max(worst["spec"], abs(sum(r["spectrum"]["dE"]) / (r["rm2.5"]["rho_a"] / 2) - 1))
print(n, "rows; worst deviations:", worst)
json.dump(dict(rows=n, worst=worst), open(L.Path(__file__).parent / "a0_units.json", "w"), indent=1)
