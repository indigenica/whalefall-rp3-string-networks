# RP³ ≃ SO(3) global-string networks relative to U(1)

[![Preprint DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23451009.svg)](https://doi.org/10.5281/zenodo.23451009)
[![reproducibility](https://github.com/indigenica/whalefall-rp3-string-networks/actions/workflows/ci.yml/badge.svg)](https://github.com/indigenica/whalefall-rp3-string-networks/actions/workflows/ci.yml)

Code, network summaries, preregistered specifications, independent audit and paper for

> J. Bautista, *Density and Goldstone emission of an RP³ ≃ SO(3) global-string network relative to U(1), and a
> massive-bath systematic of the emission estimator*, preprint (2026), Zenodo, [doi:10.5281/zenodo.23451009](https://doi.org/10.5281/zenodo.23451009).
> Source in `paper/prd_network/`.

We evolve networks of global ℤ₂ strings with vacuum manifold RP³ (a real symmetric 4×4 order parameter) and U(1)
axion-string networks in the same code, with matched Goldstone decay constant and core mass, up to 576³ (U(1) control
also at 768³), and measure:

| Quantity | Value | Where |
|---|---|---|
| String-density ratio R_ξ = ξ_RP³/ξ_U(1) | ≈ 1.6 for every protocol; observed 1.57 (from above) and 1.68 (from below) at ℓ = 5.61 | U6, U9, U10b, U11 |
| Massive-bath systematic of the kinetic emission estimator | 19 % (RP³) and 28 % (U(1)) of the sparse-start emission; > 50 % for dense starts | U9, U10a, U11, T7 |
| Emission ratio with the initial bath drained, R_Γ = Γ_RP³/Γ_U(1) | 1.56 ± 0.09 (stat) ± 0.08 (snapshot sampling) | U11, U18 |
| Spectral slope q_RP³(ℓ = 5) | 0.87 ± 0.06; the momentum-window convention moves it by up to 0.04 and a spectrum-weighted fit gives 0.94 ± 0.04, so it is quoted as compatible with 1 | U11, audit |
| Drain and time step (384³, three seeds) | the drain end and rate change R_Γ by at most 1.9 ± 1.4 % (doubling the rate lowers R_ξ by 5.6 ± 1.5 %); dτ/2 changes R_Γ by −0.10 ± 0.14 % | U24, U25 |
| Resolution / volume | R_ξ changes by −0.26 % between m_r a dx = 1 and 0.6; a box with 4.6× more Hubble volumes changes R_ξ by −0.6 ± 5.5 % and R_Γ by +0.7 ± 4.9 % | U16 |

## Reproduce

Requirements: Python ≥ 3.10 with `numpy`, `scipy`, `matplotlib` (`pip install -r requirements.txt`), and `gcc` with
OpenMP for the lattice kernel (compiled on first use into `network/_build/`). TeX Live with `latexmk` and `revtex4-2` (Debian/Ubuntu: `texlive-publishers`) only for the PDF.

**Tested environment.** CI runs on Ubuntu (GitHub `ubuntu-latest`) with Python 3.11 and the newest `numpy`, `scipy` and
`matplotlib` allowed by `requirements.txt` at the time of the run, and the system `gcc` with OpenMP. The analysis numbers
do not depend on those versions (the audit step compares them to 1e-8). Without `gcc` the kernel tests are reported as
`SKIP`, never as passed, and `make test` does not print `ALL PASSED`; CI sets `REQUIRE_ALL=1`, so there a skipped test fails.

```bash
make test      # kernel tests, incl. the reproduction of a stored 128^3 U2 run: identical string counts, rho_a to 1e-4 (~1 min)
make audit     # independent re-implementation of every estimator over the stored networks (audit/u4_u11/)
make numbers   # every number quoted in main.tex, recomputed from the stored results
make figures   # paper figures from the stored network summaries
make paper     # paper/prd_network/main.pdf
make validate  # kernel validations V1-V4 (linear spectrum, string detector, energy conservation, estimator)
```

`make test audit numbers figures` runs in a few minutes on a laptop and is what CI runs. The audit rewrites
`audit/u4_u11/a1_network.json`; `git diff` on that file should be empty.

**Production runs.** `make network-u2` and `make network-u3` reproduce the 256³ networks (a few CPU hours, ~1.5 GB).
The 512–768³ networks of U4b–U11, U16 and U18 need ~30–50 GB of RAM and tens of cores; the scripts are
`network/u4b_produce.py`, `network/u16_produce.py` and `network/u18_produce.py` (one network per call, resumable), and the
analyses `network/u*_analyze.py`. The 384³ drain and time-step tests (U24, U25) fit in 16 GB and four cores
(`network/u24_produce.py`, which skips finished networks, and `network/u25_produce.py`). Every analysis script reads only the stored summaries in `network/results/`.

## Layout

```
network/            lattice kernel (kern.c, C + OpenMP), initial conditions and observables (lattice.py),
                    validations, production and analysis scripts for each study
network/results/    one JSON summary per network (runs/*.json: snapshots with ell, xi, Goldstone energy and spectra)
                    and one *_RESULT.json per study with its preregistered verdict
audit/u4_u11/       independent audit: estimators re-implemented from the specifications only (a_lib.py)
audit/u2/, audit/u4_prelaunch/   audit of the first networks; pre-launch review of the large-network code
docs/               the dated preregistrations (preregistros/ESPEC_U*.md) and reports (informes/INFORME_U*.md)
paper/prd_network/  main.tex, refs.bib, figure script and figures
results/s3_sectors.json   tree-level widths of the massive modes (used by network/t7_drain_fade.py)
EXPORT_PROVENANCE.json    source commit and sha256 of every file (see below)
```

`docs/README.md` maps the internal study codes (U2 … U25) to the sections of the paper.

## Preregistration and audit

Every analysis that gives a verdict was written as a dated specification (estimators, windows, thresholds, verdict
logic) before its networks were produced; later changes are dated amendments at the end of the same file, never
rewrites. The specifications and reports are the original records and are in Spanish. The audit in `audit/u4_u11/`
re-implemented the estimators without reading the analysis code and reproduces the density, emission and drain
numbers of the production networks to ≤ 4×10⁻³ relative (R_Γ itself to 2×10⁻⁵); the resolution (U16), sampling (U18),
drain (U24) and time-step (U25) tests were not audited. The audit was carried out by a separate session of the AI
assistant used for the code (Claude Code), independent of the audited code but not of the tool (paper, §5).

## Provenance

This repository is a one-way export of the larger project [`indigenica/whalefall`](https://github.com/indigenica/whalefall),
a toy cosmology in which this network feeds a dark-matter component. Nothing here depends on that model.
`EXPORT_PROVENANCE.json` records the whalefall commit of each export and the hash of every file. Changes are made
in whalefall and exported again; files are not edited here.

## Citation and licence

Please cite the paper: J. Bautista (2026), Zenodo, [doi:10.5281/zenodo.23451009](https://doi.org/10.5281/zenodo.23451009) (also in `CITATION.cff`, and under "Cite this repository" on GitHub).

Code: GPL-3.0-or-later. Paper, figures, data and documentation: CC BY-NC 4.0. Research and teaching use is free; commercial use of the text, figures or data, or use of the code in closed-source software, needs the author's permission. Details in `LICENSE.md`.
