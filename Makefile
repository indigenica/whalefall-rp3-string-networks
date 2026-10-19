# whalefall-rp3-string-networks: code, data and paper for the RP^3 vs U(1) global-string network paper.
# Generated from indigenica/whalefall (see EXPORT_PROVENANCE.json). Do not edit here: edit whalefall and re-export.
PY ?= python3
export PYTHONDONTWRITEBYTECODE = 1
export OMP_NUM_THREADS ?= 2

.PHONY: all test validate audit numbers figures paper clean network-u2 network-u3

all: test audit numbers figures ## everything that runs in a few minutes on a laptop

test:                           ## kernel tests: linear spectrum, Z2 detection, reproduction of a stored U2 run (string counts exact, rho_a to 1e-4), drain/time-grid identity (~1 min; gcc + OpenMP)
	$(PY) -B network/test_network.py

validate:                       ## kernel validations V1-V4 of the U2 specification (~1-2 min); rewrites network/results/u2/validation.json
	cd network && $(PY) -B validate.py

audit:                          ## independent re-implementation of the estimators over the stored networks; compares with the published numbers
	cd audit/u4_u11 && $(PY) -B a0_units.py && $(PY) -B a1_network.py

numbers:                        ## every number quoted in main.tex, recomputed from the stored results (exit 1 on any mismatch)
	$(PY) -B paper/prd_network/paper_numbers.py

figures:                        ## paper figures from the stored network summaries (fig_bracket.pdf, fig_bath.pdf)
	$(PY) -B paper/prd_network/make_figures.py

paper: figures                  ## PDF (needs latexmk and a TeX distribution)
	cd paper/prd_network && latexmk -pdf -interaction=nonstopmode main.tex

# ---- production (hours of CPU; the 512-768^3 networks need a large machine, see README) ----
network-u2:                     ## the 20 networks of U2 (256^3; ~3 h with 2 cores) + analysis
	cd network && for S in D A_U1 A_RP3a E B C; do $(PY) -B produce.py $$S; done && $(PY) -B analyze.py

network-u3:                     ## U3b core-preparation networks (~2 h) + gates
	cd network && for S in F1_U1 F2 F1_RP3; do $(PY) -B produce_u3.py $$S; done && $(PY) -B analyze_u3.py

clean:
	rm -rf network/_build paper/prd_network/*.aux paper/prd_network/*.log paper/prd_network/*.fls paper/prd_network/*.fdb_latexmk paper/prd_network/*.out paper/prd_network/*.bbl paper/prd_network/*.blg
