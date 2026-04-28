PYTHON ?= python
OUT ?= outputs

.PHONY: install install-dev test lint format run robustness effect-sign pc-sensitivity projection clean

install:
	$(PYTHON) -m pip install -r requirements.txt

install-dev:
	$(PYTHON) -m pip install -r requirements-dev.txt

test:
	$(PYTHON) -m compileall .
	$(PYTHON) -m pytest -q

lint:
	$(PYTHON) -m ruff check .

format:
	$(PYTHON) -m black *.py tests

run:
	$(PYTHON) sdi_pipeline.py --out $(OUT)

robustness:
	$(PYTHON) sdi_robustness_checks.py --out $(OUT) --sdi_col SDI --traits headsmut_greenhouse_avg,headsmut_highest_score --n_pcs 3 --n_perm 2000

effect-sign:
	$(PYTHON) effect_sign_enrichment.py --out $(OUT) --disease_pattern "GWAS_headsmut*.csv" --top_pcts "0.5,1,2,5" --use_fdr

pc-sensitivity:
	$(PYTHON) gwas_pc_sensitivity.py --out $(OUT) --trait SDI --pc_list "0,1,2,3,4,6,8,10" --include_race "0,1"

projection:
	$(PYTHON) genotype_to_spectral_projection.py --input $(OUT)/Supplementary_Data_S1.xlsx --sheet S1_Phenotypes --out $(OUT)/projection

clean:
	rm -rf __pycache__ .pytest_cache .ruff_cache *.egg-info build dist
