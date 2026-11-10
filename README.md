# Does Wealth Burn Carbon? — Week 2: Statistical Analysis and Hypothesis Testing Framework

CO₂ emissions per person vs economic development, 1990–2022. This repository is the *execution phase* of the Week 1 project: five research questions, seven pre-specified hypotheses (H1–H6 confirmatory with Holm correction, H7 exploratory), a data-integrity test suite, a full testing framework in Python, figures, and a report.

**Order of work (see `git log`):** the analysis plan (`docs/analysis_plan.md`) was committed *before any data were opened or test was run*; code, results and report were committed afterwards.

## Data (downloaded by `src/acquire.py`, cached in `data/raw/`)

* Our World in Data CO₂ dataset — `https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv`
* World Bank Indicators API v2 — `NY.GDP.PCAP.KD` (primary GDP per capita, constant 2015 US$), `NY.GDP.PCAP.PP.KD` (PPP, sensitivity), `SP.POP.TOTL`, `SP.URB.TOTL.IN.ZS`, `NV.IND.TOTL.ZS`, plus the country-metadata endpoint (income group, region)

No API key needed. `data/raw/manifest.json` records each file's URL, download date (UTC), size and SHA-256.

## Repository structure

```
README.md, requirements.txt      overview; pinned package versions
docs/analysis_plan.md            pre-specification (first commit)
docs/time_log.csv                planned hours (+ blank column for your actual hours)
src/acquire.py                   download + cache (never overwrites raw data)
src/clean.py                     Week 1 cleaning rules, analysis tables, sample flow, cleaning logs
src/tests.py                     hypothesis tests H1-H7, robustness, power, panel sensitivity
src/make_figures.py              11 figures (PNG, 200 dpi, colour-blind-safe)
tests/test_data_integrity.py     25 pytest checks
notebooks/02_hypothesis_testing.ipynb   runs the whole pipeline and displays results
data/raw, data/processed         raw downloads; cleaned tables and logs
results/                         CSV tables of every test, check and sample flow
figures/                         PNG figures
report/report_content.md         the report (Markdown, ready to paste into Word)
report/report_notes.md           private notes: figure map, checks, troubleshooting log, viva Q&A
report/_build/                   template + script that splices real code/tables into the report
```

## How to reproduce

```
python -m venv .venv
.venv\Scripts\activate                 # Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
python src/acquire.py                  # ~1 min; cached on later runs
python src/clean.py
python src/tests.py                    # ~4-5 min (10,000 bootstraps, 100,000 Monte Carlo permutations)
python src/make_figures.py
pytest                                 # 25 integrity tests
jupyter nbconvert --to notebook --execute --inplace notebooks/02_hypothesis_testing.ipynb
python report/_build/build_report.py   # optional: regenerate report_content.md from results
```

Seeds are fixed (42), so repeated runs give identical result tables. Tested with Python 3.13.3 on Windows 11.

## Where to find things

* Every test result: `results/hypothesis_results.csv` (hypothesis, variant, test, n, statistic, raw p, Holm p, BH p, effect size, 95% CI, decision)
* Assumption checks: `results/assumption_checks.csv`; post hoc: `results/h2_posthoc.csv`; H4 slopes/coefficients/VIF; H6 residuals; H7 coefficients; power: `results/power_*.csv`
* Who entered each test: `results/sample_flow.csv`, `data/processed/exclusions_by_test.csv`
* Every dropped/changed data row: `data/processed/cleaning_log_rows.csv`
