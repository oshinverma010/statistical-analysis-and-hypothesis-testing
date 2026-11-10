"""Data-integrity tests. Run with:  pytest -v   (after src/acquire.py and src/clean.py).

Expected counts are regression guards for the cached raw files recorded in data/raw/manifest.json
(download date 2026-11-05). If you deliberately refresh the raw data, re-check and update EXPECTED.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RAW, PROC, RES = ROOT / "data" / "raw", ROOT / "data" / "processed", ROOT / "results"

EXPECTED = {
    "owid_raw_rows": 50_411,
    "owid_rows_1990_2022": 8_375,
    "owid_country_rows": 7_227,
    "owid_noniso_entities": 35,
    "wb_meta_rows": 296,
    "wb_aggregates": 79,
    "wb_countries": 217,
    "panel_rows": 7_557,
    "join_both": 6_831, "join_owid_only": 396, "join_wb_only": 330,
    "cross_section_rows": 229,
    "n_H1": 198, "n_H2": 205, "n_H3": 146, "n_H4": 198, "n_H5": 74, "n_H6": 184, "n_H7": 194,
}


def _load_hyp_module():
    spec = importlib.util.spec_from_file_location("hyp_tests", ROOT / "src" / "tests.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def panel():
    return pd.read_csv(PROC / "panel_1990_2022.csv")


@pytest.fixture(scope="session")
def cs():
    return pd.read_csv(PROC / "cross_section_2022.csv")


@pytest.fixture(scope="session")
def paired():
    return pd.read_csv(PROC / "paired_1990_2022.csv")


@pytest.fixture(scope="session")
def samples():
    return _load_hyp_module().load_samples()


@pytest.fixture(scope="session")
def wb_meta():
    return pd.read_csv(RAW / "wb_country_metadata.csv", keep_default_na=False)


@pytest.fixture(scope="session")
def owid_raw():
    return pd.read_csv(RAW / "owid-co2-data.csv")


# ---------------------------------------------------------------- keys / duplicates
def test_no_duplicate_iso3_year_in_panel(panel):
    assert not panel.duplicated(["iso3", "year"]).any()


def test_cross_section_and_paired_one_row_per_country(cs, paired):
    assert cs["iso3"].is_unique and paired["iso3"].is_unique


def test_years_within_window(panel):
    assert panel["year"].between(1990, 2022).all()


# ---------------------------------------------------------------- aggregates
def test_no_world_bank_aggregates_in_country_tables(panel, cs, paired, wb_meta):
    agg = set(wb_meta.loc[wb_meta["region"] == "Aggregates", "iso3"])
    for df in (panel, cs, paired):
        assert not set(df["iso3"]) & agg


def test_no_owid_aggregates_in_country_tables(panel, cs, owid_raw):
    noiso_names = set(owid_raw.loc[owid_raw["iso_code"].isna(), "country"]) - {"Kosovo"}
    assert len(noiso_names) == EXPECTED["owid_noniso_entities"]
    for df in (panel, cs):
        assert df["iso3"].notna().all()
        assert df["iso3"].str.fullmatch(r"[A-Z]{3}").all()  # ISO3 only; no 'OWID_*' pseudo codes, no names
        assert not df["iso3"].str.startswith("OWID").any()


def test_kosovo_crosswalk_applied(panel):
    k = panel[(panel["iso3"] == "XKX") & panel["co2_pc"].notna()]
    assert len(k) > 0, "Kosovo (XKX) CO2 rows should be joined via the manual crosswalk"
    assert k["gdp_pc"].notna().any(), "Kosovo should also have World Bank GDP after the join"


def test_extremes_retained(cs):
    assert "QAT" in set(cs.loc[cs["co2_pc"].notna(), "iso3"]), "Qatar is a legitimate extreme and must be kept"


# ---------------------------------------------------------------- log transforms
@pytest.mark.parametrize("hyp", ["H1", "H2", "H3", "H4", "H7"])
def test_log_inputs_are_finite_and_cross_section_values_positive(samples, hyp):
    d = samples[hyp]
    assert np.isfinite(d["lg_co2"]).all()
    assert (10 ** d["lg_co2"] > 0).all()
    if hyp in ("H1", "H4", "H7"):
        assert np.isfinite(d["lg_gdp"]).all()


def test_non_positive_values_never_reach_log(samples):
    h5 = samples["H5"]
    assert np.isfinite(h5["log_ratio_co2"]).all()
    assert (h5["co2_pc_1990"] > 0).all() and (h5["co2_pc_2022"] > 0).all()
    h6 = samples["H6"]
    for c in ["co2_pc_1990", "co2_pc_2022", "gdp_pc_1990", "gdp_pc_2022"]:
        assert (h6[c] > 0).all()


# ---------------------------------------------------------------- ranges
def test_shares_within_0_100(panel, cs):
    for col in ["urban_pct", "industry_pct"]:
        v = panel[col].dropna()
        assert v.between(0, 100).all(), col
    for col in ["urban_pct_latest", "industry_pct_latest"]:
        assert cs[col].dropna().between(0, 100).all(), col


def test_covariate_lookback_window_respected(cs):
    for col in ["urban_pct_year", "industry_pct_year"]:
        assert cs[col].dropna().between(2018, 2022).all()


# ---------------------------------------------------------------- no imputation
def test_outcome_and_main_predictor_not_imputed(panel, owid_raw):
    """Non-null counts of CO2 pc and GDP pc in the panel must equal the raw non-null counts (no filling)."""
    raw = owid_raw[(owid_raw["year"].between(1990, 2022)) & owid_raw["iso_code"].notna() | ((owid_raw["country"] == "Kosovo") & owid_raw["year"].between(1990, 2022))]
    assert panel["co2_pc"].notna().sum() == raw["co2_per_capita"].notna().sum()
    wb = pd.read_csv(RAW / "wb_NY.GDP.PCAP.KD.csv", keep_default_na=False)
    wb["value"] = pd.to_numeric(wb["value"].replace("", np.nan))
    meta = pd.read_csv(RAW / "wb_country_metadata.csv", keep_default_na=False)
    wb = wb[wb["iso3"].isin(set(meta.loc[meta["region"] != "Aggregates", "iso3"]))]
    assert panel["gdp_pc"].notna().sum() == wb["value"].notna().sum()


def test_panel_values_equal_raw_for_spot_check(panel, owid_raw):
    usa = owid_raw[(owid_raw["iso_code"] == "USA") & (owid_raw["year"] == 2019)]["co2_per_capita"].iloc[0]
    assert panel[(panel["iso3"] == "USA") & (panel["year"] == 2019)]["co2_pc"].iloc[0] == pytest.approx(usa)


# ---------------------------------------------------------------- counts & join report
def test_expected_row_counts(panel, cs, owid_raw, wb_meta):
    assert len(owid_raw) == EXPECTED["owid_raw_rows"]
    assert len(wb_meta) == EXPECTED["wb_meta_rows"]
    assert (wb_meta["region"] == "Aggregates").sum() == EXPECTED["wb_aggregates"]
    assert (wb_meta["region"] != "Aggregates").sum() == EXPECTED["wb_countries"]
    assert len(panel) == EXPECTED["panel_rows"]
    assert len(cs) == EXPECTED["cross_section_rows"]


def test_join_count_reported_and_consistent(panel):
    steps = pd.read_csv(RES / "cleaning_log.csv")
    note = steps.loc[steps["step"] == "join", "note"].iloc[0]
    for key in ["both", "owid_only", "wb_only"]:
        assert f"{key}=" in note
    exp = EXPECTED["join_both"] + EXPECTED["join_owid_only"] + EXPECTED["join_wb_only"]
    assert exp == len(panel) == EXPECTED["panel_rows"]
    assert f"both={EXPECTED['join_both']:,}" in note
    assert f"owid_only={EXPECTED['join_owid_only']:,}" in note
    assert f"wb_only={EXPECTED['join_wb_only']:,}" in note


def test_sample_sizes_per_hypothesis(samples):
    for h in ["H1", "H2", "H3", "H4", "H5", "H6", "H7"]:
        assert len(samples[h]) == EXPECTED[f"n_{h}"], h


def test_sample_flow_file_matches_samples(samples):
    flow = pd.read_csv(RES / "sample_flow.csv")
    last = flow.groupby("hypothesis").tail(1).set_index("hypothesis")["n_remaining"]
    for h in ["H1", "H2", "H3", "H4", "H5", "H6", "H7"]:
        assert last[h] == len(samples[h]), h


# ---------------------------------------------------------------- group definitions
def test_income_groups_are_the_four_standard_ones(cs):
    assert set(cs["income_group"].dropna()) == {"Low income", "Lower middle income", "Upper middle income", "High income"}


def test_h5_only_high_income(samples):
    assert (samples["H5"]["income_group"].astype(str) == "High income").all()


def test_decoupled_flag_definition(samples):
    d = samples["H6"]
    expect = ((d["gdp_pc_2022"] > d["gdp_pc_1990"]) & (d["co2_pc_2022"] < d["co2_pc_1990"])).astype(float)
    assert (d["decoupled"] == expect).all()


# ---------------------------------------------------------------- reproducibility
def test_raw_files_match_manifest_checksums():
    manifest = json.loads((RAW / "manifest.json").read_text())
    assert len(manifest) == 7
    for name, info in manifest.items():
        h = hashlib.sha256((RAW / name).read_bytes()).hexdigest()
        assert h == info["sha256"], f"{name} changed since download"
