"""Clean and merge OWID + World Bank data (Week 1 rules reused).

Rules
-----
1. Separate non-country aggregates: OWID rows with no ISO3 code, World Bank rows whose region is "Aggregates".
2. Join on ISO3 + year only, never on country names. The only name-based step is a tiny manual crosswalk
   for Kosovo (OWID has no ISO code for it; World Bank uses XKX).
3. Never impute the outcome (CO2 per capita) or the main predictor (GDP per capita).
   Only the exploratory covariates (industry share, urbanisation) may use the latest value in 2018-2022.
4. Keep legitimate extremes (e.g. Qatar); log10 scales are used downstream.
5. Restrict to 1990-2022.
6. Every dropped / changed row is written to data/processed/cleaning_log_rows.csv, step totals to results/cleaning_log.csv.

Outputs (data/processed): panel_1990_2022.csv, cross_section_2022.csv, paired_1990_2022.csv,
exclusions_by_test.csv, cleaning_log_rows.csv ; results: cleaning_log.csv, sample_flow.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"

Y0, Y1 = 1990, 2022
COVARIATE_WINDOW = (2018, 2022)  # allowed look-back for exploratory covariates only
KOSOVO_CROSSWALK = {"Kosovo": "XKX"}  # the single manual name -> ISO3 crosswalk
INCOME_ORDER = ["Low income", "Lower middle income", "Upper middle income", "High income"]
OIL_GAS_EXPORTERS = ["QAT", "KWT", "ARE", "SAU", "BHR", "OMN", "BRN", "TTO", "LBY", "IRQ", "IRN", "DZA",
                     "AZE", "KAZ", "TKM", "NGA", "AGO", "GNQ", "GAB", "COG", "VEN", "RUS", "NOR"]  # fixed in analysis_plan.md
WB_COLS = {
    "NY.GDP.PCAP.KD": "gdp_pc",
    "NY.GDP.PCAP.PP.KD": "gdp_pc_ppp",
    "SP.POP.TOTL": "population_wb",
    "SP.URB.TOTL.IN.ZS": "urban_pct",
    "NV.IND.TOTL.ZS": "industry_pct",
}


class CleaningLog:
    """Collects step-level counts and row-level drop/change records."""

    def __init__(self) -> None:
        self.steps: list[dict] = []
        self.rows: list[pd.DataFrame] = []

    def step(self, step: str, rule: str, before: int, after: int, note: str = "") -> None:
        self.steps.append({"step": step, "rule": rule, "rows_before": before, "rows_after": after,
                           "rows_dropped_or_changed": before - after, "note": note})
        print(f"[{step}] {rule}: {before:,} -> {after:,}  {note}")

    def rows_dropped(self, df: pd.DataFrame, source: str, step: str, reason: str, key_cols: list[str]) -> None:
        if len(df) == 0:
            return
        out = df[key_cols].copy()
        out.columns = ["key1", "key2"] if len(key_cols) == 2 else ["key1"]
        if "key2" not in out:
            out["key2"] = ""
        out.insert(0, "source", source)
        out.insert(1, "step", step)
        out["reason"] = reason
        self.rows.append(out)

    def frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        rows = pd.concat(self.rows, ignore_index=True) if self.rows else pd.DataFrame()
        return pd.DataFrame(self.steps), rows


# --------------------------------------------------------------------------- loading
def load_wb_metadata(log: CleaningLog) -> tuple[pd.DataFrame, set[str]]:
    m = pd.read_csv(RAW / "wb_country_metadata.csv", keep_default_na=False)
    agg = m[m["region"] == "Aggregates"]
    countries = m[m["region"] != "Aggregates"].copy()
    log.step("wb_meta", "drop World Bank aggregates (region == 'Aggregates')", len(m), len(countries))
    log.rows_dropped(agg, "wb_metadata", "wb_meta", "aggregate", ["iso3", "name"])
    countries["income_group"] = pd.Categorical(countries["income_group"], INCOME_ORDER, ordered=True)
    return countries[["iso3", "name", "region", "income_group"]].reset_index(drop=True), set(agg["iso3"])


def load_wb_indicators(countries: pd.DataFrame, agg_codes: set[str], log: CleaningLog) -> pd.DataFrame:
    frames = []
    for code, col in WB_COLS.items():
        d = pd.read_csv(RAW / f"wb_{code}.csv", keep_default_na=False)
        d["value"] = pd.to_numeric(d["value"].replace("", np.nan))
        n0 = len(d)
        blank = d[d["iso3"] == ""]  # aggregates (income-level groups) come back with blank ISO3
        aggs = d[d["iso3"].isin(agg_codes)]
        keep = d[d["iso3"].isin(set(countries["iso3"]))]
        log.rows_dropped(blank, f"wb_{code}", "wb_ind", "aggregate (blank ISO3)", ["wb_name", "year"])
        log.rows_dropped(aggs, f"wb_{code}", "wb_ind", "aggregate (region == Aggregates)", ["iso3", "year"])
        other = d[~d.index.isin(blank.index) & ~d.index.isin(aggs.index) & ~d.index.isin(keep.index)]
        log.rows_dropped(other, f"wb_{code}", "wb_ind", "ISO3 not in country metadata", ["iso3", "year"])
        log.step(f"wb_ind:{code}", "keep country rows only", n0, len(keep),
                 f"{int(keep['value'].isna().sum()):,} null values kept as missing (not imputed)")
        if keep.duplicated(["iso3", "year"]).any():
            raise ValueError(f"duplicate (iso3, year) keys in World Bank {code}")
        frames.append(keep[["iso3", "year", "value"]].rename(columns={"value": col}))
    out = frames[0]
    for f in frames[1:]:
        out = out.merge(f, on=["iso3", "year"], how="outer", validate="one_to_one")
    return out


def load_owid(log: CleaningLog) -> pd.DataFrame:
    o = pd.read_csv(RAW / "owid-co2-data.csv", usecols=["country", "iso_code", "year", "population", "co2", "co2_per_capita"])
    n0 = len(o)
    o = o[(o["year"] >= Y0) & (o["year"] <= Y1)].copy()
    log.step("owid_years", f"restrict to {Y0}-{Y1}", n0, len(o))
    # Kosovo crosswalk (the only name-based rule)
    kos = o["iso_code"].isna() & o["country"].isin(KOSOVO_CROSSWALK)
    log.rows_dropped(o[kos], "owid", "owid_crosswalk", "Kosovo: iso_code NaN -> XKX (changed)", ["country", "year"])
    o.loc[kos, "iso_code"] = o.loc[kos, "country"].map(KOSOVO_CROSSWALK)
    log.step("owid_crosswalk", "manual crosswalk Kosovo -> XKX", int(kos.sum()), 0, "rows changed (not dropped)")
    noiso = o["iso_code"].isna()
    log.rows_dropped(o[noiso], "owid", "owid_agg", "aggregate / non-country (no ISO3)", ["country", "year"])
    o2 = o[~noiso].rename(columns={"iso_code": "iso3", "co2": "co2_total_mt", "co2_per_capita": "co2_pc",
                                   "population": "population_owid"})
    log.step("owid_agg", "drop OWID rows with no ISO3 (aggregates, regions, shipping, aviation...)", len(o), len(o2),
             f"{o.loc[noiso, 'country'].nunique()} non-country entities")
    if o2.duplicated(["iso3", "year"]).any():
        raise ValueError("duplicate (iso3, year) keys in OWID")
    return o2.drop(columns="country").reset_index(drop=True)


# --------------------------------------------------------------------------- build tables
def build_panel(log: CleaningLog) -> tuple[pd.DataFrame, pd.DataFrame]:
    countries, agg_codes = load_wb_metadata(log)
    wb = load_wb_indicators(countries, agg_codes, log)
    owid = load_owid(log)

    merged = owid.merge(wb, on=["iso3", "year"], how="outer", indicator="_join", validate="one_to_one")
    join_counts = merged["_join"].value_counts().to_dict()
    log.step("join", "outer join OWID x World Bank on (iso3, year)", len(owid) + len(wb), len(merged),
             f"both={join_counts.get('both', 0):,}, owid_only={join_counts.get('left_only', 0):,}, "
             f"wb_only={join_counts.get('right_only', 0):,}")
    owid_only = merged[merged["_join"] == "left_only"].drop_duplicates("iso3")
    log.rows_dropped(owid_only.assign(src="owid_only"), "join", "join", "OWID country without World Bank record (no income group / GDP)",
                     ["iso3", "src"])
    # attach metadata (income group, region)
    panel = merged.merge(countries.rename(columns={"name": "wb_name"}), on="iso3", how="left", validate="many_to_one")

    # ---- covariate range checks (shares must be within 0-100) ----
    for col in ["urban_pct", "industry_pct"]:
        bad = panel[col].notna() & ~panel[col].between(0, 100)
        log.rows_dropped(panel[bad], "panel", "range_check", f"{col} outside 0-100 -> set to NaN", ["iso3", "year"])
        log.step(f"range:{col}", "shares outside 0-100 set to NaN", int(panel[col].notna().sum()),
                 int(panel[col].notna().sum() - bad.sum()))
        panel.loc[bad, col] = np.nan

    # ---- non-positive values cannot be log-transformed: flag, do not alter ----
    for col in ["co2_pc", "gdp_pc", "gdp_pc_ppp"]:
        nonpos = panel[col].notna() & (panel[col] <= 0)
        log.rows_dropped(panel[nonpos], "panel", "log_guard", f"{col} <= 0: excluded from log-based analyses (value kept)", ["iso3", "year"])
        log.step(f"logguard:{col}", "non-positive values excluded from log transforms", int(panel[col].notna().sum()),
                 int(panel[col].notna().sum() - nonpos.sum()))

    panel["population"] = panel["population_wb"].combine_first(panel["population_owid"])
    panel["oil_gas_exporter"] = panel["iso3"].isin(OIL_GAS_EXPORTERS)
    panel = panel.drop(columns="_join").sort_values(["iso3", "year"]).reset_index(drop=True)
    return panel, countries


def pick_covariate(panel: pd.DataFrame, col: str) -> pd.DataFrame:
    """Latest observation in COVARIATE_WINDOW per country (exploratory covariates only). Returns iso3, value, year_used."""
    w = panel[panel["year"].between(*COVARIATE_WINDOW) & panel[col].notna()].sort_values("year")
    last = w.groupby("iso3").tail(1)
    return last[["iso3", col, "year"]].rename(columns={col: f"{col}_latest", "year": f"{col}_year"})


def cross_section(panel: pd.DataFrame) -> pd.DataFrame:
    base_cols = ["iso3", "wb_name", "region", "income_group", "population", "oil_gas_exporter",
                 "co2_pc", "co2_total_mt", "gdp_pc", "gdp_pc_ppp", "urban_pct", "industry_pct"]
    cs = panel[panel["year"] == Y1][base_cols].copy()
    # countries known to either source but absent in 2022 rows (e.g. Kosovo has no WB row issues) are rebuilt from the panel's country list
    allc = panel[["iso3", "wb_name", "region", "income_group", "oil_gas_exporter"]].drop_duplicates("iso3")
    cs = allc.merge(cs.drop(columns=["wb_name", "region", "income_group", "oil_gas_exporter"]), on="iso3", how="left")
    for col in ["urban_pct", "industry_pct"]:
        cs = cs.merge(pick_covariate(panel, col), on="iso3", how="left")
    cs["population"] = cs.groupby("iso3")["population"].transform("first")
    # population fallback for the "small country" filter: the latest available population in the window (descriptive attribute, not analysed)
    pop_latest = panel[panel["population"].notna()].sort_values("year").groupby("iso3")["population"].last()
    cs["population"] = cs["population"].fillna(cs["iso3"].map(pop_latest))
    cs["small_pop"] = cs["population"] < 1_000_000
    return cs.sort_values("iso3").reset_index(drop=True)


def paired_table(panel: pd.DataFrame, cs: pd.DataFrame) -> pd.DataFrame:
    cols = ["co2_pc", "gdp_pc", "gdp_pc_ppp"]
    a = panel[panel["year"] == Y0].set_index("iso3")[cols].add_suffix(f"_{Y0}")
    b = panel[panel["year"] == Y1].set_index("iso3")[cols].add_suffix(f"_{Y1}")
    p = cs[["iso3", "wb_name", "region", "income_group", "population", "oil_gas_exporter", "small_pop"]].merge(
        a.reset_index(), on="iso3", how="left").merge(b.reset_index(), on="iso3", how="left")
    pos = lambda c: p[c].where(p[c] > 0)  # noqa: E731  (non-positive -> NaN for ratios)
    p["log_ratio_co2"] = np.log(pos(f"co2_pc_{Y1}") / pos(f"co2_pc_{Y0}"))
    p["log_ratio_gdp"] = np.log(pos(f"gdp_pc_{Y1}") / pos(f"gdp_pc_{Y0}"))
    p["pct_change_co2"] = 100 * (np.exp(p["log_ratio_co2"]) - 1)
    ok = p["log_ratio_co2"].notna() & p["log_ratio_gdp"].notna()
    p["decoupled"] = np.where(ok, ((p["log_ratio_gdp"] > 0) & (p["log_ratio_co2"] < 0)).astype(float), np.nan)
    pos2 = lambda c: p[c].where(p[c] > 0)  # noqa: E731
    ok2 = (pos2(f"gdp_pc_ppp_{Y0}").notna() & pos2(f"gdp_pc_ppp_{Y1}").notna() & p["log_ratio_co2"].notna())
    p["decoupled_ppp"] = np.where(ok2, ((np.log(pos2(f"gdp_pc_ppp_{Y1}") / pos2(f"gdp_pc_ppp_{Y0}")) > 0) &
                                        (p["log_ratio_co2"] < 0)).astype(float), np.nan)
    return p


# --------------------------------------------------------------------------- sample flow per hypothesis
def sample_flow(cs: pd.DataFrame, pr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sequential exclusion table: for each test, which countries drop out first and why."""
    pos = lambda s: s.notna() & (s > 0)  # noqa: E731
    pr = pr.set_index("iso3")
    cs = cs.set_index("iso3")
    spec = {
        "H1": [("no CO2 pc in 2022 (or <= 0)", pos(cs["co2_pc"])), ("no GDP pc in 2022 (or <= 0)", pos(cs["gdp_pc"]))],
        "H2": [("no CO2 pc in 2022 (or <= 0)", pos(cs["co2_pc"])), ("no World Bank income group", cs["income_group"].notna())],
        "H3": [("no CO2 pc in 2022 (or <= 0)", pos(cs["co2_pc"])), ("no World Bank income group", cs["income_group"].notna()),
               ("income group not High / Lower-middle / Low", cs["income_group"].isin(["High income", "Lower middle income", "Low income"]))],
        "H4": [("no CO2 pc in 2022 (or <= 0)", pos(cs["co2_pc"])), ("no GDP pc in 2022 (or <= 0)", pos(cs["gdp_pc"])),
               ("no World Bank income group", cs["income_group"].notna())],
        "H5": [("not a high-income country", cs["income_group"].eq("High income")),
               ("no CO2 pc in 1990 (or <= 0)", pos(pr["co2_pc_1990"]).reindex(cs.index, fill_value=False)),
               ("no CO2 pc in 2022 (or <= 0)", pos(cs["co2_pc"]))],
        "H6": [("no World Bank income group", cs["income_group"].notna()),
               ("no CO2 pc in 1990 (or <= 0)", pos(pr["co2_pc_1990"]).reindex(cs.index, fill_value=False)),
               ("no CO2 pc in 2022 (or <= 0)", pos(cs["co2_pc"])),
               ("no GDP pc in 1990 (or <= 0)", pos(pr["gdp_pc_1990"]).reindex(cs.index, fill_value=False)),
               ("no GDP pc in 2022 (or <= 0)", pos(cs["gdp_pc"]))],
        "H7": [("no CO2 pc in 2022 (or <= 0)", pos(cs["co2_pc"])), ("no GDP pc in 2022 (or <= 0)", pos(cs["gdp_pc"])),
               ("no urbanisation value in 2018-2022", cs["urban_pct_latest"].notna()),
               ("no industry share in 2018-2022", cs["industry_pct_latest"].notna())],
    }
    flow_rows, excl_rows = [], []
    for h, rules in spec.items():
        alive = pd.Series(True, index=cs.index)
        flow_rows.append({"hypothesis": h, "stage": "all country entities (WB country list U OWID countries)", "n_remaining": int(alive.sum()), "n_excluded": 0})
        for reason, ok in rules:
            ok = ok.reindex(cs.index).fillna(False).astype(bool)
            drop = alive & ~ok
            for iso in cs.index[drop]:
                excl_rows.append({"hypothesis": h, "iso3": iso, "name": cs.at[iso, "wb_name"], "reason": reason})
            alive &= ok
            flow_rows.append({"hypothesis": h, "stage": reason, "n_remaining": int(alive.sum()), "n_excluded": int(drop.sum())})
    return pd.DataFrame(flow_rows), pd.DataFrame(excl_rows)


def main() -> None:
    PROC.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    log = CleaningLog()
    panel, _ = build_panel(log)
    cs = cross_section(panel)
    pr = paired_table(panel, cs)
    flow, excl = sample_flow(cs, pr)
    steps, rows = log.frames()

    # names for OWID-only entities (no WB name): fall back to ISO3 in exclusion table
    excl["name"] = excl["name"].fillna(excl["iso3"])
    panel.to_csv(PROC / "panel_1990_2022.csv", index=False)
    cs.to_csv(PROC / "cross_section_2022.csv", index=False)
    pr.to_csv(PROC / "paired_1990_2022.csv", index=False)
    excl.to_csv(PROC / "exclusions_by_test.csv", index=False)
    rows.to_csv(PROC / "cleaning_log_rows.csv", index=False)
    steps.to_csv(RES / "cleaning_log.csv", index=False)
    flow.to_csv(RES / "sample_flow.csv", index=False)
    print(f"\npanel {panel.shape}, cross-section {cs.shape}, paired {pr.shape}, row-level log {len(rows):,} rows")
    print(flow.groupby("hypothesis").tail(1)[["hypothesis", "n_remaining"]].to_string(index=False))


if __name__ == "__main__":
    main()
