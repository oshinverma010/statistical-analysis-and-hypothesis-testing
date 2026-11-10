"""Hypothesis tests H1-H7 exactly as pre-specified in docs/analysis_plan.md.

Run:  python src/tests.py     (writes results/*.csv)

Every number in the report comes from this module. Seeds are fixed (SEED) so results are reproducible.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scikit_posthocs as sp
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.stats.power import TTestIndPower, TTestPower

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
SEED = 42
B = 10_000  # bootstrap / permutation resamples
ALPHA = 0.05
INCOME_ORDER = ["Low income", "Lower middle income", "Upper middle income", "High income"]
HIGH, LMIC_LOW = "High income", ["Lower middle income", "Low income"]
WILCOXON_ARE = 0.955  # asymptotic relative efficiency of Wilcoxon vs t under normality (power approximation)

ASSUMPTIONS: list[dict] = []  # collected assumption checks
RESULTS: list[dict] = []  # collected test results


# =========================================================================== data
def load_samples() -> dict[str, pd.DataFrame]:
    """Analysis-ready samples for each hypothesis (only variables the test needs must be observed)."""
    cs = pd.read_csv(PROC / "cross_section_2022.csv")
    pr = pd.read_csv(PROC / "paired_1990_2022.csv")
    for d in (cs, pr):
        d["income_group"] = pd.Categorical(d["income_group"], INCOME_ORDER, ordered=True)
    pos = lambda s: s.where(s > 0)  # noqa: E731  non-positive -> NaN so log10 is never taken of <= 0
    cs["lg_co2"] = np.log10(pos(cs["co2_pc"]))
    cs["lg_gdp"] = np.log10(pos(cs["gdp_pc"]))
    cs["lg_gdp_ppp"] = np.log10(pos(cs["gdp_pc_ppp"]))
    s = {}
    s["H1"] = cs.dropna(subset=["lg_co2", "lg_gdp"]).reset_index(drop=True)
    s["H2"] = cs.dropna(subset=["lg_co2", "income_group"]).reset_index(drop=True)
    s["H3"] = s["H2"][s["H2"]["income_group"].isin([HIGH] + LMIC_LOW)].reset_index(drop=True)
    s["H4"] = s["H1"].dropna(subset=["income_group"]).reset_index(drop=True)
    p5 = pr[pr["income_group"] == HIGH]
    s["H5"] = p5.dropna(subset=["log_ratio_co2"]).reset_index(drop=True)
    s["H6"] = pr.dropna(subset=["income_group", "decoupled"]).reset_index(drop=True)
    s["H7"] = s["H1"].dropna(subset=["urban_pct_latest", "industry_pct_latest"]).reset_index(drop=True)
    s["cs"], s["pr"] = cs, pr
    return s


def variants(df: pd.DataFrame, ppp_col: str | None = None) -> dict[str, pd.DataFrame]:
    """Robustness subsets pre-specified in the plan."""
    out = {"excluding population < 1 million": df[~df["small_pop"].astype(bool)],
           "excluding oil/gas exporters": df[~df["oil_gas_exporter"].astype(bool)]}
    return out


# =========================================================================== helpers
def note(msg: str) -> None:
    print("   " + msg)


def add_result(h: str, variant: str, test: str, n: int, stat: float, p: float, es_name: str, es: float,
               lo: float, hi: float, primary: bool = False, direction_ok: bool | None = None, comment: str = "") -> None:
    RESULTS.append(dict(hypothesis=h, variant=variant, test=test, n=n, statistic=stat, p_raw=p, effect_size_name=es_name,
                        effect_size=es, ci_low=lo, ci_high=hi, primary=primary, direction_ok=direction_ok, comment=comment))
    note(f"{h} [{variant}] {test}: n={n}, stat={stat:.4g}, p={p:.4g}, {es_name}={es:.3f} [{lo:.3f}, {hi:.3f}] {comment}")


def check(h: str, name: str, stat: float, p: float, verdict: str, detail: str = "") -> None:
    ASSUMPTIONS.append(dict(hypothesis=h, check=name, statistic=stat, p_value=p, verdict=verdict, detail=detail))
    note(f"CHECK {h} {name}: stat={stat:.4g}, p={p:.4g} -> {verdict} {detail}")


def shapiro_check(h: str, label: str, x) -> float:
    x = np.asarray(x, float)
    w, p = stats.shapiro(x)
    check(h, f"Shapiro-Wilk: {label}", w, p, "normal (p>0.05)" if p > ALPHA else "NOT normal (p<=0.05)", f"n={len(x)}")
    return p


def boot_ci(n: int, stat_fn, rng: np.random.Generator, reps: int = B) -> tuple[float, float]:
    """Percentile bootstrap 95% CI. stat_fn receives an index array of resampled rows."""
    vals = np.empty(reps)
    for i in range(reps):
        vals[i] = stat_fn(rng.integers(0, n, n))
    vals = vals[np.isfinite(vals)]
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def cohens_d(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    sp_ = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / sp_)


def rank_biserial_mw(a, b) -> float:
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2 * u / (len(a) * len(b)) - 1)


def holm_adjust(p: np.ndarray) -> np.ndarray:
    return multipletests(p, method="holm")[1]


def decision_text(p_adj: float, h1: str) -> str:
    return f"Reject H0 (supports H1: {h1})" if p_adj < ALPHA else f"Fail to reject H0 (no evidence of: {h1})"


# =========================================================================== H1
def spearman_stat(x, y) -> float:
    return stats.spearmanr(x, y).statistic


def run_h1(S: dict) -> None:
    print("\n=== H1: GDP pc vs CO2 pc, 2022 (Spearman, one-sided rho > 0) ===")
    d = S["H1"]
    shapiro_check("H1", "log10 GDP pc", d["lg_gdp"])
    shapiro_check("H1", "log10 CO2 pc", d["lg_co2"])
    shapiro_check("H1", "raw GDP pc (before log)", d["gdp_pc"])
    shapiro_check("H1", "raw CO2 pc (before log)", d["co2_pc"])
    sk = stats.skew(d["co2_pc"]), stats.skew(d["lg_co2"])
    check("H1", "skewness raw vs log10 CO2 pc", sk[0], float("nan"), "log transform reduces skew" if abs(sk[1]) < abs(sk[0]) else "no gain",
          f"skew raw={sk[0]:.2f}, skew log10={sk[1]:.2f}")
    z = np.abs(stats.zscore(d["lg_co2"]))
    check("H1", "outliers |z|>3 on log10 CO2 pc", float((z > 3).sum()), float("nan"), "none" if (z > 3).sum() == 0 else "present (kept, rank test robust)",
          ", ".join(d.loc[z > 3, "iso3"]))

    def run(df, x, label, primary=False):
        rho_res = stats.spearmanr(df[x], df["lg_co2"], alternative="greater")
        rng = np.random.default_rng(SEED)
        xv, yv = df[x].to_numpy(), df["lg_co2"].to_numpy()
        lo, hi = boot_ci(len(df), lambda i: spearman_stat(xv[i], yv[i]), rng)
        add_result("H1", label, "Spearman rho (one-sided, greater)", len(df), rho_res.statistic, rho_res.pvalue, "Spearman rho",
                   rho_res.statistic, lo, hi, primary, rho_res.statistic > 0)
        return rho_res

    run(d, "lg_gdp", "primary", True)
    # Pearson on log10 + Fisher z CI
    r, p_r = stats.pearsonr(d["lg_gdp"], d["lg_co2"], alternative="greater")
    ci = stats.pearsonr(d["lg_gdp"], d["lg_co2"]).confidence_interval(0.95)
    add_result("H1", "companion: Pearson on log10 values", "Pearson r (one-sided, greater)", len(d), r, p_r, "Pearson r", r, ci.low, ci.high,
               False, r > 0)
    # permutation test
    rng = np.random.default_rng(SEED)
    obs = spearman_stat(d["lg_gdp"], d["lg_co2"])
    ry = d["lg_co2"].to_numpy()
    perm = np.array([spearman_stat(d["lg_gdp"], rng.permutation(ry)) for _ in range(B)])
    p_perm = (np.sum(perm >= obs) + 1) / (B + 1)
    add_result("H1", "robustness: permutation test", "Spearman permutation (B=10,000)", len(d), obs, p_perm, "Spearman rho", obs, np.nan, np.nan,
               False, obs > 0)
    run(d.dropna(subset=["lg_gdp_ppp"]), "lg_gdp_ppp", "robustness: PPP GDP")
    for lab, sub in variants(d).items():
        run(sub, "lg_gdp", f"robustness: {lab}")
    # Spearman rho by quintile is not needed; monotonic check via LOWESS in figures


# =========================================================================== H2
def run_h2(S: dict) -> None:
    print("\n=== H2: CO2 pc differs across four income groups (Kruskal-Wallis) ===")
    d = S["H2"]
    groups = [d.loc[d["income_group"] == g, "co2_pc"].to_numpy() for g in INCOME_ORDER]
    lgroups = [d.loc[d["income_group"] == g, "lg_co2"].to_numpy() for g in INCOME_ORDER]
    pd.DataFrame({"income_group": INCOME_ORDER, "n": [len(g) for g in groups],
                  "median_co2_pc": [np.median(g) for g in groups],
                  "geom_mean_co2_pc": [10 ** np.mean(l) for l in lgroups],
                  "mean_log10": [np.mean(l) for l in lgroups], "sd_log10": [np.std(l, ddof=1) for l in lgroups]}
                 ).to_csv(RES / "h2_group_descriptives.csv", index=False)
    for g, a, l in zip(INCOME_ORDER, groups, lgroups):
        shapiro_check("H2", f"raw CO2 pc, {g}", a)
        shapiro_check("H2", f"log10 CO2 pc, {g}", l)
    lev_w, lev_p = stats.levene(*lgroups, center="median")
    check("H2", "Levene/Brown-Forsythe on log10 CO2 pc", lev_w, lev_p, "equal variances" if lev_p > ALPHA else "UNEQUAL variances -> use Welch ANOVA")
    lev_w2, lev_p2 = stats.levene(*groups, center="median")
    check("H2", "Levene/Brown-Forsythe on raw CO2 pc", lev_w2, lev_p2, "equal variances" if lev_p2 > ALPHA else "UNEQUAL variances")

    n = len(d)

    def eps2(idx):
        gs = [d["co2_pc"].to_numpy()[idx][d["income_group"].to_numpy()[idx] == g] for g in INCOME_ORDER]
        if min(len(g) for g in gs) < 2:
            return np.nan
        return stats.kruskal(*gs).statistic / (len(idx) - 1)

    def run_kw(df, label, primary=False):
        gs = [df.loc[df["income_group"] == g, "co2_pc"].to_numpy() for g in INCOME_ORDER]
        H, p = stats.kruskal(*gs)
        e = H / (len(df) - 1)
        arr_g, arr_v = df["income_group"].to_numpy(), df["co2_pc"].to_numpy()

        def f(idx):
            gg = [arr_v[idx][arr_g[idx] == g] for g in INCOME_ORDER]
            return stats.kruskal(*gg).statistic / (len(idx) - 1) if min(len(x) for x in gg) > 1 else np.nan

        lo, hi = boot_ci(len(df), f, np.random.default_rng(SEED))
        add_result("H2", label, "Kruskal-Wallis H (df=3)", len(df), H, p, "epsilon-squared", e, lo, hi, primary, True)
        return H, p

    run_kw(d, "primary", True)

    # parametric comparison
    f_, p_f = stats.f_oneway(*lgroups)
    ss_between = sum(len(l) * (l.mean() - d["lg_co2"].mean()) ** 2 for l in lgroups)
    eta2 = ss_between / ((d["lg_co2"] - d["lg_co2"].mean()) ** 2).sum()
    add_result("H2", "comparison: one-way ANOVA on log10", "ANOVA F (classical)", n, f_, p_f, "eta-squared", eta2, np.nan, np.nan, False, True)
    import pingouin as pg
    w = pg.welch_anova(dv="lg_co2", between="income_group", data=d.assign(income_group=d["income_group"].astype(str)))
    add_result("H2", "comparison: Welch ANOVA on log10", "Welch ANOVA F", n, float(w["F"].iloc[0]), float(w["p_unc"].iloc[0]), "partial eta-squared",
               float(w["np2"].iloc[0]), np.nan, np.nan, False, True)

    # permutation Kruskal-Wallis
    rng = np.random.default_rng(SEED)
    obs = stats.kruskal(*groups).statistic
    gv = d["income_group"].to_numpy()
    vals = d["co2_pc"].to_numpy()
    cnt = 0
    for _ in range(B):
        perm = rng.permutation(gv)
        cnt += stats.kruskal(*[vals[perm == g] for g in INCOME_ORDER]).statistic >= obs
    add_result("H2", "robustness: permutation Kruskal-Wallis", "KW permutation (B=10,000)", n, obs, (cnt + 1) / (B + 1), "epsilon-squared",
               obs / (n - 1), np.nan, np.nan, False, True)
    for lab, sub in variants(d).items():
        run_kw(sub, f"robustness: {lab}")

    # post hoc
    dunn = sp.posthoc_dunn(d.assign(income_group=d["income_group"].astype(str)), val_col="co2_pc", group_col="income_group", p_adjust="holm")
    rows = []
    for i, a in enumerate(INCOME_ORDER):
        for b in INCOME_ORDER[i + 1:]:
            xa, xb = d.loc[d["income_group"] == a, "co2_pc"], d.loc[d["income_group"] == b, "co2_pc"]
            u = stats.mannwhitneyu(xa, xb, alternative="two-sided")
            rows.append(dict(group_a=a, group_b=b, n_a=len(xa), n_b=len(xb), median_a=xa.median(), median_b=xb.median(),
                             dunn_p_holm=dunn.loc[a, b], mw_p_raw=u.pvalue, rank_biserial=2 * u.statistic / (len(xa) * len(xb)) - 1))
    ph = pd.DataFrame(rows)
    ph["mw_p_holm"] = holm_adjust(ph["mw_p_raw"].to_numpy())
    ph.to_csv(RES / "h2_posthoc.csv", index=False)
    note("post hoc (Dunn-Holm / MW-Holm):\n" + ph.round(4).to_string(index=False))


# =========================================================================== H3
def run_h3(S: dict) -> None:
    print("\n=== H3: High income vs Lower-middle + Low (Welch t on log10 CO2 pc, one-sided) ===")
    d = S["H3"]
    a = d.loc[d["income_group"] == HIGH, "lg_co2"].to_numpy()
    b = d.loc[d["income_group"].isin(LMIC_LOW), "lg_co2"].to_numpy()
    shapiro_check("H3", "log10 CO2 pc, High", a)
    shapiro_check("H3", "log10 CO2 pc, Lower-middle + Low", b)
    lw, lp = stats.levene(a, b, center="median")
    check("H3", "Levene (Brown-Forsythe) High vs LM+Low", lw, lp, "equal variances" if lp > ALPHA else "UNEQUAL variances -> Welch justified",
          f"SD ratio={a.std(ddof=1) / b.std(ddof=1):.2f}")
    note(f"n_High={len(a)}, n_LM+Low={len(b)}")

    def welch(a, b, label, primary=False, raw=None):
        t = stats.ttest_ind(a, b, equal_var=False, alternative="greater")
        diff = a.mean() - b.mean()
        va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
        se = np.sqrt(va + vb)
        df_ = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
        tc = stats.t.ppf(0.975, df_)
        rng = np.random.default_rng(SEED)
        bs = [cohens_d(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(B)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        add_result("H3", label, f"Welch t-test (df={df_:.1f}, one-sided greater)", len(a) + len(b), t.statistic, t.pvalue, "Cohen's d", cohens_d(a, b),
                   lo, hi, primary, diff > 0,
                   f"mean diff (log10)={diff:.3f} [{diff - tc * se:.3f}, {diff + tc * se:.3f}]; geometric-mean ratio={10 ** diff:.2f}x")
        return diff

    welch(a, b, "primary", True)
    t_s = stats.ttest_ind(a, b, equal_var=True, alternative="greater")
    add_result("H3", "comparison: Student t (equal variance)", "Student t (one-sided greater)", len(a) + len(b), t_s.statistic, t_s.pvalue,
               "Cohen's d", cohens_d(a, b), np.nan, np.nan, False, True)
    mw = stats.mannwhitneyu(a, b, alternative="greater")
    rng = np.random.default_rng(SEED)
    bs = [2 * stats.mannwhitneyu(rng.choice(a, len(a)), rng.choice(b, len(b)), alternative="two-sided").statistic / (len(a) * len(b)) - 1
          for _ in range(2000)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    add_result("H3", "robustness: Mann-Whitney U", "Mann-Whitney U (one-sided greater)", len(a) + len(b), mw.statistic, mw.pvalue,
               "rank-biserial r", 2 * mw.statistic / (len(a) * len(b)) - 1, lo, hi, False, True)
    # permutation of group labels
    rng = np.random.default_rng(SEED)
    allv = np.concatenate([a, b])
    obs = a.mean() - b.mean()
    cnt = 0
    for _ in range(B):
        pm = rng.permutation(allv)
        cnt += (pm[:len(a)].mean() - pm[len(a):].mean()) >= obs
    add_result("H3", "robustness: permutation of mean difference", "permutation test (B=10,000)", len(allv), obs, (cnt + 1) / (B + 1),
               "mean diff (log10)", obs, np.nan, np.nan, False, True)
    for lab, sub in variants(d).items():
        sa = sub.loc[sub["income_group"] == HIGH, "lg_co2"].to_numpy()
        sb = sub.loc[sub["income_group"].isin(LMIC_LOW), "lg_co2"].to_numpy()
        welch(sa, sb, f"robustness: {lab}")


# =========================================================================== H4
def _fit_full(df, xcol="lg_gdp"):
    return ols(f"lg_co2 ~ {xcol} * C(income_group, Treatment('{HIGH}'))", data=df.assign(income_group=df["income_group"].astype(str))).fit()


def _interaction_F(y, X_full, X_red):
    sse_f = float(((y - X_full @ np.linalg.lstsq(X_full, y, rcond=None)[0]) ** 2).sum())
    sse_r = float(((y - X_red @ np.linalg.lstsq(X_red, y, rcond=None)[0]) ** 2).sum())
    q = X_full.shape[1] - X_red.shape[1]
    return ((sse_r - sse_f) / q) / (sse_f / (len(y) - X_full.shape[1]))


def run_h4(S: dict) -> None:
    print("\n=== H4: slope heterogeneity (interaction F-test, HC3) ===")
    d = S["H4"].copy()
    d["income_group"] = d["income_group"].astype(str)
    grp = "C(income_group, Treatment('High income'))"

    def analyse(df, xcol, label, primary=False, diagnostics=False):
        df = df.copy()
        df["income_group"] = df["income_group"].astype(str)
        m_full = ols(f"lg_co2 ~ {xcol} * {grp}", data=df).fit()
        m_main = ols(f"lg_co2 ~ {xcol} + {grp}", data=df).fit()
        classical = anova_lm(m_full, typ=2)
        robust = anova_lm(m_full, typ=2, robust="hc3")
        inter = [i for i in classical.index if ":" in i][0]
        F_c, p_c = classical.loc[inter, "F"], classical.loc[inter, "PR(>F)"]
        F_r, p_r = robust.loc[inter, "F"], robust.loc[inter, "PR(>F)"]
        dr2 = m_full.rsquared - m_main.rsquared
        f2 = dr2 / (1 - m_full.rsquared)
        y = df["lg_co2"].to_numpy()

        def f(idx):
            sub = df.iloc[idx]
            try:
                a = ols(f"lg_co2 ~ {xcol} * {grp}", data=sub).fit().rsquared
                b = ols(f"lg_co2 ~ {xcol} + {grp}", data=sub).fit().rsquared
                return a - b
            except Exception:
                return np.nan

        lo, hi = boot_ci(len(df), f, np.random.default_rng(SEED), reps=2000)
        add_result("H4", label, "Interaction F-test, HC3 Wald-type (anova_lm robust='hc3')", len(df), F_r, p_r, "delta R-squared (interaction block)", dr2, lo, hi,
                   primary, True, f"classical F={F_c:.3f}, p={p_c:.4g}; Cohen f2={f2:.3f}; R2 full={m_full.rsquared:.3f}")
        return m_full, m_main, F_r

    m_full, m_main, F_r = analyse(d, "lg_gdp", "primary", True)
    coef = pd.DataFrame({"coef": m_full.params, "se_HC3": m_full.get_robustcov_results("HC3").bse,
                         "p_HC3": m_full.get_robustcov_results("HC3").pvalues})
    ci = m_full.get_robustcov_results("HC3").conf_int()
    coef["ci_low"], coef["ci_high"] = ci[:, 0], ci[:, 1]
    coef.round(5).to_csv(RES / "h4_coefficients_hc3.csv")
    # per-group slopes
    rows = []
    for g in INCOME_ORDER:
        sub = d[d["income_group"] == g]
        mg = ols("lg_co2 ~ lg_gdp", data=sub).fit(cov_type="HC3")
        c = mg.conf_int().loc["lg_gdp"]
        rows.append(dict(income_group=g, n=len(sub), slope=mg.params["lg_gdp"], se_HC3=mg.bse["lg_gdp"], ci_low=c[0], ci_high=c[1],
                         r2=mg.rsquared, intercept=mg.params["Intercept"]))
    pd.DataFrame(rows).to_csv(RES / "h4_slopes_by_group.csv", index=False)
    note("per-group slopes:\n" + pd.DataFrame(rows).round(3).to_string(index=False))

    # diagnostics
    bp = het_breuschpagan(m_full.resid, m_full.model.exog)
    check("H4", "Breusch-Pagan (full model)", bp[0], bp[1], "homoscedastic" if bp[1] > ALPHA else "HETEROSCEDASTIC -> HC3 SEs justified")
    shapiro_check("H4", "residuals of full model", m_full.resid)
    cd = m_full.get_influence().cooks_distance[0]
    thr = 4 / len(d)
    flagged = d.loc[cd > thr, "iso3"].tolist()
    check("H4", "Cook's distance > 4/n", float((cd > thr).sum()), thr, f"{len(flagged)} influential points", ", ".join(flagged))
    X = pd.DataFrame(m_main.model.exog, columns=m_main.model.exog_names)
    vifs = {c: variance_inflation_factor(X.values, i) for i, c in enumerate(X.columns) if c != "Intercept"}
    pd.Series(vifs, name="VIF").round(3).to_csv(RES / "h4_vif_main_effects.csv")
    check("H4", "max VIF (main-effects model)", max(vifs.values()), float("nan"), "OK (<5)" if max(vifs.values()) < 5 else "elevated (>=5)",
          "interaction-model VIFs are inflated by construction and not used")
    # robustness
    analyse(d.dropna(subset=["lg_gdp_ppp"]), "lg_gdp_ppp", "robustness: PPP GDP")
    for lab, sub in variants(d).items():
        analyse(sub, "lg_gdp", f"robustness: {lab}")
    analyse(d.loc[cd <= thr], "lg_gdp", "robustness: Cook's-distance outliers removed")
    # permutation tests
    rng = np.random.default_rng(SEED)
    y = d["lg_co2"].to_numpy()
    x = d["lg_gdp"].to_numpy()
    codes = pd.get_dummies(d["income_group"], drop_first=True).to_numpy(float)
    ones = np.ones(len(d))

    def design(gcodes):
        return np.column_stack([ones, x, gcodes]), np.column_stack([ones, x, gcodes, gcodes * x[:, None]])

    Xr, Xf = design(codes)
    obsF = _interaction_F(y, Xf, Xr)
    cnt = 0
    for _ in range(B):
        perm = codes[rng.permutation(len(d))]
        Xr_p, Xf_p = design(perm)
        cnt += _interaction_F(y, Xf_p, Xr_p) >= obsF
    add_result("H4", "robustness: income-group label permutation (as planned)", "permutation of group labels (B=10,000)", len(d), obsF, (cnt + 1) / (B + 1),
               "delta R-squared", m_full.rsquared - m_main.rsquared, np.nan, np.nan, False, True,
               "tests stronger null (group irrelevant altogether)")
    # Freedman-Lane (added): permute residuals of the reduced model -> tests the interaction null only
    beta_r = np.linalg.lstsq(Xr, y, rcond=None)[0]
    fit_r, res_r = Xr @ beta_r, y - Xr @ beta_r
    cnt = 0
    for _ in range(B):
        ys = fit_r + rng.permutation(res_r)
        cnt += _interaction_F(ys, Xf, Xr) >= obsF
    add_result("H4", "robustness: Freedman-Lane residual permutation (added)", "Freedman-Lane permutation (B=10,000)", len(d), obsF, (cnt + 1) / (B + 1),
               "delta R-squared", m_full.rsquared - m_main.rsquared, np.nan, np.nan, False, True, "added after planning; tests the interaction null exactly")


# =========================================================================== H5
def run_h5(S: dict) -> None:
    print("\n=== H5: high-income countries, CO2 pc 1990 -> 2022 (Wilcoxon signed-rank, one-sided less) ===")
    d = S["H5"]
    shapiro_p = shapiro_check("H5", "log ratio ln(CO2pc_2022/CO2pc_1990)", d["log_ratio_co2"])
    note(f"skewness of log ratio: {stats.skew(d['log_ratio_co2']):.3f}; n zero differences: {(d['log_ratio_co2'] == 0).sum()}")

    def run(df, label, primary=False):
        x = df["log_ratio_co2"].to_numpy()
        w = stats.wilcoxon(x, alternative="less")
        ranks = stats.rankdata(np.abs(x))
        rbc = (ranks[x > 0].sum() - ranks[x < 0].sum()) / ranks.sum()  # matched-pairs rank-biserial (negative = decreases)
        rng = np.random.default_rng(SEED)

        def f(idx):
            xx = x[idx]
            xx = xx[xx != 0]
            if len(xx) < 3:
                return np.nan
            rk = stats.rankdata(np.abs(xx))
            return (rk[xx > 0].sum() - rk[xx < 0].sum()) / rk.sum()

        lo, hi = boot_ci(len(x), f, rng)
        med_pct = 100 * (np.exp(np.median(x)) - 1)
        mlo, mhi = boot_ci(len(x), lambda i: 100 * (np.exp(np.median(x[i])) - 1), np.random.default_rng(SEED))
        add_result("H5", label, "Wilcoxon signed-rank (one-sided, less)", len(x), w.statistic, w.pvalue, "matched-pairs rank-biserial r", rbc, lo, hi,
                   primary, rbc < 0, f"median change={med_pct:.1f}% [{mlo:.1f}, {mhi:.1f}]; {(x < 0).sum()} of {len(x)} countries decreased")
        return x

    x = run(d, "primary", True)
    t = stats.ttest_1samp(x, 0, alternative="less")
    dz = x.mean() / x.std(ddof=1)
    lo, hi = boot_ci(len(x), lambda i: x[i].mean() / x[i].std(ddof=1), np.random.default_rng(SEED))
    add_result("H5", "companion: paired t-test on log ratios" + ("" if shapiro_p > ALPHA else " (normality VIOLATED - shown for transparency only)"),
               "paired t-test (one-sided, less)", len(x), t.statistic, t.pvalue, "Cohen's d_z", dz, lo, hi, False, dz < 0)
    nneg = int((x < 0).sum())
    bt = stats.binomtest(nneg, len(x), 0.5, alternative="greater")
    bt_ci = stats.binomtest(nneg, len(x), 0.5, alternative="two-sided").proportion_ci(0.95)  # two-sided CI for the proportion
    add_result("H5", "robustness: sign test", "sign test (one-sided)", len(x), nneg, bt.pvalue, "proportion decreasing", nneg / len(x),
               bt_ci.low, bt_ci.high, False, nneg / len(x) > 0.5)
    from scipy.stats.mstats import winsorize
    xw = np.asarray(winsorize(x, limits=(0.10, 0.10)))
    tw = stats.ttest_1samp(xw, 0, alternative="less")
    dzw = xw.mean() / xw.std(ddof=1)
    add_result("H5", "robustness: 10%-winsorised paired t-test", "paired t-test on winsorised log ratios", len(xw), tw.statistic, tw.pvalue,
               "Cohen's d_z", dzw, np.nan, np.nan, False, dzw < 0)
    for lab, sub in variants(d).items():
        run(sub, f"robustness: {lab}")


# =========================================================================== H6
def _chi2_from_table(tab: np.ndarray) -> tuple[float, float]:
    tab = tab[tab.sum(1) > 0][:, tab.sum(0) > 0]
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return np.nan, np.nan
    exp = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / tab.sum()
    chi2 = ((tab - exp) ** 2 / exp).sum()
    return float(chi2), float(np.sqrt(chi2 / (tab.sum() * (min(tab.shape) - 1))))


def run_h6(S: dict) -> None:
    print("\n=== H6: absolute decoupling independent of income group (chi-square / Monte Carlo) ===")
    d = S["H6"]

    def table(df, col="decoupled"):
        return pd.crosstab(df["income_group"].astype(str), df[col].astype(int)).reindex(INCOME_ORDER).fillna(0).rename(columns={0: "Not decoupled", 1: "Decoupled"})

    def analyse(df, col, label, primary=False):
        t = table(df, col)
        for c in ["Not decoupled", "Decoupled"]:
            if c not in t:
                t[c] = 0
        chi2, p_chi, dof, exp = stats.chi2_contingency(t.to_numpy(), correction=False)
        min_exp = exp.min()
        n_lt5 = int((exp < 5).sum())
        mc_p = np.nan
        if min_exp < 5 or primary:
            rng = np.random.default_rng(SEED)
            g = pd.Categorical(df["income_group"].astype(str), INCOME_ORDER).codes
            y = df[col].astype(int).to_numpy()
            obs, _ = _chi2_from_table(t.to_numpy())
            cnt = 0
            for _ in range(100_000):
                yy = rng.permutation(y)
                tab = np.bincount(g * 2 + yy, minlength=8).reshape(4, 2)
                cnt += _chi2_from_table(tab)[0] >= obs - 1e-9
            mc_p = (cnt + 1) / 100_001
        V = np.sqrt(chi2 / (t.to_numpy().sum() * 1))
        gcodes = pd.Categorical(df["income_group"].astype(str), INCOME_ORDER).codes
        yy = df[col].astype(int).to_numpy()

        def fV(idx):
            tab = np.bincount(gcodes[idx] * 2 + yy[idx], minlength=8).reshape(4, 2)
            return _chi2_from_table(tab)[1]

        lo, hi = boot_ci(len(df), fV, np.random.default_rng(SEED))
        use_mc = min_exp < 5
        p_used = mc_p if use_mc else p_chi
        add_result("H6", label, f"Chi-square test of independence (df={dof})" + (" -> Monte Carlo p used" if use_mc else ""), len(df), chi2, p_used,
                   "Cramer's V", V, lo, hi, primary, True,
                   f"chi-square p={p_chi:.4g}; Monte Carlo p={mc_p:.4g}; min expected={min_exp:.2f}; cells with expected<5: {n_lt5} of {exp.size}")
        return t, exp

    t, exp = analyse(d, "decoupled", "primary", True)
    exp_df = pd.DataFrame(exp, index=t.index, columns=t.columns)
    check("H6", "minimum expected count", float(exp.min()), float("nan"), "all >= 5 -> chi-square valid" if exp.min() >= 5 else "some < 5 -> Monte Carlo (Fisher-Freeman-Halton analogue) is primary",
          f"{int((exp < 5).sum())} of {exp.size} cells below 5")
    n = t.to_numpy().sum()
    rowp, colp = t.sum(axis=1).to_numpy()[:, None] / n, t.sum(axis=0).to_numpy()[None, :] / n
    adj = (t.to_numpy() - exp) / np.sqrt(exp * (1 - rowp) * (1 - colp))
    long = []
    for i, g in enumerate(t.index):
        for j, c in enumerate(t.columns):
            long.append(dict(income_group=g, outcome=c, observed=int(t.iat[i, j]), expected=exp[i, j], std_adj_residual=adj[i, j],
                             row_pct_decoupled=100 * t.iat[i, 1] / t.iloc[i].sum() if c == "Decoupled" else np.nan))
    pd.DataFrame(long).to_csv(RES / "h6_contingency_and_residuals.csv", index=False)
    note("contingency table:\n" + t.to_string() + "\nexpected:\n" + exp_df.round(2).to_string() + "\nstd. adj. residuals:\n" +
         pd.DataFrame(adj, index=t.index, columns=t.columns).round(2).to_string())
    # 2x2 Fisher: High vs rest
    hi_mask = d["income_group"].astype(str) == HIGH
    t2 = np.array([[(hi_mask & (d["decoupled"] == 1)).sum(), (hi_mask & (d["decoupled"] == 0)).sum()],
                   [(~hi_mask & (d["decoupled"] == 1)).sum(), (~hi_mask & (d["decoupled"] == 0)).sum()]])
    orr, pf = stats.fisher_exact(t2)
    lor = np.log(orr) if orr not in (0, np.inf) else np.nan
    se = np.sqrt((1 / (t2 + 0.5)).sum())
    add_result("H6", "check: 2x2 Fisher exact (High vs not High)", "Fisher exact (two-sided)", int(t2.sum()), orr, pf, "odds ratio", orr,
               float(np.exp(np.log((t2[0, 0] + .5) * (t2[1, 1] + .5) / ((t2[0, 1] + .5) * (t2[1, 0] + .5))) - 1.96 * se)),
               float(np.exp(np.log((t2[0, 0] + .5) * (t2[1, 1] + .5) / ((t2[0, 1] + .5) * (t2[1, 0] + .5))) + 1.96 * se)), False, True,
               "CI: Haldane-corrected log-OR (Woolf)")
    dd = d.dropna(subset=["decoupled_ppp"])
    analyse(dd, "decoupled_ppp", "robustness: PPP-based GDP condition")
    for lab, sub in variants(d).items():
        analyse(sub, "decoupled", f"robustness: {lab}")


# =========================================================================== H7 (exploratory)
def run_h7(S: dict) -> None:
    print("\n=== H7 (EXPLORATORY): log CO2 pc ~ log GDP pc + industry share + urbanisation (HC3) ===")
    d = S["H7"]

    def fit(df, xcol, label):
        f = f"lg_co2 ~ {xcol} + industry_pct_latest + urban_pct_latest"
        m = ols(f, data=df).fit(cov_type="HC3")
        mc = ols(f, data=df).fit()
        rows = []
        for v in ["industry_pct_latest", "urban_pct_latest"]:
            red = ols(f.replace(f" + {v}", ""), data=df).fit()
            pr2 = (red.ssr - mc.ssr) / red.ssr
            rows.append(dict(variant=label, term=v, coef=m.params[v], se_HC3=m.bse[v], p_HC3=m.pvalues[v], ci_low=m.conf_int().loc[v, 0],
                             ci_high=m.conf_int().loc[v, 1], partial_r2=pr2, n=len(df)))
        rows.append(dict(variant=label, term=xcol, coef=m.params[xcol], se_HC3=m.bse[xcol], p_HC3=m.pvalues[xcol], ci_low=m.conf_int().loc[xcol, 0],
                         ci_high=m.conf_int().loc[xcol, 1], partial_r2=(ols(f"lg_co2 ~ industry_pct_latest + urban_pct_latest", data=df).fit().ssr - mc.ssr) /
                         ols(f"lg_co2 ~ industry_pct_latest + urban_pct_latest", data=df).fit().ssr, n=len(df)))
        wald = m.f_test("industry_pct_latest = 0, urban_pct_latest = 0")
        r2_base = ols(f"lg_co2 ~ {xcol}", data=df).fit().rsquared
        return rows, float(wald.fvalue), float(wald.pvalue), mc.rsquared, r2_base, m, mc

    rows, F, p, r2, r2b, m, mc = fit(d, "lg_gdp", "primary (exploratory)")
    X = pd.DataFrame(mc.model.exog, columns=mc.model.exog_names)
    vif = {c: variance_inflation_factor(X.values, i) for i, c in enumerate(X.columns) if c != "Intercept"}
    check("H7", "max VIF", max(vif.values()), float("nan"), "OK (<5)" if max(vif.values()) < 5 else "elevated", str({k: round(v, 2) for k, v in vif.items()}))
    bp = het_breuschpagan(mc.resid, mc.model.exog)
    check("H7", "Breusch-Pagan", bp[0], bp[1], "homoscedastic" if bp[1] > ALPHA else "HETEROSCEDASTIC -> HC3 used")
    shapiro_check("H7", "residuals", mc.resid)
    allrows = list(rows)
    ppp = d.dropna(subset=["lg_gdp_ppp"])
    r2_, F2, p2, *_ = fit(ppp, "lg_gdp_ppp", "robustness: PPP GDP")
    allrows += r2_
    yr = d[(d["industry_pct_year"] == 2022) & (d["urban_pct_year"] == 2022)]
    r3_, F3, p3, *_ = fit(yr, "lg_gdp", "robustness: 2022-only covariates (no look-back)")
    allrows += r3_
    pd.DataFrame(allrows).to_csv(RES / "h7_coefficients.csv", index=False)
    note("H7 coefficients:\n" + pd.DataFrame(allrows).round(4).to_string(index=False))
    # joint test row (exploratory: NOT in Holm family)
    def dr2(idx):
        sub = d.iloc[idx]
        return (ols("lg_co2 ~ lg_gdp + industry_pct_latest + urban_pct_latest", data=sub).fit().rsquared
                - ols("lg_co2 ~ lg_gdp", data=sub).fit().rsquared)

    lo7, hi7 = boot_ci(len(d), dr2, np.random.default_rng(SEED), reps=2000)
    add_result("H7", "primary (exploratory)", "Joint HC3 Wald F: industry share & urbanisation (exploratory)", len(d), F, p, "increment in R-squared",
               r2 - r2b, lo7, hi7, False, None, f"R2 full={r2:.3f} vs GDP only={r2b:.3f}; EXPLORATORY, unadjusted p")
    add_result("H7", "robustness: PPP GDP", "Joint HC3 Wald F (exploratory)", len(ppp), F2, p2, "see h7_coefficients.csv", np.nan, np.nan, np.nan, False, None, "exploratory")
    add_result("H7", "robustness: 2022-only covariates", "Joint HC3 Wald F (exploratory)", len(yr), F3, p3, "see h7_coefficients.csv", np.nan, np.nan, np.nan, False, None, "exploratory")


# =========================================================================== power
def run_power(S: dict) -> None:
    print("\n=== Power analysis (H3, H5) ===")
    n_high = int((S["H3"]["income_group"] == HIGH).sum())
    n_low = len(S["H3"]) - n_high
    n5 = len(S["H5"])
    rows = []
    mde3 = TTestIndPower().solve_power(effect_size=None, nobs1=n_high, alpha=ALPHA, power=0.8, ratio=n_low / n_high, alternative="larger")
    mde5 = TTestPower().solve_power(effect_size=None, nobs=n5, alpha=ALPHA, power=0.8, alternative="larger")
    mde5w = TTestPower().solve_power(effect_size=None, nobs=n5 * WILCOXON_ARE, alpha=ALPHA, power=0.8, alternative="larger")
    rows += [dict(hypothesis="H3", n=f"{n_high} vs {n_low}", test="Welch/independent t (power via pooled-variance t, one-sided)", mde_d_80pct_power=mde3),
             dict(hypothesis="H5", n=str(n5), test="paired t (one-sided)", mde_d_80pct_power=mde5),
             dict(hypothesis="H5", n=f"{n5} (effective {n5 * WILCOXON_ARE:.1f})", test="Wilcoxon (t-power with ARE=0.955 approximation)", mde_d_80pct_power=mde5w)]
    for e in (0.2, 0.5, 0.8):
        rows.append(dict(hypothesis="H3", n=f"{n_high} vs {n_low}", test=f"power at d={e}", mde_d_80pct_power=np.nan,
                         power=TTestIndPower().power(e, nobs1=n_high, alpha=ALPHA, ratio=n_low / n_high, alternative="larger")))
        rows.append(dict(hypothesis="H5", n=str(n5), test=f"power at d_z={e} (paired t)", mde_d_80pct_power=np.nan,
                         power=TTestPower().power(e, nobs=n5, alpha=ALPHA, alternative="larger")))
    note(f"H3 MDE d={mde3:.3f}; H5 MDE d_z={mde5:.3f} (paired t), {mde5w:.3f} (Wilcoxon approx.)")
    pd.DataFrame(rows).to_csv(RES / "power_mde.csv", index=False)
    grid = np.linspace(0.05, 1.6, 120)
    curve = pd.DataFrame({"effect_size": grid,
                          "power_H3": [TTestIndPower().power(e, nobs1=n_high, alpha=ALPHA, ratio=n_low / n_high, alternative="larger") for e in grid],
                          "power_H5": [TTestPower().power(e, nobs=n5, alpha=ALPHA, alternative="larger") for e in grid],
                          "power_H5_wilcoxon_approx": [TTestPower().power(e, nobs=n5 * WILCOXON_ARE, alpha=ALPHA, alternative="larger") for e in grid]})
    curve.to_csv(RES / "power_curve.csv", index=False)


# =========================================================================== panel sensitivity (optional)
def run_panel(S: dict) -> None:
    print("\n=== Sensitivity: pooled country-year panel, year effects, country-clustered SEs ===")
    p = pd.read_csv(PROC / "panel_1990_2022.csv")
    p = p[(p["co2_pc"] > 0) & (p["gdp_pc"] > 0) & p["income_group"].notna()].copy()
    p["lg_co2"], p["lg_gdp"] = np.log10(p["co2_pc"]), np.log10(p["gdp_pc"])
    m = smf.ols("lg_co2 ~ lg_gdp + C(year)", data=p).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(p["iso3"])[0]})
    ci = m.conf_int().loc["lg_gdp"]
    mi = smf.ols(f"lg_co2 ~ lg_gdp * C(income_group, Treatment('{HIGH}')) + C(year)", data=p).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(p["iso3"])[0]})
    names = [n for n in mi.params.index if "lg_gdp:" in n]
    R = np.zeros((len(names), len(mi.params)))
    for i, nm in enumerate(names):
        R[i, list(mi.params.index).index(nm)] = 1
    w = mi.wald_test(R, scalar=True)
    pd.DataFrame([dict(model="log10 CO2 pc ~ log10 GDP pc + year effects (cluster SE)", n_country_years=len(p), n_countries=p["iso3"].nunique(),
                       slope=m.params["lg_gdp"], se=m.bse["lg_gdp"], ci_low=ci[0], ci_high=ci[1], p=m.pvalues["lg_gdp"]),
                  dict(model="income-group x slope interaction (joint Wald, cluster SE)", n_country_years=len(p), n_countries=p["iso3"].nunique(),
                       slope=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan, p=float(w.pvalue))]).to_csv(RES / "panel_sensitivity.csv", index=False)
    note(f"pooled slope={m.params['lg_gdp']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] (n={len(p)} country-years, {p['iso3'].nunique()} countries); interaction Wald p={float(w.pvalue):.3g}")


# =========================================================================== driver
def finalize() -> pd.DataFrame:
    df = pd.DataFrame(RESULTS)
    prim = df[(df["primary"]) & df["hypothesis"].isin(["H1", "H2", "H3", "H4", "H5", "H6"])].copy()
    prim = prim.set_index("hypothesis").loc[["H1", "H2", "H3", "H4", "H5", "H6"]]
    adj = holm_adjust(prim["p_raw"].to_numpy())
    bh = multipletests(prim["p_raw"].to_numpy(), method="fdr_bh")[1]
    h1txt = {"H1": "GDP pc and CO2 pc positively associated (rho>0)", "H2": "CO2 pc differs across income groups",
             "H3": "High-income mean log10 CO2 pc exceeds Lower-middle+Low", "H4": "income-emissions slope differs by income group",
             "H5": "CO2 pc of high-income countries fell 1990-2022", "H6": "decoupling depends on income group"}
    df["p_holm"] = np.nan
    df["p_bh"] = np.nan
    df["decision"] = ""
    for h, a, b in zip(prim.index, adj, bh):
        i = df.index[(df["hypothesis"] == h) & df["primary"]][0]
        df.loc[i, ["p_holm", "p_bh"]] = [a, b]
        df.loc[i, "decision"] = decision_text(a, h1txt[h])
    prim_p = df[df["primary"]].set_index("hypothesis")["p_raw"]
    for i, r in df[~df["primary"]].iterrows():
        if r["hypothesis"] in prim_p.index and r["variant"].startswith("robustness"):
            ok = r["direction_ok"] if pd.notna(r["direction_ok"]) else None
            df.loc[i, "decision"] = f"direction {'agrees' if ok else 'DISAGREES'} with primary; raw p {'<' if r['p_raw'] < ALPHA else '>='} 0.05"
        elif r["hypothesis"] == "H7":
            df.loc[i, "decision"] = "exploratory (unadjusted)"
    cols = ["hypothesis", "variant", "test", "n", "statistic", "p_raw", "p_holm", "p_bh", "effect_size_name", "effect_size", "ci_low", "ci_high",
            "decision", "comment", "primary", "direction_ok"]
    return df[cols]


def run_all() -> pd.DataFrame:
    warnings.filterwarnings("ignore", category=FutureWarning)
    RES.mkdir(exist_ok=True)
    RESULTS.clear()
    ASSUMPTIONS.clear()
    S = load_samples()
    print({k: len(v) for k, v in S.items() if k not in ("cs", "pr")})
    for fn in (run_h1, run_h2, run_h3, run_h4, run_h5, run_h6, run_h7, run_power, run_panel):
        fn(S)
    out = finalize()
    out.to_csv(RES / "hypothesis_results.csv", index=False)
    pd.DataFrame(ASSUMPTIONS).to_csv(RES / "assumption_checks.csv", index=False)
    print("\n=== PRIMARY RESULTS (Holm over H1-H6) ===")
    print(out[out["primary"]][["hypothesis", "test", "n", "statistic", "p_raw", "p_holm", "effect_size_name", "effect_size", "ci_low", "ci_high", "decision"]]
          .round(4).to_string(index=False))
    return out


if __name__ == "__main__":
    run_all()
