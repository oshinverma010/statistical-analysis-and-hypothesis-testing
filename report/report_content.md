# Does Wealth Burn Carbon? A Statistical Analysis and Hypothesis-Testing Framework for CO₂ Emissions and Economic Development, 1990–2022

*Week 2 Task – Statistical Analysis and Hypothesis Testing Framework (execution phase of the Week 1 project)*

---

## 1. Introduction

### 1.1 Context

Week 1 of this project framed one policy-relevant question: does getting richer inevitably mean emitting more carbon dioxide? Climate targets and development goals are often presented as competing priorities. If emissions per person simply track income, growth in poorer countries threatens the global carbon budget; if wealthy countries have begun to *decouple* growth from emissions, there is room for a different story. Week 1 produced the theme, data sources and cleaning rules. This report is the execution phase: a framework of research questions, testable hypotheses, an experimental design and statistical tests run on real public data.

Statistics is treated here as a chain of justified decisions rather than a menu of tests. Every step, from which countries are counted to which p-value enters the multiplicity correction, was written in an analysis plan (`docs/analysis_plan.md`) and committed to version control *before any test was run*. The results in Section 6 are therefore confirmatory for H1–H6 and explicitly exploratory for H7.

![Figure 1.1: The overall experimental framework, from research questions to interpretation. Colours group the phases (blue = question and hypotheses, green = data, amber = assumptions and test choice, red = testing and correction, pink = interpretation).](../figures/fig_framework_flowchart.png)

### 1.2 Research questions

| ID | Research question | Why it matters |
|----|-------------------|----------------|
| RQ1 | Is income associated with emissions per person across countries? | The baseline "wealth burns carbon" claim. |
| RQ2 | Do income groups differ in emissions per person? | Tests whether the association is visible at the level of policy-relevant country categories. |
| RQ3 | Does the income–emissions slope differ by income group? | A flatter slope at high income would be a first sign of saturation or decoupling. |
| RQ4 | Have wealthy countries decoupled growth from emissions? | The temporal question: have emissions fallen while GDP per person rose? |
| RQ5 | (Exploratory) After controlling for income, do industrial share and urbanisation explain emissions per person? | Probes structural drivers beyond income. |

### 1.3 Scope

The analysis covers countries (not regions, not individuals), the years 1990–2022, and one pollutant measure: territorial CO₂ from fossil fuels and industry, excluding land-use change. Statements are therefore about *countries* and *production-based* emissions. They are not statements about individual consumers, about embodied emissions in traded goods, or about other greenhouse gases. These limits are carried through to Section 9.

---

## 2. Hypothesis Formulation

### 2.1 Notation and conventions

Let *y* = log₁₀(CO₂ emissions per capita, tonnes per person) and *x* = log₁₀(GDP per capita, constant 2015 US$). Logarithms are used because both variables are positive and span about three orders of magnitude across countries (Section 4.2). Tests that are rank-based (Spearman, Kruskal–Wallis, Wilcoxon) are unaffected by the transformation; tests on means or regression slopes use the logged variables.

### 2.2 Hypothesis table

| ID | RQ | H0 | H1 | Direction | Test (pre-specified primary) | Unit of analysis |
|----|----|----|----|-----------|------------------------------|------------------|
| H1 | RQ1 | Spearman ρ between GDP pc and CO₂ pc = 0 | ρ > 0 | one-sided | Spearman rank correlation | Country, 2022 cross-section |
| H2 | RQ2 | CO₂ pc distribution is the same in all four World Bank income groups | At least one group differs | two-sided | Kruskal–Wallis, Dunn post hoc (Holm) | Country, 2022, grouped by income group |
| H3 | RQ2 | Mean log₁₀ CO₂ pc: High = (Lower-middle + Low) | High > (Lower-middle + Low) | one-sided | Welch's t-test on log₁₀ values | Country, 2022 |
| H4 | RQ3 | All income-group × log-GDP interaction coefficients = 0 | At least one ≠ 0 | F-test (upper tail) | OLS with interaction, HC3-robust F-test | Country, 2022 |
| H5 | RQ4 | Median of ln(CO₂pc₂₀₂₂ / CO₂pc₁₉₉₀) among high-income countries = 0 | Median < 0 | one-sided | Wilcoxon signed-rank | Country, paired 1990 → 2022 |
| H6 | RQ4 | Absolute decoupling is independent of income group | Not independent | two-sided | Chi-square test of independence (Monte Carlo if any expected count < 5) | Country, 1990 → 2022 |
| H7 | RQ5 | Industry share and urbanisation have zero partial coefficients given log GDP pc | At least one ≠ 0 | two-sided | Multiple OLS, HC3 SEs (**exploratory**) | Country, 2022 |

### 2.3 Why each hypothesis is formulated this way

**H1 (one-sided Spearman).** The "wealth burns carbon" hypothesis predicts a positive association, so a one-sided alternative is justified and more powerful. Spearman is preferred to Pearson as the primary test because it assumes neither linearity nor bivariate normality; it only asks whether richer countries tend to rank higher in emissions. Pearson on log₁₀ values is a companion, outside the confirmatory family.

**H2 and H3 (two group questions).** H2 is the omnibus question (any difference among four groups); H3 is the sharper, directional question that policy debates actually pose: do high-income countries emit more than the poorer half of the world? An omnibus test does not say *which* groups differ (hence the post hoc), whereas a two-group directional contrast is more powerful and easier to communicate. Lower-middle and Low income are pooled in H3 because the Low group is small (n = 25) and the pooled contrast reads as "rich vs poor".

**H4 (interaction).** Fitting four regressions and comparing slopes by eye is not a test; the interaction term in one model is. The null is that a single elasticity of emissions to income describes all four groups. HC3 robust covariance is specified because residual variance commonly differs between groups.

**H5 (paired, one-sided).** The same country measured twice is a paired design; treating 1990 and 2022 as independent groups would waste the pairing. The ratio 2022/1990 is analysed on the natural-log scale so that halving and doubling are symmetric around zero. Restricting H5 to high-income countries asks exactly "have the wealthy decoupled?".

**H6 (two categorical variables).** H5 says whether emissions fell on average among the wealthy; H6 asks whether the *probability* of absolute decoupling (GDP per capita up *and* CO₂ per capita down) differs across income groups. This is a contingency-table question, so a chi-square test of independence is natural, with Monte Carlo or Fisher-type exact methods as the fallback when expected counts are small.

**H7 (exploratory).** Industry share and urbanisation are plausible structural drivers, but Week 1 gave no firm prediction of their size or sign. Labelling H7 exploratory keeps it out of the confirmatory family and avoids inflating the number of "findings".

### 2.4 Decision rules (fixed in advance)

* α = 0.05; two-sided unless a direction is stated (H1, H3, H5 are one-sided).
* Confirmatory family H1–H6, one primary p-value each, **Holm–Bonferroni** correction across the six. H7 is reported separately with unadjusted p-values.
* Every result is reported with *n*, effect size and 95% confidence interval alongside the p-value.
* A non-significant result is worded "no evidence of …", never "no effect". Likewise a p-value is never reported as "the probability that H0 is true".

---

## 3. Experimental Design and Data

### 3.1 Sources

| Source | Content | Access |
|--------|---------|--------|
| Our World in Data (OWID) CO₂ dataset | CO₂ per capita (fossil fuels + industry, excluding land-use change), total CO₂, population | `https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv` (HTTP 200 verified; 14.4 MB) |
| World Bank Indicators API v2 | `NY.GDP.PCAP.KD` (GDP pc, constant 2015 US$), `NY.GDP.PCAP.PP.KD` (PPP), `SP.POP.TOTL`, `SP.URB.TOTL.IN.ZS`, `NV.IND.TOTL.ZS` | `https://api.worldbank.org/v2/country/all/indicator/<code>?date=1990:2022` |
| World Bank country-metadata endpoint | Income group and region for each economy | `https://api.worldbank.org/v2/country` |

All URLs responded on the first attempt, so no substitute source was needed. Downloads were made on 2026-11-05; each file's SHA-256 is stored in `data/raw/manifest.json`.

### 3.2 Sampling frame and unit of analysis

The sampling frame is the set of 229 country-level entities that appear in either source after aggregates are removed: the 217 economies on the World Bank country list plus 12 OWID-only territories (for example Taiwan, Anguilla and the Cook Islands) that have no World Bank record. This is a census of available countries rather than a random sample. The inferential question is therefore whether the *observed* pattern across countries is stronger than chance association would produce, not whether a sample generalises to a larger population of countries (see Section 9).

The unit of analysis is **the country**. Country-years are deliberately *not* treated as independent observations: successive years of one country are strongly autocorrelated, and using 6,000+ country-years would suggest false precision. H1–H4 and H7 use the 2022 cross-section; H5 and H6 use one paired 1990-to-2022 change per country; a pooled panel with clustered standard errors is only a sensitivity check.

### 3.3 Cleaning rules (reused from Week 1)

1. **Separate non-country aggregates.** OWID rows without an ISO code (regions, "World", international shipping and aviation, and income-group aggregates; 35 entities) are removed. World Bank rows whose region is "Aggregates" (79 of 296 metadata rows) are removed. Five World Bank income-level aggregate series come back with a blank ISO3 in the indicator files and are removed by the same rule.
2. **Join on ISO3 + year only, never on names.** A manual crosswalk is applied only for Kosovo, whose OWID rows carry no ISO code while the World Bank uses `XKX`; 33 rows are re-keyed.
3. **Never impute the outcome or the main predictor.** CO₂ per capita and GDP per capita are used only where observed. The exploratory covariates (industry share, urbanisation) may use the latest value in 2018–2022 when 2022 is missing (H7 only; a 2022-only sensitivity is also run).
4. **Constant-2015-US$ GDP per capita is primary; PPP GDP is a sensitivity check.**
5. **Keep legitimate extremes** (for example Qatar at 37.9 t CO₂ per person) and handle skew with log₁₀ scales rather than deleting observations.
6. **Restrict to 1990–2022** (OWID extends to 2024).
7. **Log every dropped or changed row** (`data/processed/cleaning_log_rows.csv`, 9,367 rows; step totals in `results/cleaning_log.csv`).

**Cleaning log (step totals from `results/cleaning_log.csv`)**

| Step | Rule | Rows before | Rows after | Note |
|---|---|---|---|---|
| wb_meta | drop World Bank aggregates (region == 'Aggregates') | 296 | 217 |  |
| wb_ind:NY.GDP.PCAP.KD | keep country rows only | 8,778 | 7,161 | 459 null values kept as missing (not imputed) |
| wb_ind:NY.GDP.PCAP.PP.KD | keep country rows only | 8,778 | 7,161 | 761 null values kept as missing (not imputed) |
| wb_ind:SP.POP.TOTL | keep country rows only | 8,778 | 7,161 | 0 null values kept as missing (not imputed) |
| wb_ind:SP.URB.TOTL.IN.ZS | keep country rows only | 8,778 | 7,161 | 0 null values kept as missing (not imputed) |
| wb_ind:NV.IND.TOTL.ZS | keep country rows only | 8,778 | 7,161 | 1,119 null values kept as missing (not imputed) |
| owid_years | restrict to 1990-2022 | 50,411 | 8,375 |  |
| owid_crosswalk | manual crosswalk Kosovo -> XKX | 33 | 0 | rows changed (not dropped) |
| owid_agg | drop OWID rows with no ISO3 (aggregates, regions, shipping, aviation...) | 8,375 | 7,227 | 35 non-country entities |
| join | outer join OWID x World Bank on (iso3, year) | 14,388 | 7,557 | both=6,831, owid_only=396, wb_only=330 |
| range:urban_pct | shares outside 0-100 set to NaN | 7,161 | 7,161 |  |
| range:industry_pct | shares outside 0-100 set to NaN | 6,042 | 6,042 |  |
| logguard:co2_pc | non-positive values excluded from log transforms | 7,045 | 7,035 |  |
| logguard:gdp_pc | non-positive values excluded from log transforms | 6,702 | 6,702 |  |
| logguard:gdp_pc_ppp | non-positive values excluded from log transforms | 6,400 | 6,400 |  |

The join of OWID and World Bank country-years produced 7,557 country-year rows: 6,831 present in both sources, 396 only in OWID (territories without a World Bank record) and 330 only in the World Bank file (for example Puerto Rico and Guam, which OWID does not report). Ten country-year rows with CO₂ per capita ≤ 0 exist (Timor-Leste 1994–1997 and Kosovo 1994–1999, all early-record zeros); they are retained in the panel but excluded wherever a logarithm is taken. None occurs in 2022 or 1990 for any analysed country.

### 3.4 Who enters each test, and why others do not

| Test | Countries entering | Exclusions (real counts, applied sequentially) |
|------|--------------------|-----------------------------------------------|
| H1 | **198** | 229 entities → 15 with no 2022 CO₂ per capita (mostly small territories such as Puerto Rico, Guam, the Cayman Islands, Monaco, San Marino) → 16 with no 2022 GDP per capita (e.g. Eritrea, North Korea, South Sudan, Yemen, Liechtenstein, New Caledonia and nine OWID-only territories) |
| H2 | **205** | 229 → 15 no 2022 CO₂ → 9 with no World Bank income group (Taiwan and eight small OWID-only territories) |
| H3 | **146** (74 High vs 72 Lower-middle + Low) | 205 → 59 Upper-middle-income countries are outside the contrast |
| H4 | **198** | same as H1; all have an income group (21 Low, 47 Lower-middle, 59 Upper-middle, 71 High) |
| H5 | **74** | 86 high-income entities → 12 with no 1990 CO₂ record (micro-territories such as Gibraltar, Guam, Monaco, Puerto Rico) |
| H6 | **184** | 217 entities with an income group → 18 with no 1990 CO₂ (including Kosovo, whose series starts in 1994, and Timor-Leste, Namibia, Eritrea) → 14 with no 1990 GDP (e.g. Afghanistan, Serbia, Montenegro, Djibouti) → 1 with no 2022 GDP (Yemen) |
| H7 | **194** | as H1 → 4 lacking any industry-share value in 2018–2022 (Nauru, Somalia, Tuvalu, Venezuela) |

Sample sizes differ across tests because each test uses only the variables it needs; a country is never dropped from a test for a variable that test does not use. Robustness subsets additionally exclude 66 entities with population < 1 million, and 23 pre-listed oil and gas exporters (list fixed in the analysis plan before the data were examined).

### 3.5 Integrity checks

A pytest suite (`tests/test_data_integrity.py`, 25 tests) verifies: no duplicate (iso3, year) keys; no World Bank or OWID aggregates in country tables and only three-letter ISO codes; Kosovo joined through the crosswalk; Qatar retained; every value entering a logarithm is positive; urbanisation and industry shares lie within 0–100 and look-back years within 2018–2022; non-null counts of CO₂ and GDP in the cleaned panel equal the raw non-null counts (proof that nothing was imputed); expected row counts at each stage; the join count reported in the cleaning log; per-hypothesis sample sizes; the definition of the decoupling flag; and raw-file checksums against the manifest. All 25 tests pass (Section 7).

---

## 4. Statistical Testing Framework

### 4.1 Logic of test selection

For each hypothesis the primary test was fixed in the plan from the *structure* of the question (two continuous variables, k independent groups, paired measures, moderated slope, two categorical variables), and then each test's assumptions were checked against the real data. The rule adopted in the plan is that a violated assumption never causes a switch away from a pre-specified primary test; it triggers a documented fallback only for the *parametric comparison* test. This protects against choosing the test that gives the smaller p-value.

![Figure 4.1: Test-selection decision tree showing the branches actually taken in this project. P-values on the boxes come from the assumption checks in Section 7.](../figures/fig_test_selection_tree.png)

### 4.2 Assumption checks and what they showed

Normality is assessed with Shapiro–Wilk (Shapiro and Wilk, 1965) together with QQ plots, because with n ≈ 200 Shapiro–Wilk flags trivial departures and the plot shows whether the departure matters. Equal variances are assessed with the median-centred Levene (Brown–Forsythe) test. Heteroscedasticity of regression residuals uses Breusch–Pagan (Breusch and Pagan, 1979), influence uses Cook's distance (threshold 4/n), and collinearity uses variance-inflation factors.

![Figure 4.2: Histograms and normal QQ plots, raw versus log-transformed, for CO₂ per capita, GDP per capita and the high-income CO₂ ratio 2022/1990. The log transform removes most of the right skew (skewness 2.58 → −0.58 for CO₂ pc) but not all departures from normality.](../figures/fig_assumption_checks.png)

**Assumption checks (from `results/assumption_checks.csv`)**

| ID | Check | Statistic | p | Verdict |
|---|---|---|---|---|
| H1 | Shapiro-Wilk: log10 GDP pc | 0.980 | 0.0069 | NOT normal (p<=0.05) |
| H1 | Shapiro-Wilk: log10 CO2 pc | 0.966 | 9.83e-05 | NOT normal (p<=0.05) |
| H1 | Shapiro-Wilk: raw GDP pc (before log) | 0.713 | 3.43e-18 | NOT normal (p<=0.05) |
| H1 | Shapiro-Wilk: raw CO2 pc (before log) | 0.733 | 1.37e-17 | NOT normal (p<=0.05) |
| H1 | skewness raw vs log10 CO2 pc | 2.583 | – | log transform reduces skew |
| H1 | outliers |z|>3 on log10 CO2 pc | 0.000 | – | none |
| H2 | Shapiro-Wilk: raw CO2 pc, Low income | 0.483 | 2.56e-08 | NOT normal (p<=0.05) |
| H2 | Shapiro-Wilk: log10 CO2 pc, Low income | 0.880 | 0.0069 | NOT normal (p<=0.05) |
| H2 | Shapiro-Wilk: raw CO2 pc, Lower middle income | 0.810 | 2.62e-06 | NOT normal (p<=0.05) |
| H2 | Shapiro-Wilk: log10 CO2 pc, Lower middle income | 0.963 | 0.1353 | normal (p>0.05) |
| H2 | Shapiro-Wilk: raw CO2 pc, Upper middle income | 0.808 | 2.49e-07 | NOT normal (p<=0.05) |
| H2 | Shapiro-Wilk: log10 CO2 pc, Upper middle income | 0.981 | 0.4946 | normal (p>0.05) |
| H2 | Shapiro-Wilk: raw CO2 pc, High income | 0.790 | 6.71e-09 | NOT normal (p<=0.05) |
| H2 | Shapiro-Wilk: log10 CO2 pc, High income | 0.982 | 0.3648 | normal (p>0.05) |
| H2 | Levene/Brown-Forsythe on log10 CO2 pc | 0.258 | 0.8557 | equal variances |
| H2 | Levene/Brown-Forsythe on raw CO2 pc | 14.970 | 7.82e-09 | UNEQUAL variances |
| H3 | Shapiro-Wilk: log10 CO2 pc, High | 0.982 | 0.3648 | normal (p>0.05) |
| H3 | Shapiro-Wilk: log10 CO2 pc, Lower-middle + Low | 0.982 | 0.3971 | normal (p>0.05) |
| H3 | Levene (Brown-Forsythe) High vs LM+Low | 14.891 | 0.0002 | UNEQUAL variances -> Welch justified |
| H4 | Breusch-Pagan (full model) | 2.836 | 0.8998 | homoscedastic |
| H4 | Shapiro-Wilk: residuals of full model | 0.974 | 0.0010 | NOT normal (p<=0.05) |
| H4 | Cook's distance > 4/n | 8.000 | 0.0202 | 8 influential points |
| H4 | max VIF (main-effects model) | 8.321 | – | elevated (>=5) |
| H5 | Shapiro-Wilk: log ratio ln(CO2pc_2022/CO2pc_1990) | 0.966 | 0.0460 | NOT normal (p<=0.05) |
| H6 | minimum expected count | 6.413 | – | all >= 5 -> chi-square valid |
| H7 | max VIF | 1.841 | – | OK (<5) |
| H7 | Breusch-Pagan | 2.090 | 0.5540 | homoscedastic |
| H7 | Shapiro-Wilk: residuals | 0.995 | 0.7726 | normal (p>0.05) |

Key observations, all from the table above:

* **Raw CO₂ and GDP per capita are strongly right-skewed** (skewness 2.58 and 2.16; Shapiro–Wilk p < 0.0001), which justifies both the log₁₀ transform and a rank-based primary test for H1 and H2.
* **After log₁₀, the group distributions are close to normal** for Lower-middle, Upper-middle and High income (p = 0.14, 0.49, 0.36) and for the two H3 groups (p = 0.36, 0.40). The Low-income group still fails (p = 0.007, n = 25), so the omnibus parametric test is only a comparison for H2.
* **Variances are equal across the four groups on the log scale** (Levene p = 0.86) but **unequal for the H3 contrast** (p = 0.0002; SD ratio High : Lower-middle+Low = 0.63), which is the data-based reason for Welch's test over Student's.
* **The H4 model has non-normal residuals** (Shapiro p = 0.001), no evidence of heteroscedasticity by Breusch–Pagan (p = 0.90) and **a maximum VIF of 8.3**. The high VIF is expected: log GDP and the income-group dummies are mechanically related because the World Bank defines income groups from income. It is not a data-entry problem and does not bias the F-test, but it inflates the standard errors of individual coefficients.
* **The log ratios for H5 are marginally non-normal** (p = 0.046, skewness 0.53), so the paired t-test is shown only as a companion.
* **All expected counts in the H6 table are at least 5** (minimum 6.41), so the ordinary chi-square approximation is valid; a Monte Carlo p-value is still computed as a check.

### 4.3 The tests in plain language

* **Spearman rank correlation (H1)** correlates ranks: when one country is richer than another, is it also usually the higher emitter? It is invariant to monotone transforms, so it is robust to skew and to outliers such as Qatar.
* **Kruskal–Wallis (H2)** (Kruskal and Wallis, 1952†) ranks all countries together and asks whether the average rank differs across groups more than chance would allow. Because it only says *that* groups differ, **Dunn's post hoc test** (Dunn, 1964) compares pairs, with **Holm** adjustment over six pairs; Mann–Whitney is a cross-check.
* **Welch's t-test (H3)** (Welch, 1947) compares two means without assuming equal variances by adjusting the degrees of freedom. With unequal group variances (Levene p = 0.0002) Student's test would be miscalibrated. The log₁₀ scale makes the mean difference interpretable: a difference of 1 on the log₁₀ scale means a ten-fold ratio of geometric means.
* **OLS with an interaction (H4)** lets each income group have its own intercept and slope. The F-test on the block of interaction coefficients asks whether letting slopes differ improves the fit more than chance. The robust (HC3) covariance (MacKinnon–White family; White, 1980) keeps the test valid if residual variance differs across groups.
* **Wilcoxon signed-rank (H5)** (Wilcoxon, 1945) ranks the absolute within-country changes and asks whether positive or negative changes carry the larger ranks. It uses the sign and size of each change but needs only a roughly symmetric distribution of log-ratios, not normality.
* **Chi-square test of independence (H6)** compares the observed counts of decoupled/not-decoupled countries in each income group with the counts expected if income group and decoupling were unrelated.

### 4.4 Effect sizes and confidence intervals

| Hypothesis | Effect size | Rule of thumb used for interpretation (Cohen, 1988†; Tomczak and Tomczak, 2014†) | CI method |
|------------|-------------|----------------------------------------|-----------|
| H1 | Spearman ρ | ≈ 0.1 small, 0.3 medium, 0.5 large | Percentile bootstrap over countries, B = 10,000 |
| H2 | ε² (= H/(n−1)) | ≈ 0.01 small, 0.08 medium, 0.26 large | Bootstrap |
| H3 | Cohen's d (pooled SD); mean difference on log₁₀ scale | 0.2 / 0.5 / 0.8 | Stratified bootstrap for d; Welch CI for the difference |
| H4 | ΔR² of the interaction block (Cohen's f²) | f² ≈ 0.02 / 0.15 / 0.35 | Bootstrap (2,000) |
| H5 | Matched-pairs rank-biserial r; median % change | like r | Bootstrap |
| H6 | Cramér's V; standardised residuals | like r for 2 columns | Bootstrap |

The bootstrap (Efron, 1979†) resamples countries with replacement (seed 42); its intervals are two-sided 95% intervals even where the test is one-sided.

### 4.5 Multiple testing

With six confirmatory tests, the chance of at least one false rejection at α = 0.05 would be about 26% if all nulls were true and the tests independent. Holm's step-down procedure (Holm, 1979†) controls the family-wise error rate at 5% without assuming independence and is uniformly more powerful than plain Bonferroni. Benjamini–Hochberg (1995), which controls the false discovery rate, is shown as a secondary column. Holm is the *decision* rule because there are few claims and each is a headline finding, where family-wise control is appropriate.

### 4.6 Power analysis

With the real sample sizes, `statsmodels` gives the **minimum detectable effect at 80% power** (α = 0.05, one-sided): **d = 0.414** for H3 (74 vs 72 countries) and **d_z = 0.292** for H5 (n = 74; 0.299 after a Wilcoxon-efficiency adjustment of 0.955). H3 has 91% power to detect a medium effect (d = 0.5), while H5 has 52% power to detect a small effect (d_z = 0.2) and 99.6% power for d_z = 0.5. Power matters for interpreting H5 later: the observed paired effect (|d_z| = 0.25) is *below* the minimum detectable effect, so a non-significant result would have been uninformative, and the significant result should be read as a borderline-sized effect, not a large one.

![Figure 4.3: Power curves for H3 and H5 at the real sample sizes. The observed H3 effect (d = 3.06) lies far to the right of the plotted range.](../figures/fig_power_curve.png)

---

## 5. Step-by-Step Methodology

Each step shows real code from `src/` and explains what it does and why. The code uses pandas, SciPy (Virtanen et al., 2020†), statsmodels (Seabold and Perktold, 2010), scikit-posthocs and pingouin (Vallat, 2018).

### Step 1 – Acquire data without overwriting

```python
# from src/acquire.py
def fetch_owid(manifest: dict) -> Path:
    dest = RAW / "owid-co2-data.csv"
    if dest.exists():
        print(f"cached: {dest.name}")
        return dest
    print("downloading OWID CO2 CSV ...")
    r = _get(OWID_URL)
    dest.write_bytes(r.content)
    _record(manifest, dest, OWID_URL)
```

An existing file is never overwritten, so the analysed bytes cannot change silently; the manifest records the download date, size and SHA-256. World Bank calls go through a retry-with-backoff helper; in this project the first attempt always succeeded, so the retry path was not exercised.

### Step 2 – Remove aggregates and apply the Kosovo crosswalk

```python
# from src/clean.py
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
```

Rows with no ISO code are regions or non-country series (35 entities) and must not enter a country-level test. Kosovo is rescued by an explicit, logged crosswalk; no other name-based matching is allowed. The final check raises an error if a duplicated key survives.

### Step 3 – Join on ISO3 + year and report the join

```python
# from src/clean.py
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
```

`validate="one_to_one"` makes pandas fail loudly if the key is not unique. An *outer* join is used so that unmatched rows are counted and reported rather than silently dropped (`both = 6,831`, `owid_only = 396`, `wb_only = 330`).

### Step 4 – Guard logarithms

```python
# from src/tests.py
pos = lambda s: s.where(s > 0)  # noqa: E731  non-positive -> NaN so log10 is never taken of <= 0
cs["lg_co2"] = np.log10(pos(cs["co2_pc"]))
cs["lg_gdp"] = np.log10(pos(cs["gdp_pc"]))
cs["lg_gdp_ppp"] = np.log10(pos(cs["gdp_pc_ppp"]))
s = {}
s["H1"] = cs.dropna(subset=["lg_co2", "lg_gdp"]).reset_index(drop=True)
```

Non-positive values become missing *only inside the log calculation*, never in the stored data; each test's sample is the countries observed on the variables that test needs.

### Step 5 – Bootstrap helper (fixed seed)

```python
# from src/tests.py
def boot_ci(n: int, stat_fn, rng: np.random.Generator, reps: int = B) -> tuple[float, float]:
    """Percentile bootstrap 95% CI. stat_fn receives an index array of resampled rows."""
    vals = np.empty(reps)
    for i in range(reps):
        vals[i] = stat_fn(rng.integers(0, n, n))
    vals = vals[np.isfinite(vals)]
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
```

One helper provides every percentile bootstrap interval, so the procedure is identical across hypotheses; each call gets its own `np.random.default_rng(42)`, so results do not depend on execution order.

### Step 6 – H1: Spearman with bootstrap CI

```python
# from src/tests.py
rho_res = stats.spearmanr(df[x], df["lg_co2"], alternative="greater")
rng = np.random.default_rng(SEED)
xv, yv = df[x].to_numpy(), df["lg_co2"].to_numpy()
lo, hi = boot_ci(len(df), lambda i: spearman_stat(xv[i], yv[i]), rng)
add_result("H1", label, "Spearman rho (one-sided, greater)", len(df), rho_res.statistic, rho_res.pvalue, "Spearman rho",
           rho_res.statistic, lo, hi, primary, rho_res.statistic > 0)
```

`alternative="greater"` implements the one-sided H1. The interval comes from the bootstrap, not from a normal approximation, because the sampling distribution of ρ is skewed near 1.

### Step 7 – H3: Welch's t-test, mean-difference CI and Cohen's d

```python
# from src/tests.py
t = stats.ttest_ind(a, b, equal_var=False, alternative="greater")
diff = a.mean() - b.mean()
va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
se = np.sqrt(va + vb)
df_ = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
tc = stats.t.ppf(0.975, df_)
rng = np.random.default_rng(SEED)
bs = [cohens_d(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(B)]
lo, hi = np.percentile(bs, [2.5, 97.5])
```

`equal_var=False` selects Welch. The Welch–Satterthwaite degrees of freedom are computed explicitly so the confidence interval for the mean difference matches the test. Cohen's d is bootstrapped by resampling *within* each group so group sizes stay fixed.

### Step 8 – H4: interaction F-test with robust covariance

```python
# from src/tests.py
m_full = ols(f"lg_co2 ~ {xcol} * {grp}", data=df).fit()
m_main = ols(f"lg_co2 ~ {xcol} + {grp}", data=df).fit()
classical = anova_lm(m_full, typ=2)
robust = anova_lm(m_full, typ=2, robust="hc3")
inter = [i for i in classical.index if ":" in i][0]
F_c, p_c = classical.loc[inter, "F"], classical.loc[inter, "PR(>F)"]
F_r, p_r = robust.loc[inter, "F"], robust.loc[inter, "PR(>F)"]
```

Both the classical and the HC3-robust F-tests are computed; the robust one is the pre-specified primary p-value.

### Step 9 – H5: Wilcoxon signed-rank on paired log ratios

```python
# from src/tests.py
w = stats.wilcoxon(x, alternative="less")
ranks = stats.rankdata(np.abs(x))
rbc = (ranks[x > 0].sum() - ranks[x < 0].sum()) / ranks.sum()  # matched-pairs rank-biserial (negative = decreases)
```

The matched-pairs rank-biserial correlation is computed from the ranks of the absolute changes; negative values mean that decreases dominate.

### Step 10 – H6: chi-square with Monte Carlo check

```python
# from src/tests.py
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
```

`correction=False` because Yates' correction applies only to 2 × 2 tables. The Monte Carlo p-value permutes the decoupling labels 100,000 times (an analogue of the Fisher–Freeman–Halton exact test) and adds one to numerator and denominator so a p-value of exactly zero is never reported.

### Step 11 – Holm correction across the confirmatory family

```python
# from src/tests.py
prim = df[(df["primary"]) & df["hypothesis"].isin(["H1", "H2", "H3", "H4", "H5", "H6"])].copy()
prim = prim.set_index("hypothesis").loc[["H1", "H2", "H3", "H4", "H5", "H6"]]
adj = holm_adjust(prim["p_raw"].to_numpy())
bh = multipletests(prim["p_raw"].to_numpy(), method="fdr_bh")[1]
```

### Step 12 – Power analysis at the real n

```python
# from src/tests.py
mde3 = TTestIndPower().solve_power(effect_size=None, nobs1=n_high, alpha=ALPHA, power=0.8, ratio=n_low / n_high, alternative="larger")
mde5 = TTestPower().solve_power(effect_size=None, nobs=n5, alpha=ALPHA, power=0.8, alternative="larger")
mde5w = TTestPower().solve_power(effect_size=None, nobs=n5 * WILCOXON_ARE, alpha=ALPHA, power=0.8, alternative="larger")
```

`solve_power` with `effect_size=None` returns the smallest standardised effect that the design can detect with 80% power.

---

## 6. Results

### 6.1 Summary of confirmatory results

| ID | Primary test | n | Statistic | Raw p | Holm-adj. p | Effect size [95% CI] | Decision (α = 0.05) |
|---|---|---|---|---|---|---|---|
| H1 | Spearman rho (one-sided, greater) | 198 | 0.855 | 4.98e-58 | 2.99e-57 | Spearman rho = 0.855 [0.802, 0.893] | Reject H0 |
| H2 | Kruskal-Wallis H (df=3) | 205 | 146 | 1.65e-31 | 6.60e-31 | epsilon-squared = 0.717 [0.644, 0.782] | Reject H0 |
| H3 | Welch t-test (df=119.3, one-sided greater) | 146 | 18.355 | 7.45e-37 | 3.73e-36 | Cohen's d = 3.057 [2.680, 3.550] | Reject H0 |
| H4 | Interaction F-test, HC3 Wald-type (anova_lm robust='hc3') | 198 | 3.691 | 0.0129 | 0.0155 | delta R-squared (interaction block) = 0.010 [0.003, 0.026] | Reject H0 |
| H5 | Wilcoxon signed-rank (one-sided, less) | 74 | 938 | 0.0077 | 0.0155 | matched-pairs rank-biserial r = -0.324 [-0.568, -0.067] | Reject H0 |
| H6 | Chi-square test of independence (df=3) | 184 | 52.064 | 2.90e-11 | 8.71e-11 | Cramer's V = 0.532 [0.413, 0.655] | Reject H0 |
| H7* | Joint HC3 Wald F: industry share & urbanisation (exploratory) | 194 | 30.817 | 2.56e-12 | not adjusted | increment in R-squared = 0.050 [0.027, 0.085] | exploratory |

*\*H7 is exploratory and outside the Holm family.*

*Notes.* H1, H3 and H5 are one-sided tests as pre-specified; H2, H4 and H6 are omnibus tests of a two-sided or upper-tail alternative (see Section 2). "Raw p" values below 10⁻⁴ are shown in scientific notation; values this extreme are not meaningful beyond "very small", and the effect size and confidence interval are the informative quantities. The H5 statistic is the Wilcoxon sum of positive-change ranks (V). For H6 the p-value shown is the ordinary chi-square p-value because all expected counts exceed 5; the Monte Carlo p-value is 1 × 10⁻⁵, which is the smallest value that 100,000 permutations can resolve. Effect-size CIs are 95% percentile bootstrap intervals.

**All six confirmatory null hypotheses are rejected after Holm correction.** The strength of evidence is, however, very different across hypotheses: H1, H2, H3 and H6 have effect sizes in the "large" range with confidence intervals far from zero, whereas H4 and H5 are small-to-moderate effects whose Holm-adjusted p-values (0.0155 each) are far closer to the 0.05 threshold.

### 6.2 H1 – Income and emissions per person

Across 198 countries, GDP per capita and CO₂ per capita are strongly and positively associated: Spearman ρ = 0.855 (95% bootstrap CI 0.802–0.893); the Pearson correlation of the log₁₀ values is almost identical (r = 0.852, 95% CI 0.808–0.886), which indicates that the relationship is close to linear on the log–log scale (OLS slope 0.83). In ordinary terms: a ten-fold difference in income per person goes with roughly a 6.8-fold difference in emissions per person (10^0.83 ≈ 6.8). The expected outcome in the plan (ρ ≈ 0.8 or higher) was met.

![Figure 6.1: H1 – log–log scatter of CO₂ per capita against GDP per capita (2022) with OLS fit, LOWESS smooth and Spearman ρ. Only extreme countries are labelled. LOWESS flattens at high income, previewing H4.](../figures/fig_h1_scatter.png)

The LOWESS curve flattens above roughly 20,000 US$ (Figure 6.1), a visual hint of the slope heterogeneity tested in H4. Qatar, Kuwait, Brunei, Mongolia and Venezuela lie well above the fitted line; Macao and Switzerland lie well below it.

*Robustness.* The result is stable with PPP GDP (ρ = 0.869, n = 195), without countries under one million people (ρ = 0.862, n = 155), without oil and gas exporters (ρ = 0.884, n = 175) and under a permutation test (p = 1 × 10⁻⁴, the smallest value 10,000 permutations can resolve).

### 6.3 H2 and H3 – Income groups

Median CO₂ per capita in 2022 rises at each step of the income ladder: 0.16 t (Low, n = 25), 0.76 t (Lower-middle, n = 47), 3.03 t (Upper-middle, n = 59) and 6.86 t (High, n = 74). Kruskal–Wallis rejects equality of the four distributions (H = 146.3, df = 3, Holm-adjusted p = 6.6 × 10⁻³¹; ε² = 0.717, 95% CI 0.644–0.782). **All six Dunn post hoc comparisons are significant after Holm adjustment**; the closest pair is Low vs Lower-middle (adjusted p = 0.014, rank-biserial r = −0.84).

![Figure 6.2: H2 – CO₂ per capita by World Bank income group (log scale), with individual countries, box summaries and the Kruskal–Wallis result.](../figures/fig_h2_boxplot_by_income.png)

The parametric comparison agrees (ANOVA on log₁₀ values F = 205.7, η² = 0.754; Welch ANOVA F = 176.6), and excluding small countries (ε² = 0.740) or oil/gas exporters (ε² = 0.742) changes nothing.

For **H3** the geometric-mean CO₂ per capita is 7.14 t in high-income countries against a pooled lower group that is 14.1 times smaller: the mean difference in log₁₀ CO₂ per capita is 1.149 (95% CI 1.025–1.273). Welch's t = 18.36 on 119.3 degrees of freedom (Holm-adjusted p = 3.7 × 10⁻³⁶); Cohen's d = 3.06 (95% CI 2.68–3.55). This is an extraordinarily large standardised effect, reflecting how tightly income group and emissions are tied (and, in part, how income groups are defined; Section 9). Mann–Whitney agrees (rank-biserial r = 0.983, 95% CI 0.963–0.997), as do a permutation test and the two exclusion variants (d = 2.90 and 3.07).

### 6.4 H4 – Do slopes differ by income group?

The interaction F-test with HC3 covariance gives F = 3.69 (raw p = 0.0129; Holm-adjusted p = 0.0155), so H0 is rejected at the 5% family-wise level. The interaction block adds only ΔR² = 0.010 (95% CI 0.003–0.026; Cohen's f² = 0.051): income group changes the *slope* detectably but modestly. The group-specific elasticities show where:

| Income group | n | Slope of log₁₀ CO₂ pc on log₁₀ GDP pc | 95% CI (HC3) |
|--------------|---|--------------------------------------:|--------------|
| Low | 21 | 0.874 | 0.225 – 1.524 |
| Lower-middle | 47 | 0.841 | 0.476 – 1.206 |
| Upper-middle | 59 | 0.694 | 0.251 – 1.136 |
| High | 71 | 0.220 | −0.007 – 0.448 |

![Figure 6.3: H4 – Within-group fits of log CO₂ per capita on log GDP per capita; the slope in high-income countries is markedly flatter.](../figures/fig_h4_slopes_by_income.png)

In poorer countries a doubling of income is associated with a near-proportional rise in emissions (slope ≈ 0.8–0.9), while among high-income countries income level has almost no relationship with emissions (0.22, CI includes 0). That agrees with the expectation in the plan (flatter slope at high income), but **this is the least robust of the six confirmatory results**:

* Excluding countries with population < 1 million raises the p-value to **0.074** (n = 155).
* The planned label-permutation test gives p = **0.052**. It shuffles income-group labels, so it examines a stronger null (group is irrelevant altogether) than the interaction null; a Freedman–Lane (1983†) residual permutation, added after planning because it tests exactly the interaction null, gives p = 0.024.
* Removing the 8 influential countries (Cook's D > 4/n: Afghanistan, Costa Rica, Kazakhstan, Mongolia, Qatar, Syria, Ukraine, Venezuela) *strengthens* the result (p = 1.2 × 10⁻⁵, n = 190); PPP GDP (p = 0.023) and excluding exporters (p = 0.007) agree.
* The VIF of 8.3 and the non-normal residuals do not undermine the F-test, because HC3 and permutation tests do not rely on them, but they reflect the circularity discussed in Section 9.

An honest summary is therefore: *there is evidence that the slope is flatter at high income, but the evidence is sensitive to which countries are included, and the size of the interaction is small in terms of variance explained.* A pooled panel with year effects and country-clustered standard errors (6,447 country-years, 203 countries) gives a pooled elasticity of 0.97 (95% CI 0.89–1.05) and a joint interaction Wald test p = 4 × 10⁻⁵, consistent with heterogeneity; it is a sensitivity analysis, not a replacement for the primary test.

### 6.5 H5 – Have high-income countries cut emissions per person?

Among the 74 high-income countries with CO₂ in both years, emissions per person fell in 50 (68%) and rose in 24. The median change is **−15.9%** (95% bootstrap CI −24.8% to −9.6%). The Wilcoxon signed-rank test rejects "no decrease" (V = 938, raw p = 0.0077, Holm-adjusted p = 0.0155) with a matched-pairs rank-biserial correlation of −0.32 (95% CI −0.57 to −0.07): a moderate effect, with a confidence interval that approaches zero.

![Figure 6.4: H5 – CO₂ per person in 1990 (open circle) and 2022 (filled) for 74 high-income countries, sorted by change. Blue = fell, vermilion = rose.](../figures/fig_h5_paired_change.png)

The picture is heterogeneous. The largest proportional falls are in Luxembourg, Estonia, Nauru and Curaçao, with further large reductions in the United Kingdom, Lithuania and Denmark; the largest proportional increases are in Guyana, Oman, Panama, Seychelles and South Korea, and several other Gulf, Caribbean and Latin American economies also rose (Figure 6.4). The expected outcome in the plan (a slight decline with substantial heterogeneity) matches what was observed.

*Assumptions and robustness.* The log ratios are marginally non-normal (Shapiro p = 0.046), so the paired t-test is only a companion (t = −2.14, one-sided p = 0.018, d_z = −0.25, 95% CI −0.52 to −0.02); it is *less* significant than the Wilcoxon test because a few large increases inflate the variance. The sign test (p = 0.0017; proportion decreasing 0.68, 95% CI 0.56–0.78), a 10%-winsorised paired t-test (p = 0.004), excluding small countries (n = 49, median −20.4%, p = 0.0018) and excluding exporters (n = 64, median −17.0%, p = 0.0015) all agree. As noted in Section 4.6, the observed effect is just below the minimum the design could detect with 80% power, so the evidence is real but not overwhelming.

### 6.6 H6 – Is decoupling related to income group?

A country is classified as *absolutely decoupled* if its GDP per capita rose and its CO₂ per capita fell between 1990 and 2022. Decoupling occurred in 1 of 20 low-income countries (5.0%), 4 of 43 lower-middle-income (9.3%), 11 of 54 upper-middle-income (20.4%) and 43 of 67 high-income countries (64.2%). The chi-square test rejects independence (χ² = 52.06, df = 3, p = 2.9 × 10⁻¹¹; Holm-adjusted p = 8.7 × 10⁻¹¹; Monte Carlo p = 1 × 10⁻⁵) with a large effect (Cramér's V = 0.532, 95% CI 0.413–0.655). The standardised residuals show where independence fails: high-income countries decouple far more often than expected (+7.06) and the other three groups less often (−2.75, −3.65, −2.19 for Low, Lower-middle, Upper-middle).

![Figure 6.5: H6 – Observed counts (expected under independence in brackets) of decoupled and non-decoupled countries by income group; colour shows the standardised residual.](../figures/fig_h6_decoupling_heatmap.png)

A 2 × 2 Fisher exact test (High vs not High) gives an odds ratio of 11.3 (95% CI 5.3–22.4). PPP-based GDP in the decoupling definition (V = 0.560), excluding small countries (V = 0.602) and excluding exporters (V = 0.598) all strengthen the association. Note that 36% of high-income countries did *not* decouple.

### 6.7 H7 – Exploratory: structural drivers

*This analysis is exploratory and its p-values are not corrected for multiplicity.* Adding industry share of GDP and urbanisation to log GDP per capita (n = 194) raises R² from 0.730 to 0.780 (ΔR² = 0.050; bootstrap 95% CI 0.027–0.085); the joint HC3 Wald test gives F = 30.8 (p = 2.6 × 10⁻¹²). The contributions are very unequal. Industry share has a clear positive coefficient (0.0095 log₁₀ units per percentage point, 95% CI 0.0067–0.0124; partial R² = 0.173), meaning that each additional percentage point of industry in GDP is associated with about 2.2% more CO₂ per person at the same income. Urbanisation shows no evidence of an association once income and industry are controlled (coefficient 0.0019, 95% CI −0.0006 to 0.0044; partial R² = 0.013). Income remains the dominant term (coefficient 0.781, partial R² = 0.616). Diagnostics are clean (maximum VIF 1.84, Breusch–Pagan p = 0.55, residual Shapiro p = 0.77), and the conclusion holds with PPP GDP (industry partial R² = 0.089; urbanisation p = 0.42) and with 2022-only covariates (n = 191; industry partial R² = 0.180).

### 6.8 All effect sizes together

![Figure 6.6: Effect sizes with 95% bootstrap confidence intervals for each hypothesis (each on its own scale, because the metrics differ).](../figures/fig_effect_size_forest.png)

### 6.9 Robustness summary

| ID | Robustness variant | n | Statistic | Raw p | Effect size [95% CI] | Direction vs primary |
|---|---|---|---|---|---|---|
| H1 | permutation test | 198 | 0.855 | 1.00e-04 | 0.855 | agrees |
| H1 | PPP GDP | 195 | 0.869 | 4.64e-61 | 0.869 [0.814, 0.907] | agrees |
| H1 | excluding population < 1 million | 155 | 0.862 | 2.28e-47 | 0.862 [0.805, 0.900] | agrees |
| H1 | excluding oil/gas exporters | 175 | 0.884 | 2.64e-59 | 0.884 [0.834, 0.918] | agrees |
| H2 | permutation Kruskal-Wallis | 205 | 146 | 1.00e-04 | 0.717 | agrees |
| H2 | excluding population < 1 million | 159 | 117 | 3.51e-25 | 0.740 [0.661, 0.804] | agrees |
| H2 | excluding oil/gas exporters | 182 | 134 | 6.29e-29 | 0.742 [0.667, 0.807] | agrees |
| H3 | Mann-Whitney U | 146 | 5.28e+03 | 5.71e-25 | 0.983 [0.963, 0.997] | agrees |
| H3 | permutation of mean difference | 146 | 1.149 | 1.00e-04 | 1.149 | agrees |
| H3 | excluding population < 1 million | 115 | 16.416 | 1.12e-31 | 2.900 [2.508, 3.430] | agrees |
| H3 | excluding oil/gas exporters | 132 | 17.971 | 1.59e-33 | 3.074 [2.668, 3.630] | agrees |
| H4 | PPP GDP | 195 | 3.262 | 0.0227 | 0.011 [0.002, 0.029] | agrees |
| H4 | excluding population < 1 million | 155 | 2.355 | 0.0744 | 0.009 [0.002, 0.028] | agrees |
| H4 | excluding oil/gas exporters | 175 | 4.130 | 0.0074 | 0.011 [0.003, 0.028] | agrees |
| H4 | Cook's-distance outliers removed | 190 | 9.105 | 1.20e-05 | 0.018 [0.007, 0.036] | agrees |
| H4 | income-group label permutation (as planned) | 198 | 3.228 | 0.0517 | 0.010 | agrees |
| H4 | Freedman-Lane residual permutation (added) | 198 | 3.228 | 0.0244 | 0.010 | agrees |
| H5 | sign test | 74 | 50.000 | 0.0017 | 0.676 [0.557, 0.780] | agrees |
| H5 | 10%-winsorised paired t-test | 74 | -2.710 | 0.0042 | -0.315 | agrees |
| H5 | excluding population < 1 million | 49 | 325 | 0.0018 | -0.469 [-0.757, -0.153] | agrees |
| H5 | excluding oil/gas exporters | 64 | 595 | 0.0015 | -0.428 [-0.681, -0.161] | agrees |
| H6 | PPP-based GDP condition | 180 | 56.532 | 3.23e-12 | 0.560 [0.446, 0.679] | agrees |
| H6 | excluding population < 1 million | 148 | 53.607 | 1.36e-11 | 0.602 [0.477, 0.730] | agrees |
| H6 | excluding oil/gas exporters | 161 | 57.634 | 1.88e-12 | 0.598 [0.476, 0.722] | agrees |

### 6.10 Expected versus observed outcomes

| Hypothesis | Expected in the plan | Observed |
|------------|----------------------|----------|
| H1 | ρ around 0.8 or higher | ρ = 0.855 [0.802, 0.893] ✔ |
| H2 | Clear monotone ordering, large effect | All six pairs differ; ε² = 0.72 ✔ |
| H3 | Very large effect | d = 3.06 ✔ (larger than most would expect) |
| H4 | Some evidence of slope difference; least certain | Evidence present (Holm p = 0.0155) but sensitive to country inclusion; small ΔR² ✔ (as flagged) |
| H5 | Slightly negative median, heterogeneous, moderate effect | Median −15.9%, 50/74 fell; r = −0.32 ✔ |
| H6 | Dependence, decoupling concentrated in high income | 64% vs 5–20%; V = 0.53 ✔ |
| H7 | Modest added explanatory power | ΔR² = 0.05; industry yes, urbanisation no (partly ✔) |

### 6.11 Deviations from the plan

Three additions and one clarification were made after the plan was committed; none changes a primary test or decision. (1) A **Freedman–Lane permutation** was added for H4 after the planned label permutation turned out to test a stronger null. Both are reported. (2) The **Monte Carlo p-value for H6** was computed although the chi-square conditions were met, as an additional check. (3) The **sign-test confidence interval** was changed from one-sided to the two-sided 95% interval after the first run returned an upper limit of 1.0 by construction. (4) The plan's H2 robustness list mentions a PPP variant; PPP is irrelevant to H2 because it involves CO₂ only, so only the exclusion and permutation variants were run.

---

## 7. Quality Control and Validation

**Automated data checks.** The 25 pytest checks described in Section 3.5 all pass on the final, cleanly rebuilt data (`pytest` output: *25 passed*). They are regression guards: if a future raw download changes any count or checksum, the suite fails and signals that the analysis must be re-examined before results are reused.

**Assumption checks.** Every test is preceded by printed assumption checks (Section 4.2 and `results/assumption_checks.csv`): Shapiro–Wilk and QQ plots, Levene, Breusch–Pagan, Cook's distance, VIF, outlier counts, and minimum expected cell counts. Violations found: raw-scale non-normality everywhere; non-normal log₁₀ CO₂ in the Low-income group; unequal variances for H3; non-normal regression residuals for H4; borderline non-normal log ratios for H5; VIF 8.3 in the H4 main-effects model. The fallback or protection used in each case (rank tests, Welch, HC3, permutation tests, companion-only status of the paired t-test) is documented in Sections 4 and 6.

**Sensitivity analyses.** Every confirmatory hypothesis has at least three robustness variants chosen in advance: an alternative GDP price base (PPP), exclusion of small-population countries, exclusion of oil/gas exporters, permutation tests, winsorising or outlier removal. Only H4 shows a meaningful change (p = 0.074 when small countries are excluded).

**Cross-checking approaches.** Parametric and non-parametric tests were run side by side for H2, H3 and H5, and the pooled panel with clustered standard errors was run as an independent check of H1/H4 conclusions.

**Reproducibility.** Seeds are fixed (42) in every bootstrap and permutation; package versions are pinned in `requirements.txt`; raw data are cached and checksummed; the notebook and scripts run end to end from an empty `data/` folder; and the final results were checked to be *identical* between the first run and the clean rebuild (Appendix A).

---

## 8. Expected Challenges and Solutions

| # | Challenge | Type | Solution applied | Outcome |
|---|-----------|------|------------------|---------|
| 1 | Aggregates mixed with countries (OWID regions, World, shipping; World Bank income and region groups) | Anticipated | Drop OWID rows without ISO3 and World Bank rows with region "Aggregates" | 35 + 79 entities separated; pytest confirms none remain |
| 2 | Country names differ between sources | Anticipated | Join only on ISO3 + year; single crosswalk for Kosovo | Kosovo joined (33 rows); test asserts it |
| 3 | OWID appeared to contain 7,656 duplicate (iso_code, year) keys | Real | Cause: all aggregates share the missing ISO code; removing them first leaves zero duplicates | Guard in code and test |
| 4 | World Bank indicator files return five aggregate series with blank ISO3, which pandas reads as missing | Real | Read with `keep_default_na=False`; classify aggregates from the metadata endpoint | Rows removed and logged |
| 5 | Non-positive CO₂ per capita values (Timor-Leste 1994–97, Kosovo 1994–99) make logarithms undefined | Real | Excluded from log transforms only; stored values untouched | 10 rows; none in analysed years |
| 6 | Missing data for key variables (e.g. GDP missing for 16 entities in 2022; no 1990 data for several new states) | Anticipated | No imputation of outcome or main predictor; per-test sample sizes and exclusion lists reported | Samples from 74 to 205 |
| 7 | Skewed, heavy-tailed distributions | Anticipated | log₁₀ transforms plus rank-based primary tests | Skewness 2.58 → −0.58 |
| 8 | Unequal variances (H3) | Real | Welch instead of Student; Brown–Forsythe test | Levene p = 0.0002 |
| 9 | Normality violated in the Low-income group (H2), in H4 residuals, and marginally in H5 | Real | Kruskal–Wallis and Wilcoxon were primary by design; HC3 and permutation tests for H4 | Parametric tests shown only as comparisons |
| 10 | Collinearity between log GDP and income-group dummies (VIF 8.3) | Real | Reported; not used to drop terms because the F-test of the block is unaffected | Documented as a consequence of circular income groups |
| 11 | Planned label-permutation test for H4 tested a stronger null than the interaction null (p = 0.052) | Real | Added Freedman–Lane residual permutation (p = 0.024), reporting both | Evidence for H4 judged "present but sensitive" |
| 12 | Small expected counts in the H6 table | Anticipated | Planned Monte Carlo/Fisher fallback; not needed (minimum expected 6.41) but run as a check | Same conclusion |
| 13 | One-sided sign-test interval returned an upper bound of 1.0 | Real | Two-sided CI requested explicitly | Corrected |
| 14 | pandas deprecation warnings on `sum(1)` | Real | Switched to `sum(axis=1)` | Silent run |
| 15 | Module name `src/tests.py` clashes with the `tests/` folder | Real | pytest loads the module by file path with `importlib` | No import ambiguity |
| 16 | Long run time (about 4 min 17 s) from 10,000-fold bootstraps and 100,000 Monte Carlo permutations in Python loops | Real | Accepted; deterministic seeds make the cost one-off | Whole pipeline reproducible |
| 17 | Overlapping boxes and labels in the first versions of the decision-tree and assumption figures | Real | Widened layouts, reserved headroom for annotation boxes, re-rendered and inspected each figure | Fixed |
| 18 | Multiple testing across six hypotheses | Anticipated | Holm correction, plus BH as secondary | All six remain significant |
| 19 | Non-independence and ecological inference | Anticipated | Country as unit, never country-year; panel with clustered SEs as sensitivity; explicit caution in interpretation | See Section 9 |

---

## 9. Limitations and Threats to Validity

1. **Circularity of income grouping.** The World Bank assigns countries to income groups on the basis of income per person. Comparing groups (H2, H3) on *emissions* is not circular, but H4 and H6 use groups as a moderator or classifier of a relationship in which income is already a variable, so the group effects partly restate the income effect. The unusually high VIF in H4 (8.3) and the very large effect sizes in H3 (d = 3.06) are symptoms of this overlap. Group membership is also taken from the *current* World Bank classification (the 2026 download), not from the classification in 1990 or 2022 — for example Russia is currently classified as high income and Venezuela as lower-middle income — so each group mixes countries that moved between categories over the period.
2. **Ecological fallacy.** All statements are about countries. A country-level association between income and emissions (Robinson, 1950) says nothing about whether richer *individuals* emit more, and the within-group slope for high-income countries (≈ 0.2) does not mean that a richer household in such a country has no extra emissions.
3. **Correlation is not causation.** Income, industry structure, energy mix, climate, geography and policy are intertwined. H1–H7 describe association and decoupling in time; none identifies causes. The observation that high-income countries decoupled more often is not evidence that wealth *causes* decoupling.
4. **Non-independence of countries.** Countries share trade, technology, fuel markets and regional policy, so the effective sample size is smaller than the nominal one and p-values are optimistic. The sensitivity analyses reduce, but do not remove, this concern.
5. **Measurement and price base.** CO₂ is territorial (fossil fuels and industry, excluding land-use change): emissions embodied in imports are attributed to the producer, which flatters wealthy importers and penalises exporters. GDP per person in constant 2015 US$ and in PPP use different price bases and give different elasticities (in the exploratory H7 model the income coefficient is 0.78 with constant-US$ GDP and 1.03 with PPP GDP; H4's p-value moves from 0.013 to 0.023); results are therefore reported for both.
6. **Endpoint and definition sensitivity (H5, H6).** Comparing two years (1990 and 2022) is sensitive to the starting year: economic collapses and restructuring around 1990 affect the baseline of some countries, and the binary decoupling flag ignores the magnitude of changes and the path in between.
7. **Missing data are not random.** Countries lacking 2022 GDP or 1990 CO₂ are mostly small territories or states with disrupted statistics, so results apply to countries with reasonably complete statistical systems.
8. **Researcher degrees of freedom.** Holm covers the six confirmatory p-values; the robustness variants and H7 are uncorrected and serve to judge consistency, not to hunt for significance. Fixing tests, directions and subsets in advance is the main protection against selective reporting.
9. **Statistical versus substantive significance.** With effect sizes this large the p-values for H1, H2, H3 and H6 are essentially zero; effect sizes and intervals carry the interpretation.
10. **Limited power in H5.** The minimum detectable effect (d_z = 0.29) exceeds the observed paired effect (|d_z| = 0.25 under the t-test), so a real effect of this size could easily have been missed.

---

## 10. Time Allocation and Project Plan

The work was planned as 32 hours across eight phases, with an explicit troubleshooting buffer (12.5%) for the kinds of problems recorded in Section 8.

| Phase | Planned hours | Main deliverable | Actual hours |
|-------|--------------:|------------------|--------------|
| Planning and hypothesis design | 4 | `docs/analysis_plan.md`, committed before any test | *(author to enter)* |
| Data acquisition | 3 | `src/acquire.py`, raw data and manifest | *(author to enter)* |
| Cleaning and quality control | 5 | `src/clean.py`, cleaning log, pytest suite | *(author to enter)* |
| Assumption checks and testing | 7 | `src/tests.py`, `results/*.csv` | *(author to enter)* |
| Troubleshooting buffer | 4 | Fixes recorded in Section 8 | *(author to enter)* |
| Figures | 2 | 11 PNG figures | *(author to enter)* |
| Report writing | 6 | This report | *(author to enter)* |
| Review and viva preparation | 1 | `report_notes.md`, final checks | *(author to enter)* |
| **Total** | **32** | | |

![Figure 10.1: Planned allocation of the 32 project hours. The actual hours column of `docs/time_log.csv` is left for the author to fill; the figure shows actual bars automatically once values are entered.](../figures/fig_timeline_gantt.png)

For orientation, the measured computing times in this project were about 53 seconds for the data downloads and about 4 minutes 17 seconds for the complete hypothesis-testing run; the cleaning step and the pytest suite each complete in seconds. Computing time was therefore a small part of the budget; most of the time is consumed by design decisions, assumption checking, interpretation and documentation.

---

## 11. Conclusion

The framework was designed before it was run, and it was run as designed. The data support a clear, if unsurprising, answer to the first three research questions: across countries in 2022, richer countries emit far more CO₂ per person (ρ = 0.855), emissions rise steadily through the four income groups (ε² = 0.72), and high-income countries emit about 14 times more per person, in geometric-mean terms, than lower-middle- and low-income countries combined.

The more interesting results are the nuanced ones. The elasticity of emissions to income is flatter in high-income countries (about 0.2, against 0.7–0.9 elsewhere), but the evidence for this interaction is sensitive to which countries are included and explains only about 1% of the variance. Among high-income countries emissions per person fell in two thirds of cases, with a median reduction of about 16%, and absolute decoupling is dramatically more common in high-income countries (64%) than elsewhere (5–20%). Yet over a third of high-income countries did not decouple, and the evidence for a decline, while significant, is a moderate effect that the design only just had the power to detect. In the exploratory analysis, industrial share of GDP, not urbanisation, adds explanatory power beyond income.

The answer to "Does wealth burn carbon?" is therefore: across countries, yes, strongly; but among the wealthiest the link is weaker, and a majority of them have reduced emissions per person since 1990, which is consistent with decoupling being possible. The data cannot say *why* decoupling happened (policy, structural change or relocated production) nor whether it is fast enough for climate targets; those questions need consumption-based emissions and causal designs.

---

## 12. References

**Data sources (URLs fetched and verified on 2026-11-05)**

* Our World in Data. *CO₂ and Greenhouse Gas Emissions* data (`owid-co2-data.csv` and `owid-co2-codebook.csv`). GitHub repository `owid/co2-data`. https://github.com/owid/co2-data (file: https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv). The codebook lists the Global Carbon Budget (2025) as the source of the CO₂ series.
* World Bank. *Indicators API v2* (indicators NY.GDP.PCAP.KD, NY.GDP.PCAP.PP.KD, SP.POP.TOTL, SP.URB.TOTL.IN.ZS, NV.IND.TOTL.ZS) and *country metadata*. https://api.worldbank.org/v2/

**Methodological and software references**

* Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society, Series B*, 57(1), 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
* Breusch, T. S. and Pagan, A. R. (1979). A simple test for heteroscedasticity and random coefficient variation. *Econometrica*, 47(5), 1287–1294. doi:10.2307/1911963
* Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics*, 6(3), 241–252. doi:10.1080/00401706.1964.10490181
* Robinson, W. S. (1950). Ecological correlations and the behavior of individuals. *American Sociological Review*, 15(3), 351–357. doi:10.2307/2087176
* Seabold, S. and Perktold, J. (2010). Statsmodels: econometric and statistical modeling with Python. *Proceedings of the 9th Python in Science Conference*, 92–96. doi:10.25080/majora-92bf1922-011
* Shapiro, S. S. and Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). *Biometrika*, 52(3/4), 591–611. doi:10.2307/2333709
* Vallat, R. (2018). Pingouin: statistics in Python. *Journal of Open Source Software*, 3(31), 1026. doi:10.21105/joss.01026
* Welch, B. L. (1947). The generalization of "Student's" problem when several different population variances are involved. *Biometrika*, 34(1/2), 28–35. doi:10.2307/2332510
* White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. *Econometrica*, 48(4), 817–838. doi:10.2307/1912934
* Wilcoxon, F. (1945). Individual comparisons by ranking methods. *Biometrics Bulletin*, 1(6), 80–83. doi:10.2307/3001968

**Well-known references (†) – details cited from memory; verify before submitting (see `report_notes.md`)**

* † Efron, B. (1979). Bootstrap methods: another look at the jackknife. *The Annals of Statistics*, 7(1), 1–26.
* † Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics*, 6(2), 65–70.
* † Kruskal, W. H. and Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association*, 47(260), 583–621.
* † Freedman, D. and Lane, D. (1983). A nonstochastic interpretation of reported significance levels. *Journal of Business & Economic Statistics*, 1(4), 292–298.
* † Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.). Lawrence Erlbaum.
* † Tomczak, M. and Tomczak, E. (2014). The need to report effect size estimates revisited. *Trends in Sport Sciences*, 21(1), 19–25.
* † Virtanen, P. et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.

---

## Appendix A. Repository Structure and How to Reproduce

```
README.md                      overview and run instructions
requirements.txt               pinned package versions
data/raw/                      cached downloads + manifest.json (date, URL, SHA-256)
data/processed/                panel, cross-section, paired table, exclusions, row-level cleaning log
docs/analysis_plan.md          pre-specified plan (first git commit)
docs/time_log.csv              planned (and, once entered, actual) hours
src/acquire.py                 download and cache raw data
src/clean.py                   apply Week 1 rules, build analysis tables, sample flow
src/tests.py                   hypothesis tests H1-H7, power, panel sensitivity
src/make_figures.py            all figures
tests/test_data_integrity.py   25 pytest integrity checks
notebooks/02_hypothesis_testing.ipynb   executes the whole pipeline
results/                       CSV tables of every test and check
figures/                       PNG figures (200 dpi)
report/                        this report and the author's notes
```

To reproduce (Windows PowerShell or any shell):

```
python -m venv .venv
.venv\Scripts\activate              # (Linux/macOS: source .venv/bin/activate)
pip install -r requirements.txt
python src/acquire.py               # downloads once; cached afterwards
python src/clean.py
python src/tests.py                 # about 4 min; writes results/*.csv
python src/make_figures.py
pytest
jupyter nbconvert --to notebook --execute --inplace notebooks/02_hypothesis_testing.ipynb
```

## Appendix B. Full Results Table

Every test that was run, in the order of the hypotheses. "Primary" rows carry the Holm- and BH-adjusted p-values; other rows are companion, comparison, check or robustness results. Permutation p-values equal to 1 × 10⁻⁴ are the smallest value 10,000 permutations can resolve. The same data are in `results/hypothesis_results.csv`.

| ID | Variant | Test | n | Statistic | Raw p | Holm p | BH p | Effect size [95% CI] | Decision / note |
|---|---|---|---|---|---|---|---|---|---|
| H1 | primary | Spearman rho (one-sided, greater) | 198 | 0.855 | 4.98e-58 | 2.99e-57 | 2.99e-57 | Spearman rho: 0.855 [0.802, 0.893] | Reject H0 (supports H1: GDP pc and CO2 pc positively associated (rho>0)) |
| H1 | companion: Pearson on log10 values | Pearson r (one-sided, greater) | 198 | 0.852 | 3.06e-57 | – | – | Pearson r: 0.852 [0.808, 0.886] |  |
| H1 | robustness: permutation test | Spearman permutation (B=10,000) | 198 | 0.855 | 1.00e-04 | – | – | Spearman rho: 0.855 | direction agrees with primary; raw p < 0.05 |
| H1 | robustness: PPP GDP | Spearman rho (one-sided, greater) | 195 | 0.869 | 4.64e-61 | – | – | Spearman rho: 0.869 [0.814, 0.907] | direction agrees with primary; raw p < 0.05 |
| H1 | robustness: excluding population < 1 million | Spearman rho (one-sided, greater) | 155 | 0.862 | 2.28e-47 | – | – | Spearman rho: 0.862 [0.805, 0.900] | direction agrees with primary; raw p < 0.05 |
| H1 | robustness: excluding oil/gas exporters | Spearman rho (one-sided, greater) | 175 | 0.884 | 2.64e-59 | – | – | Spearman rho: 0.884 [0.834, 0.918] | direction agrees with primary; raw p < 0.05 |
| H2 | primary | Kruskal-Wallis H (df=3) | 205 | 146 | 1.65e-31 | 6.60e-31 | 3.30e-31 | epsilon-squared: 0.717 [0.644, 0.782] | Reject H0 (supports H1: CO2 pc differs across income groups) |
| H2 | comparison: one-way ANOVA on log10 | ANOVA F (classical) | 205 | 206 | 5.31e-61 | – | – | eta-squared: 0.754 |  |
| H2 | comparison: Welch ANOVA on log10 | Welch ANOVA F | 205 | 177 | 2.00e-35 | – | – | partial eta-squared: 0.754 |  |
| H2 | robustness: permutation Kruskal-Wallis | KW permutation (B=10,000) | 205 | 146 | 1.00e-04 | – | – | epsilon-squared: 0.717 | direction agrees with primary; raw p < 0.05 |
| H2 | robustness: excluding population < 1 million | Kruskal-Wallis H (df=3) | 159 | 117 | 3.51e-25 | – | – | epsilon-squared: 0.740 [0.661, 0.804] | direction agrees with primary; raw p < 0.05 |
| H2 | robustness: excluding oil/gas exporters | Kruskal-Wallis H (df=3) | 182 | 134 | 6.29e-29 | – | – | epsilon-squared: 0.742 [0.667, 0.807] | direction agrees with primary; raw p < 0.05 |
| H3 | primary | Welch t-test (df=119.3, one-sided greater) | 146 | 18.355 | 7.45e-37 | 3.73e-36 | 2.24e-36 | Cohen's d: 3.057 [2.680, 3.550] | Reject H0 (supports H1: High-income mean log10 CO2 pc exceeds Lower-middle+Low) |
| H3 | comparison: Student t (equal variance) | Student t (one-sided greater) | 146 | 18.464 | 4.26e-40 | – | – | Cohen's d: 3.057 |  |
| H3 | robustness: Mann-Whitney U | Mann-Whitney U (one-sided greater) | 146 | 5.28e+03 | 5.71e-25 | – | – | rank-biserial r: 0.983 [0.963, 0.997] | direction agrees with primary; raw p < 0.05 |
| H3 | robustness: permutation of mean difference | permutation test (B=10,000) | 146 | 1.149 | 1.00e-04 | – | – | mean diff (log10): 1.149 | direction agrees with primary; raw p < 0.05 |
| H3 | robustness: excluding population < 1 million | Welch t-test (df=110.2, one-sided greater) | 115 | 16.416 | 1.12e-31 | – | – | Cohen's d: 2.900 [2.508, 3.430] | direction agrees with primary; raw p < 0.05 |
| H3 | robustness: excluding oil/gas exporters | Welch t-test (df=100.9, one-sided greater) | 132 | 17.971 | 1.59e-33 | – | – | Cohen's d: 3.074 [2.668, 3.630] | direction agrees with primary; raw p < 0.05 |
| H4 | primary | Interaction F-test, HC3 Wald-type (anova_lm robust='hc3') | 198 | 3.691 | 0.0129 | 0.0155 | 0.0129 | delta R-squared (interaction block): 0.010 [0.003, 0.026] | Reject H0 (supports H1: income-emissions slope differs by income group) |
| H4 | robustness: PPP GDP | Interaction F-test, HC3 Wald-type (anova_lm robust='hc3') | 195 | 3.262 | 0.0227 | – | – | delta R-squared (interaction block): 0.011 [0.002, 0.029] | direction agrees with primary; raw p < 0.05 |
| H4 | robustness: excluding population < 1 million | Interaction F-test, HC3 Wald-type (anova_lm robust='hc3') | 155 | 2.355 | 0.0744 | – | – | delta R-squared (interaction block): 0.009 [0.002, 0.028] | direction agrees with primary; raw p >= 0.05 |
| H4 | robustness: excluding oil/gas exporters | Interaction F-test, HC3 Wald-type (anova_lm robust='hc3') | 175 | 4.130 | 0.0074 | – | – | delta R-squared (interaction block): 0.011 [0.003, 0.028] | direction agrees with primary; raw p < 0.05 |
| H4 | robustness: Cook's-distance outliers removed | Interaction F-test, HC3 Wald-type (anova_lm robust='hc3') | 190 | 9.105 | 1.20e-05 | – | – | delta R-squared (interaction block): 0.018 [0.007, 0.036] | direction agrees with primary; raw p < 0.05 |
| H4 | robustness: income-group label permutation (as planned) | permutation of group labels (B=10,000) | 198 | 3.228 | 0.0517 | – | – | delta R-squared: 0.010 | direction agrees with primary; raw p >= 0.05 |
| H4 | robustness: Freedman-Lane residual permutation (added) | Freedman-Lane permutation (B=10,000) | 198 | 3.228 | 0.0244 | – | – | delta R-squared: 0.010 | direction agrees with primary; raw p < 0.05 |
| H5 | primary | Wilcoxon signed-rank (one-sided, less) | 74 | 938 | 0.0077 | 0.0155 | 0.0093 | matched-pairs rank-biserial r: -0.324 [-0.568, -0.067] | Reject H0 (supports H1: CO2 pc of high-income countries fell 1990-2022) |
| H5 | companion: paired t-test on log ratios (normality VIOLATED - shown for transparency only) | paired t-test (one-sided, less) | 74 | -2.136 | 0.0180 | – | – | Cohen's d_z: -0.248 [-0.516, -0.015] |  |
| H5 | robustness: sign test | sign test (one-sided) | 74 | 50.000 | 0.0017 | – | – | proportion decreasing: 0.676 [0.557, 0.780] | direction agrees with primary; raw p < 0.05 |
| H5 | robustness: 10%-winsorised paired t-test | paired t-test on winsorised log ratios | 74 | -2.710 | 0.0042 | – | – | Cohen's d_z: -0.315 | direction agrees with primary; raw p < 0.05 |
| H5 | robustness: excluding population < 1 million | Wilcoxon signed-rank (one-sided, less) | 49 | 325 | 0.0018 | – | – | matched-pairs rank-biserial r: -0.469 [-0.757, -0.153] | direction agrees with primary; raw p < 0.05 |
| H5 | robustness: excluding oil/gas exporters | Wilcoxon signed-rank (one-sided, less) | 64 | 595 | 0.0015 | – | – | matched-pairs rank-biserial r: -0.428 [-0.681, -0.161] | direction agrees with primary; raw p < 0.05 |
| H6 | primary | Chi-square test of independence (df=3) | 184 | 52.064 | 2.90e-11 | 8.71e-11 | 4.35e-11 | Cramer's V: 0.532 [0.413, 0.655] | Reject H0 (supports H1: decoupling depends on income group) |
| H6 | check: 2x2 Fisher exact (High vs not High) | Fisher exact (two-sided) | 184 | 11.310 | 2.60e-12 | – | – | odds ratio: 11.310 [5.326, 22.398] |  |
| H6 | robustness: PPP-based GDP condition | Chi-square test of independence (df=3) | 180 | 56.532 | 3.23e-12 | – | – | Cramer's V: 0.560 [0.446, 0.679] | direction agrees with primary; raw p < 0.05 |
| H6 | robustness: excluding population < 1 million | Chi-square test of independence (df=3) | 148 | 53.607 | 1.36e-11 | – | – | Cramer's V: 0.602 [0.477, 0.730] | direction agrees with primary; raw p < 0.05 |
| H6 | robustness: excluding oil/gas exporters | Chi-square test of independence (df=3) | 161 | 57.634 | 1.88e-12 | – | – | Cramer's V: 0.598 [0.476, 0.722] | direction agrees with primary; raw p < 0.05 |
| H7 | primary (exploratory) | Joint HC3 Wald F: industry share & urbanisation (exploratory) | 194 | 30.817 | 2.56e-12 | – | – | increment in R-squared: 0.050 [0.027, 0.085] | exploratory (unadjusted) |
| H7 | robustness: PPP GDP | Joint HC3 Wald F (exploratory) | 192 | 8.860 | 0.0002 | – | – | see h7_coefficients.csv: – | exploratory (unadjusted) |
| H7 | robustness: 2022-only covariates | Joint HC3 Wald F (exploratory) | 191 | 30.842 | 2.66e-12 | – | – | see h7_coefficients.csv: – | exploratory (unadjusted) |
