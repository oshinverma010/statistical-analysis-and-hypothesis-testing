# Analysis Plan (pre-specification) — "Does Wealth Burn Carbon?" Week 2

**Status:** written and committed *before* any hypothesis test was run and before any data file was opened.
The only things done before this commit were: creating the project skeleton, installing packages, and confirming that the
two data URLs respond (HTTP 200 / JSON header). No summary statistic of the data had been looked at.

Project: CO2 emissions vs economic development, 1990-2022 (continuation of the Week 1 project).

## 1. Research questions

| ID | Research question |
|----|-------------------|
| RQ1 | Is income associated with emissions per person across countries? |
| RQ2 | Do income groups differ in emissions per person? |
| RQ3 | Does the income-emissions slope differ by income group? |
| RQ4 | Have wealthy countries decoupled growth from emissions? |
| RQ5 | (exploratory) After controlling for income, do industrial share and urbanisation explain emissions per person? |

## 2. Data and definitions (fixed in advance)

* **CO2 per capita** = OWID `co2_per_capita` (production-based CO2, tonnes per person). Source: Our World in Data CO2 dataset.
* **GDP per capita (primary)** = World Bank `NY.GDP.PCAP.KD` (constant 2015 US$). **Sensitivity:** `NY.GDP.PCAP.PP.KD` (PPP, constant international $).
* **Other covariates:** `SP.POP.TOTL`, `SP.URB.TOTL.IN.ZS` (urban % of population), `NV.IND.TOTL.ZS` (industry value added, % of GDP).
* **Income group / region:** World Bank country-metadata endpoint (current classification; it is *not* year-specific — discussed as a limitation).
  Four groups: Low, Lower-middle, Upper-middle, High. Countries "Not classified" are excluded from group-based tests.
* **Country definition:** OWID rows with a non-missing ISO3 code that does not start with `OWID_` (Kosovo, `OWID_KOS`, is mapped by a manual crosswalk to `XKX`);
  World Bank rows whose region is "Aggregates" are excluded. Join on ISO3 + year only, never on names.
* **Window:** 1990-2022 inclusive. **Reference years:** 2022 (cross-section) and 1990 vs 2022 (paired changes).
* **Missing data:** the outcome (CO2 pc) and the main predictor (GDP pc) are never imputed. A country enters a test only if all the variables
  *that test needs* are observed in the required year(s).
  Only for the exploratory H7 covariates (industry share, urbanisation) the latest observation within 2018-2022 may be used when 2022 is missing
  (documented in the cleaning log; complete-case 2022-only is the sensitivity).
* **Extremes:** legitimate extreme values (e.g. Qatar) are retained; log10 transformations are used.
* **Unit of analysis:** the country (never the country-year). H1-H4, H7: the 2022 country cross-section. H5-H6: per-country change 1990 -> 2022.
* **Randomness:** all bootstrap / permutation procedures use a fixed seed (42); B = 10,000 resamples unless noted.

## 3. Hypotheses

Notation: `y = log10(CO2 per capita)`, `x = log10(GDP per capita)`. alpha = 0.05.

| ID | RQ | H0 | H1 | Direction | Primary test (pre-specified) | Unit / sample | Why this test |
|----|----|----|----|-----------|------------------------------|---------------|---------------|
| H1 | RQ1 | Spearman rho = 0 between GDP pc and CO2 pc | rho > 0 | one-sided (greater) | Spearman rank correlation, 2022 | country, 2022 | Monotone association without assuming linearity or normality; distributions are heavy-tailed. |
| H2 | RQ2 | CO2 pc (2022) has the same distribution in all four income groups | at least one group differs | two-sided | Kruskal-Wallis; Dunn post hoc (Holm) | country, 2022, by group | k = 4 independent groups, skewed data; rank-based, robust to outliers. |
| H3 | RQ2 | mean log10 CO2 pc (High) = mean (Lower-middle + Low) | High > Lower-middle + Low | one-sided (greater) | Welch's t-test on log10 values | country, 2022 | Two independent groups, unequal sizes and probably unequal variances: Welch does not assume equal variance. |
| H4 | RQ3 | all income-group x log-GDP interaction coefficients = 0 | at least one interaction != 0 | F-test (upper tail) | OLS `y ~ x * C(income_group)`, `anova_lm` (type II) with HC3 robust SEs for coefficients | country, 2022 | Tests slope heterogeneity directly in one model rather than four separate regressions. |
| H5 | RQ4 | median of paired log ratio ln(CO2pc_2022 / CO2pc_1990) = 0 among high-income countries | median < 0 | one-sided (less) | Wilcoxon signed-rank on paired log ratios | high-income countries with CO2 in both years | Paired design; log ratios may be non-normal; Wilcoxon needs only symmetry. |
| H6 | RQ4 | absolute decoupling (GDP pc up AND CO2 pc down, 1990-2022) is independent of income group | not independent | two-sided | Chi-square test of independence | countries with GDP and CO2 in both years and an income group | Two categorical variables. |
| H7 | RQ5 | partial coefficients of industry share and urbanisation on y (controlling for x) = 0 | at least one != 0 | two-sided | Multiple OLS with HC3 SEs; partial R^2 | country, 2022 | **Exploratory.** |

### Test-selection rules (decided now, applied mechanically)

* **H1:** Spearman is primary regardless of normality. Pearson on log10 values is reported alongside (it is not part of the multiplicity family).
  95% CI by percentile bootstrap over countries. Robustness: PPP GDP; excluding countries with population < 1 million; excluding the pre-listed oil/gas exporters;
  permutation test for rho.
* **H2:** Kruskal-Wallis is primary (pre-specified; the parametric test is a comparison, not a replacement). Assumption checks are printed: Shapiro-Wilk per group on raw and log10 values,
  Levene (median-centred, Brown-Forsythe). Parametric comparison: one-way ANOVA on log10 values, or Welch ANOVA if Levene p < 0.05. Post hoc: Dunn with Holm adjustment (six pairs),
  with Mann-Whitney + Holm as a cross-check. Effect size: epsilon-squared, bootstrap 95% CI. Robustness: permutation Kruskal-Wallis; PPP-independent (CO2 only), so exclusions of small countries and oil/gas exporters.
* **H3:** Welch's t-test on log10 CO2 pc is primary. Mann-Whitney U is the robustness check. Cohen's d (pooled-SD version) and Hedges' g with bootstrap 95% CI. Also the CI of the mean difference (Welch-Satterthwaite).
  Robustness: permutation test of the mean difference; excluding small countries; excluding oil/gas exporters.
* **H4:** Model `y ~ x * C(income_group)`, reference group = High income. Diagnostics: Breusch-Pagan, Cook's distance (flag > 4/n), VIF (on the main-effects model; interaction VIFs are inflated by construction and are reported only as such),
  Shapiro-Wilk on residuals. The F-test for the interaction block is computed with the HC3 covariance (Wald F) as the primary p-value, and the classical F (anova_lm) is reported alongside. Robustness: PPP GDP; excluding small countries; removing Cook's-distance outliers;
  permutation of the income-group labels (interaction F). Effect size: increase in R^2 from adding the interaction block; per-group slopes with 95% CI.
* **H5:** Wilcoxon signed-rank (alternative = "less") on ln(CO2pc_2022 / CO2pc_1990). A paired t-test on the same log ratios is reported only if Shapiro-Wilk on the differences has p > 0.05; otherwise it is reported as "not appropriate" but still shown for transparency.
  Effect size: matched-pairs rank-biserial correlation and Cohen's d_z with bootstrap CI; the median % change with bootstrap CI. Robustness: sign test; excluding small countries; excluding oil/gas exporters; 10%-winsorised paired t-test.
* **H6:** Decoupling indicator (1/0): GDP pc (constant 2015 US$) 2022 > 1990 AND CO2 pc 2022 < 1990. Chi-square test of independence (no continuity correction for r x c).
  **If any expected cell count < 5**, the primary p-value is the Monte Carlo permutation p-value (100,000 permutations; a Fisher-Freeman-Halton analogue) and the chi-square p-value is reported as secondary; a 2x2 Fisher exact test (High vs non-High) is added as a check.
  Effect size: Cramer's V (bootstrap CI); standardised (adjusted) residuals per cell. Robustness: PPP-based GDP condition; excluding small countries.
* **Normality rule used for the *comparison* parametric tests:** Shapiro-Wilk p > 0.05 *and* an acceptable QQ plot. A violated assumption never causes a switch away from a pre-specified primary test; it causes a documented fallback only for the *parametric comparison*.
  Because of n > 30 per group in most groups the t-test and ANOVA are expected to be approximately robust (central limit theorem), which is why a Welch test on log10 values is acceptable as a primary test for H3.

### Expected outcomes (written before seeing data)

* H1: strong positive monotone association (rho around 0.8 or higher) — emissions rise with income across countries.
* H2: significant differences; clear ordering Low < Lower-middle < Upper-middle < High, with large effect size.
* H3: High-income mean log10 CO2 pc clearly greater than that of Lower-middle + Low (very large effect).
* H4: some evidence that the slope is not the same across groups (flatter slope in High-income countries) — least certain hypothesis.
* H5: median change slightly negative (decline) among high-income countries, but with substantial heterogeneity; moderate effect.
* H6: dependence between income group and decoupling is expected: absolute decoupling mostly among high-income countries; poorer countries mostly "coupled" (both up).
* H7: industry share and urbanisation add only modest explanatory power beyond income.

## 4. Decision rules and multiplicity

* **alpha = 0.05.** Two-sided unless a direction is stated in the table above (H1, H3, H5 are one-sided).
* **Confirmatory family = H1-H6.** One *primary* p-value per hypothesis (H1 Spearman; H2 Kruskal-Wallis; H3 Welch; H4 interaction Wald F with HC3; H5 Wilcoxon; H6 chi-square / Monte Carlo as defined above).
  **Holm-Bonferroni correction** over these six p-values (family-wise error rate controlled at 0.05). Benjamini-Hochberg adjusted values are shown as a secondary column only.
* **H7 is exploratory** and is reported separately, with unadjusted p-values, clearly labelled; it is not part of the Holm family.
  Post hoc comparisons within H2 use their own Holm adjustment (six pairwise comparisons) and are not part of the family of six.
* Every result is reported with: n, test statistic, raw p, Holm-adjusted p (for H1-H6), an effect size, and a 95% confidence interval.
* **Wording rule:** a non-significant result is reported as "no evidence of ..." and never as "no effect"; a p-value is never described as the probability that H0 is true.
* **Decision:** reject H0 for a confirmatory hypothesis iff the Holm-adjusted p < 0.05. Robustness variants are judged on direction and effect size agreement with the primary result, not by re-applying a threshold.
* All hypotheses stay in the report whatever the result; no hypothesis, test or exclusion is added or removed after seeing results. Any deviation will be listed explicitly under "Deviations from the plan" in the report.

## 5. Power analysis (planned)

* H3 and H5: using the *real* sample sizes after cleaning, `statsmodels` power classes compute the minimum detectable standardised effect (Cohen's d / d_z) at 80% power, alpha = 0.05, one-sided,
  and a power curve against effect size. For H5 the Wilcoxon test's power is approximated by the paired t-test power with an asymptotic relative efficiency correction of 0.955 (documented approximation).

## 6. Threats to validity acknowledged in advance

* World Bank income group is itself defined from gross national income per capita, so comparing groups on GDP is partly circular. H2/H3 compare groups on *emissions*, not GDP; H4 and H6 use group as a moderator/classifier, which is acceptable but inherits this overlap.
* Countries are not independent draws (trade, shared energy systems, spatial dependence), so p-values are somewhat optimistic.
* Country-level associations are not individual-level (ecological fallacy) and not causal.
* Production-based emissions ignore emissions embodied in trade; constant-US$ GDP and PPP GDP have different price bases.
* Country-years are *not* treated as independent observations. Optional sensitivity: pooled panel OLS 1990-2022 with year effects and country-clustered SEs.

## 7. Pre-listed oil/gas exporter set for robustness (ISO3)

QAT, KWT, ARE, SAU, BHR, OMN, BRN, TTO, LBY, IRQ, IRN, DZA, AZE, KAZ, TKM, NGA, AGO, GNQ, GAB, COG, VEN, RUS, NOR.
(Chosen from general knowledge of large hydrocarbon-export dependence, fixed before looking at the data; it is a sensitivity filter, not a classification claim.)
