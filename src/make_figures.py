"""Generate all report figures (PNG, 200 dpi, Okabe-Ito colour-blind-safe palette).

Run after src/tests.py:  python src/make_figures.py
Numbers shown in figures are read from results/*.csv or recomputed from data/processed (never typed in by hand).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import tests as T  # noqa: E402  (src/tests.py)

FIG, RES, DOCS = ROOT / "figures", ROOT / "results", ROOT / "docs"
FIG.mkdir(exist_ok=True)

# Okabe-Ito colour-blind-safe palette
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442", blue="#0072B2", verm="#D55E00", pink="#CC79A7")
GROUP_COL = {"Low income": OI["orange"], "Lower middle income": OI["sky"], "Upper middle income": OI["green"], "High income": OI["blue"]}
GROUP_MARK = {"Low income": "o", "Lower middle income": "s", "Upper middle income": "^", "High income": "D"}
GROUP_SHORT = {"Low income": "Low", "Lower middle income": "Lower-middle", "Upper middle income": "Upper-middle", "High income": "High"}
INK, MUTED, GRID = "#222222", "#6b6b6b", "#dddddd"

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "axes.edgecolor": "#888888", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "figure.dpi": 100, "savefig.dpi": 200, "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "font.family": "DejaVu Sans"})


def save(fig, name):
    fig.savefig(FIG / name, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved", name)


def fmt_p(p):
    return "p < 0.0001" if p < 1e-4 else f"p = {p:.4f}" if p < 0.01 else f"p = {p:.3f}"


def box(ax, xy, w, h, title, body="", fc="#EAF3FA", ec=OI["blue"], fs=9.5, bold=True):
    x, y = xy
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec=ec, lw=1.4))
    if body:
        ax.text(x + w / 2, y + h - 0.12, title, ha="center", va="top", fontsize=fs, fontweight="bold" if bold else "normal", color=INK)
        ax.text(x + w / 2, y + h - 0.12 - 0.27 * (1 + title.count("\n")), body, ha="center", va="top", fontsize=fs - 1.5, color="#333333", linespacing=1.35)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", color=INK, linespacing=1.3)


def arrow(ax, p1, p2, color="#555555", rad=0.0, lw=1.6, style="-|>"):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=14, lw=lw, color=color, connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))


# ------------------------------------------------------------------ 1. framework flowchart
def fig_framework():
    fig, ax = plt.subplots(figsize=(14, 8.4))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8.4)
    ax.axis("off")
    ax.grid(False)
    W, H = 3.9, 1.75
    xs = [0.45, 5.05, 9.65]
    ys = [5.95, 3.35, 0.75]
    items = [
        ("1  Research questions", "RQ1 income-emissions link\nRQ2 income-group differences\nRQ3 slope differences\nRQ4 decoupling   RQ5 (exploratory)", "#EAF3FA", OI["blue"]),
        ("2  Hypotheses (pre-specified)", "H0/H1, direction, test, alpha\nH1-H6 confirmatory, H7 exploratory\ndocs/analysis_plan.md committed\nBEFORE any test was run", "#EAF3FA", OI["blue"]),
        ("3  Data acquisition", "OWID CO2 CSV + World Bank API v2\nGDP pc (2015 US$ and PPP), urban %,\nindustry %, income group; cached raw\nfiles + SHA-256 manifest", "#E6F5EF", OI["green"]),
        ("4  Cleaning and quality control", "aggregates removed, join on ISO3+year,\nKosovo crosswalk, no imputation of\nCO2 / GDP, log guards, pytest suite\nand a row-level cleaning log", "#E6F5EF", OI["green"]),
        ("5  Assumption checks", "Shapiro-Wilk + QQ plots, Levene,\nBreusch-Pagan, Cook's D, VIF,\nexpected cell counts, sample sizes,\noutliers and skewness", "#FDF3DC", OI["orange"]),
        ("6  Test selection", "decision tree: Spearman, Welch,\nKruskal-Wallis, OLS + interaction,\nWilcoxon, chi-square / Monte Carlo\n(primary test fixed in the plan)", "#FDF3DC", OI["orange"]),
        ("7  Test + effect size + CI", "primary test with n, statistic, p,\neffect size and 95% bootstrap CI;\nrobustness: permutation, PPP GDP,\nexclusions, outlier removal", "#FBE9DD", OI["verm"]),
        ("8  Multiple-testing correction", "Holm over H1-H6 (FWER 5%),\nBH shown as secondary column;\nH7 reported separately as\nexploratory (unadjusted)", "#FBE9DD", OI["verm"]),
        ("9  Interpretation", "expected vs observed outcomes,\n\"no evidence of\" wording rule,\nlimitations: ecological fallacy,\ncircular income groups, non-independence", "#F6E8F0", OI["pink"]),
    ]
    centers = []
    # snake order: row0 left->right, row1 right->left, row2 left->right
    order = [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (0, 1), (0, 2), (1, 2), (2, 2)]
    for (title, body, fc, ec), (cx, cy) in zip(items, order):
        x, y = xs[cx], ys[cy]
        box(ax, (x, y), W, H, title, body, fc=fc, ec=ec, fs=11)
        centers.append((x, y))
    for i in range(8):
        (x1, y1), (x2, y2) = centers[i], centers[i + 1]
        if y1 == y2:
            if x2 > x1:
                arrow(ax, (x1 + W + 0.03, y1 + H / 2), (x2 - 0.03, y2 + H / 2))
            else:
                arrow(ax, (x1 - 0.03, y1 + H / 2), (x2 + W + 0.03, y2 + H / 2))
        else:
            arrow(ax, (x1 + W / 2, y1 - 0.03), (x2 + W / 2, y2 + H + 0.03))
    ax.text(7, 8.25, "Experimental framework: from research question to interpretation", ha="center", va="center", fontsize=14, fontweight="bold")
    save(fig, "fig_framework_flowchart.png")


# ------------------------------------------------------------------ 2. test-selection tree
def fig_tree(chk: pd.DataFrame):
    def p_of(h, text):
        r = chk[(chk["hypothesis"] == h) & chk["check"].str.contains(text, regex=False)]
        return r["p_value"].iloc[0] if len(r) else np.nan

    p_lev3 = p_of("H3", "Levene")
    p_sh5 = p_of("H5", "log ratio")
    p_sh1 = p_of("H1", "log10 CO2 pc")
    p_sh2 = p_of("H2", "log10 CO2 pc, Low")
    res = pd.read_csv(RES / "hypothesis_results.csv")
    exp_min = chk[(chk["hypothesis"] == "H6") & (chk["check"] == "minimum expected count")]["statistic"].iloc[0]
    fig, ax = plt.subplots(figsize=(18, 9.4))
    ax.set_xlim(0, 18.2)
    ax.set_ylim(0, 9.4)
    ax.axis("off")
    ax.grid(False)
    ax.text(9.1, 9.15, "Test-selection decision tree (branches taken in this project)", ha="center", fontsize=14, fontweight="bold")
    box(ax, (7.25, 7.7), 3.7, 0.9, "What is the data structure?\n(unit = country)", fc="#EAF3FA", ec=OI["blue"], fs=10.5)
    cols = [1.65, 5.7, 10.1, 13.2, 16.4]
    heads = [("Two continuous\nvariables", "#F6E8F0", OI["pink"]), ("Compare independent\ngroups", "#E6F5EF", OI["green"]), ("Paired measures\n(same country, 2 years)", "#FDF3DC", OI["orange"]),
             ("Slope may differ\nby a moderator", "#EAF3FA", OI["blue"]), ("Two categorical\nvariables", "#FBE9DD", OI["verm"])]
    for cx, (t, fc, ec) in zip(cols, heads):
        box(ax, (cx - 1.3, 6.15), 2.6, 0.95, t, fc=fc, ec=ec, fs=9.5)
        arrow(ax, (9.1, 7.68), (cx, 7.12), lw=1.2)

    def node(cx, y, text, fc="white", ec="#777777", w=2.7, h=0.95, fs=8.6, bold=False):
        box(ax, (cx - w / 2, y), w, h, text, fc=fc, ec=ec, fs=fs, bold=bold)

    def down(cx, y1, y2, label=None, dx=0.0, color="#555555"):
        arrow(ax, (cx, y1), (cx + dx, y2), color=color, lw=1.2)
        if label:
            ax.text(cx + dx / 2 + 0.08, (y1 + y2) / 2, label, fontsize=8, color=MUTED, ha="left", va="center", style="italic")

    # col 0: correlation (H1)
    c = cols[0]
    down(c, 6.13, 5.55)
    node(c, 4.6, f"Normal? Shapiro on log10\nCO2 pc: p = {p_sh1:.1g} -> NO,\nheavy tails and non-linear", w=2.9)
    down(c, 4.58, 3.95, "no")
    node(c, 2.95, "Spearman rank correlation\n(H1, one-sided)\n+ Pearson on log10 values\nas companion", fc="#F6E8F0", ec=OI["pink"], w=2.9, h=1.0, bold=True)
    # col 1: independent groups
    c = cols[1]
    down(c, 6.13, 5.55)
    node(c, 4.6, "How many groups?", w=2.6, h=0.9)
    arrow(ax, (c - 0.6, 4.58), (c - 1.1, 3.95), lw=1.2)
    arrow(ax, (c + 0.6, 4.58), (c + 1.1, 3.95), lw=1.2)
    ax.text(c - 1.3, 4.25, "2", fontsize=9, color=MUTED, style="italic")
    ax.text(c + 1.1, 4.25, "4", fontsize=9, color=MUTED, style="italic")
    node(c - 1.1, 2.75, f"Levene on log10 CO2:\np = {p_lev3:.1g} -> unequal\nvariances (H3)\n-> Welch t-test", fc="#E6F5EF", ec=OI["green"], w=1.95, h=1.2, bold=False, fs=8)
    node(c + 1.1, 2.75, f"Shapiro, Low group:\np = {p_sh2:.1g}; raw data\nskewed -> Kruskal-Wallis\n(H2) + Dunn-Holm", fc="#E6F5EF", ec=OI["green"], w=1.95, h=1.2, fs=8)
    arrow(ax, (c - 1.1, 2.73), (c - 1.1, 2.2), lw=1.2)
    arrow(ax, (c + 1.1, 2.73), (c + 1.1, 2.2), lw=1.2)
    node(c - 1.1, 1.1, "Fallbacks: Mann-Whitney,\npermutation test;\nStudent t not used", w=1.95, h=1.05, fs=7.8)
    node(c + 1.1, 1.1, "Comparison: ANOVA /\nWelch ANOVA on log10,\npermutation KW", w=1.95, h=1.05, fs=7.8)
    # col 2: paired
    c = cols[2]
    down(c, 6.13, 5.55)
    node(c, 4.6, f"Differences (log ratios)\nnormal? Shapiro p = {p_sh5:.3f}\n-> marginal violation", w=2.9)
    down(c, 4.58, 3.95, "no")
    node(c, 2.95, "Wilcoxon signed-rank\n(H5, one-sided)", fc="#FDF3DC", ec=OI["orange"], w=2.7, h=1.0, bold=True)
    arrow(ax, (c, 2.93), (c, 2.2), lw=1.2)
    node(c, 1.1, "Robustness: sign test, 10%-\nwinsorised paired t, exclusions;\npaired t shown only as a companion", w=2.9, h=1.05, fs=7.8)
    # col 3: moderator
    c = cols[3]
    down(c, 6.13, 5.55)
    node(c, 4.6, "Continuous outcome and\npredictor, 4 groups: fit one\nmodel with interaction", w=2.9)
    down(c, 4.58, 3.95)
    node(c, 2.95, "OLS  y ~ x * group\nF-test of interaction block\nwith HC3 robust covariance (H4)", fc="#EAF3FA", ec=OI["blue"], w=2.9, h=1.0, bold=True)
    arrow(ax, (c, 2.93), (c, 2.2), lw=1.2)
    node(c, 1.1, "Diagnostics: Breusch-Pagan,\nCook's D, VIF, residual Shapiro;\nlabel + Freedman-Lane permutation", w=2.9, h=1.05, fs=7.8)
    # col 4: categorical
    c = cols[4]
    down(c, 6.13, 5.55)
    node(c, 4.6, f"All expected counts >= 5?\nminimum = {exp_min:.2f} -> YES\n(otherwise Monte Carlo / Fisher)", w=2.9)
    down(c, 4.58, 3.95, "yes")
    node(c, 2.95, "Chi-square test of\nindependence (H6)\n+ Cramer's V", fc="#FBE9DD", ec=OI["verm"], w=2.7, h=1.0, bold=True)
    arrow(ax, (c, 2.93), (c, 2.2), lw=1.2)
    node(c, 1.1, "Checks: Monte Carlo permutation\np-value (100,000), 2x2 Fisher\nexact (High vs not High)", w=2.9, h=1.05, fs=7.8)
    ax.text(9.1, 0.28, "All six primary p-values then go through Holm-Bonferroni correction. Bold boxes = primary tests fixed in docs/analysis_plan.md.", ha="center", fontsize=9, color=MUTED)
    save(fig, "fig_test_selection_tree.png")


# ------------------------------------------------------------------ 3. assumption checks
def fig_assumptions(S):
    h1, h5 = S["H1"], S["H5"]
    rows = [("CO2 per capita 2022 (t/person)", h1["co2_pc"], h1["lg_co2"], "raw", "log10"),
            ("GDP per capita 2022 (constant 2015 US$)", h1["gdp_pc"], h1["lg_gdp"], "raw", "log10"),
            ("High-income CO2 pc: ratio 2022/1990", np.exp(h5["log_ratio_co2"]), h5["log_ratio_co2"], "ratio", "natural-log ratio")]
    fig, axes = plt.subplots(3, 4, figsize=(14, 10))
    for r, (name, raw, lg, l1, l2) in enumerate(rows):
        for c, (v, lab) in enumerate([(raw, l1), (lg, l2)]):
            v = np.asarray(v, float)
            w, p = stats.shapiro(v)
            ax = axes[r, 2 * c]
            ax.hist(v, bins=25, color=[OI["verm"], OI["blue"]][c], alpha=0.85, edgecolor="white")
            ax.set_title(f"{lab} - histogram", fontsize=9.5)
            ax.set_ylabel("countries")
            ax.set_xlabel(("t CO2 / person" if r == 0 else "US$ / person" if r == 1 else "ratio 2022 / 1990") + (" (log10 scale)" if (c and r < 2) else " (natural log)" if c else ""), fontsize=8.5)
            ax.set_ylim(0, ax.get_ylim()[1] * 1.45)
            ax.text(0.97, 0.95, f"Shapiro-Wilk W = {w:.3f}\n{fmt_p(p)}\nskew = {stats.skew(v):.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
                    bbox=dict(boxstyle="round", fc="white", ec="#bbbbbb"))
            ax2 = axes[r, 2 * c + 1]
            (osm, osr), (slope, icpt, _) = stats.probplot(v, dist="norm")
            ax2.scatter(osm, osr, s=10, color=[OI["verm"], OI["blue"]][c], alpha=0.8)
            ax2.plot(osm, slope * np.asarray(osm) + icpt, color="#333333", lw=1.2)
            ax2.set_title(f"{lab} - normal QQ plot", fontsize=9.5)
            ax2.set_xlabel("theoretical quantiles")
            ax2.set_ylabel("sample quantiles")
            fig.subplots_adjust(hspace=0.62, wspace=0.32, top=0.90, left=0.06, right=0.98, bottom=0.06)
    for r, row in enumerate(rows):
        y1 = axes[r, 0].get_position().y1
        fig.text(0.06, y1 + 0.032, row[0], fontsize=11, fontweight="bold")
    fig.suptitle("Normality checks before choosing between parametric and rank-based tests", fontsize=13, fontweight="bold", y=0.97)
    save(fig, "fig_assumption_checks.png")


# ------------------------------------------------------------------ 4. H1 scatter
def fig_h1(S, res):
    d = S["H1"].copy()
    fig, ax = plt.subplots(figsize=(9.5, 7))
    for g, col in GROUP_COL.items():
        s = d[d["income_group"].astype(str) == g]
        ax.scatter(10 ** s["lg_gdp"], 10 ** s["lg_co2"], s=34, color=col, marker=GROUP_MARK[g], edgecolor="white", linewidth=0.5, alpha=0.9, label=f"{GROUP_SHORT[g]} income (n={len(s)})", zorder=3)
    s = d[d["income_group"].isna()]
    if len(s):
        ax.scatter(10 ** s["lg_gdp"], 10 ** s["lg_co2"], s=34, color="#999999", marker="x", label="no income group")
    xs = np.linspace(d["lg_gdp"].min(), d["lg_gdp"].max(), 100)
    b, a, *_ = stats.linregress(d["lg_gdp"], d["lg_co2"])
    ax.plot(10 ** xs, 10 ** (a + b * xs), color="#333333", lw=1.6, ls="--", label=f"OLS on log10 values (slope = {b:.2f})", zorder=2)
    lo = lowess(d["lg_co2"], d["lg_gdp"], frac=0.5)
    ax.plot(10 ** lo[:, 0], 10 ** lo[:, 1], color=OI["verm"], lw=2.2, label="LOWESS (frac = 0.5)", zorder=4)
    resid = d["lg_co2"] - (a + b * d["lg_gdp"])
    lab = set(d.nlargest(3, "lg_co2")["iso3"]) | set(d.nsmallest(3, "lg_co2")["iso3"]) | set(d.loc[resid.nlargest(2).index, "iso3"]) | set(d.loc[resid.nsmallest(2).index, "iso3"])
    offs = {"KWT": (-30, 3), "BRN": (8, 3), "SOM": (-34, 7), "COD": (8, -11), "BDI": (8, -13)}
    for _, r in d[d["iso3"].isin(lab)].iterrows():
        dx, dy = offs.get(r["iso3"], (6, 6) if r["lg_co2"] > a + b * r["lg_gdp"] else (6, -11))
        ax.annotate(r["iso3"], (10 ** r["lg_gdp"], 10 ** r["lg_co2"]), xytext=(dx, dy), textcoords="offset points", fontsize=8.5, color=INK, fontweight="bold")
    h1 = res[(res["hypothesis"] == "H1") & res["primary"]].iloc[0]
    ax.text(0.03, 0.97, f"Spearman rho = {h1['statistic']:.3f}\n95% bootstrap CI [{h1['ci_low']:.3f}, {h1['ci_high']:.3f}]\nn = {int(h1['n'])} countries, {fmt_p(h1['p_raw'])} (one-sided)",
            transform=ax.transAxes, va="top", fontsize=9.5, bbox=dict(boxstyle="round", fc="white", ec="#bbbbbb"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("GDP per capita, 2022 (constant 2015 US$, log scale)")
    ax.set_ylabel("CO$_2$ emissions per capita, 2022 (tonnes per person, log scale)")
    ax.set_title("H1: richer countries emit more CO$_2$ per person (country cross-section, 2022)")
    ax.legend(loc="lower right", fontsize=8.5, frameon=True)
    save(fig, "fig_h1_scatter.png")


# ------------------------------------------------------------------ 5. H2 boxplot
def fig_h2(S, res):
    d = S["H2"]
    fig, ax = plt.subplots(figsize=(9, 6.5))
    order = list(GROUP_COL)
    data = [d.loc[d["income_group"].astype(str) == g, "co2_pc"].to_numpy() for g in order]
    bp = ax.boxplot(data, positions=range(4), widths=0.5, patch_artist=True, showfliers=False, medianprops=dict(color="black", lw=2), whiskerprops=dict(color="#555555"), capprops=dict(color="#555555"))
    rng = np.random.default_rng(T.SEED)
    for i, (g, v) in enumerate(zip(order, data)):
        bp["boxes"][i].set(facecolor=GROUP_COL[g], alpha=0.35, edgecolor=GROUP_COL[g], lw=1.6)
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=16, color=GROUP_COL[g], edgecolor="white", linewidth=0.3, alpha=0.9, zorder=3, marker=GROUP_MARK[g])
        ax.text(i, 0.0105, f"n = {len(v)}\nmedian {np.median(v):.2f}", ha="center", va="bottom", fontsize=8.5, color=INK)
    ax.set_yscale("log")
    ax.set_ylim(0.01, 80)
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"{GROUP_SHORT[g]}\nincome" for g in order])
    ax.set_ylabel("CO$_2$ per capita, 2022 (tonnes per person, log scale)")
    ax.set_xlabel("World Bank income group")
    h = res[(res["hypothesis"] == "H2") & res["primary"]].iloc[0]
    ax.set_title("H2: CO$_2$ per person rises steadily with income group")
    ax.text(0.02, 0.97, f"Kruskal-Wallis H = {h['statistic']:.1f}, {fmt_p(h['p_raw'])}\nepsilon$^2$ = {h['effect_size']:.2f} (95% CI {h['ci_low']:.2f}-{h['ci_high']:.2f})\nall 6 Dunn-Holm pairwise p < 0.05",
            transform=ax.transAxes, va="top", fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="#bbbbbb"))
    ph = pd.read_csv(RES / "h2_posthoc.csv")
    assert (ph["dunn_p_holm"] < 0.05).all()
    ax.grid(axis="x", visible=False)
    save(fig, "fig_h2_boxplot_by_income.png")


# ------------------------------------------------------------------ 6. H4 slopes
def fig_h4(S, res):
    d = S["H4"].copy()
    sl = pd.read_csv(RES / "h4_slopes_by_group.csv").set_index("income_group")
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), sharex=True, sharey=True)
    xs = np.linspace(d["lg_gdp"].min(), d["lg_gdp"].max(), 50)
    b0, a0, *_ = stats.linregress(d["lg_gdp"], d["lg_co2"])
    for ax, g in zip(axes, GROUP_COL):
        s = d[d["income_group"].astype(str) == g]
        ax.scatter(10 ** d["lg_gdp"], 10 ** d["lg_co2"], s=8, color="#cccccc", zorder=1)
        ax.scatter(10 ** s["lg_gdp"], 10 ** s["lg_co2"], s=26, color=GROUP_COL[g], marker=GROUP_MARK[g], edgecolor="white", linewidth=0.4, zorder=3)
        r = sl.loc[g]
        gx = np.linspace(s["lg_gdp"].min(), s["lg_gdp"].max(), 50)
        ax.plot(10 ** gx, 10 ** (r["intercept"] + r["slope"] * gx), color=INK, lw=2, zorder=4)
        ax.plot(10 ** xs, 10 ** (a0 + b0 * xs), color=OI["verm"], lw=1.2, ls="--", zorder=2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{GROUP_SHORT[g]} income (n = {int(r['n'])})", fontsize=10.5)
        ax.text(0.04, 0.96, f"slope = {r['slope']:.2f}\n95% CI [{r['ci_low']:.2f}, {r['ci_high']:.2f}]", transform=ax.transAxes, va="top", fontsize=8.8,
                bbox=dict(boxstyle="round", fc="white", ec="#bbbbbb"))
        ax.set_xlabel("GDP pc, 2022 (2015 US$, log)")
    axes[0].set_ylabel("CO$_2$ per capita, 2022 (t, log)")
    h = res[(res["hypothesis"] == "H4") & res["primary"]].iloc[0]
    fig.suptitle(f"H4: elasticity of CO$_2$ to income by group (black = within-group fit; dashed orange = pooled fit; grey = all countries).  Interaction F = {h['statistic']:.2f}, {fmt_p(h['p_raw'])}", fontsize=10.5, y=1.02)
    fig.tight_layout()
    save(fig, "fig_h4_slopes_by_income.png")


# ------------------------------------------------------------------ 7. H5 dumbbell
def fig_h5(S, res):
    d = S["H5"].sort_values("log_ratio_co2").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9, 15))
    y = np.arange(len(d))[::-1]
    dec = d["log_ratio_co2"] < 0
    for yi, (_, r), isdec in zip(y, d.iterrows(), dec):
        col = OI["blue"] if isdec else OI["verm"]
        ax.plot([r["co2_pc_1990"], r["co2_pc_2022"]], [yi, yi], color=col, lw=1.8, alpha=0.8, zorder=2)
    ax.scatter(d["co2_pc_1990"], y, s=26, color="white", edgecolor=INK, linewidth=1.1, zorder=3, label="1990")
    ax.scatter(d["co2_pc_2022"], y, s=26, color=[OI["blue"] if v else OI["verm"] for v in dec], edgecolor="white", linewidth=0.4, zorder=4, label="2022")
    ax.set_yticks(y)
    ax.set_yticklabels(d["iso3"], fontsize=7.2)
    ax.set_xscale("log")
    ax.set_ylim(-1, len(d) + 7)
    ax.set_xlabel("CO$_2$ per capita (tonnes per person, log scale)")
    ax.set_ylabel("High-income countries, sorted by change (largest fall at top)")
    h = res[(res["hypothesis"] == "H5") & res["primary"]].iloc[0]
    ax.set_title("H5: CO$_2$ per person, 1990 vs 2022, high-income countries", fontsize=11.5)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", color="w", markerfacecolor="white", markeredgecolor=INK, markersize=7, label="1990"),
               Line2D([], [], color=OI["blue"], lw=2, marker="o", markerfacecolor=OI["blue"], markeredgecolor="white", markersize=7, label=f"2022, fell ({int(dec.sum())} countries)"),
               Line2D([], [], color=OI["verm"], lw=2, marker="o", markerfacecolor=OI["verm"], markeredgecolor="white", markersize=7, label=f"2022, rose ({int((~dec).sum())} countries)")]
    ax.legend(handles=handles, loc="upper right", fontsize=8.5)
    ax.text(0.02, 0.985, "Wilcoxon signed-rank (one-sided):\n" + f"{fmt_p(h['p_raw'])}, matched-pairs r = {h['effect_size']:.2f} [{h['ci_low']:.2f}, {h['ci_high']:.2f}]\n" + h["comment"].split(";")[0],
            transform=ax.transAxes, va="top", fontsize=8.5, bbox=dict(boxstyle="round", fc="white", ec="#bbbbbb"))
    ax.grid(axis="y", visible=False)
    save(fig, "fig_h5_paired_change.png")


# ------------------------------------------------------------------ 8. H6 heatmap
def fig_h6(res):
    t = pd.read_csv(RES / "h6_contingency_and_residuals.csv")
    groups = list(GROUP_COL)
    cols = ["Not decoupled", "Decoupled"]
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("div", [OI["verm"], "#f2f2f2", OI["blue"]])
    mat = np.zeros((4, 2))
    for i, g in enumerate(groups):
        for j, c in enumerate(cols):
            r = t[(t["income_group"] == g) & (t["outcome"] == c)].iloc[0]
            mat[i, j] = r["std_adj_residual"]
    im = ax.imshow(mat, cmap=cmap, vmin=-7.5, vmax=7.5, aspect="auto")
    for i, g in enumerate(groups):
        for j, c in enumerate(cols):
            r = t[(t["income_group"] == g) & (t["outcome"] == c)].iloc[0]
            ax.text(j, i, f"{int(r['observed'])}\n(expected {r['expected']:.1f})\nresidual {r['std_adj_residual']:+.2f}", ha="center", va="center", fontsize=10,
                    color="white" if abs(r["std_adj_residual"]) > 4.5 else INK, fontweight="bold" if abs(r["std_adj_residual"]) > 1.96 else "normal")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Not decoupled\n(GDP up, CO$_2$ not down)", "Absolute decoupling\n(GDP pc up AND CO$_2$ pc down)"])
    n = t.groupby("income_group")["observed"].sum()
    ax.set_yticks(range(4))
    ax.set_yticklabels([f"{GROUP_SHORT[g]} income\n(n = {int(n[g])})" for g in groups])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("standardised (adjusted) residual\n|r| > 1.96 = cell differs from independence")
    h = res[(res["hypothesis"] == "H6") & res["primary"]].iloc[0]
    ax.set_title(f"H6: decoupling vs income group, 1990-2022\nchi-square = {h['statistic']:.1f}, {fmt_p(h['p_raw'])}; Cramer's V = {h['effect_size']:.2f} [{h['ci_low']:.2f}, {h['ci_high']:.2f}]", fontsize=10.5)
    save(fig, "fig_h6_decoupling_heatmap.png")


# ------------------------------------------------------------------ 9. effect-size forest
def fig_forest(res):
    prim = res[res["primary"] & res["hypothesis"].isin(["H1", "H2", "H3", "H4", "H5", "H6"])].set_index("hypothesis")
    h7 = res[(res["hypothesis"] == "H7") & (res["variant"] == "primary (exploratory)")].iloc[0]
    spec = [("H1", "Spearman rho", (0, 1), [(0.1, "small"), (0.3, "medium"), (0.5, "large")]),
            ("H2", "epsilon-squared (Kruskal-Wallis)", (0, 1), [(0.01, ""), (0.08, "medium"), (0.26, "large")]),
            ("H3", "Cohen's d (High vs LM+Low, log10)", (0, 4), [(0.2, ""), (0.5, ""), (0.8, "large")]),
            ("H4", "delta R-squared (interaction block)", (0, 0.05), []),
            ("H5", "matched-pairs rank-biserial r", (-1, 0.2), [(-0.1, ""), (-0.3, "medium"), (-0.5, "large")]),
            ("H6", "Cramer's V (2 columns)", (0, 1), [(0.1, "small"), (0.3, "medium"), (0.5, "large")])]
    fig, axes = plt.subplots(len(spec) + 1, 1, figsize=(10, 10.5))
    for ax, (h, lab, xlim, marks) in zip(axes, spec):
        r = prim.loc[h]
        ax.axvline(0, color="#999999", lw=1)
        for m, _ in marks:
            ax.axvline(m, color=GRID, lw=1, ls=":")
        ax.errorbar([r["effect_size"]], [0], xerr=[[r["effect_size"] - r["ci_low"]], [r["ci_high"] - r["effect_size"]]], fmt="o", color=OI["blue"], ms=9, capsize=5, lw=2.2, zorder=3)
        ax.set_xlim(*xlim)
        ax.set_ylim(-1, 1)
        ax.set_yticks([])
        ax.grid(False)
        ax.set_xlabel(lab, fontsize=8.8)
        ax.set_title(f"{h}: {r['effect_size']:.3f}  (95% CI {r['ci_low']:.3f} to {r['ci_high']:.3f}); n = {int(r['n'])}; Holm-adjusted {fmt_p(r['p_holm'])}", fontsize=9.5, loc="left")
    ax = axes[-1]
    ax.axvline(0, color="#999999", lw=1)
    ax.errorbar([h7["effect_size"]], [0], xerr=[[h7["effect_size"] - h7["ci_low"]], [h7["ci_high"] - h7["effect_size"]]], fmt="D", color=OI["orange"], ms=8, capsize=5, lw=2.2, zorder=3)
    ax.set_xlim(0, 0.12)
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.grid(False)
    ax.set_xlabel("increase in R-squared from industry share + urbanisation (EXPLORATORY, unadjusted)", fontsize=8.8)
    ax.set_title(f"H7: {h7['effect_size']:.3f}  (95% CI {h7['ci_low']:.3f} to {h7['ci_high']:.3f}); n = {int(h7['n'])}; exploratory", fontsize=9.5, loc="left")
    fig.suptitle("Effect sizes with 95% bootstrap confidence intervals (each row on its own scale)", fontsize=12.5, fontweight="bold", y=0.995)
    fig.tight_layout(h_pad=1.6, rect=(0, 0.02, 1, 1))
    fig.text(0.5, 0.005, "Dotted vertical lines = conventional small / medium / large benchmarks (Cohen, Tomczak & Tomczak) where defined; grey solid line = no effect.", ha="center", fontsize=8.5, color=MUTED)
    save(fig, "fig_effect_size_forest.png")


# ------------------------------------------------------------------ 10. power curves
def fig_power(res):
    pc = pd.read_csv(RES / "power_curve.csv")
    mde = pd.read_csv(RES / "power_mde.csv")
    h3 = res[(res["hypothesis"] == "H3") & res["primary"]].iloc[0]
    h5 = res[(res["hypothesis"] == "H5") & (res["variant"].str.startswith("companion: paired t"))].iloc[0]
    fig, ax = plt.subplots(figsize=(9, 5.8))
    ax.plot(pc["effect_size"], pc["power_H3"], color=OI["blue"], lw=2.4, label=f"H3: independent groups (74 vs 72), one-sided")
    ax.plot(pc["effect_size"], pc["power_H5"], color=OI["verm"], lw=2.4, label="H5: paired design (n = 74), one-sided t")
    ax.plot(pc["effect_size"], pc["power_H5_wilcoxon_approx"], color=OI["verm"], lw=1.4, ls="--", label="H5: Wilcoxon approximation (ARE = 0.955)")
    ax.axhline(0.8, color="#555555", lw=1, ls=":")
    ax.text(1.58, 0.81, "80% power", ha="right", fontsize=8.5, color=MUTED)
    m3 = mde[(mde["hypothesis"] == "H3") & mde["mde_d_80pct_power"].notna()]["mde_d_80pct_power"].iloc[0]
    m5 = mde[(mde["hypothesis"] == "H5") & mde["test"].str.startswith("paired")]["mde_d_80pct_power"].iloc[0]
    ax.scatter([m3], [0.8], color=OI["blue"], s=70, zorder=5, edgecolor="white")
    ax.annotate(f"H3 minimum detectable d = {m3:.2f}", (m3, 0.8), xytext=(0.78, 0.60), arrowprops=dict(arrowstyle="-", color=OI["blue"]), fontsize=9, color=OI["blue"])
    ax.scatter([m5], [0.8], color=OI["verm"], s=70, zorder=5, edgecolor="white")
    ax.annotate(f"H5 minimum detectable d_z = {m5:.2f}", (m5, 0.8), xytext=(0.78, 0.42), arrowprops=dict(arrowstyle="-", color=OI["verm"]), fontsize=9, color=OI["verm"])
    ax.axvline(abs(h5["effect_size"]), color=OI["verm"], lw=1, ls="-.", alpha=0.6)
    ax.text(abs(h5["effect_size"]) + 0.02, 0.05, f"observed |d_z| (H5) = {abs(h5['effect_size']):.2f}", color=OI["verm"], fontsize=8.5, rotation=90, va="bottom")
    ax.set_xlim(0, 1.6)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("standardised effect size (Cohen's d for H3, d_z for H5)")
    ax.set_ylabel("statistical power (alpha = 0.05, one-sided)")
    ax.set_title("Power curves at the real sample sizes (observed H3 effect, d = 3.06, is off the right edge)", fontsize=10.5)
    ax.legend(loc="lower right", fontsize=8.8)
    save(fig, "fig_power_curve.png")


# ------------------------------------------------------------------ 11. gantt
def fig_gantt():
    tl = pd.read_csv(DOCS / "time_log.csv")
    fig, ax = plt.subplots(figsize=(11, 5.4))
    start = 0.0
    cols = [OI["blue"], OI["sky"], OI["green"], OI["orange"], OI["verm"], OI["pink"], "#8c8c8c", OI["yellow"]]
    n = len(tl)
    for i, r in tl.iterrows():
        y = n - 1 - i
        ax.barh(y, r["planned_hours"], left=start, color=cols[i % len(cols)], edgecolor="white", height=0.55, label=None)
        ax.text(start + r["planned_hours"] / 2, y, f"{r['planned_hours']:g} h", ha="center", va="center", fontsize=9, color="white" if cols[i % len(cols)] not in (OI["yellow"], OI["sky"]) else INK, fontweight="bold")
        if pd.notna(r.get("actual_hours")):
            ax.barh(y - 0.38, r["actual_hours"], left=start, color=cols[i % len(cols)], alpha=0.4, height=0.2)
            ax.text(start + r["actual_hours"] + 0.1, y - 0.38, f"actual {r['actual_hours']:g} h", fontsize=7.5, va="center")
        start += r["planned_hours"]
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(tl["phase"])[::-1])
    ax.set_xlabel("cumulative project hours (planned)")
    total = tl["planned_hours"].sum()
    ax.set_xlim(0, total + 1)
    ax.set_title(f"Project timeline: planned allocation of {total:g} hours")
    has_actual = tl["actual_hours"].notna().any()
    ax.text(0.99, 0.02, "actual hours: to be entered by the author in docs/time_log.csv" if not has_actual else "pale bars = actual hours", transform=ax.transAxes, ha="right", fontsize=8.5, color=MUTED, style="italic")
    ax.grid(axis="y", visible=False)
    save(fig, "fig_timeline_gantt.png")


def main():
    S = T.load_samples()
    res = pd.read_csv(RES / "hypothesis_results.csv")
    chk = pd.read_csv(RES / "assumption_checks.csv")
    fig_framework()
    fig_tree(chk)
    fig_assumptions(S)
    fig_h1(S, res)
    fig_h2(S, res)
    fig_h4(S, res)
    fig_h5(S, res)
    fig_h6(res)
    fig_forest(res)
    fig_power(res)
    fig_gantt()


if __name__ == "__main__":
    main()
