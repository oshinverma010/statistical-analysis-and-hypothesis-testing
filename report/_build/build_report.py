"""Expand @@TOKENS@@ in report/_build/report_template.md into report/report_content.md.

Code snippets are copied verbatim from src/ and tables are rendered from results/*.csv, so the report cannot drift from the code or the data.
"""
import re, textwrap
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"


def fp(p):
    if pd.isna(p): return "–"
    return f"{p:.2e}" if p < 1e-4 else f"{p:.4f}"


def fn(x, d=3):
    return "–" if pd.isna(x) else f"{x:.{d}g}" if abs(x) >= 100 else f"{x:.{d}f}"


def eff(r):
    if pd.isna(r["effect_size"]): return "–"
    s = f"{r['effect_size']:.3f}"
    if pd.notna(r["ci_low"]): s += f" [{r['ci_low']:.3f}, {r['ci_high']:.3f}]"
    return s


def md(rows, header):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def snip(spec):
    path, start, end, lang = spec.split("|")
    lines = (ROOT / path).read_text(encoding="utf-8").splitlines()
    i = next(k for k, l in enumerate(lines) if start in l)
    j = next(k for k in range(i, len(lines)) if end in lines[k])
    block = textwrap.dedent("\n".join(lines[i:j + 1]))
    return f"```{lang}\n# from {path}\n{block}\n```"


def primary_table():
    r = pd.read_csv(RES / "hypothesis_results.csv")
    rows = []
    for h in ["H1", "H2", "H3", "H4", "H5", "H6"]:
        x = r[(r.hypothesis == h) & r.primary].iloc[0]
        dec = x["decision"].split(" (")[0]
        rows.append([h, x["test"], int(x["n"]), fn(x["statistic"]), fp(x["p_raw"]), fp(x["p_holm"]), f"{x['effect_size_name']} = {eff(x)}", dec])
    x = r[(r.hypothesis == "H7") & (r.variant == "primary (exploratory)")].iloc[0]
    rows.append(["H7*", x["test"], int(x["n"]), fn(x["statistic"]), fp(x["p_raw"]), "not adjusted", f"{x['effect_size_name']} = {eff(x)}", "exploratory"])
    return md(rows, ["ID", "Primary test", "n", "Statistic", "Raw p", "Holm-adj. p", "Effect size [95% CI]", "Decision (α = 0.05)"]) + "\n\n*\\*H7 is exploratory and outside the Holm family.*"


def robust_table():
    r = pd.read_csv(RES / "hypothesis_results.csv")
    r = r[r.variant.str.startswith("robustness") & r.hypothesis.isin(["H1", "H2", "H3", "H4", "H5", "H6"])]
    rows = [[x.hypothesis, x.variant.replace("robustness: ", ""), int(x.n), fn(x.statistic), fp(x.p_raw), eff(x), "agrees" if x.direction_ok else "DISAGREES"] for _, x in r.iterrows()]
    return md(rows, ["ID", "Robustness variant", "n", "Statistic", "Raw p", "Effect size [95% CI]", "Direction vs primary"])


def full_table():
    r = pd.read_csv(RES / "hypothesis_results.csv")
    rows = []
    for _, x in r.iterrows():
        rows.append([x.hypothesis, x.variant, x.test.replace("|", "/"), int(x.n), fn(x.statistic), fp(x.p_raw), fp(x.p_holm), fp(x.p_bh),
                     f"{x.effect_size_name}: {eff(x)}", (x.decision if isinstance(x.decision, str) else "")])
    return md(rows, ["ID", "Variant", "Test", "n", "Statistic", "Raw p", "Holm p", "BH p", "Effect size [95% CI]", "Decision / note"])


def clean_table():
    c = pd.read_csv(RES / "cleaning_log.csv").fillna("")
    rows = [[x.step, x.rule, f"{x.rows_before:,}", f"{x.rows_after:,}", x.note] for x in c.itertuples()]
    return "**Cleaning log (step totals from `results/cleaning_log.csv`)**\n\n" + md(rows, ["Step", "Rule", "Rows before", "Rows after", "Note"])


def assumption_table():
    a = pd.read_csv(RES / "assumption_checks.csv")
    rows = [[x.hypothesis, x.check, fn(x.statistic, 3), fp(x.p_value) if pd.notna(x.p_value) else "–", x.verdict] for x in a.itertuples()]
    return "**Assumption checks (from `results/assumption_checks.csv`)**\n\n" + md(rows, ["ID", "Check", "Statistic", "p", "Verdict"])


text = (ROOT / "report" / "_build" / "report_template.md").read_text(encoding="utf-8")
text = re.sub(r"@@SNIP (.+?)@@", lambda m: snip(m.group(1)), text)
for tok, fn_ in {"PRIMARY_TABLE": primary_table, "ROBUST_TABLE": robust_table, "FULL_TABLE": full_table, "CLEAN_LOG_TABLE": clean_table,
                 "ASSUMPTION_TABLE": assumption_table}.items():
    text = text.replace(f"@@{tok}@@", fn_())
assert "@@" not in text, re.findall(r"@@.{0,40}", text)
(ROOT / "report" / "report_content.md").write_text(text, encoding="utf-8")
# word counts: prose only (outside code fences and tables)
inside, prose = False, 0
for line in text.splitlines():
    if line.startswith("```"): inside = not inside; continue
    if inside or line.startswith("|") or line.startswith("!["): continue
    prose += len(re.findall(r"\w+", line))
print("prose words:", prose, "| total words incl. tables/code:", len(re.findall(r"\w+", text)))
