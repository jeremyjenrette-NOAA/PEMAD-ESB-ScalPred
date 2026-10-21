#!/usr/bin/env python3
"""Manuscript figures for the leave-one-year-out tests (run from s/).

  Fig. S5  tierb2026_effort.pdf  - 2026 calibration accuracy vs annotation effort
  Fig. S6  tiera_loyo_bias.pdf   - total bias by held-out year (GAMs refit without it)

Reads the CSVs written by tierb/tierb_effort.R and tiera/tiera_loyo.R; writes PDF + PNG
to ../figures/revision/ and copies them to the Overleaf ms_figures folder if it exists.
"""
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("../figures/revision"); OUT.mkdir(parents=True, exist_ok=True)
OVERLEAF = Path.home() / "Library/CloudStorage/Dropbox/Apps/Overleaf/NEFSC-scallop-detection/ms_figures"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"     # categorical slots 1-3
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "pdf.fonttype": 42})

def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    if OVERLEAF.exists():
        for ext in ("pdf", "png"):
            shutil.copy(OUT / f"{name}.{ext}", OVERLEAF / f"{name}.{ext}")
    plt.close(fig)

# ---------------------------------------------------------------- Fig. S5: effort
s = pd.read_csv("../data/processed/tierb2026/effort_summary.csv")
series = [("Two-stage", "full", "Two-stage, transferred forms", BLUE, "o"),
          ("Two-stage", "simple", "Two-stage, simplest forms", ORANGE, "s"),
          ("YOLO F1-threshold count", "full", "YOLOv12 F$_1$-threshold count", AQUA, "^")]
fr = sorted(s.fraction.unique()); xpos = {f: i for i, f in enumerate(fr)}
xlab = [f"{f*100:.0f}%\n1 in {50/(0.6*f):,.0f}" for f in fr]
fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0), sharex=True)
RMSE_CAP = {"GB": 0.65, "MAB": 0.22}
for j, reg in enumerate(["GB", "MAB"]):
    for i, (metric, ylab) in enumerate([("pct_bias_total", "Holdout total bias (%)"),
                                         ("rmse", "Holdout RMSE (scallops per image)")]):
        ax = axes[i, j]; style(ax)
        if metric == "pct_bias_total": ax.axhline(0, color=INK2, lw=0.8, ls="--")
        for k, (meth, var, lab, col, mk) in enumerate(series):
            d = s[(s.region == reg) & (s.method == meth) & (s.variant == var)].sort_values("fraction")
            x = np.array([xpos[f] for f in d.fraction]) + (k - 1) * 0.12
            med, lo, hi = d[f"{metric}_med"].values, d[f"{metric}_min"].values, d[f"{metric}_max"].values
            if metric == "rmse":
                cap = RMSE_CAP[reg]
                for xx, h in zip(x, hi):
                    if h > cap:
                        ax.annotate(f"max {h:.1f}", (xx, cap), xytext=(0, -2), textcoords="offset points",
                                    ha="center", va="top", fontsize=6, color=INK2)
                hi = np.minimum(hi, cap)
            ax.vlines(x, lo, hi, color=col, lw=1.2)
            ax.plot(x, med, color=col, lw=1.5, marker=mk, ms=4.5, mec="white", mew=0.8, label=lab)
        if metric == "rmse": ax.set_ylim(top=RMSE_CAP[reg])
        if j == 0: ax.set_ylabel(ylab)
        if i == 0:
            nsc = s[(s.region == reg) & (s.method == "Two-stage") & (s.variant == "full")].sort_values("fraction")
            ax.set_title(f"{'Georges Bank' if reg == 'GB' else 'Mid-Atlantic Bight'}", fontsize=9, color=INK)
        ax.set_xticks(range(len(fr))); ax.set_xticklabels(xlab)
fig.supxlabel("Share of 2026 Calibration-fitting used to fit the GAMs (and the equivalent annotation rate, frames per survey frames)",
              fontsize=8, color=INK, y=0.035)
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.06, 1, 1))
save(fig, "tierb2026_effort")

# ---------------------------------------------------------------- Fig. S6: Tier A
m = pd.read_csv("../data/processed/tiera/tiera_holdout_metrics.csv")
b = pd.read_csv("../data/processed/tiera/tiera_holdout_bootstrap.csv").dropna(subset=["region"])
meths = [("Two-stage, selected", "Two-stage calibration", BLUE, "o"),
         ("YOLOv12 F1", "YOLOv12 F$_1$-threshold count", AQUA, "^"),
         ("Cascade F1", "Cascade R-CNN F$_1$-threshold count", ORANGE, "s")]
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.3), gridspec_kw={"width_ratios": [4, 3]})
for ax, reg in zip(axes, ["GB", "MAB"]):
    style(ax); ax.axhline(0, color=INK2, lw=0.8, ls="--")
    yrs = sorted(m[(m.region == reg)].held_year.unique())
    for k, (meth, lab, col, mk) in enumerate(meths):
        for t, yr in enumerate(yrs):
            x = t + (k - 1) * 0.22
            lo = m[(m.fit_type == "Leave-year-out") & (m.region == reg) & (m.held_year == yr) & (m.method == meth)]
            po = m[(m.fit_type == "Pooled (all years)") & (m.region == reg) & (m.held_year == yr) & (m.method == meth)]
            ci = b[(b.region == reg) & (b.held_year == yr) & (b.method == meth)]
            if len(ci): ax.vlines(x, ci.bias_lo.iloc[0], ci.bias_hi.iloc[0], color=col, lw=1.2)
            ax.plot(x, lo.pct_bias_total.iloc[0], marker=mk, color=col, ms=5.5, mec=col,
                    ls="none", label=lab if t == 0 else None)
            ax.plot(x + 0.09, po.pct_bias_total.iloc[0], marker=mk, mfc="white", mec=col, mew=1.1, ms=4.5,
                    ls="none", label=None)
    ntrue = m[(m.region == reg) & (m.fit_type == "Leave-year-out") & (m.method == "YOLOv12 F1")].set_index("held_year").total_true
    ax.set_xticks(range(len(yrs)))
    ax.set_xticklabels([f"{y}\n{int(ntrue[y]):,} scallops" for y in yrs])
    ax.set_title("Georges Bank" if reg == "GB" else "Mid-Atlantic Bight", fontsize=9, color=INK)
    ax.set_xlabel("Held-out survey year")
axes[0].set_ylabel("Holdout test total bias (%)")
ylim = (min(a.get_ylim()[0] for a in axes), max(a.get_ylim()[1] for a in axes))
for a in axes: a.set_ylim(ylim)
h, l = axes[0].get_legend_handles_labels()
h.append(plt.Line2D([], [], marker="o", color=INK2, ls="none", ms=5)); l.append("GAMs refit without that year (95% CI)")
h.append(plt.Line2D([], [], marker="o", mfc="white", mec=INK2, ls="none", ms=4.5)); l.append("GAMs fit on all years")
fig.legend(h, l, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12))
fig.tight_layout(rect=(0, 0.06, 1, 1))
save(fig, "tiera_loyo_bias")
print("Wrote", *sorted(p.name for p in OUT.iterdir()))
