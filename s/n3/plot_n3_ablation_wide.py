#!/usr/bin/env python3
"""Landscape (16:9) version of the ablation figure (Fig. 4) for presentation slides.
Same data and colours as plot_n3_ablation.py, transposed: workflows run along the x-axis,
grouped by calibration step; rows are RMSE and total bias, columns are the two regions.
Reads ../data/processed/n3/n3_holdout_metrics.csv (n3_ablation.R).
Run from s/:  python3 n3/plot_n3_ablation_wide.py
Writes ../figures/presentation/2226_n3_ablation_wide.pdf|png (and copies the PNG to
~/saltnoaa/presentations/figures if that folder exists)."""
import shutil
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = pd.read_csv("../data/processed/n3/n3_holdout_metrics.csv")
OUT = Path("../figures/presentation")
PRES = Path.home() / "saltnoaa/presentations/figures"
TIERS = ["F1 threshold", "Stage 1 only", "Stage 2 on raw counts", "Two-stage, one detector",
         "Two-stage, two detectors", "Two-stage, three detectors"]
TIER_LAB = {"F1 threshold": "$F_1$\nthreshold", "Stage 1 only": "Stage 1\nonly",
            "Stage 2 on raw counts": "Stage 2,\nraw counts", "Two-stage, one detector": "Two-stage\n1 detector",
            "Two-stage, two detectors": "Two-stage\n2 detectors", "Two-stage, three detectors": "Two-\nstage\n3 det."}
SHORT = {"YOLOv12": "YOLO", "Cascade R-CNN": "Cascade", "RT-DETR-L": "RT-DETR",
         "YOLOv12 + Cascade R-CNN": "YOLO + Cascade", "YOLOv12 + RT-DETR-L": "YOLO + RT-DETR",
         "Cascade R-CNN + RT-DETR-L": "Cascade + RT-DETR", "YOLOv12 + Cascade R-CNN + RT-DETR-L": "All three"}
ORDER = ["YOLOv12", "Cascade R-CNN", "RT-DETR-L", "YOLOv12 + Cascade R-CNN", "YOLOv12 + RT-DETR-L",
         "Cascade R-CNN + RT-DETR-L", "YOLOv12 + Cascade R-CNN + RT-DETR-L"]
COL = {"F1 threshold": "#8c8c8c", "Stage 1 only": "#52514e", "Stage 2 on raw counts": "#d6457a",
       "Two-stage, one detector": "#1baf7a", "Two-stage, two detectors": "#2a78d6", "Two-stage, three detectors": "#6250d6"}
REF = ("Two-stage, one detector", "YOLOv12")          # recommended workflow, outlined and banded

M["tier"] = M.method.str.split(": ").str[0]
M["set"] = M.method.str.split(": ").str[1]

# x positions: one slot per detector set, a gap between calibration steps
cols, x, xt, xl, spans = [], 0.0, [], [], []
for t in TIERS:
    sets = [s for s in ORDER if ((M.tier == t) & (M.set == s)).any()]
    x0 = x
    for s in sets:
        cols.append((t, s, x)); xt.append(x); xl.append(SHORT[s]); x += 1
    spans.append((t, x0 - 0.45, x - 0.55))
    x += 1.0

plt.rcParams.update({"font.size": 12, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(2, 2, figsize=(13.33, 7.5), sharex=True)
rows = [("rmse", "rmse_lo", "rmse_hi", "RMSE (scallops per image)"),
        ("pct_bias_total", "bias_lo", "bias_hi", "Total bias (%)")]
for j, (reg, title) in enumerate([("GB", "Georges Bank"), ("MAB", "Mid-Atlantic Bight")]):
    for i, (met, lo, hi, ylab) in enumerate(rows):
        ax = axes[i, j]
        ref = M[(M.region == reg) & (M.tier == REF[0]) & (M.set == REF[1])].iloc[0]
        ax.axhspan(ref[lo], ref[hi], color=COL[REF[0]], alpha=0.10, lw=0)
        if met == "pct_bias_total": ax.axhline(0, color="#52514e", lw=0.9, ls="--")
        for t, s, xx in cols:
            r = M[(M.region == reg) & (M.tier == t) & (M.set == s)]
            if not len(r): continue
            r = r.iloc[0]; c = COL[t]; is_ref = (t, s) == REF
            ax.vlines(xx, r[lo], r[hi], color=c, lw=2.2)
            ax.plot(xx, r[met], "o", color=c, ms=9 if is_ref else 7, mec="black" if is_ref else c, mew=1.4)
        ax.grid(axis="y", color="#e6e5e0", lw=0.7); ax.set_axisbelow(True)
        if j == 0: ax.set_ylabel(ylab)
        if i == 0:
            ax.set_title(title, fontsize=14, pad=44)
            ymax = ax.get_ylim()[1]
            for t, a, b in spans:          # calibration-step labels above the top row
                ax.annotate(TIER_LAB[t], xy=((a + b) / 2, 1.0), xycoords=("data", "axes fraction"),
                            xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                            fontsize=10.5, color=COL[t], fontweight="bold", linespacing=0.95)
        for t, a, b in spans[1:]:           # faint dividers between calibration steps
            ax.axvline(a - 0.5, color="#d9d8d2", lw=0.8)
for ax in axes[1]:
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=10, rotation=55, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.6, x - 0.6)
fig.text(0.5, 0.005, "Points: Holdout test estimate; bars: 95% transect-block bootstrap interval; "
         "outlined point and green band: recommended two-stage YOLOv12", ha="center", fontsize=10, color="#52514e")
fig.tight_layout(rect=(0, 0.03, 1, 1), h_pad=1.2, w_pad=2.0)
OUT.mkdir(parents=True, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"2226_n3_ablation_wide.{ext}", dpi=300, bbox_inches="tight")
if PRES.exists(): shutil.copy(OUT / "2226_n3_ablation_wide.png", PRES / "2226_n3_ablation_wide.png")
print("Wrote", OUT / "2226_n3_ablation_wide.pdf", "and .png")
