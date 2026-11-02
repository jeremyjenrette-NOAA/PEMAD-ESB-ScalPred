#!/usr/bin/env python3
"""Ablation figure (main text): Holdout RMSE and total bias for every calibration step and
detector combination, with short labels. Reads ../data/processed/n3/n3_holdout_metrics.csv
(n3_ablation.R). Run from s/:  python3 n3/plot_n3_ablation.py
Writes ../figures/ms_figures/2226_n3_ablation.pdf|png (and copies to Overleaf if present)."""
import shutil
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = pd.read_csv("../data/processed/n3/n3_holdout_metrics.csv")
FIG = Path("../figures/ms_figures")
OVL = Path.home() / "Library/CloudStorage/Dropbox/Apps/Overleaf/NEFSC-scallop-detection/ms_figures"
TIERS = ["F1 threshold", "Stage 1 only", "Stage 2 on raw counts", "Two-stage, one detector",
         "Two-stage, two detectors", "Two-stage, three detectors"]
TIER_LAB = {"F1 threshold": "$F_1$ threshold", "Stage 1 only": "Stage 1 only",
            "Stage 2 on raw counts": "Stage 2 on raw counts", "Two-stage, one detector": "Two-stage, one detector",
            "Two-stage, two detectors": "Two-stage, two detectors", "Two-stage, three detectors": "Two-stage, three detectors"}
SHORT = {"YOLOv12": "YOLOv12", "Cascade R-CNN": "Cascade", "RT-DETR-L": "RT-DETR",
         "YOLOv12 + Cascade R-CNN": "YOLO + Cascade", "YOLOv12 + RT-DETR-L": "YOLO + RT-DETR",
         "Cascade R-CNN + RT-DETR-L": "Cascade + RT-DETR", "YOLOv12 + Cascade R-CNN + RT-DETR-L": "All three"}
ORDER = ["YOLOv12", "Cascade R-CNN", "RT-DETR-L", "YOLOv12 + Cascade R-CNN", "YOLOv12 + RT-DETR-L",
         "Cascade R-CNN + RT-DETR-L", "YOLOv12 + Cascade R-CNN + RT-DETR-L"]
COL = {"F1 threshold": "#8c8c8c", "Stage 1 only": "#52514e", "Stage 2 on raw counts": "#d6457a",
       "Two-stage, one detector": "#1baf7a", "Two-stage, two detectors": "#2a78d6", "Two-stage, three detectors": "#6250d6"}
REF = ("Two-stage, one detector", "YOLOv12")

M["tier"] = M.method.str.split(": ").str[0]
M["set"] = M.method.str.split(": ").str[1]
rows, y, yt, ylab, heads = [], 0, [], [], []
for t in TIERS:
    sets = [s for s in ORDER if ((M.tier == t) & (M.set == s)).any()]
    heads.append((y - 0.15, TIER_LAB[t]))
    y += 0.55
    for s in sets:
        rows.append((t, s, y)); yt.append(y); ylab.append(SHORT[s]); y += 1
    y += 0.35
plt.rcParams.update({"font.size": 8.5, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 4, figsize=(7.2, 6.3), sharey=True)
panels = [("GB", "rmse", "Georges Bank\nRMSE (scallops per image)"), ("GB", "pct_bias_total", "Georges Bank\ntotal bias (%)"),
          ("MAB", "rmse", "Mid-Atlantic Bight\nRMSE (scallops per image)"), ("MAB", "pct_bias_total", "Mid-Atlantic Bight\ntotal bias (%)")]
for ax, (reg, met, title) in zip(axes, panels):
    lo, hi = ("rmse_lo", "rmse_hi") if met == "rmse" else ("bias_lo", "bias_hi")
    if met == "pct_bias_total": ax.axvline(0, color="#52514e", lw=0.8, ls="--")
    for t, s, yy in rows:
        r = M[(M.region == reg) & (M.tier == t) & (M.set == s)]
        if not len(r): continue
        r = r.iloc[0]; c = COL[t]
        ax.hlines(yy, r[lo], r[hi], color=c, lw=1.4)
        ax.plot(r[met], yy, "o", color=c, ms=5 if (t, s) == REF else 4, mec="black" if (t, s) == REF else c, mew=0.9)
    ax.set_title(title, fontsize=8.5)
    ax.grid(axis="x", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    ref = M[(M.region == reg) & (M.tier == REF[0]) & (M.set == REF[1])].iloc[0]
    ax.axvspan(ref[lo], ref[hi], color="#1baf7a", alpha=0.08, lw=0)
axes[0].set_yticks(yt); axes[0].set_yticklabels(ylab)
for yy, h in heads:
    axes[0].annotate(h, xy=(0, yy), xycoords=("axes fraction", "data"), xytext=(-4, 0), textcoords="offset points",
                     ha="right", va="center", fontsize=8.5, fontweight="bold", annotation_clip=False)
axes[0].tick_params(axis="y", length=0)
axes[0].set_ylim(y - 0.2, -0.7)
fig.tight_layout(w_pad=0.6)
FIG.mkdir(parents=True, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(FIG / f"2226_n3_ablation.{ext}", dpi=300, bbox_inches="tight")
    if OVL.exists(): shutil.copy(FIG / f"2226_n3_ablation.{ext}", OVL / f"2226_n3_ablation.{ext}")
print("Wrote 2226_n3_ablation.pdf/png")
