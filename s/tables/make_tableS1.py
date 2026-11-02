#!/usr/bin/env python3
"""Table S1 (image-level, covariate and spatial error of the three plotted workflows, by region)
from error.R outputs. Run from s/:  python3 tables/make_tableS1.py  -> tables/tableS1_error_diagnostics.tex"""
import csv
from pathlib import Path
P = Path("../data/processed")
rd = lambda f: list(csv.DictReader(open(P / f)))
M = ["Two-stage YOLOv12", "YOLOv12 F1", "Cascade F1"]
R = ["GB", "MAB"]
img = {(r["region"], r["Method"]): r for r in rd("holdout_image_error_summary_2226.csv")}
env = {(r["region"], r["Method"], r["Covariate"]): r for r in rd("holdout_covariate_stability_2226.csv")}
sp  = {(r["region"], r["Method"]): r for r in rd("holdout_spatial_stability_2226.csv")}
COV = [("Altitude (m)", "Vehicle altitude (m)"), ("Bottom depth (m)", "Bottom depth (m)"),
       ("Field of view (m²)", "Image field of view (m$^2$)"), ("Temperature (°C)", "Water temperature ($^\\circ$C)"),
       ("Backscatter (NTU)", "Backscatter (NTU)")]
f = lambda x, d=2: f"{float(x):.{d}f}".replace("-", "$-$")
def row(label, get, d=2):
    return " & ".join([label] + [f(get(r, m), d) for r in R for m in M]) + r" \\"
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Holdout test error of two-stage \gls{yolo}v12 and the two $F_1$ threshold workflows across images, environmental covariates and space, by region.}",
     r"\label{tab:error_diagnostics}", r"\resizebox{\textwidth}{!}{%", r"\begin{tabular}{lcccccc}", r"\toprule",
     r" & \multicolumn{3}{c}{\textbf{Georges Bank}} & \multicolumn{3}{c}{\textbf{Mid-Atlantic Bight}} \\",
     r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
     r"\textbf{Metric} & \textbf{Two-stage} & \textbf{YOLOv12 $F_1$} & \textbf{Cascade $F_1$} & \textbf{Two-stage} & \textbf{YOLOv12 $F_1$} & \textbf{Cascade $F_1$} \\ \midrule",
     r"\textbf{Image level} & & & & & & \\",
     row(r"\hspace{1em} SD of error", lambda r, m: img[(r, m)]["SD"]),
     row(r"\hspace{1em} MAE", lambda r, m: img[(r, m)]["MAE"]),
     row(r"\hspace{1em} Largest image error", lambda r, m: img[(r, m)]["Max_Error"], 1),
     r"\midrule", r"\textbf{Largest covariate-bin bias} & & & & & & \\"]
L += [row(r"\hspace{1em} " + lab, lambda r, m, c=c: env[(r, m, c)]["Peak_Bin_Bias"]) for c, lab in COV]
L += [r"\midrule", r"\textbf{Spatial grid cells}$^{\text{a}}$ & & & & & & \\",
      row(r"\hspace{1em} MAE of cell means", lambda r, m: sp[(r, m)]["Grid_MAE"]),
      row(r"\hspace{1em} Largest cell bias", lambda r, m: sp[(r, m)]["Max_Cell_Bias"]),
      r"\bottomrule", r"\end{tabular}%", "}", r"\vspace{2mm}", r"\begin{flushleft}", r"\footnotesize",
      r"Errors are predicted minus manual scallops per image. Two-stage: two-stage calibration of \gls{yolo}v12. Largest covariate-bin bias: the mean error of the covariate decile (within region) furthest from zero (Fig.~5). "
      r"$^{\text{a}}$ $0.1^\circ \times 0.1^\circ$ cells with at least three Holdout test images.",
      r"\end{flushleft}", r"\end{table}"]
Path("tables").mkdir(exist_ok=True)
Path("tables/tableS1_error_diagnostics.tex").write_text("\n".join(L) + "\n")
print("wrote tables/tableS1_error_diagnostics.tex")
