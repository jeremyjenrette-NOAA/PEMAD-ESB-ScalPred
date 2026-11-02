#!/usr/bin/env python3
"""Build Table 3 (ablation, main text) and Tables S4-S5 (image-level selection under both
rules; concurvity of the submitted forms) from s/n3/n3_ablation.R outputs.
Run from s/:  python3 n3/make_n3_tables.py
Reads  ../data/processed/n3/*.csv and ../data/processed/n3_onese/*.csv
Writes tables/table3_ablation.tex, tables/tableS4_imagelevel_selection.tex, tables/tableS5_concurvity.tex
"""
import csv, sys
from pathlib import Path
base = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../data/processed")
out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("tables")
out.mkdir(exist_ok=True)
rd = lambda p: list(csv.DictReader(open(p)))
H  = {(r["region"], r["method"]): r for r in rd(base / "n3/n3_holdout_metrics.csv")}
H1 = {(r["region"], r["method"]): r for r in rd(base / "n3_onese/n3_holdout_metrics.csv")}
CV  = rd(base / "n3/n3_cv_all.csv"); CV1 = rd(base / "n3_onese/n3_cv_all.csv")
CO  = {(r["region"], r["method"]): r for r in rd(base / "n3/n3_concurvity.csv")}
CO1 = {(r["region"], r["method"]): r for r in rd(base / "n3_onese/n3_concurvity.csv")}
M = "−"
def num(x, d): s = f"{float(x):.{d}f}"; return s.replace("-", "$-$")
def ci(r, k, lo, hi, d): return f"{num(r[k], d)} ({num(r[lo], d)}, {num(r[hi], d)})"
REG = [("GB", "Georges Bank"), ("MAB", "Mid-Atlantic Bight")]

# ---------------- Table 3 ----------------
ROWS = [("$F_1$ threshold", "YOLOv12", "F1 threshold: YOLOv12", 1),
        ("Stage 1 only", "YOLOv12", "Stage 1 only: YOLOv12", 1),
        ("Stage 2 on raw counts", "YOLOv12", "Stage 2 on raw counts: YOLOv12", 1),
        ("Two-stage", "YOLOv12", "Two-stage, one detector: YOLOv12", 1),
        ("Two-stage", "YOLOv12 + Cascade R-CNN", "Two-stage, two detectors: YOLOv12 + Cascade R-CNN", 2),
        ("Two-stage", "All three detectors", "Two-stage, three detectors: YOLOv12 + Cascade R-CNN + RT-DETR-L", 3)]
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Holdout test accuracy of each calibration step and detector combination.}",
     r"\label{tab:ablation}", r"\resizebox{\textwidth}{!}{%", r"\begin{tabular}{llccccc}", r"\hline",
     r" & & & \multicolumn{2}{c}{\textbf{Georges Bank}} & \multicolumn{2}{c}{\textbf{Mid-Atlantic Bight}} \\",
     r"\textbf{Workflow} & \textbf{Detectors} & \textbf{Passes} & \textbf{RMSE} & \textbf{Total bias (\%)} & \textbf{RMSE} & \textbf{Total bias (\%)} \\ \hline"]
for lab, det, meth, p in ROWS:
    cells = []
    for r, _ in REG:
        h = H[(r, meth)]
        cells += [ci(h, "rmse", "rmse_lo", "rmse_hi", 2), ci(h, "pct_bias_total", "bias_lo", "bias_hi", 1)]
    row = [lab, det, str(p)] + cells
    if meth == "Two-stage, one detector: YOLOv12": row = [rf"\textbf{{{c}}}" for c in row]
    L.append(" & ".join(row) + r" \\")
    if lab == "Stage 2 on raw counts": L.append(r"\hline")
L += [r"\hline", r"\end{tabular}%", "}", r"\vspace{2mm}", r"\begin{flushleft}", r"\footnotesize",
      r"RMSE: scallops per image. Total bias: (estimated $-$ manual total)/manual total. Brackets: 95\% transect-block bootstrap intervals ($B = 1000$). "
      r"Passes: detector inference passes per image. Stage 1 only sums the detection-level \gls{gam} probabilities; stage 2 on raw counts fits the image-level \gls{gam} to the uncalibrated detection counts. "
      r"Image-level forms were selected by lowest block cross-validated RMSE on Calibration-fitting (Supplementary Section 7). Bold: recommended workflow. All detectors and combinations are shown in Fig.~S10.",
      r"\end{flushleft}", r"\end{table}"]
(out / "table3_ablation.tex").write_text("\n".join(L) + "\n")

# ---------------- Table S4: selection under both rules ----------------
FORM = {"s_raw": r"$s(r)$", "te_raw_low": r"$te(r, f)$", "s_level": r"$s(\ell)$", "te_lowconf": r"$te(\ell, f)$",
        "te_box": r"$te(\ell, b)$", "crowd": r"$s(\ell) + s(o) + s(n)$", "structure": r"$s(\ell) + s(b) + s(f) + s(o) + s(n)$",
        "level": r"$s(L)$", "level_contrast": r"$s(L) + s(D)$", "te_level_contr": r"$te(L, D)$",
        "te_lc_crowd": r"$te(L, D) + s(o) + s(n)$", "lc_lowconf": r"$s(L) + s(D) + s(f)$",
        "te_L_H1": r"$te(L, H_1) + s(H_2)$", "lc_crowd": r"$s(L) + s(H_1) + s(H_2) + s(o) + s(n)$"}
def form3(t, f):
    if t == "f":
        return {"level": r"$s(L)$", "level_contrast": r"$s(L) + s(H_1) + s(H_2)$",
                "lc_lowconf": r"$s(L) + s(H_1) + s(H_2) + s(f)$"}.get(f, FORM[f])
    return FORM[f]
TL = {"c": "Stage 2 on raw counts", "d": "Two-stage, one detector", "e": "Two-stage, two detectors", "f": "Two-stage, three detectors"}
SH = {"YOLOv12": "Y", "Cascade R-CNN": "C", "RT-DETR-L": "R"}
def short(s): return " + ".join(SH[x] for x in s.split(" + "))
def pick(rows, key):
    return {(r["tier"], r["set"], r["region"]): r for r in rows if r[key] == "TRUE"}
P0, P1 = pick(CV, "best"), pick(CV1, "pick")
ncand = {}
for r in CV: ncand.setdefault((r["tier"], r["set"], r["region"]), [0, 0]); ncand[(r["tier"], r["set"], r["region"])][0] += 1; ncand[(r["tier"], r["set"], r["region"])][1] += r["within_1se"] == "TRUE"
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Image-level model selected for each detector set by lowest cross-validated RMSE (primary) and by the one-standard-error rule (sensitivity), with concurvity and Holdout test accuracy.}",
     r"\label{tab:imagelevel_selection}", r"\resizebox{\textwidth}{!}{%", r"\begin{tabular}{lllcccccccccc}", r"\hline",
     r" & & & & \multicolumn{5}{c}{\textbf{Lowest CV RMSE (primary)}} & \multicolumn{4}{c}{\textbf{One-SE rule}} \\",
     r"\textbf{Workflow} & \textbf{Set} & \textbf{Region} & \textbf{Within 1 SE} & \textbf{Form} & \textbf{CV RMSE (SE)} & \textbf{Conc.} & \textbf{RMSE} & \textbf{Bias (\%)} & \textbf{Form} & \textbf{Conc.} & \textbf{RMSE} & \textbf{Bias (\%)} \\ \hline"]
prev = None
for (t, s, reg), a in P0.items():
    b = P1[(t, s, reg)]; meth = f"{TL[t]}: {s}"
    h0, h1 = H[(reg, meth)], H1[(reg, meth)]
    c0, c1 = float(CO[(reg, meth)]["worst_concurvity"]), float(CO1[(reg, meth)]["worst_concurvity"])
    nc = ncand[(t, s, reg)]
    if prev and prev != t: L.append(r"\hline")
    prev = t
    L.append(" & ".join([TL[t] if reg == "GB" else "", short(s) if reg == "GB" else "", reg, f"{nc[1]}/{nc[0]}",
        form3(t, a["form"]), f"{float(a['cv_rmse']):.2f} ({float(a['cv_rmse_se']):.2f})", f"{c0:.2f}",
        num(h0["rmse"], 2), num(h0["pct_bias_total"], 1),
        form3(t, b["form"]) if b["form"] != a["form"] else "same", f"{c1:.2f}", num(h1["rmse"], 2), num(h1["pct_bias_total"], 1)]) + r" \\")
L += [r"\hline", r"\end{tabular}%", "}", r"\vspace{2mm}", r"\begin{flushleft}", r"\footnotesize",
      r"Set: Y, \gls{yolo}v12; C, Cascade R-CNN; R, RT-DETR-L. Within 1 SE: candidate forms whose five-fold block cross-validated RMSE was within one standard error of the lowest. "
      r"$r = \log(1 + \text{raw detections})$; $\ell = \log(1 + \sum p_{\text{detection}})$; $f$, fraction of detections below the regional $F_1$ threshold; $b$, median box size (square root of area, mm); $o$, fraction of boxes overlapping another box; $n$, median nearest-neighbour distance between box centres divided by box size. "
      r"For several detectors, $L$ is the mean of their $\ell$, $D = \ell_Y - \ell_C$ (or the analogous difference for the pair), and $H_1 = \ell_Y - \ell_C$, $H_2 = \ell_R - (\ell_Y + \ell_C)/2$; crowding terms use \gls{yolo}v12 boxes. "
      r"Conc.: worst-case concurvity of any smooth term (0, none; 1, a term is an exact function of the others). RMSE (scallops per image) and total bias are on the Holdout test. "
      r"One-SE rule: the simplest listed form within one standard error of the lowest CV RMSE, as used for the detection-level models.",
      r"\end{flushleft}", r"\end{table}"]
(out / "tableS4_imagelevel_selection.tex").write_text("\n".join(L) + "\n")

# ---------------- Table S5: concurvity of the submitted forms ----------------
OLD = [("M10_YoloLogDiff", r"$s(\ell_Y, |\Delta|)$"), ("M11_Additive", r"$s(\ell_Y) + s(\ell_C) + s(|\Delta|)$"),
       ("M12_Tensor", r"$te(\ell_Y, \ell_C)$"), ("M13_TensorDiff", r"$ti(\ell_Y, \ell_C) + |\Delta|$"),
       ("M14_Champion", r"$s(\ell_Y, |\Delta|) + s(\ell_C, |\Delta|)$"), ("C10_TensorCrowd", r"$te(\ell_Y, \ell_C) + s(o) + s(n)$")]
L = [r"\begin{table}[htbp]", r"\centering",
     r"\caption{Worst-case concurvity of the two-detector image-level forms in the original submission, refitted on Calibration-fitting.}",
     r"\label{tab:old_concurvity}", r"\begin{tabular}{llcc}", r"\hline",
     r"\textbf{Model} & \textbf{Form} & \textbf{Georges Bank} & \textbf{Mid-Atlantic Bight} \\ \hline"]
for k, f in OLD:
    m = f"Submitted form {k}"
    L.append(f"{k.split('_')[0]} & {f} & {float(CO[('GB', m)]['worst_concurvity']):.2f} & {float(CO[('MAB', m)]['worst_concurvity']):.2f}" + r" \\")
L += [r"\hline", r"\end{tabular}", r"\vspace{2mm}", r"\begin{flushleft}", r"\footnotesize",
      r"$|\Delta| = \log(1 + |\sum p_C - \sum p_Y|)$, the inter-detector disagreement term of the submitted models, which is a function of $\ell_Y$ and $\ell_C$ and drops the sign of the disagreement. M14 was the submitted selected model (SM14). "
      r"Single-term forms (M10, M12) have no concurvity by construction; they are shown for completeness.",
      r"\end{flushleft}", r"\end{table}"]
(out / "tableS5_concurvity.tex").write_text("\n".join(L) + "\n")
print("wrote", *sorted(p.name for p in out.glob("table*")))
