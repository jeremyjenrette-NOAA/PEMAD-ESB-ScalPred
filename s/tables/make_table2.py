#!/usr/bin/env python3
"""Table 2 (detection-level GAM selection) from the block-fold CV results, plus an
intercept-only row per region (R1.10/R1.11). Run from s/:  python tables/make_table2.py
Writes tables/table2_detgam_cv.tex and tables/table2_detgam_cv.csv."""
import re
import numpy as np
import pandas as pd

FORMS = {
    "M00": "Intercept only (regional base rate)",
    "M01": r"\textit{s}(Conf)",
    "M02": r"\textit{s}(Conf) + \textit{s}(Lat, Lon)",
    "M03": r"M02 + \textit{s}(Depth)",
    "M04": r"M03 + \textit{s}(Alt) + \textit{s}(FOV)",
    "M05": r"M03 + \textit{s}(Alt) + \textit{s}(Box)",
    "M06": r"M03 + \textit{s}(Backscatter)",
    "M07": r"\textit{te}(Conf, Depth) + \textit{s}(Lat, Lon)",
    "M08": r"\textit{te}(Conf, Alt) + \textit{s}(Lat, Lon) + \textit{s}(Depth)",
    "M09": r"\textit{te}(Conf, Alt) + \textit{te}(Conf, FOV) + \textit{s}(Lat, Lon)",
    "M10": r"\textit{te}(Conf, Box) + \textit{te}(Conf, Alt) + \textit{s}(Lat, Lon)",
    "M11": r"\textit{te}(Conf, Alt) + \textit{te}(Conf, FOV) + \textit{te}(Conf, Box) + \textit{s}(Lat, Lon)",
    "M12": r"M11 + \textit{s}(Depth)",
    "M13": r"M11 + \textit{s}(Backscatter)",
    "M14": r"M11 + \textit{s}(Temperature)",
    "M15": r"M11 + \textit{te}(Conf, Depth)",
}
DET = {"YOLOv12": ("../data/processed/YOLOv12_formula_cv_2226.csv",
                   "../data/raw_model/2226data_yolo12n_slurm_598585/autotest2226_yolo12n.csv"),
       "Cascade R-CNN": ("../data/processed/CascadeR-CNN_formula_cv_2226.csv",
                         "../data/raw_model/2226data_cascade_slurm_598695/autotest2226_viame_cascade.csv")}

sid = lambda s: re.sub(r"\.[A-Za-z0-9]+$", "", s)
split = pd.read_csv("../data/raw/dataset_split_2226.csv", usecols=["imagename", "dataset", "longitude"])
excl = pd.read_csv("../data/raw/annotation_exclusions_2226.csv")
excl = set(excl.imagename[excl.status.str.strip().str.lower() == "confirmed"].map(sid))
split["image_id"] = split.imagename.map(sid)
split = split[~split.image_id.isin(excl)].drop_duplicates("image_id")
split["region"] = np.where(split.longitude >= -71, "GB", "MAB")
train = split[split.dataset == "test_GAM_train"][["image_id", "region"]]

rows = []
for det, (cvf, autof) in DET.items():
    cv = pd.read_csv(cvf)
    cv["mid"] = cv.candidate.str.slice(0, 3)
    a = pd.read_csv(autof)
    a = a[a.Spname.str.lower() == "scallop"].assign(image_id=lambda x: x.Imagename.map(sid))
    a = a.merge(train, on="image_id")
    a["y"] = a.truedetect.astype(str).str.lower() == "true"
    for reg, g in a.groupby("region"):
        c = cv[cv.region == reg]
        n_check = int(c.n_det.iloc[0])
        assert len(g) == n_check, f"{det} {reg}: {len(g)} detections vs {n_check} in CV file"
        p = g.y.mean()
        null_ll = -(p * np.log(p) + (1 - p) * np.log(1 - p))
        per_img = g.groupby("image_id").agg(n=("y", "size"), tp=("y", "sum"))
        null_rmse = float(np.sqrt(np.mean((per_img.tp - p * per_img.n) ** 2)))
        best_ll = c.log_loss.min()
        conf_ll = c.loc[c.mid == "M01", "log_loss"].iloc[0]
        for _, r in c.iterrows():
            rows.append(dict(detector=det, region=reg, model=r.mid, d_ll=1e3 * r.d_logloss, se=1e3 * r.se,
                             rmse=r.rmse_tp, within_1se=bool(r.within_1se), chosen=bool(r.chosen),
                             edf=r.complexity))
        rows.append(dict(detector=det, region=reg, model="M00", d_ll=1e3 * (null_ll - best_ll), se=np.nan,
                         rmse=null_rmse, within_1se=False, chosen=False, edf=1.0,
                         conf_share=(null_ll - conf_ll) / (null_ll - best_ll)))
t = pd.DataFrame(rows)
t.to_csv("tables/table2_detgam_cv.csv", index=False)
share = t[t.model == "M00"].set_index(["detector", "region"]).conf_share
print("Share of the null-to-best CV log-loss reduction achieved by confidence alone (M01):")
print((100 * share).round(1).to_string())

def cell(r):
    if r.model == "M00":
        s = f"{r.d_ll:.1f}"
    elif r.d_ll == 0:
        s = "0"
    else:
        s = f"{r.d_ll:.1f} ({r.se:.1f})"
    s += f" & {r.rmse:.2f}"
    if r.chosen:
        s = " & ".join(rf"\textbf{{{x.strip()}}}" for x in s.split("&"))
    elif r.within_1se:
        s = s.replace(" & ", r"$^{\dagger}$ & ", 1)
    return s

cols = [("YOLOv12", "GB"), ("YOLOv12", "MAB"), ("Cascade R-CNN", "GB"), ("Cascade R-CNN", "MAB")]
lines = [r"\begin{table}[htbp]", r"\centering",
         r"\caption{Detection-level \gls{gam} selection by five-fold block cross-validation on the Calibration-fitting subset.}",
         r"\label{tab:det_gam_selection}", r"\resizebox{\textwidth}{!}{%", r"\begin{tabular}{ll" + "cc" * 4 + "}", r"\toprule",
         r" & & \multicolumn{4}{c}{\textbf{YOLOv12}} & \multicolumn{4}{c}{\textbf{Cascade R-CNN}} \\",
         r"\cmidrule(lr){3-6} \cmidrule(lr){7-10}",
         r" & & \multicolumn{2}{c}{GB} & \multicolumn{2}{c}{MAB} & \multicolumn{2}{c}{GB} & \multicolumn{2}{c}{MAB} \\",
         r"\cmidrule(lr){3-4} \cmidrule(lr){5-6} \cmidrule(lr){7-8} \cmidrule(lr){9-10}",
         r"\textbf{Model} & \textbf{Formula} & $\Delta$LL & RMSE & $\Delta$LL & RMSE & $\Delta$LL & RMSE & $\Delta$LL & RMSE \\",
         r"\midrule"]
for m in FORMS:
    cells = [cell(t[(t.detector == d) & (t.region == rg) & (t.model == m)].iloc[0]) for d, rg in cols]
    lines.append(f"{m} & {FORMS[m]} & " + " & ".join(cells) + r" \\")
    if m in ("M00", "M06", "M10"):
        lines.append(r"\midrule")
lines += [r"\bottomrule", r"\end{tabular}%", r"}", r"\vspace{2mm}", r"\begin{flushleft}", r"\footnotesize",
          r"$\Delta$LL: out-of-fold log loss minus the lowest log loss in that region ($\times 10^{-3}$; block standard error in brackets). "
          r"RMSE: out-of-fold error of the summed detection probabilities against true-positive counts per image. "
          r"Bold: selected model (simplest model within one standard error of the lowest log loss); $\dagger$: also within one standard error. "
          r"Models were fitted separately for each region and detector; M00 is the intercept-only (regional base-rate) model. "
          r"For RT-DETR-L (Fig.~4; Table~S5), the same procedure selected M15 on Georges Bank and M10 in the Mid-Atlantic Bight.",
          r"\end{flushleft}", r"\end{table}"]
open("tables/table2_detgam_cv.tex", "w").write("\n".join(lines) + "\n")
print("Wrote tables/table2_detgam_cv.tex")
