# ======================================================================
# make_tableS7_covariates.R - Table S7 (R1.9): every calibration-model input,
# with its unit, source, stage and range by region (median [5th-95th percentile])
# Run from s/:  Rscript tables/make_tableS7_covariates.R
# Needs: ../data/processed/eval_2226.rds, ../data/processed/n3/n3_results.rds
# Writes tables/tableS7_covariates.tex
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
models <- readRDS("../data/processed/eval_2226.rds")
y <- models[[intersect(c("YOLO", "YOLOv12"), names(models))[1]]]
img <- y$img %>% mutate(region = as.character(region))
det <- y$det %>% mutate(region = as.character(region))
n3i <- readRDS("../data/processed/n3/n3_results.rds")$img %>% mutate(region = as.character(region))
# Stage-2 shares and spacing are only informative when an image has >= 2 detections
# (with one detection overlap = 0, spacing = cap 5, low-confidence share 0 or 1), so their
# ranges are given over those images; box size over images with >= 1 detection.
n3_1 <- n3i %>% filter(n_det_all_Y >= 1); n3_2 <- n3i %>% filter(n_det_all_Y >= 2)
share2 <- n3i %>% group_by(region) %>% summarise(p1 = mean(n_det_all_Y >= 1), p2 = mean(n_det_all_Y >= 2))
pc <- function(r, v) sprintf("%.0f\\%%", 100 * share2[[v]][share2$region == r])

clean <- function(v) { v <- as.numeric(v); v[v <= -90] <- NA; v }
rng <- function(v, d) { v <- clean(v); q <- quantile(v, c(.5, .05, .95), na.rm = TRUE)
  f <- function(x) formatC(x, format = "f", digits = d, big.mark = ",")
  sprintf("%s [%s, %s]", f(q[1]), f(q[2]), f(q[3])) }
row <- function(df, col, d) sapply(c("GB", "MAB"), function(r) rng(df[[col]][df$region == r], d))

# name, symbol/term, unit, source, stage, data frame, column, digits
V <- list(
  list("Detection confidence", "$c$", "0--1", "Detector output", "1", det, "conf", 2),
  list("Box area", "box area", "px$^2$", "Detector output", "1", det, "boxsize", 0),
  list("Vehicle altitude", "altitude", "m", "HabCam altimeter", "1", img, "altitude", 2),
  list("Field of view", "FOV", "m$^2$", "Altitude, roll and pitch", "1", img, "field_of_view_sq_meter", 2),
  list("Bottom depth", "depth", "m", "Altimeter + CTD depth", "1", img, "bottom_depth", 0),
  list("Water temperature", "temperature", "$^\\circ$C", "CTD", "1", img, "ctd_temperature_celsius", 1),
  list("Optical backscatter", "backscatter", "NTU", "Fluorometer", "1", img, "fluorometer_backscatter_ntu", 0),
  list("Latitude", "latitude", "$^\\circ$N", "Vehicle navigation", "1", img, "latitude", 2),
  list("Longitude", "longitude", "$^\\circ$E", "Vehicle navigation", "1", img, "longitude", 2),
  list("Stage-1 count", "$\\hat{N}$ ($\\ell = \\log(1 + \\hat{N})$)", "scallops", "Sum of stage-1 probabilities", "2", n3i, "pred_Y", 1),
  list("Low-confidence share", "$f$", "0--1", "Share of detections below the $F_1$ threshold$^{\\text{b}}$", "2", n3_2, "lowconf_frac_Y", 2),
  list("Median box size", "$b$", "mm", "Square root of box area, scaled by mm per pixel$^{\\text{a}}$", "2", n3_1, "med_box_mm_Y", 0),
  list("Overlap share", "$o$", "0--1", "Share of boxes overlapping another box$^{\\text{b}}$", "2", n3_2, "overlap_frac_Y", 2),
  list("Nearest-neighbour spacing", "$n$", "ratio", "Median centre distance divided by box size (capped at 5)$^{\\text{b}}$", "2", n3_2, "nn_ratio_Y", 2))

body <- vapply(V, function(v) {
  r <- row(v[[6]], v[[7]], v[[8]])
  sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", v[[1]], v[[2]], v[[3]], v[[4]], v[[5]], r[["GB"]], r[["MAB"]])
}, character(1))
body <- append(body, "\\midrule", after = 9)
L <- c("\\begin{table}[htbp]", "\\centering",
  "\\caption{Inputs to the calibration models, with their units, sources and ranges by region.}",
  "\\label{tab:covariates}", "\\resizebox{\\textwidth}{!}{%", "\\begin{tabular}{lllllcc}", "\\toprule",
  "\\textbf{Variable} & \\textbf{Term} & \\textbf{Unit} & \\textbf{Source} & \\textbf{Stage} & \\textbf{Georges Bank} & \\textbf{Mid-Atlantic Bight} \\\\",
  "\\midrule", body, "\\bottomrule", "\\end{tabular}%", "}", "\\vspace{2mm}", "\\begin{flushleft}", "\\footnotesize",
  paste0("Ranges: median [5th, 95th percentile] over Calibration-set images (stage-1 environmental inputs and stage-2 features) or detections (confidence and box area), for \\gls{yolo}v12. ",
         "Sensor values recorded as missing ($\\le -90$) were replaced by the cruise and region median before model fitting and are excluded here. ",
         "Stage-2 features are computed from \\gls{yolo}v12 boxes in each image; images without detections take 0 for shares and 5 for spacing. ",
         sprintf("$^{\\text{a}}$Images with at least one detection (Georges Bank %s, Mid-Atlantic Bight %s of Calibration-set images). ", pc("GB", "p1"), pc("MAB", "p1")),
         sprintf("$^{\\text{b}}$Images with at least two detections (%s and %s), because with a single detection the overlap share is 0, the spacing takes its cap and the low-confidence share is 0 or 1. ", pc("GB", "p2"), pc("MAB", "p2")),
         "Distributions are shown in Fig.~\\ref{fig:final_gams}."),
  "\\end{flushleft}", "\\end{table}")
dir.create("tables", showWarnings = FALSE)
writeLines(L, "tables/tableS7_covariates.tex")
cat("Wrote tables/tableS7_covariates.tex\n"); cat(body, sep = "\n")
