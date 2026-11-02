# ======================================================================
# performance.R - detection performance on the Calibration set (Table 1, Fig. S2)
# for YOLOv12-nano, Cascade R-CNN and RT-DETR-L, by region.
# Run from s/:  Rscript performance.R
# Needs: ../data/processed/eval_2226.rds (makedat.R)
#        ../data/raw_model/2226data_rtdetr_l_gcp_vit_20260924_190833/autotest2226_rtdetr_l.csv
# Writes ../data/processed/detection_performance_2226.csv, tables/table1_detection.tex,
#        ../figures/ms_figures/2226_model_performance_summary.pdf|png
# ======================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(stringr); library(janitor) })
source("./procfunc.R")

models <- readRDS("../data/processed/eval_2226.rds")
names(models)[names(models) == "YOLO"] <- "YOLOv12"
names(models)[names(models) == "Cascade"] <- "Cascade R-CNN"

# RT-DETR-L: same Calibration-set images and matching as the other two (as in n3_ablation.R)
rt_csv <- "../data/raw_model/2226data_rtdetr_l_gcp_vit_20260924_190833/autotest2226_rtdetr_l.csv"
if (file.exists(rt_csv)) {
  img_meta <- models$YOLOv12$img %>% select(-any_of(c("n_auto", "man_density", "auto_density")))
  at <- read.csv(rt_csv) %>% clean_names() %>%
    mutate(image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$"),
           truedetect = tolower(as.character(truedetect)) == "true", spname = tolower(spname)) %>%
    filter(spname == "scallop", image_id %in% img_meta$image_id)
  models[["RT-DETR-L"]] <- list(
    det = at %>% select(-any_of("imagename")) %>% left_join(img_meta %>% select(image_id, region), by = "image_id"),
    img = img_meta)
}
DET <- c("YOLOv12", "Cascade R-CNN", "RT-DETR-L")
models <- models[intersect(DET, names(models))]

out <- evaluate_pr_models(models, stratify_region = TRUE, conf_grid = seq(0, 1, by = 0.0005))
pr_all <- out$pr_all %>%
  mutate(detector = str_remove(model, " (GB|MAB)$"), region = str_extract(model, "(GB|MAB)$"),
         detector = factor(detector, levels = DET),
         region = recode(region, GB = "Georges Bank", MAB = "Mid-Atlantic Bight"))
best <- pr_all %>% group_by(model) %>% filter(f1 == max(f1, na.rm = TRUE)) %>% slice_max(conf, n = 1) %>% ungroup() %>%
  left_join(out$map_summary %>% rename(map010 = Average_Precision), by = "model") %>%
  arrange(detector, region)
write.csv(best %>% select(detector, region, conf, f1, precision, recall, map010, TP, FP, FN),
          "../data/processed/detection_performance_2226.csv", row.names = FALSE)
print(best %>% select(detector, region, conf, f1, precision, recall, map010, TP, FP, FN), width = Inf)

# ---------------- Table 1 ----------------
dir.create("tables", showWarnings = FALSE)
lab <- c("YOLOv12" = "YOLOv12-\\textit{nano}", "Cascade R-CNN" = "Cascade R-CNN", "RT-DETR-L" = "RT-DETR-L")
n3f <- function(x) sprintf("%.3f", x); big <- function(x) format(x, big.mark = ",", trim = TRUE)
rows <- character(0)
for (d in levels(best$detector)) {
  b <- best %>% filter(detector == d); if (!nrow(b)) next
  if (length(rows)) rows <- c(rows, "\\midrule")
  for (i in seq_len(nrow(b))) rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
    if (i == 1) lab[[d]] else "", if (b$region[i] == "Georges Bank") "GB" else "MAB",
    n3f(b$f1[i]), n3f(b$precision[i]), n3f(b$recall[i]), n3f(b$map010[i]), big(b$TP[i]), big(b$FP[i]), big(b$FN[i])))
}
n_img <- nrow(models$YOLOv12$img); n_sc <- sum(models$YOLOv12$img$n_annotations)
writeLines(c("\\begin{table}[htbp]", "\\centering",
  "\\caption{Detection performance at the maximum-$F_1$ confidence threshold, by model and region, on the Calibration set.}",
  "\\label{tab:detection_performance_metrics}", "\\resizebox{\\textwidth}{!}{%", "\\begin{tabular}{llccccrrr}", "\\toprule",
  "\\textbf{Model} & \\textbf{Region} & \\textbf{$F_1$ Score} & \\textbf{Precision} & \\textbf{Recall} & \\textbf{$\\text{mAP}_{0.10}$} & \\textbf{TP} & \\textbf{FP} & \\textbf{FN} \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}%", "}", "\\vspace{2mm}", "\\begin{flushleft}", "\\footnotesize",
  sprintf("The Calibration set comprises %s images and %s scallops. $\\text{mAP}_{0.10}$: mean Average Precision at IoU 0.1. TP, FP and FN: true positive, false positive and (image-level) false negative counts at the maximum-$F_1$ threshold. GB: Georges Bank; MAB: Mid-Atlantic Bight.",
          big(n_img), big(n_sc)),
  "\\end{flushleft}", "\\end{table}"), "tables/table1_detection.tex")

# ---------------- Figure (supplement): PR and F1 curves ----------------
COL <- c("YOLOv12" = "#CB181D", "Cascade R-CNN" = "#2171B5", "RT-DETR-L" = "#238B45")
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
p_pr <- ggplot(pr_all %>% filter(!is.na(precision), !is.na(recall)),
               aes(recall, precision, colour = detector, linetype = region)) +
  geom_path(linewidth = 0.8) + scale_colour_manual(values = COL) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Recall", y = "Precision", colour = NULL, linetype = NULL) + th
p_f1 <- ggplot(pr_all %>% filter(!is.na(f1)), aes(conf, f1, colour = detector, linetype = region)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = best, aes(conf, f1), size = 2, show.legend = FALSE) +
  scale_colour_manual(values = COL) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Confidence threshold", y = expression(F[1] ~ "score"), colour = NULL, linetype = NULL) + th
fig <- (p_pr | p_f1) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
ggsave("../figures/ms_figures/2226_model_performance_summary.pdf", fig, width = 9, height = 4.6)
ggsave("../figures/ms_figures/2226_model_performance_summary.png", fig, width = 9, height = 4.6, dpi = 300, bg = "white")
cat("Wrote Table 1 (tables/table1_detection.tex) and the PR/F1 figure\n")
