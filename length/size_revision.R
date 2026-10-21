# ==============================================================================
# size_revision.R - shell-height validation for the revision (N1, N2, N14)
#
#   N1  matched-pair accuracy: predicted vs manual shell height, OLS slope/R2,
#       Bland-Altman bias and limits, MAE and bias by size bin and region
#   N2  < 40 mm vs >= 40 mm: per-class counts and size-frequency overlap (OVL)
#   N14 counts per size bin (0-40, 40-80, 80-120, 120-180 mm), manual vs automated,
#       and the 2 x 2 confusion matrix across the 40 mm cutoff
#
# Pairs come from match_size_pairs.py, which repeats the evaluator's Hungarian
# match. (measuring.R joined detections to annotations on `index`, which is the
# detection's rank in the image, so most multi-scallop pairs were wrong; it also
# filtered its "F1 threshold" curve on pred_p instead of confidence.)
#
# Workflows for counts and size-frequency (Holdout test, every detection, not only TPs):
#   F1 threshold     detections with conf >= the regional F1 threshold
#   Calibrated       every detection weighted by its detection-GAM probability pred_p
#   Manual           every annotation (<= 180 mm)
#
# Run from length/:  python match_size_pairs.py   (once; writes size_pairs_2226.csv)
#                    Rscript size_revision.R
# Outputs: ../data/processed/size_revision/*.csv, ../figures/ms_figures/2226_size_*.pdf|png
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(purrr) })

MAX_SH_MM <- 180
SIZE_CUT  <- 40
BIN_EDGES <- c(0, 40, 80, 120, 180)
FREQ_BIN  <- 5
OUT_DIR <- "../data/processed/size_revision"; FIG_DIR <- "../figures/ms_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
MODEL_LEVELS  <- c("YOLOv12", "Cascade R-CNN", "RT-DETR-L")
REGION_LEVELS <- c("GB", "MAB")
size_bin <- function(x) cut(x, BIN_EDGES, include.lowest = TRUE, right = FALSE,
                           labels = c("0-40", "40-80", "80-120", "120-180"))
size_cls <- function(x) factor(if_else(x < SIZE_CUT, "< 40 mm", ">= 40 mm"), levels = c("< 40 mm", ">= 40 mm"))
se <- function(x) sd(x) / sqrt(length(x))

# ------------------------------------------------------------------------------
# Data
# ------------------------------------------------------------------------------
pairs <- read.csv("../data/processed/size_pairs_2226.csv", stringsAsFactors = FALSE) %>%
  filter(true_mm <= MAX_SH_MM) %>%
  mutate(model = factor(model, levels = MODEL_LEVELS), region = factor(region, levels = REGION_LEVELS),
         err = pred_mm - true_mm, mean_mm = (pred_mm + true_mm) / 2,
         bin = size_bin(true_mm), cls_true = size_cls(true_mm), cls_pred = size_cls(pred_mm))

split <- read.csv("../data/raw/dataset_split_2226.csv", stringsAsFactors = FALSE)
hold  <- split$imagename[split$dataset == "test_GAM_test"]
gt <- read.csv("../data/raw/groundtruth2226.csv", stringsAsFactors = FALSE) %>%
  filter(imagename %in% hold, tolower(trimws(label)) == "scallop") %>%
  mutate(true_px = sapply(str_extract_all(geometry_text, "[-+]?[0-9]*\\.?[0-9]+"), function(p) {
           p <- as.numeric(p); if (length(p) >= 4) sqrt((p[3] - p[1])^2 + (p[4] - p[2])^2) else NA }),
         true_mm = true_px * millimeter_per_pixel,
         region = factor(if_else(longitude >= -71, "GB", "MAB"), levels = REGION_LEVELS))
n_over <- sum(gt$true_mm > MAX_SH_MM, na.rm = TRUE)
gt <- gt %>% filter(!is.na(true_mm), true_mm <= MAX_SH_MM)
cat(sprintf("Holdout annotations used: %d (%d over %d mm excluded)\n", nrow(gt), n_over, MAX_SH_MM))

# All Holdout detections with calibrated probabilities (YOLOv12 and Cascade R-CNN)
load_det <- function(f, model) {
  readRDS(f)$det_eval %>%
    filter(dataset == "test_GAM_test") %>%
    transmute(model = model, region = factor(as.character(region), levels = REGION_LEVELS),
              detectid, conf, pred_p, regional_f1_thresh, truedetect,
              pred_mm = ((abs(b_rx - t_lx) + abs(b_ry - t_ly)) / 2) * millimeter_per_pixel)
}
dets <- bind_rows(load_det("../data/processed/YOLOv12_deteval_2226.rds", "YOLOv12"),
                  load_det("../data/processed/CascadeR-CNN_deteval_2226.rds", "Cascade R-CNN")) %>%
  mutate(model = factor(model, levels = MODEL_LEVELS)) %>%
  filter(!is.na(pred_mm))
print(dets %>% group_by(model, region) %>% summarise(f1_thresh = first(regional_f1_thresh), n = n(), .groups = "drop"))

# ------------------------------------------------------------------------------
# N1: matched-pair accuracy
# ------------------------------------------------------------------------------
n1_fit <- function(d) {
  fit <- lm(pred_mm ~ true_mm, data = d); ba <- lm(err ~ mean_mm, data = d)
  tibble(n = nrow(d), slope = coef(fit)[2], intercept = coef(fit)[1], r2 = summary(fit)$r.squared,
         bias = mean(d$err), bias_se = se(d$err), sd_diff = sd(d$err),
         loa_lo = mean(d$err) - 1.96 * sd(d$err), loa_hi = mean(d$err) + 1.96 * sd(d$err),
         mae = mean(abs(d$err)), mae_se = se(abs(d$err)),
         prop_bias_slope = coef(ba)[2], prop_bias_p = summary(ba)$coefficients[2, 4])
}
n1_overall <- pairs %>% group_by(model) %>% group_modify(~ n1_fit(.x)) %>% ungroup() %>% mutate(region = "All", .after = model)
n1_region  <- pairs %>% group_by(model, region) %>% group_modify(~ n1_fit(.x)) %>% ungroup()
n1 <- bind_rows(n1_overall, n1_region %>% mutate(region = as.character(region)))
write.csv(n1, file.path(OUT_DIR, "N1_shell_height_agreement.csv"), row.names = FALSE)
n1_bins <- pairs %>% group_by(model, region, bin) %>%
  summarise(n = n(), bias = mean(err), bias_se = se(err), mae = mean(abs(err)), mae_se = se(abs(err)),
            .groups = "drop")
write.csv(n1_bins, file.path(OUT_DIR, "N1_error_by_size_bin.csv"), row.names = FALSE)
cat("\nN1 agreement (matched pairs):\n"); print(n1 %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = Inf)
cat("\nN1 error by size bin:\n"); print(n1_bins %>% mutate(across(where(is.numeric), ~ round(.x, 2))), n = Inf)

lab <- n1_region %>% mutate(txt = sprintf("n = %s\nslope %.2f, R² %.2f\nbias %.1f mm, MAE %.1f mm",
                                         format(n, big.mark = ","), slope, r2, bias, mae))
p_sc <- ggplot(pairs, aes(true_mm, pred_mm)) +
  geom_point(alpha = 0.08, size = 0.5, colour = "#2a78d6") +
  geom_abline(linetype = 2, colour = "grey40") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#eb6834", linewidth = 0.6) +
  geom_vline(xintercept = SIZE_CUT, linetype = 3, colour = "grey50") +
  geom_text(data = lab, aes(x = 5, y = 175, label = txt), hjust = 0, vjust = 1, size = 2.6, colour = "grey15") +
  facet_grid(region ~ model) + coord_equal(xlim = c(0, 180), ylim = c(0, 180)) +
  labs(x = "Manual shell height (mm)", y = "Predicted shell height (mm)") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIG_DIR, "2226_size_scatter.pdf"), p_sc, width = 7.5, height = 5.4)
ggsave(file.path(FIG_DIR, "2226_size_scatter.png"), p_sc, width = 7.5, height = 5.4, dpi = 300, bg = "white")

p_ba <- ggplot(pairs, aes(mean_mm, err)) +
  geom_point(alpha = 0.08, size = 0.5, colour = "#2a78d6") +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(data = n1_region, aes(yintercept = bias), colour = "#eb6834") +
  geom_hline(data = n1_region, aes(yintercept = loa_lo), colour = "#eb6834", linetype = 2) +
  geom_hline(data = n1_region, aes(yintercept = loa_hi), colour = "#eb6834", linetype = 2) +
  facet_grid(region ~ model) + coord_cartesian(xlim = c(0, 180), ylim = c(-40, 40)) +
  labs(x = "Mean of predicted and manual shell height (mm)", y = "Predicted − manual (mm)") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIG_DIR, "2226_size_bland_altman.pdf"), p_ba, width = 7.5, height = 4.6)
ggsave(file.path(FIG_DIR, "2226_size_bland_altman.png"), p_ba, width = 7.5, height = 4.6, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# N14 (and N2 counts): counts per size bin and per 40 mm class, all detections
# ------------------------------------------------------------------------------
auto_counts <- function(binfun) {
  bind_rows(
    dets %>% filter(conf >= regional_f1_thresh) %>% mutate(k = binfun(pred_mm)) %>%
      count(model, region, k, name = "count") %>% mutate(workflow = "F1 threshold"),
    dets %>% mutate(k = binfun(pred_mm)) %>% group_by(model, region, k) %>%
      summarise(count = sum(pred_p, na.rm = TRUE), .groups = "drop") %>% mutate(workflow = "Calibrated (sum of p)")
  ) %>% filter(!is.na(k))
}
man_counts <- function(binfun) gt %>% mutate(k = binfun(true_mm)) %>% count(region, k, name = "manual")
count_table <- function(binfun) {
  auto_counts(binfun) %>%
    left_join(man_counts(binfun), by = c("region", "k")) %>%
    mutate(pct_diff = 100 * (count - manual) / manual) %>%
    arrange(model, region, workflow, k)
}
n14_bins <- count_table(size_bin) %>% rename(size_bin = k)
n2_cls   <- count_table(size_cls) %>% rename(size_class = k)
write.csv(n14_bins, file.path(OUT_DIR, "N14_counts_by_size_bin.csv"), row.names = FALSE)
write.csv(n2_cls,   file.path(OUT_DIR, "N2_counts_by_40mm_class.csv"), row.names = FALSE)
cat("\nN14 counts by size bin (all Holdout detections):\n")
print(n14_bins %>% mutate(count = round(count), pct_diff = round(pct_diff, 1)), n = Inf)
cat("\nN2 counts by 40 mm class:\n")
print(n2_cls %>% mutate(count = round(count), pct_diff = round(pct_diff, 1)), n = Inf)

# 2 x 2 confusion across the 40 mm cutoff (matched pairs)
conf40 <- pairs %>% count(model, region, cls_true, cls_pred) %>%
  group_by(model, region, cls_true) %>% mutate(pct_of_true = 100 * n / sum(n)) %>% ungroup()
write.csv(conf40, file.path(OUT_DIR, "N14_confusion_40mm.csv"), row.names = FALSE)
cat("\n40 mm confusion (matched pairs; rows = manual class):\n"); print(conf40 %>% mutate(pct_of_true = round(pct_of_true, 1)), n = Inf)

# N2: matched-pair accuracy and size-frequency overlap within each 40 mm class
ovl <- function(a, b, wa = NULL, wb = NULL) {
  br <- seq(0, MAX_SH_MM, by = FREQ_BIN)
  h <- function(x, w) { if (is.null(w)) w <- rep(1, length(x)); v <- tapply(w, cut(x, br, include.lowest = TRUE), sum, na.rm = TRUE); v[is.na(v)] <- 0; v / sum(v) }
  sum(pmin(h(a, wa), h(b, wb)))
}
n2_ovl <- expand_grid(model = c("YOLOv12", "Cascade R-CNN"), region = REGION_LEVELS,
                      size_class = c("All", "< 40 mm", ">= 40 mm")) %>%
  pmap_dfr(function(model, region, size_class) {
    keep <- function(x) if (size_class == "All") rep(TRUE, length(x)) else as.character(size_cls(x)) == size_class
    m  <- gt$true_mm[gt$region == region]; m <- m[keep(m)]
    dd <- dets[dets$model == model & dets$region == region & dets$pred_mm <= MAX_SH_MM, ]
    d1 <- dd[dd$conf >= dd$regional_f1_thresh, ]; d1 <- d1[keep(d1$pred_mm), ]
    d2 <- dd[keep(dd$pred_mm), ]
    tibble(model, region, size_class, OVL_F1 = ovl(m, d1$pred_mm), OVL_calibrated = ovl(m, d2$pred_mm, wb = d2$pred_p))
  })
n2_err <- pairs %>% filter(model != "RT-DETR-L") %>% group_by(model, region, size_class = as.character(cls_true)) %>%
  summarise(n_pairs = n(), bias = mean(err), mae = mean(abs(err)), .groups = "drop")
n2 <- n2_ovl %>% left_join(n2_err, by = c("model", "region", "size_class"))
write.csv(n2, file.path(OUT_DIR, "N2_size_class_metrics.csv"), row.names = FALSE)
cat("\nN2 size-class metrics:\n"); print(n2 %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf, width = Inf)

# ------------------------------------------------------------------------------
# Replacement for the size-frequency figure: region x model, 5 mm bins,
# manual vs F1-threshold vs calibrated (all detections), 40 mm cutoff marked
# ------------------------------------------------------------------------------
freq <- function(x, w = NULL) {
  br <- seq(0, MAX_SH_MM, by = FREQ_BIN); if (is.null(w)) w <- rep(1, length(x))
  v <- tapply(w, cut(x, br, include.lowest = TRUE), sum, na.rm = TRUE); v[is.na(v)] <- 0
  tibble(mid = br[-1] - FREQ_BIN / 2, pct = 100 * as.numeric(v) / sum(v))
}
fd <- bind_rows(
  expand_grid(model = c("YOLOv12", "Cascade R-CNN"), region = REGION_LEVELS) %>%
    pmap_dfr(function(model, region) bind_rows(
      freq(gt$true_mm[gt$region == region]) %>% mutate(workflow = "Manual"),
      { d <- dets[dets$model == model & dets$region == region & dets$pred_mm <= MAX_SH_MM, ]
        bind_rows(freq(d$pred_mm[d$conf >= d$regional_f1_thresh]) %>% mutate(workflow = "F1 threshold"),
                  freq(d$pred_mm, d$pred_p) %>% mutate(workflow = "Calibrated (sum of p)")) }
    ) %>% mutate(model = model, region = region))
) %>% mutate(model = factor(model, levels = MODEL_LEVELS), region = factor(region, levels = REGION_LEVELS),
             workflow = factor(workflow, levels = c("Manual", "F1 threshold", "Calibrated (sum of p)")))
olab <- n2 %>% filter(size_class == "All") %>%
  mutate(model = factor(model, levels = MODEL_LEVELS), region = factor(region, levels = REGION_LEVELS),
         txt = sprintf("OVL F1 %.3f\nOVL calibrated %.3f", OVL_F1, OVL_calibrated))
p_fd <- ggplot(fd, aes(mid, pct)) +
  geom_area(data = filter(fd, workflow == "Manual"), fill = "grey85", colour = "grey35", linewidth = 0.4) +
  geom_line(data = filter(fd, workflow != "Manual"), aes(colour = workflow), linewidth = 0.6) +
  geom_vline(xintercept = SIZE_CUT, linetype = 3, colour = "grey40") +
  geom_text(data = olab, aes(x = 178, y = Inf, label = txt), hjust = 1, vjust = 1.3, size = 2.6, inherit.aes = FALSE) +
  scale_colour_manual(values = c("F1 threshold" = "#2a78d6", "Calibrated (sum of p)" = "#eb6834"),
                      labels = c(expression(F[1]~"threshold"), "Calibrated (sum of p)")) +
  facet_grid(region ~ model, scales = "free_y") +
  labs(x = "Shell height (mm, 5 mm bins)", y = "Share of scallops (%)", colour = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(FIG_DIR, "2226_size_frequency_region.pdf"), p_fd, width = 7.5, height = 5)
ggsave(file.path(FIG_DIR, "2226_size_frequency_region.png"), p_fd, width = 7.5, height = 5, dpi = 300, bg = "white")
cat("\nDone. Tables in", OUT_DIR, "; figures 2226_size_*.pdf in", FIG_DIR, "\n")
