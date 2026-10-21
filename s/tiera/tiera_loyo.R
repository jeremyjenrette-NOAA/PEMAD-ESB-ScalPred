# ======================================================================
# tiera_loyo.R - Tier A: leave-one-year-out test of the GAM calibration
#
# Question: does the error-environment relationship learned by the GAMs hold
# in a survey year they never saw? Detectors are the current ones (their
# Detector-training set includes every year), so only the GAMs hold the year
# out. Tier B (s/tierb/) is the stricter test that also holds 2026 out of the
# detector.
#
# For each held-out year Y in 2022, 2023, 2024, 2026:
#   * stage 1 (detection GAMs, YOLOv12 and Cascade R-CNN) and stage 2
#     (image-level Tweedie GAMs) are refit on Calibration-fitting images from
#     the OTHER years (stage 1 cross-fitted on the existing block folds, so
#     stage-2 training inputs are out-of-fold);
#   * F1 thresholds are tuned on the same training images;
#   * everything is scored on Year Y's Holdout test images.
# The "pooled" fit (all years' Calibration-fitting, same fixed forms) is
# scored on the same images, so each year has a with/without comparison.
#
# Model forms are FIXED from the main CV selection (no re-selection):
#   stage 1  YOLOv12  GB M15_All_te    MAB M01_ConfOnly
#            Cascade  GB M12_PB_Depth  MAB M10_teBoxsize
#            (random effects dropped: a held-out year's cruises are unseen,
#             so cruise terms would be integrated out at prediction anyway)
#   stage 2  one detector  GB C07_teLowConf  MAB C08_Structure (best YOLO-only)
#            two detectors GB M06_MeanLog    MAB C10_TensorCrowd (best forms using both,
#                                            no |difference| terms - concurvity)
#            selected      GB C07_teLowConf  MAB C10_TensorCrowd (main analysis)
# RT-DETR-L joins as a third detector once its stage-1 form is chosen (N3).
#
# Run from s/:  Rscript tiera/tiera_loyo.R      (about 1.5-2.5 h on a laptop)
# Needs: ../data/processed/eval_2226.rds (makedat.R), ../data/raw/dataset_split_2226.csv
# Outputs: ../data/processed/tiera/*.csv|rds, ../figures/tiera/*.png
# ======================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(janitor)
  library(mgcv); library(ggplot2); library(purrr)
})
DET_ENGINE <- "gam"
source("./gamfunc.R")   # sources cvfunc.R

OUT_DIR <- "../data/processed/tiera"; FIG_DIR <- "../figures/tiera"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
N_FOLDS <- 5
N_BOOT  <- as.integer(Sys.getenv("TIERA_BOOT", "500"))
SEED    <- 2226
HOLD_YEARS <- c(2022, 2023, 2024, 2026)

log_file <- file.path(OUT_DIR, "tiera_run_log.txt"); cat("", file = log_file)
say <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE) }
say("Tier A leave-one-year-out run ", format(Sys.time()))

# ----------------------------------------------------------------------
# Fixed forms
# ----------------------------------------------------------------------
S1 <- list(
  YOLOv12 = list(
    GB  = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) +
            te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + te(conf, bottom_depth, k = 5),
    MAB = truedetect ~ s(conf, k = 5)),
  `Cascade R-CNN` = list(
    GB  = truedetect ~ te(conf, altitude, k = 5) + s(latitude, longitude, k = 7) +
            te(conf, field_of_view_sq_meter, k = 7) + te(conf, boxsize, k = 7) + s(bottom_depth, k = 5),
    MAB = truedetect ~ te(conf, boxsize, k = 5) + te(conf, altitude, k = 5) + s(latitude, longitude, k = 7))
)
S2_FORMS <- list(
  C07_teLowConf   = n_annotations ~ te(log_yolo_pred, lowconf_frac, k = c(5, 4)),
  C08_Structure   = n_annotations ~ s(log_yolo_pred) + s(med_box_mm, k = 5) + s(lowconf_frac, k = 5) +
                      s(overlap_frac, k = 5) + s(nn_ratio, k = 5),
  M06_MeanLog     = n_annotations ~ s(log_pred_mean),
  C10_TensorCrowd = n_annotations ~ te(log_yolo_pred, log_cascade_pred) + s(overlap_frac, k = 5) +
                      s(nn_ratio, k = 5)
)
S2_VARIANTS <- list(
  "Two-stage, one detector"  = c(GB = "C07_teLowConf", MAB = "C08_Structure"),
  "Two-stage, two detectors" = c(GB = "M06_MeanLog",   MAB = "C10_TensorCrowd"),
  "Two-stage, selected"      = c(GB = "C07_teLowConf", MAB = "C10_TensorCrowd")
)

# ----------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------
models <- readRDS("../data/processed/eval_2226.rds")
names(models)[names(models) == "YOLO"] <- "YOLOv12"
names(models)[names(models) == "Cascade"] <- "Cascade R-CNN"
meta <- read.csv("../data/raw/dataset_split_2226.csv") %>% clean_names() %>%
  mutate(image_id = str_remove(imagename, "\\.[A-Za-z0-9]+$")) %>%
  distinct(image_id, .keep_all = TRUE)
for (d in names(models)) {
  if (!"year" %in% names(models[[d]]$img)) stop("eval_2226.rds has no year column; re-run makedat.R")
  if (!"fold" %in% names(models[[d]]$det)) stop("eval_2226.rds has no fold column; re-run makedat.R")
}
counts <- models$YOLOv12$img %>% count(year, region, dataset, wt = NULL, name = "images")
print(counts)

f1_threshold <- function(det, img) {
  n_gt <- sum(img$n_annotations)
  map_dfr(seq(0.10, 0.95, by = 0.01), function(t) {
    k <- det$conf >= t; tp <- sum(det$truedetect & k); fp <- sum(!det$truedetect & k)
    p <- tp / max(tp + fp, 1); r <- tp / max(n_gt, 1)
    tibble(conf = t, f1 = ifelse(p + r > 0, 2 * p * r / (p + r), 0))
  }) %>% filter(f1 == max(f1)) %>% slice_max(conf, n = 1) %>% pull(conf)
}

# Stage 1 for one detector on one training/test configuration
run_stage1 <- function(d, hold_year) {
  m <- models[[d]]
  keep_img <- function(x) {
    if (is.na(hold_year)) x %>% filter(dataset %in% c("test_GAM_train", "test_GAM_test"))
    else x %>% filter((dataset == "test_GAM_train" & year != hold_year) |
                      (dataset == "test_GAM_test"  & year == hold_year))
  }
  img <- keep_img(m$img); det <- keep_img(m$det)
  f1 <- list(GB = NA_real_, MAB = NA_real_)
  for (r in c("GB", "MAB"))
    f1[[r]] <- f1_threshold(det %>% filter(region == r, dataset == "test_GAM_train"),
                            img %>% filter(region == r, dataset == "test_GAM_train"))
  train <- det_fit_rows(det %>% filter(dataset == "test_GAM_train"))
  gams <- lapply(c(GB = "GB", MAB = "MAB"), function(r) {
    fit <- try_fit_det_gam(S1[[d]][[r]], train %>% filter(region == r), engine = DET_ENGINE,
                           label = paste(d, r, "final"))
    if (!fit_ok(fit)) stop("Stage-1 fit unusable: ", d, " ", r, " (held-out year ", hold_year, ")")
    fit
  })
  cf <- crossfit_detection_gams(det, fixed_formulas = S1[[d]], K = N_FOLDS,
                                nested_selection = FALSE, engine = DET_ENGINE)
  ev <- suppressWarnings(test_calibration_gams(gams, cf$det, img, f1, test_pred_source = "full"))
  ev$f1 <- f1
  ev
}

# Image-level predictions on the held-out images for every method
run_config <- function(hold_year) {
  lab <- if (is.na(hold_year)) "pooled" else as.character(hold_year)
  say("\n=== Configuration: ", lab, " ===")
  t0 <- Sys.time()
  ev <- lapply(setNames(names(S1), names(S1)), run_stage1, hold_year = hold_year)
  img <- build_synergy_dataset(ev$YOLOv12, ev$`Cascade R-CNN`, meta) %>%
    mutate(log_ratio = log((pred_yolo + 1) / (pred_cascade + 1)))
  img <- add_crowding_features(img, image_detection_features(ev$YOLOv12$det_eval))
  f1n <- function(e, nm) e$img_eval %>% select(image_id, !!nm := predicted_f1_number)
  img <- img %>% left_join(f1n(ev$YOLOv12, "yolo_f1"), by = "image_id") %>%
    left_join(f1n(ev$`Cascade R-CNN`, "cascade_f1"), by = "image_id")
  out <- img %>% filter(dataset == "test_GAM_test") %>%
    select(image_id, region, year, transect_block_id, n_annotations, field_of_view_sq_meter,
           yolo_f1, cascade_f1, pred_yolo, pred_cascade)
  for (v in names(S2_VARIANTS)) {
    ids <- S2_VARIANTS[[v]]
    g <- tryCatch(fit_synergy_gams(img, S2_FORMS[[ids[["GB"]]]], S2_FORMS[[ids[["MAB"]]]]),
                  error = function(e) { say("  !! ", v, " failed: ", conditionMessage(e)); NULL })
    if (is.null(g)) { out[[v]] <- NA_real_; next }
    pr <- test_synergy_gams(g, img)$img_eval %>% select(image_id, pred_synergy)
    out <- out %>% left_join(pr %>% rename(!!v := pred_synergy), by = "image_id")
  }
  thr <- bind_rows(lapply(names(ev), function(d) tibble(detector = d, region = c("GB", "MAB"),
                                                         f1_threshold = unlist(ev[[d]]$f1))))
  say(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(pred = out %>% mutate(fit = lab, .before = 1), thr = thr %>% mutate(fit = lab))
}

configs <- c(list(run_config(NA)), lapply(HOLD_YEARS, run_config))
preds <- bind_rows(lapply(configs, `[[`, "pred"))
thrs  <- bind_rows(lapply(configs, `[[`, "thr"))
saveRDS(preds, file.path(OUT_DIR, "tiera_holdout_predictions.rds"))
write.csv(thrs, file.path(OUT_DIR, "tiera_f1_thresholds.csv"), row.names = FALSE)

# ----------------------------------------------------------------------
# Metrics: each held-out year scored with its leave-year-out fit and with
# the pooled fit on the same images
# ----------------------------------------------------------------------
METHODS <- c("YOLOv12 F1" = "yolo_f1", "Cascade F1" = "cascade_f1",
             "YOLOv12 det-GAM" = "pred_yolo", "Cascade det-GAM" = "pred_cascade",
             setNames(names(S2_VARIANTS), names(S2_VARIANTS)))
pooled <- preds %>% filter(fit == "pooled")
eval_tab <- bind_rows(
  preds %>% filter(fit != "pooled") %>% mutate(fit_type = "Leave-year-out", held_year = as.integer(fit)),
  bind_rows(lapply(HOLD_YEARS, function(y) pooled %>% filter(year == y) %>%
                     mutate(fit_type = "Pooled (all years)", held_year = y)))
)
metric_row <- function(y, yhat) tibble(rmse = sqrt(mean((y - yhat)^2)), mae = mean(abs(y - yhat)),
                                       r2_1to1 = r2_1to1(y, yhat),
                                       pct_bias_total = 100 * (sum(yhat) - sum(y)) / sum(y),
                                       calibration_slope = unname(coef(lm(y ~ yhat))[2]),
                                       total_true = sum(y), total_pred = sum(yhat), n_images = length(y))
metrics <- eval_tab %>% group_by(fit_type, held_year, region) %>%
  group_modify(function(d, k) map_dfr(names(METHODS), function(m)
    if (all(is.na(d[[METHODS[[m]]]]))) NULL else
      metric_row(d$n_annotations, d[[METHODS[[m]]]]) %>% mutate(method = m, .before = 1))) %>%
  ungroup()
write.csv(metrics, file.path(OUT_DIR, "tiera_holdout_metrics.csv"), row.names = FALSE)

# Block-bootstrap intervals for the leave-year-out fits, per held-out year
boot <- bind_rows(lapply(HOLD_YEARS, function(y) {
  d <- eval_tab %>% filter(fit_type == "Leave-year-out", held_year == y)
  if (!nrow(d)) return(NULL)
  cols <- METHODS[sapply(METHODS, function(c) !all(is.na(d[[c]])))]
  block_bootstrap_metrics(d, cols, block = "transect_block_id", ref = "Two-stage, selected",
                          B = N_BOOT, seed = SEED) %>% mutate(held_year = y, .before = 1)
}))
write.csv(boot, file.path(OUT_DIR, "tiera_holdout_bootstrap.csv"), row.names = FALSE)

# Cost of holding the year out: leave-year-out minus pooled, same images
cost <- metrics %>% select(fit_type, held_year, region, method, rmse, pct_bias_total) %>%
  pivot_wider(names_from = fit_type, values_from = c(rmse, pct_bias_total)) %>%
  mutate(d_rmse = `rmse_Leave-year-out` - `rmse_Pooled (all years)`,
         d_abs_bias = abs(`pct_bias_total_Leave-year-out`) - abs(`pct_bias_total_Pooled (all years)`))
write.csv(cost, file.path(OUT_DIR, "tiera_cost_of_holdout.csv"), row.names = FALSE)

say("\nLeave-year-out holdout metrics:")
print(metrics %>% filter(fit_type == "Leave-year-out") %>%
        select(held_year, region, method, rmse, r2_1to1, pct_bias_total, calibration_slope, n_images),
      n = Inf, width = Inf)
say("\nCost of holding the year out (leave-year-out minus pooled):")
print(cost %>% select(held_year, region, method, d_rmse, d_abs_bias), n = Inf, width = Inf)

# ----------------------------------------------------------------------
# Figure: total bias by held-out year, leave-year-out vs pooled
# ----------------------------------------------------------------------
pdat <- metrics %>% filter(method %in% c("YOLOv12 F1", "Cascade F1", "Two-stage, one detector",
                                         "Two-stage, two detectors", "Two-stage, selected")) %>%
  left_join(boot %>% select(held_year, region, method, bias_lo, bias_hi), by = c("held_year", "region", "method")) %>%
  mutate(bias_lo = if_else(fit_type == "Leave-year-out", bias_lo, NA_real_),
         bias_hi = if_else(fit_type == "Leave-year-out", bias_hi, NA_real_),
         held_year = factor(held_year))
p <- ggplot(pdat, aes(held_year, pct_bias_total, colour = method, shape = fit_type)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = bias_lo, ymax = bias_hi), width = 0, position = position_dodge(0.6), na.rm = TRUE) +
  geom_point(position = position_dodge(0.6), size = 2) +
  facet_wrap(~ region) +
  labs(x = "Held-out survey year", y = "Holdout test total bias (%)", colour = NULL, shape = NULL,
       title = "Tier A: GAMs refit without the held-out year") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(FIG_DIR, "tiera_bias_by_year.png"), p, width = 8, height = 4.8, dpi = 300, bg = "white")
say("\nDone. Outputs in ", OUT_DIR, " and ", FIG_DIR)
