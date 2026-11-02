# ======================================================================
# tierb_effort_cv.R - Tier B annotation effort, with a ROTATING holdout
#
# tierb_effort.R thins Calibration-fitting and scores every refit on ONE fixed
# 2026 Holdout test, so a quirk of that split sits in every point of Fig. S5.
# Here the 2026 annotated images are cut into K_OUT folds of whole 20-image
# transect chunks. For each fold: thin the other folds to an annotation rate,
# refit both GAM stages on them (model forms transferred from the main
# analysis, exactly as in tierb_effort.R), predict the held-out fold. Out-of-fold
# predictions are pooled over all folds, so every 2026 image is tested once and
# each (rate, draw) gets one total bias / RMSE over ALL of 2026, with a
# transect-block bootstrap interval, plus the fold-to-fold spread.
#
# Annotation rates match Fig. S5 (1 annotated frame per 1,667 / 833 / 333 / 167 /
# 83 survey frames) plus all training folds (about 1 in 56).
#
# Run from s/ (same folder layout as the repo):
#   Rscript tierb/tierb_effort_cv.R
# Options (environment variables):
#   EFFORT_CV_CORES=4      parallel workers (forked; macOS/Linux). Default 1.
#   EFFORT_CV_FOLDS=10     outer folds. Default 10.
#   EFFORT_CV_REPS=3       thinning draws per rate below 100%. Default 3.
#   EFFORT_CV_VARIANT=full "full" (transferred forms) or "simple". Default full.
#   EFFORT_CV_BOOT=500     bootstrap replicates for pooled intervals. Default 500.
# Each finished fit is cached in OUT/cache/, so an interrupted run resumes.
# About 160 fits; roughly 2 min per fit on one core.
#
# Needs (relative to s/):
#   tierb/tierb_common.R, gamfunc.R, cvfunc.R
#   ../data/raw/dataset_split_tierb2026.csv
#   ../data/raw/allimg_2026.RData
#   ../data/raw_model/2224scallop_yolo12n_slurm_598580/tierb2026/autotest_tierb2026_yolo12n.csv
#   R packages: dplyr tidyr stringr janitor mgcv ggplot2 purrr (parallel is base R)
# Writes ../data/processed/tierb2026/effort_cv/:
#   effort_cv_predictions.csv  out-of-fold image predictions (every rate x draw)
#   effort_cv_fold_metrics.csv metrics per held-out fold
#   effort_cv_pooled.csv       pooled metrics per rate x draw, with 95% block-bootstrap CI
#   effort_cv_summary.csv      per rate: median and range over draws, fold SD
#   effort_cv_fits.csv         forms actually used (fallbacks) and training sizes per fit
#   tierb_effort_cv_log.txt
# and ../figures/tierb2026/tierb2026_effort_cv.pdf|png
# ======================================================================
LOG_NAME <- "tierb_effort_cv_log.txt"
source("./tierb/tierb_common.R")
suppressPackageStartupMessages(library(parallel))

K_OUT   <- as.integer(Sys.getenv("EFFORT_CV_FOLDS", "10"))
N_REPS  <- as.integer(Sys.getenv("EFFORT_CV_REPS", "3"))
N_CORES <- as.integer(Sys.getenv("EFFORT_CV_CORES", "1"))
VARIANT <- Sys.getenv("EFFORT_CV_VARIANT", "full")
B_POOL  <- as.integer(Sys.getenv("EFFORT_CV_BOOT", "500"))
SURVEY_STEP <- 50                                     # survey frames per annotated 2026 frame
TARGET_RATES <- c(1667, 833, 333, 167, 83)            # 1 annotated frame per N survey frames (Fig. S5)

CV_DIR <- file.path(OUT_DIR, "effort_cv"); CACHE <- file.path(CV_DIR, "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
regions <- c("GB", "MAB")
S1_SIMPLE <- truedetect ~ s(conf, k = 5)
VARIANTS <- list(
  full   = list(s1 = S1_FORMS, s1_id = S1_IDS, s2 = S2_TRANSFER),
  simple = list(s1 = list(GB = S1_SIMPLE, MAB = S1_SIMPLE), s1_id = c(GB = "M01_ConfOnly", MAB = "M01_ConfOnly"),
                s2 = c(GB = "M04_YoloLog", MAB = "M04_YoloLog")))
stopifnot(VARIANT %in% names(VARIANTS)); V <- VARIANTS[[VARIANT]]
fit_tw <- function(form, d) tryCatch(gam(form, family = tw(link = "log"), data = d, method = "REML"),
                                     error = function(e) NULL)

# ----------------------------------------------------------------------
# 1. All annotated 2026 primary images; outer folds of whole transect chunks,
#    balanced on scallops (chunks with scallops) and images (empty chunks)
# ----------------------------------------------------------------------
pm_all <- meta %>% filter(platform == "primary", has_meta, region %in% regions)
assign_outer <- function(d, K, seed) {
  set.seed(seed)
  b <- d %>% group_by(region, transect_block_id) %>%
    summarise(n_ann = sum(n_annotations), n_img = n(), .groups = "drop") %>%
    mutate(u = runif(n())) %>% arrange(region, desc(n_ann), u)
  b$outer_fold <- NA_integer_
  for (r in unique(b$region)) {
    load <- numeric(K); cnt <- numeric(K)
    for (i in which(b$region == r)) {
      score <- if (b$n_ann[i] > 0) load else cnt
      cand <- which(score == min(score)); k <- if (length(cand) == 1) cand else sample(cand, 1)
      b$outer_fold[i] <- k; load[k] <- load[k] + b$n_ann[i]; cnt[k] <- cnt[k] + b$n_img[i]
    }
  }
  d %>% left_join(b %>% select(region, transect_block_id, outer_fold), by = c("region", "transect_block_id"))
}
pm_all <- assign_outer(pm_all, K_OUT, SEED)
fold_tab <- pm_all %>% group_by(region, outer_fold) %>%
  summarise(images = n(), scallop_images = sum(n_annotations > 0), scallops = sum(n_annotations), .groups = "drop")
say("Outer folds (images / scallop images / scallops):"); print(fold_tab, n = Inf)
train_share <- 1 - 1 / K_OUT
FRACTIONS <- c(pmin(1, SURVEY_STEP / (train_share * TARGET_RATES)), 1)
say(sprintf("Training share %.2f; thinning fractions %s", train_share, paste(round(FRACTIONS, 3), collapse = ", ")))

# Systematic thinning of the training folds within each transect (random start); as tierb_effort.R
thin_fit <- function(pm, f, seed) {
  if (f >= 1) return(pm)
  set.seed(seed)
  pm %>% group_by(transect_group) %>% arrange(transect_position, .by_group = TRUE) %>%
    mutate(fit_rank = cumsum(dataset == "test_GAM_train") - 1, start_u = runif(1),
           keep_row = dataset != "test_GAM_train" |
             (floor((fit_rank + 1) * f + start_u) - floor(fit_rank * f + start_u)) == 1) %>%
    ungroup() %>% filter(keep_row) %>% select(-fit_rank, -start_u, -keep_row)
}

# ----------------------------------------------------------------------
# 2. One fit: fold k held out, training folds thinned to fraction f.
#    Same two-stage procedure as run_one() in tierb_effort.R, but returns the
#    held-out images' predictions instead of metrics.
# ----------------------------------------------------------------------
fit_one <- function(k, f, rep) {
  seed <- SEED + round(1000 * f) + rep + 100 * k
  d <- pm_all %>% mutate(dataset = if_else(outer_fold == k, "test_GAM_test", "test_GAM_train"))
  d <- suppressMessages(assign_folds(d, K = N_FOLDS, seed = seed, block_col = "transect_block_id"))
  d <- thin_fit(d, f, seed)
  tb <- build_tables(d); det <- add_re_factors(tb$det)
  f1_list <- list(GB = NA_real_, MAB = NA_real_)
  for (r in regions)
    f1_list[[r]] <- f1_threshold(det %>% filter(region == r, dataset == "test_GAM_train"),
                                 tb$img %>% filter(region == r, dataset == "test_GAM_train"))$conf
  det$pred_p_oof <- NA_real_; s1_gams <- list(); used1 <- character(0)
  for (r in regions) {
    form <- V$s1[[r]]; id <- V$s1_id[[r]]
    tr <- det_fit_rows(det %>% filter(region == r, dataset == "test_GAM_train"))
    fit <- try_fit_det_gam(form, tr, engine = DET_ENGINE, label = paste(r, id))
    if (!fit_ok(fit)) { form <- S1_SIMPLE; id <- "M01_ConfOnly (fallback)"
      fit <- try_fit_det_gam(form, tr, engine = DET_ENGINE, label = paste(r, id)) }
    s1_gams[[r]] <- fit; used1[r] <- id
    for (j in seq_len(N_FOLDS)) {
      in_set <- det$region == r & det$dataset == "test_GAM_train" & !is.na(det$fold)
      fi <- which(in_set & det$fold != j); ho <- which(in_set & det$fold == j)
      if (!length(ho) || !length(fi)) next
      fj <- try_fit_det_gam(form, det_fit_rows(det[fi, ]), engine = DET_ENGINE, label = sprintf("%s fold %d", r, j))
      if (!fit_ok(fj)) fj <- try_fit_det_gam(S1_SIMPLE, det_fit_rows(det[fi, ]), engine = DET_ENGINE)
      if (!is.null(fj)) det$pred_p_oof[ho] <- predict_det_p(fj, det[ho, ])
    }
  }
  s1 <- suppressWarnings(test_calibration_gams(s1_gams, det, tb$img, f1_list, test_pred_source = "full"))
  img2 <- s1$img_eval %>% mutate(pred_yolo = predicted_number, log_yolo_pred = log1p(pred_yolo))
  img2 <- add_crowding_features(img2, image_detection_features(s1$det_eval))
  img2$pred_two_stage <- NA_real_; used2 <- character(0)
  for (r in regions) {
    tr <- img2$region == r & img2$dataset == "test_GAM_train"
    te <- img2$region == r & img2$dataset == "test_GAM_test"
    id <- V$s2[[r]]; g <- fit_tw(S2_CANDIDATES[[id]], img2[tr, ])
    if (is.null(g)) { id <- "M04_YoloLog (fallback)"; g <- fit_tw(S2_CANDIDATES$M04_YoloLog, img2[tr, ]) }
    used2[r] <- id
    img2$pred_two_stage[te] <- predict(g, img2[te, ], type = "response")
  }
  fit_info <- d %>% filter(dataset == "test_GAM_train") %>% group_by(region) %>%
    summarise(fit_images = n(), fit_scallop_images = sum(n_annotations > 0), fit_scallops = sum(n_annotations),
              .groups = "drop") %>%
    mutate(s1_form = unname(used1[region]), s2_form = unname(used2[region]),
           f1_threshold = unlist(f1_list[region]))
  list(pred = img2 %>% filter(dataset == "test_GAM_test") %>%
         transmute(image_id, region, transect_block_id, n_annotations,
                   f1_count = predicted_f1_number, stage1 = predicted_number, two_stage = pred_two_stage),
       info = fit_info)
}

# ----------------------------------------------------------------------
# 3. Run all fits (cached; resumable; optional parallel)
# ----------------------------------------------------------------------
grid <- tidyr::expand_grid(fraction = FRACTIONS, rep = seq_len(N_REPS), k = seq_len(K_OUT)) %>%
  filter(fraction < 1 | rep == 1) %>%
  mutate(id = sprintf("%s_f%.4f_r%d_k%02d", VARIANT, fraction, rep, k))
say("Fits: ", nrow(grid), " (", K_OUT, " folds; ", N_REPS, " draws per rate below 100%; variant ", VARIANT,
    "; ", N_CORES, " core(s))")
t0 <- Sys.time()
run_task <- function(i) {
  g <- grid[i, ]; f_cache <- file.path(CACHE, paste0(g$id, ".rds"))
  if (file.exists(f_cache)) return(f_cache)
  res <- tryCatch(fit_one(g$k, g$fraction, g$rep),
                  error = function(e) { say("  !! ", g$id, " failed: ", conditionMessage(e)); NULL })
  if (!is.null(res)) saveRDS(res, f_cache)
  say(sprintf("  [%d/%d] %s %s (%.1f min elapsed)", i, nrow(grid), g$id, if (is.null(res)) "FAILED" else "done",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (is.null(res)) NA_character_ else f_cache
}
done <- if (N_CORES > 1) unlist(mclapply(seq_len(nrow(grid)), run_task, mc.cores = N_CORES, mc.preschedule = FALSE)) else
  vapply(seq_len(nrow(grid)), run_task, character(1))
grid$ok <- !is.na(done) & file.exists(file.path(CACHE, paste0(grid$id, ".rds")))
if (any(!grid$ok)) say("!! ", sum(!grid$ok), " fits missing; rerun the script to retry them (finished fits are kept).")

res <- lapply(seq_len(nrow(grid)), function(i) if (grid$ok[i]) readRDS(file.path(CACHE, paste0(grid$id[i], ".rds"))) else NULL)
preds <- bind_rows(lapply(seq_len(nrow(grid)), function(i) if (!is.null(res[[i]]))
  res[[i]]$pred %>% mutate(fraction = grid$fraction[i], rep = grid$rep[i], outer_fold = grid$k[i])))
fits <- bind_rows(lapply(seq_len(nrow(grid)), function(i) if (!is.null(res[[i]]))
  res[[i]]$info %>% mutate(fraction = grid$fraction[i], rep = grid$rep[i], outer_fold = grid$k[i])))
rate_of <- function(f) SURVEY_STEP / (train_share * f)
preds <- preds %>% mutate(survey_frames_per_annotated = rate_of(fraction))
write.csv(preds, file.path(CV_DIR, "effort_cv_predictions.csv"), row.names = FALSE)
write.csv(fits %>% mutate(survey_frames_per_annotated = rate_of(fraction)),
          file.path(CV_DIR, "effort_cv_fits.csv"), row.names = FALSE)

# ----------------------------------------------------------------------
# 4. Metrics: per fold, and pooled over folds (only draws with every fold present)
# ----------------------------------------------------------------------
M <- c("YOLO F1-threshold count" = "f1_count", "Stage 1 (detection GAM sum)" = "stage1", "Two-stage" = "two_stage")
fold_metrics <- preds %>% group_by(fraction, rep, outer_fold) %>%
  group_modify(function(d, key) method_metrics(d, M)) %>% ungroup() %>%
  mutate(survey_frames_per_annotated = rate_of(fraction))
write.csv(fold_metrics, file.path(CV_DIR, "effort_cv_fold_metrics.csv"), row.names = FALSE)

complete <- preds %>% distinct(fraction, rep, outer_fold) %>% count(fraction, rep) %>% filter(n == K_OUT)
pooled <- preds %>% semi_join(complete, by = c("fraction", "rep")) %>%
  group_by(fraction, rep) %>%
  group_modify(function(d, key) {
    pt <- method_metrics(d, M)
    bb <- block_bootstrap_metrics(d, M, block = "transect_block_id", ref = "Two-stage", B = B_POOL,
                                  seed = SEED + round(1000 * key$fraction) + key$rep) %>%
      select(region, method, rmse_lo, rmse_hi, bias_lo, bias_hi)
    pt %>% left_join(bb, by = c("region", "method"))
  }) %>% ungroup() %>%
  mutate(survey_frames_per_annotated = rate_of(fraction))
write.csv(pooled, file.path(CV_DIR, "effort_cv_pooled.csv"), row.names = FALSE)

fold_sd <- fold_metrics %>% group_by(fraction, region, method) %>%
  summarise(fold_sd_bias = sd(pct_bias_total), fold_sd_rmse = sd(rmse), .groups = "drop")
summ <- pooled %>% group_by(fraction, survey_frames_per_annotated, region, method) %>%
  summarise(n_draws = n(),
            across(c(rmse, pct_bias_total, r2_1to1), list(med = median, min = min, max = max)),
            bias_lo_med = median(bias_lo), bias_hi_med = median(bias_hi), .groups = "drop") %>%
  left_join(fold_sd, by = c("fraction", "region", "method")) %>%
  left_join(fits %>% group_by(fraction, region) %>%
              summarise(fit_scallop_images = median(fit_scallop_images),
                        n_fallback_s1 = sum(grepl("fallback", s1_form)),
                        n_fallback_s2 = sum(grepl("fallback", s2_form)), n_fits = n(), .groups = "drop"),
            by = c("fraction", "region")) %>%
  arrange(region, method, desc(survey_frames_per_annotated))
write.csv(summ, file.path(CV_DIR, "effort_cv_summary.csv"), row.names = FALSE)
say("\nPooled out-of-fold accuracy over all 2026 images (median over draws; range; fold SD):")
print(summ %>% filter(method %in% c("Two-stage", "YOLO F1-threshold count")) %>%
        transmute(region, method, rate = sprintf("1 in %s", format(round(survey_frames_per_annotated), big.mark = ",")),
                  fit_scallop_images, rmse = round(rmse_med, 3),
                  bias = round(pct_bias_total_med, 1),
                  bias_range = sprintf("%.1f to %.1f", pct_bias_total_min, pct_bias_total_max),
                  bias_CI = sprintf("%.1f to %.1f", bias_lo_med, bias_hi_med),
                  fold_sd_bias = round(fold_sd_bias, 1), fallbacks = n_fallback_s2), n = Inf, width = Inf)

# ----------------------------------------------------------------------
# 5. Figure: pooled RMSE and total bias by annotation rate
#    points = draws (95% block-bootstrap CI), line = median over draws
# ----------------------------------------------------------------------
lab_rate <- function(x) sprintf("1 in %s", format(round(x), big.mark = ","))
pd <- pooled %>% filter(method %in% c("Two-stage", "YOLO F1-threshold count")) %>%
  mutate(rate = factor(lab_rate(survey_frames_per_annotated),
                       levels = lab_rate(sort(unique(survey_frames_per_annotated), decreasing = TRUE))),
         region = recode(region, GB = "Georges Bank", MAB = "Mid-Atlantic Bight"),
         method = recode(method, "Two-stage" = "Two-stage calibration", "YOLO F1-threshold count" = "YOLOv12 F1-threshold count")) %>%
  select(region, method, rate, rep, rmse, rmse_lo, rmse_hi, pct_bias_total, bias_lo, bias_hi) %>%
  pivot_longer(c(rmse, pct_bias_total), names_to = "metric", values_to = "value") %>%
  mutate(lo = if_else(metric == "rmse", rmse_lo, bias_lo), hi = if_else(metric == "rmse", rmse_hi, bias_hi),
         metric = recode(metric, rmse = "RMSE (scallops per image)", pct_bias_total = "Total bias (%)"))
pdm <- pd %>% group_by(region, method, rate, metric) %>% summarise(value = median(value), .groups = "drop")
dodge <- position_dodge(width = 0.5)
p <- ggplot(pd, aes(rate, value, colour = method)) +
  geom_hline(data = tibble(metric = "Total bias (%)", y = 0), aes(yintercept = y), linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi, group = interaction(method, rep)), position = dodge,
                  size = 0.2, alpha = 0.6) +
  geom_line(data = pdm, aes(group = method), linewidth = 0.8) +
  facet_grid(metric ~ region, scales = "free_y", switch = "y") +
  scale_colour_manual(values = c("Two-stage calibration" = "#2a78d6", "YOLOv12 F1-threshold count" = "#1baf7a")) +
  labs(x = "Annotation rate used to fit the GAMs (annotated frames per survey frames)", y = NULL, colour = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", strip.placement = "outside", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(FIG_DIR, "tierb2026_effort_cv.pdf"), p, width = 7, height = 5)
ggsave(file.path(FIG_DIR, "tierb2026_effort_cv.png"), p, width = 7, height = 5, dpi = 300, bg = "white")
say("\nDone. Outputs in ", CV_DIR, " | figure ", file.path(FIG_DIR, "tierb2026_effort_cv.pdf"))
