# ======================================================================
# tierb_effort.R - Tier B annotation-effort curve
#
# How many annotated 2026 images does the two-stage calibration need?
# The 2026 survey was annotated at 1 stereo pair in 50; Calibration-fitting
# holds 60% of those (about 1 in 83 survey frames). Calibration-fitting is
# thinned systematically along each transect (random start, same rule as a
# coarser annotation frequency) to 5-100%, both GAM stages are refit on the
# thinned set, and every refit is scored on the SAME, untouched Holdout test.
#
# Variants
#   full   : model forms transferred from the main analysis
#            (stage 1 GB M15_All_te / MAB M01_ConfOnly; stage 2 GB C07 / MAB C08)
#   simple : stage 1 M01_ConfOnly, stage 2 M04_YoloLog in both regions
#            (fewest parameters; the candidate for low annotation effort)
# A form that fails to fit falls back to the simple one; the form actually
# used is recorded in every row.
#
# Run from s/:  Rscript tierb/tierb_effort.R      (about 1-1.5 h on a laptop)
# Outputs: ../data/processed/tierb2026/effort_runs.csv, effort_summary.csv,
#          ../figures/tierb2026/tierb_effort_curve.png
# ======================================================================
LOG_NAME <- "tierb_effort_log.txt"
source("./tierb/tierb_common.R")

FRACTIONS   <- c(0.05, 0.10, 0.25, 0.50, 1.00)
N_REPS      <- as.integer(Sys.getenv("TIERB_EFFORT_REPS", "5"))   # random starts per fraction < 1
SURVEY_STEP <- 50     # survey frames per annotated frame in 2026
FIT_SHARE   <- 0.6    # Calibration-fitting share of annotated frames

VARIANTS <- list(
  full   = list(s1 = S1_FORMS, s1_id = S1_IDS, s2 = S2_TRANSFER),
  simple = list(s1 = list(GB = truedetect ~ s(conf, k = 5), MAB = truedetect ~ s(conf, k = 5)),
                s1_id = c(GB = "M01_ConfOnly", MAB = "M01_ConfOnly"),
                s2 = c(GB = "M04_YoloLog", MAB = "M04_YoloLog"))
)
S1_SIMPLE <- truedetect ~ s(conf, k = 5)

fit_tw <- function(form, d) tryCatch(gam(form, family = tw(link = "log"), data = d, method = "REML"),
                                     error = function(e) NULL)

# ----------------------------------------------------------------------
# Data: primary HabCam, folds assigned once on the FULL Calibration-fitting
# set, so every thinned set keeps the same block folds.
# ----------------------------------------------------------------------
regions <- c("GB", "MAB")
pm <- meta %>% filter(platform == "primary", has_meta, region %in% regions)
pm <- assign_folds(pm, K = N_FOLDS, seed = SEED, block_col = "transect_block_id")

# Systematic thinning of Calibration-fitting within each transect (random start)
thin_fit <- function(pm, f, seed) {
  if (f >= 1) return(pm)
  set.seed(seed)
  pm %>%
    group_by(transect_group) %>%
    arrange(transect_position, .by_group = TRUE) %>%
    # (column names must not be .keep/.by etc.: those are mutate() arguments)
    mutate(fit_rank = cumsum(dataset == "test_GAM_train") - 1, start_u = runif(1),
           keep_row = dataset != "test_GAM_train" |
             (floor((fit_rank + 1) * f + start_u) - floor(fit_rank * f + start_u)) == 1) %>%
    ungroup() %>%
    filter(keep_row) %>%
    select(-fit_rank, -start_u, -keep_row)
}

run_one <- function(pm_sub, v) {
  tb <- build_tables(pm_sub)
  det <- add_re_factors(tb$det)
  f1_list <- list(GB = NA_real_, MAB = NA_real_)
  for (r in regions)
    f1_list[[r]] <- f1_threshold(det %>% filter(region == r, dataset == "test_GAM_train"),
                                 tb$img %>% filter(region == r, dataset == "test_GAM_train"))$conf

  # Stage 1: final fit + own cross-fit loop (per-region fallback)
  det$pred_p_oof <- NA_real_
  s1_gams <- list(); used1 <- character(0)
  for (r in regions) {
    form <- v$s1[[r]]; id <- v$s1_id[[r]]
    tr <- det_fit_rows(det %>% filter(region == r, dataset == "test_GAM_train"))
    fit <- try_fit_det_gam(form, tr, engine = DET_ENGINE, label = paste(r, id))
    if (!fit_ok(fit)) {
      form <- S1_SIMPLE; id <- "M01_ConfOnly (fallback)"
      fit <- try_fit_det_gam(form, tr, engine = DET_ENGINE, label = paste(r, id))
    }
    s1_gams[[r]] <- fit; used1[r] <- id
    for (k in seq_len(N_FOLDS)) {
      in_set <- det$region == r & det$dataset == "test_GAM_train" & !is.na(det$fold)
      fi <- which(in_set & det$fold != k); ho <- which(in_set & det$fold == k)
      if (!length(ho) || !length(fi)) next
      fk <- try_fit_det_gam(form, det_fit_rows(det[fi, ]), engine = DET_ENGINE, label = sprintf("%s fold %d", r, k))
      if (!fit_ok(fk)) fk <- try_fit_det_gam(S1_SIMPLE, det_fit_rows(det[fi, ]), engine = DET_ENGINE)
      if (!is.null(fk)) det$pred_p_oof[ho] <- predict_det_p(fk, det[ho, ])
    }
  }
  s1 <- suppressWarnings(test_calibration_gams(s1_gams, det, tb$img, f1_list, test_pred_source = "full"))

  # Stage 2
  img2 <- s1$img_eval %>% mutate(pred_yolo = predicted_number, log_yolo_pred = log1p(pred_yolo))
  img2 <- add_crowding_features(img2, image_detection_features(s1$det_eval))
  img2$pred_two_stage <- NA_real_; used2 <- character(0)
  for (r in regions) {
    tr <- img2$region == r & img2$dataset == "test_GAM_train"
    te <- img2$region == r & img2$dataset == "test_GAM_test"
    id <- v$s2[[r]]; g <- fit_tw(S2_CANDIDATES[[id]], img2[tr, ])
    if (is.null(g)) { id <- "M04_YoloLog (fallback)"; g <- fit_tw(S2_CANDIDATES$M04_YoloLog, img2[tr, ]) }
    used2[r] <- id
    img2$pred_two_stage[te] <- predict(g, img2[te, ], type = "response")
  }

  hold <- img2 %>% filter(dataset == "test_GAM_test")
  M <- c("YOLO F1-threshold count" = "predicted_f1_number",
         "Stage 1 (detection GAM sum)" = "predicted_number",
         "Two-stage" = "pred_two_stage")
  fit_info <- pm_sub %>% filter(dataset == "test_GAM_train") %>% group_by(region) %>%
    summarise(fit_images = n(), fit_scallop_images = sum(n_annotations > 0),
              fit_scallops = sum(n_annotations), .groups = "drop")
  method_metrics(hold, M) %>%
    left_join(method_metrics(hold %>% filter(n_annotations > 0), M) %>%
                select(region, method, rmse_scallop_images = rmse, pct_bias_scallop_images = pct_bias_total),
              by = c("region", "method")) %>%
    left_join(fit_info, by = "region") %>%
    mutate(s1_form = unname(used1[region]), s2_form = unname(used2[region]))
}

# ----------------------------------------------------------------------
# Runs
# ----------------------------------------------------------------------
grid <- expand_grid(fraction = FRACTIONS, rep = seq_len(N_REPS), variant = names(VARIANTS)) %>%
  filter(fraction < 1 | rep == 1)
say("Effort runs: ", nrow(grid), " (", N_REPS, " random starts per fraction < 1)")
t0 <- Sys.time()
runs <- map_dfr(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  seed <- SEED + round(1000 * g$fraction) + g$rep
  res <- tryCatch(run_one(thin_fit(pm, g$fraction, seed), VARIANTS[[g$variant]]),
                  error = function(e) { say("  !! run failed: ", conditionMessage(e),
                                         " [in ", deparse(conditionCall(e))[1], "]"); NULL })
  say(sprintf("  [%d/%d] %s fraction %.2f rep %d done (%.1f min elapsed)", i, nrow(grid), g$variant,
              g$fraction, g$rep, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (is.null(res)) return(tibble(fraction = g$fraction, rep = g$rep, variant = g$variant, failed = TRUE))
  res %>% mutate(fraction = g$fraction, rep = g$rep, variant = g$variant, seed = seed, failed = FALSE, .before = 1)
}) %>%
  mutate(survey_frames_per_annotated = SURVEY_STEP / (FIT_SHARE * fraction))
write.csv(runs, file.path(OUT_DIR, "effort_runs.csv"), row.names = FALSE)

summ <- runs %>% filter(!failed) %>%
  group_by(variant, fraction, survey_frames_per_annotated, region, method) %>%
  summarise(n_reps = n(),
            fit_images = median(fit_images), fit_scallop_images = median(fit_scallop_images),
            fit_scallops = median(fit_scallops),
            across(c(rmse, r2_1to1, pct_bias_total, calibration_slope, pct_bias_scallop_images),
                   list(med = median, min = min, max = max)),
            n_fallback_s1 = sum(grepl("fallback", s1_form)), n_fallback_s2 = sum(grepl("fallback", s2_form)),
            .groups = "drop") %>%
  arrange(region, method, variant, fraction)
write.csv(summ, file.path(OUT_DIR, "effort_summary.csv"), row.names = FALSE)
say("\nEffort summary (Two-stage, holdout):")
print(summ %>% filter(method == "Two-stage") %>%
        select(region, variant, fraction, survey_frames_per_annotated, fit_scallop_images,
               rmse_med, rmse_min, rmse_max, pct_bias_total_med, pct_bias_total_min, pct_bias_total_max,
               n_fallback_s1, n_fallback_s2), n = Inf, width = Inf)

# ----------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------
pd <- runs %>% filter(!failed, method %in% c("Two-stage", "YOLO F1-threshold count")) %>%
  mutate(series = ifelse(method == "Two-stage", paste0("Two-stage (", variant, " forms)"),
                         "YOLO F1-threshold count")) %>%
  distinct(region, series, fraction, rep, .keep_all = TRUE) %>%
  select(region, series, fraction, fit_scallop_images, rmse, pct_bias_total) %>%
  pivot_longer(c(rmse, pct_bias_total), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, rmse = "Holdout RMSE (scallops per image)",
                         pct_bias_total = "Holdout total bias (%)"))
pm_line <- pd %>% group_by(region, series, fraction, metric) %>%
  summarise(fit_scallop_images = median(fit_scallop_images), value = median(value), .groups = "drop")
p <- ggplot(pd, aes(fit_scallop_images, value, colour = series)) +
  geom_hline(data = tibble(metric = "Holdout total bias (%)", y = 0), aes(yintercept = y),
             linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_line(data = pm_line, linewidth = 0.7) +
  scale_x_log10() +
  facet_grid(metric ~ region, scales = "free_y", switch = "y") +
  labs(x = "Scallop-containing images in Calibration-fitting (log scale)", y = NULL, colour = NULL,
       title = "Tier B: calibration accuracy vs. annotation effort (2026 Holdout test)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", strip.placement = "outside")
ggsave(file.path(FIG_DIR, "tierb_effort_curve.png"), p, width = 8, height = 5.5, dpi = 300, bg = "white")
say("\nDone. ", nrow(runs), " rows; outputs in ", OUT_DIR)
