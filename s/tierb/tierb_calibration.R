# ======================================================================
# tierb_calibration.R - Tier B: independent 2026 survey test
#
# Detector : YOLOv12n trained on 2022-2024 only (job 598580); run on every
#            annotated 2026 image (eval_tierb_2026.sh on ARC).
# Split    : split_tierb2026.py -> Calibration-fitting (test_GAM_train) /
#            Holdout test (test_GAM_test), 60/40 within 20-image transect chunks.
# Question : does the two-stage calibration, with the model FORMS chosen on
#            2022-2026 (no re-selection), correct a detector applied to a year
#            it never saw, when the GAMs are refit on a small 2026 sample?
#
# Stage 1  : detection-level binomial GAM, formula fixed from the main YOLO
#            CV selection (GB M15_All_te, MAB M01_ConfOnly), no REs (one cruise).
#            Cross-fitted on Calibration-fitting (block folds) so stage-2
#            training inputs are out-of-fold; full fit predicts the holdout.
# Stage 2  : image-level Tweedie GAM (log link). Only one detector exists here,
#            so the form is the best YOLO-only formula in the main CV
#            (GB C07_teLowConf, MAB C08_Structure). Sensitivity: re-select
#            among YOLO-only forms by block CV on 2026 Calibration-fitting.
# Platforms: primary HabCam = main result (region-specific, full covariates).
#            Polaris AUV = separate, pooled, confidence-only stage 1 (no
#            vehicle metadata needed). Stella = descriptive only (26 scallops).
# Empties  : all kept (2026 is ~90% empty). Sensitivity: holdout metrics on
#            repeated subsamples with the main-analysis empty share (12.4%).
#
# Run from s/:  Rscript tierb/tierb_calibration.R
# Inputs (copy from ARC: output/2224scallop_yolo12n_slurm_598580/eval_tierb/tierb2026/):
#   ../data/raw_model/2224scallop_yolo12n_slurm_598580/tierb2026/autotest_tierb2026_yolo12n.csv
#   ../data/raw/dataset_split_tierb2026.csv   (split_tierb2026.py)
#   ../data/raw/allimg_2026.RData             (vehicle metadata, primary HabCam)
# Outputs: ../data/processed/tierb2026/*.csv|rds, ../figures/tierb2026/*.png
# ======================================================================
source("./tierb/tierb_common.R")   # config, data, helpers (sections 0-3)

# ======================================================================
# 4. PRIMARY HabCam: region-specific two-stage calibration
# ======================================================================
say("\n===== Primary HabCam =====")
pm <- meta %>% filter(platform == "primary", has_meta)
regions <- pm %>% filter(dataset == "test_GAM_train", n_annotations > 0) %>% count(region) %>%
  filter(n >= MIN_POS_IMAGES_FIT) %>% pull(region)
say("Regions with >= ", MIN_POS_IMAGES_FIT, " scallop images in Calibration-fitting: ", paste(regions, collapse = ", "))
pm <- pm %>% filter(region %in% regions)
print(pm %>% group_by(region, partition) %>%
        summarise(images = n(), scallop_images = sum(n_annotations > 0), scallops = sum(n_annotations), .groups = "drop"))

# --- Spatial spacing (independence argument) ---
sp <- pm %>% arrange(transect_group, transect_position) %>% group_by(transect_group) %>%
  mutate(d_prev_m = haversine_m(lag(latitude), lag(longitude), latitude, longitude)) %>% ungroup()
nn_fit <- sp %>% group_by(transect_group) %>%
  mutate(fit_lat = if_else(dataset == "test_GAM_train", latitude, NA_real_),
         fit_lon = if_else(dataset == "test_GAM_train", longitude, NA_real_),
         prev_lat = fit_lat, prev_lon = fit_lon, next_lat = fit_lat, next_lon = fit_lon) %>%
  fill(prev_lat, prev_lon, .direction = "down") %>% fill(next_lat, next_lon, .direction = "up") %>%
  ungroup() %>% filter(dataset == "test_GAM_test") %>%
  mutate(d_nn_fit_m = pmin(haversine_m(latitude, longitude, prev_lat, prev_lon),
                           haversine_m(latitude, longitude, next_lat, next_lon), na.rm = TRUE))
spacing <- tibble(
  metric = c("consecutive annotated frames: median distance (m)",
             "consecutive annotated frames: median time gap (s)",
             "holdout image -> nearest Calibration-fitting image along transect: median distance (m)",
             "holdout image -> nearest Calibration-fitting image: 10th percentile (m)"),
  value = c(median(sp$d_prev_m, na.rm = TRUE), median(sp$sec_since_prev_annotated, na.rm = TRUE),
            median(nn_fit$d_nn_fit_m, na.rm = TRUE), quantile(nn_fit$d_nn_fit_m, .1, na.rm = TRUE)))
print(spacing); write.csv(spacing, file.path(OUT_DIR, "spacing_primary.csv"), row.names = FALSE)

# --- Folds (whole chunks, within region, high-density chunks balanced) ---
pm <- assign_folds(pm, K = N_FOLDS, seed = SEED, block_col = "transect_block_id")
tb <- build_tables(pm)
tb$det <- add_re_factors(tb$det)

# --- F1 thresholds on Calibration-fitting ---
f1_tab <- map_dfr(regions, function(r)
  f1_threshold(tb$det %>% filter(region == r, dataset == "test_GAM_train"),
               tb$img %>% filter(region == r, dataset == "test_GAM_train")) %>% mutate(region = r, .before = 1))
print(f1_tab); write.csv(f1_tab, file.path(OUT_DIR, "f1_thresholds_primary.csv"), row.names = FALSE)
f1_list <- list(GB = NA_real_, MAB = NA_real_)
for (r in regions) f1_list[[r]] <- f1_tab$conf[f1_tab$region == r]

# --- Stage 1: final fits + cross-fit ---
train_det <- det_fit_rows(tb$det %>% filter(dataset == "test_GAM_train"))
say("Stage-1 fitting rows: ", nrow(train_det), " (of ", sum(tb$det$dataset == "test_GAM_train"), " Calibration-fitting detections)")
s1_gams <- list()
for (r in regions) {
  fit <- try_fit_det_gam(S1_FORMS[[r]], train_det %>% filter(region == r), engine = DET_ENGINE,
                         label = paste("stage-1 final", r))
  if (!fit_ok(fit)) {
    say("!! Stage-1 ", S1_IDS[[r]], " unusable in ", r, "; falling back to M01_ConfOnly")
    S1_FORMS[[r]] <- truedetect ~ s(conf, k = 5); S1_IDS[[r]] <- "M01_ConfOnly (fallback)"
    fit <- try_fit_det_gam(S1_FORMS[[r]], train_det %>% filter(region == r), engine = DET_ENGINE)
  }
  s1_gams[[r]] <- fit
  say(sprintf("  %s stage 1 = %s | edf %.1f | dev. expl. %.3f", r, S1_IDS[[r]], sum(fit$edf), dev_expl(fit)))
}
cf <- crossfit_detection_gams(tb$det, fixed_formulas = S1_FORMS, K = N_FOLDS,
                              nested_selection = FALSE, engine = DET_ENGINE)
write.csv(cf$selection, file.path(OUT_DIR, "stage1_crossfit_folds_primary.csv"), row.names = FALSE)
s1 <- test_calibration_gams(s1_gams, cf$det, tb$img, f1_list, test_pred_source = "full")
say("Stage-1 holdout metrics:"); print(s1$metrics, width = Inf)
write.csv(s1$leakage_check, file.path(OUT_DIR, "stage1_leakage_check_primary.csv"), row.names = FALSE)

# --- Stage 2 ---
img2 <- s1$img_eval %>%
  mutate(pred_yolo = predicted_number, log_yolo_pred = log1p(pred_yolo))
img2 <- add_crowding_features(img2, image_detection_features(s1$det_eval))

fit_tw <- function(form, d) tryCatch(gam(form, family = tw(link = "log"), data = d, method = "REML"),
                                     error = function(e) { say("  !! stage-2 fit failed: ", conditionMessage(e)); NULL })
img2$pred_two_stage <- NA_real_; img2$pred_two_stage_resel <- NA_real_
s2_gams <- list(); s2_resel <- list(); s2_cv <- list()
for (r in regions) {
  tr <- img2$region == r & img2$dataset == "test_GAM_train"
  te <- img2$region == r & img2$dataset == "test_GAM_test"
  s2_gams[[r]] <- fit_tw(S2_CANDIDATES[[S2_TRANSFER[[r]]]], img2[tr, ])
  if (!is.null(s2_gams[[r]])) img2$pred_two_stage[te] <- predict(s2_gams[[r]], img2[te, ], type = "response")
  s2_cv[[r]] <- cv_select_synergy(img2, S2_CANDIDATES, r, metric = "cv_rmse") %>% mutate(region = r)
  best <- s2_cv[[r]]$Model_ID[1]
  s2_resel[[r]] <- fit_tw(S2_CANDIDATES[[best]], img2[tr, ])
  if (!is.null(s2_resel[[r]])) img2$pred_two_stage_resel[te] <- predict(s2_resel[[r]], img2[te, ], type = "response")
  say(sprintf("  %s stage 2: transferred = %s | re-selected on 2026 = %s", r, S2_TRANSFER[[r]], best))
}
write.csv(bind_rows(s2_cv), file.path(OUT_DIR, "stage2_cv_reselection_primary.csv"), row.names = FALSE)

# --- Holdout evaluation ---
METHODS <- c("Raw YOLO count (conf >= 0.1)" = "raw_detection_number",
             "YOLO F1-threshold count"      = "predicted_f1_number",
             "Stage 1 (detection GAM sum)"  = "predicted_number",
             "Two-stage (transferred forms)" = "pred_two_stage",
             "Two-stage (re-selected forms)" = "pred_two_stage_resel")
hold <- img2 %>% filter(dataset == "test_GAM_test")
m_all <- method_metrics(hold, METHODS)
say("\nHoldout test metrics (primary, all images kept):"); print(m_all, n = Inf, width = Inf)
write.csv(m_all, file.path(OUT_DIR, "holdout_metrics_primary.csv"), row.names = FALSE)

boot <- block_bootstrap_metrics(hold, METHODS, block = "transect_block_id",
                                ref = "Two-stage (transferred forms)", B = N_BOOT, seed = SEED)
write.csv(boot, file.path(OUT_DIR, "holdout_bootstrap_primary.csv"), row.names = FALSE)

hold <- hold %>% mutate(density_class = cut(n_annotations, c(-Inf, 0, 2, 10, Inf),
                                            labels = c("0", "1-2", "3-10", ">10")))
m_dens <- method_metrics(hold, METHODS, by = c("region", "density_class"))
write.csv(m_dens, file.path(OUT_DIR, "holdout_metrics_by_density_primary.csv"), row.names = FALSE)

m_pos <- method_metrics(hold %>% filter(n_annotations > 0), METHODS)
say("\nHoldout metrics on scallop images only (n_annotations > 0):"); print(m_pos, n = Inf, width = Inf)
write.csv(m_pos, file.path(OUT_DIR, "holdout_metrics_scallop_images_primary.csv"), row.names = FALSE)

m_match <- matched_empty_metrics(hold, METHODS)
say("\nMatched-empties sensitivity (empty share = ", MAIN_EMPTY_SHARE, ", ", N_MATCH_DRAWS, " draws):")
print(m_match %>% select(region, method, rmse_med, rmse_lo, rmse_hi, pct_bias_total_med), n = Inf)
write.csv(m_match, file.path(OUT_DIR, "holdout_matched_empties_primary.csv"), row.names = FALSE)

dens <- density_table(hold, METHODS)
say("\nHoldout density per m^2 (primary):"); print(dens, n = Inf)
write.csv(dens, file.path(OUT_DIR, "holdout_density_m2_primary.csv"), row.names = FALSE)

saveRDS(list(img = img2, s1_gams = s1_gams, s2_gams = s2_gams, s2_resel = s2_resel,
             s1_forms = S1_IDS, s2_forms = S2_TRANSFER, f1 = f1_tab),
        file.path(OUT_DIR, "tierb_primary_eval.rds"))

p11 <- hold %>%
  select(region, n_annotations, all_of(unname(METHODS))) %>%
  pivot_longer(-c(region, n_annotations), names_to = "col", values_to = "pred") %>%
  mutate(method = factor(names(METHODS)[match(col, METHODS)], levels = names(METHODS))) %>%
  ggplot(aes(n_annotations, pred)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_grid(region ~ method, labeller = label_wrap_gen(18)) +
  labs(x = "Manual count per image", y = "Predicted count per image",
       title = "Tier B: 2026 Holdout test (detector trained on 2022-2024)") +
  theme_bw(base_size = 9)
ggsave(file.path(FIG_DIR, "tierb_holdout_1to1_primary.png"), p11, width = 11, height = 4.5, dpi = 300, bg = "white")

# ======================================================================
# 5. POLARIS AUV: pooled, confidence-only stage 1 (+ log-count stage 2)
# ======================================================================
say("\n===== Polaris AUV =====")
pol <- meta %>% filter(platform == "polaris")
n_pos_fit <- sum(pol$dataset == "test_GAM_train" & pol$n_annotations > 0)
if (n_pos_fit >= MIN_POS_IMAGES_FIT) {
  set.seed(SEED)
  blk <- unique(pol$transect_block_id[pol$dataset == "test_GAM_train"])
  fold_map <- tibble(transect_block_id = blk, fold = sample(rep_len(seq_len(N_FOLDS), length(blk))))
  pol <- pol %>% select(-any_of("fold")) %>% left_join(fold_map, by = "transect_block_id") %>%
    mutate(fold = if_else(dataset == "test_GAM_train", fold, NA_integer_))
  pt <- build_tables(pol)
  thr <- f1_threshold(pt$det %>% filter(dataset == "test_GAM_train"), pt$img %>% filter(dataset == "test_GAM_train"))
  f_s1 <- truedetect ~ s(conf, k = 5)
  pt$det$pred_p <- NA_real_
  trd <- pt$det$dataset == "test_GAM_train"
  for (k in seq_len(N_FOLDS)) {
    fi <- trd & pt$det$fold != k; ho <- trd & pt$det$fold == k
    if (!any(ho)) next
    g <- gam(f_s1, family = binomial(), data = pt$det[fi, ], method = "REML")
    pt$det$pred_p[ho] <- predict(g, pt$det[ho, ], type = "response")
  }
  g_full <- gam(f_s1, family = binomial(), data = pt$det[trd, ], method = "REML")
  pt$det$pred_p[!trd] <- predict(g_full, pt$det[!trd, ], type = "response")
  pimg <- pt$img %>%
    left_join(pt$det %>% group_by(image_id) %>%
                summarise(raw_detection_number = n(), predicted_f1_number = sum(conf >= thr$conf),
                          predicted_number = sum(pred_p), .groups = "drop"), by = "image_id") %>%
    mutate(across(c(raw_detection_number, predicted_f1_number, predicted_number), ~ replace_na(.x, 0)),
           pred_yolo = predicted_number, log_yolo_pred = log1p(pred_yolo), region = "Polaris")
  g2 <- fit_tw(n_annotations ~ s(log_yolo_pred), pimg %>% filter(dataset == "test_GAM_train"))
  pimg$pred_two_stage <- NA_real_
  te <- pimg$dataset == "test_GAM_test"
  if (!is.null(g2)) pimg$pred_two_stage[te] <- predict(g2, pimg[te, ], type = "response")
  PM <- METHODS[1:4]; names(PM)[4] <- "Two-stage (conf-only + M04_YoloLog)"
  pm_all <- method_metrics(pimg %>% filter(dataset == "test_GAM_test"), PM)
  print(pm_all, width = Inf)
  write.csv(pm_all %>% mutate(f1_threshold = thr$conf), file.path(OUT_DIR, "holdout_metrics_polaris.csv"), row.names = FALSE)
  write.csv(matched_empty_metrics(pimg %>% filter(dataset == "test_GAM_test"), PM), file.path(OUT_DIR, "holdout_matched_empties_polaris.csv"), row.names = FALSE)
  saveRDS(pimg, file.path(OUT_DIR, "tierb_polaris_eval.rds"))
} else say("Polaris: only ", n_pos_fit, " scallop images in Calibration-fitting; skipped.")

# ======================================================================
# 6. STELLA AUV: descriptive only
# ======================================================================
st <- meta %>% filter(platform == "stella")
st_det <- det_all %>% filter(image_id %in% st$image_id)
stella <- tibble(images = nrow(st), scallop_images = sum(st$n_annotations > 0),
                 scallops_annotated = sum(st$n_annotations), yolo_detections = nrow(st_det),
                 true_positives = sum(st_det$truedetect),
                 images_with_any_detection = n_distinct(st_det$image_id))
say("\n===== Stella AUV (descriptive) ====="); print(stella)
write.csv(stella, file.path(OUT_DIR, "stella_descriptive.csv"), row.names = FALSE)

say("\nDone. Outputs in ", OUT_DIR, " and ", FIG_DIR)
