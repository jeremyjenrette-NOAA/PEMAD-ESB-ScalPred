library(dplyr)
library(tidyr)
library(mgcv)
library(ggplot2)
library(patchwork)
source("./cvfunc.R")   # predict_det_p(), add_re_factors(), r2_1to1(), fit_det_gam()

# ======================================================================
# NEW: Compare Candidate Detection-Level GAM Structures
# ======================================================================
compare_detection_gams <- function(det_df, candidate_formulas, region_focus = "GB") {
  
  # Filter strictly for the training split in the specified region
  train_data <- det_df %>%
    filter(dataset == "test_GAM_train", region == region_focus) %>%
    drop_na(truedetect, conf, bottom_depth, latitude, longitude, altitude, 
            field_of_view_sq_meter, fluorometer_backscatter_ntu, ctd_temperature_celsius,
            boxsize, iu)
  
  cat("\nEvaluating", length(candidate_formulas), "candidate models for", region_focus, "(n =", nrow(train_data), ")...\n")
  
  results <- purrr::map_dfr(names(candidate_formulas), function(mod_name) {
    form <- candidate_formulas[[mod_name]]
    
    cat("  Fitting:", mod_name, "...\n")
    fit <- gam(form, family = binomial(), data = train_data, method = "REML")
    
    # Extract evaluation metrics
    tibble(
      Model_ID = mod_name,
      Formula = paste(deparse(form), collapse = " "),
      AIC = AIC(fit),
      Deviance_Explained = summary(fit)$dev.expl,
      UBRE_Score = fit$gcv.ubre
    )
  })
  
  # Sort by lowest AIC
  results <- results %>% arrange(AIC)
  return(results)
}

# ======================================================================
# UPDATED: Fit Regional Calibration GAMs (Now accepts dynamic formulas)
# ======================================================================
# re_terms: character vector, or list(GB = , MAB = ); see RE_SPECS in cvfunc.R.
# Defaults reproduce the original fixed-effect, mgcv::gam fits.
fit_calibration_gams <- function(det_df, formula_gb, formula_mab,
                                 re_terms = character(0), engine = "gam") {
  
  train_data <- det_fit_rows(add_re_factors(det_df) %>% filter(dataset == "test_GAM_train"))
  
  cat("  Fitting Georges Bank GAM...\n")
  m_gb <- fit_det_gam(with_re(formula_gb, re_for_region(re_terms, "GB")),
                      train_data %>% filter(region == "GB"), engine = engine)
  
  cat("  Fitting Mid-Atlantic Bight GAM...\n")
  m_mab <- fit_det_gam(with_re(formula_mab, re_for_region(re_terms, "MAB")),
                       train_data %>% filter(region == "MAB"), engine = engine)
  
  return(list(GB = m_gb, MAB = m_mab))
}
# ======================================================================
# 2. Test GAMs on Holdout Data (Dynamic F1 Thresholds)
# ======================================================================
# ======================================================================
# UPDATED: Safely Predict (Processes all data, evaluates only holdout)
# ======================================================================
# test_pred_source: "full"     -> test_GAM_test uses the full-data GAMs (deployed model)
#                   "fold_avg" -> mean of the K cross-fit fold models (same models
#                                 that produced the OOF training inputs)
test_calibration_gams <- function(gams, det_df, img_df, f1_thresholds,
                                  test_pred_source = c("full", "fold_avg")) {
  test_pred_source <- match.arg(test_pred_source)
  
  # We NO LONGER filter out the training data here! We want predictions for EVERYTHING.
  det_eval <- add_re_factors(det_df)
  img_eval_base <- img_df 
  
  # 1. Safely assign regional F1 thresholds
  det_eval$regional_f1_thresh <- ifelse(det_eval$region == "GB", f1_thresholds[["GB"]], f1_thresholds[["MAB"]])
  
  # 2. Safely predict by explicitly subsetting data (prevents silent vector recycling)
  det_eval$pred_p <- NA_real_
  
  gb_idx <- which(det_eval$region == "GB")
  if (length(gb_idx) > 0) {
    det_eval$pred_p[gb_idx] <- predict_det_p(gams$GB, det_eval[gb_idx, ])
  }
  
  mab_idx <- which(det_eval$region == "MAB")
  if (length(mab_idx) > 0) {
    det_eval$pred_p[mab_idx] <- predict_det_p(gams$MAB, det_eval[mab_idx, ])
  }
  
  # 2b. Cross-fitting: for test_GAM_train detections, replace the in-sample
  #     prediction with the out-of-fold one (from crossfit_detection_gams()).
  #     test_GAM_test detections keep the full-data GAM prediction.
  #     The in-sample value is retained for the leakage diagnostic.
  det_eval$pred_p_insample <- det_eval$pred_p
  det_eval$pred_p_full     <- det_eval$pred_p
  if (test_pred_source == "fold_avg") {
    if (!"pred_p_foldavg" %in% names(det_eval)) stop("test_pred_source = 'fold_avg' needs crossfit_detection_gams(predict_test = TRUE).")
    te_idx <- which(det_eval$dataset == "test_GAM_test")
    det_eval$pred_p[te_idx] <- det_eval$pred_p_foldavg[te_idx]
  }
  if ("pred_p_oof" %in% names(det_eval)) {
    tr_idx <- which(det_eval$dataset == "test_GAM_train")
    det_eval$pred_p[tr_idx] <- det_eval$pred_p_oof[tr_idx]
    n_miss <- sum(is.na(det_eval$pred_p_oof[tr_idx]) & !is.na(det_eval$pred_p_insample[tr_idx]))
    if (n_miss > 0) warning(n_miss, " test_GAM_train detections have an in-sample but no OOF prediction.")
  } else {
    warning("No pred_p_oof column: test_GAM_train sums are IN-SAMPLE (leaky for the synergy GAM).")
  }
  
  # 3. Aggregate predictions to the image level
  img_predictions <- det_eval %>%
    group_by(image_id) %>%
    summarise(
      raw_detection_number = n(),
      predicted_f1_number  = sum(conf >= regional_f1_thresh, na.rm = TRUE),
      predicted_number     = sum(pred_p, na.rm = TRUE),
      predicted_number_insample = sum(pred_p_insample, na.rm = TRUE),
      true_positive_sum    = sum(truedetect, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 4. Join back to master image list
  img_eval <- img_eval_base %>%
    left_join(img_predictions, by = "image_id") %>%
    mutate(
      raw_detection_number = replace_na(raw_detection_number, 0),
      predicted_f1_number  = replace_na(predicted_f1_number, 0),
      predicted_number     = replace_na(predicted_number, 0),
      predicted_number_insample = replace_na(predicted_number_insample, 0),
      true_positive_sum    = replace_na(true_positive_sum, 0),
      false_negative       = n_annotations - true_positive_sum
    )
  
  # 5. Calculate metrics STRICTLY on the holdout test set to avoid data leakage
  metrics <- img_eval %>%
    filter(dataset == "test_GAM_test") %>%  # <-- Crucial line!
    group_by(region) %>%
    summarise(
      # r2_* (lm-based) measure correlation only: a biased predictor can score
      # high. r2_1to1_* measure agreement with the 1:1 line (1 - SSE/SST) and
      # are the ones shown on the 1:1 plots.
      r2_calibrated   = summary(lm(n_annotations ~ predicted_number))$adj.r.squared,
      r2_1to1_calibrated = r2_1to1(n_annotations, predicted_number),
      rmse_calibrated = sqrt(mean((n_annotations - predicted_number)^2, na.rm = TRUE)),
      mae_calibrated  = mean(abs(n_annotations - predicted_number), na.rm = TRUE),
      pct_bias_total_calibrated = 100 * (sum(predicted_number) - sum(n_annotations)) / sum(n_annotations),
      r2_f1           = summary(lm(n_annotations ~ predicted_f1_number))$adj.r.squared,
      r2_1to1_f1      = r2_1to1(n_annotations, predicted_f1_number),
      rmse_f1         = sqrt(mean((n_annotations - predicted_f1_number)^2, na.rm = TRUE)),
      mae_f1          = mean(abs(n_annotations - predicted_f1_number), na.rm = TRUE),
      pct_bias_total_f1 = 100 * (sum(predicted_f1_number) - sum(n_annotations)) / sum(n_annotations),
      n_images        = n(),
      .groups = "drop"
    )
  
  # 6. Leakage diagnostic on test_GAM_train: in-sample vs out-of-fold sums.
  #    A large gap means the old (in-sample) synergy inputs were optimistic.
  leakage_check <- img_eval %>%
    filter(dataset == "test_GAM_train") %>%
    group_by(region) %>%
    summarise(
      r2_insample   = summary(lm(n_annotations ~ predicted_number_insample))$adj.r.squared,
      rmse_insample = sqrt(mean((n_annotations - predicted_number_insample)^2, na.rm = TRUE)),
      r2_oof        = summary(lm(n_annotations ~ predicted_number))$adj.r.squared,
      rmse_oof      = sqrt(mean((n_annotations - predicted_number)^2, na.rm = TRUE)),
      mean_abs_diff = mean(abs(predicted_number - predicted_number_insample), na.rm = TRUE),
      total_true    = sum(n_annotations),
      total_insample = sum(predicted_number_insample),
      total_oof     = sum(predicted_number),
      n_images      = n(),
      .groups = "drop"
    )
  
  return(list(img_eval = img_eval, det_eval = det_eval, metrics = metrics,
              leakage_check = leakage_check))
}
# ======================================================================
# 3. Generate Diagnostic Plots (Restored Aesthetics)
# ======================================================================
# ======================================================================
# 3. Generate Diagnostic Plots (Restored Aesthetics)
# ======================================================================
generate_gam_plots <- function(eval_res, model_name, dataset_label = "2022-2024, 2026") {
  
  # STRICTLY filter to the holdout test set to prevent visual data leakage!
  img_df <- eval_res$img_eval %>% filter(dataset == "test_GAM_test")
  det_df <- eval_res$det_eval %>% filter(dataset == "test_GAM_test")
  
  # Metrics were already filtered in test_calibration_gams, but we'll pull them normally
  metrics <- eval_res$metrics
  
  # Calculate n images and spell out region names for facet labels
  img_df <- img_df %>%
    mutate(
      region_full = case_when(
        region == "GB" ~ "Georges Bank",
        region == "MAB" ~ "Mid-Atlantic Bight",
        TRUE ~ as.character(region)
      )
    ) %>%
    group_by(region_full) %>%
    mutate(facet_label = paste0(region_full, " (n = ", n(), ")")) %>%
    ungroup()
  
  metrics_labels <- metrics %>%
    mutate(
      region_full = case_when(region == "GB" ~ "Georges Bank", region == "MAB" ~ "Mid-Atlantic Bight", TRUE ~ region),
      label = sprintf("R²(1:1) = %.2f\nRMSE = %.2f", r2_1to1_calibrated, rmse_calibrated)
    ) %>%
    left_join(img_df %>% distinct(region_full, facet_label), by = "region_full")
  
  # A. Zoomed 1-to-1 Calibration Fit Plot
  p_zoomed <- ggplot(img_df, aes(x = predicted_number, y = n_annotations)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(color = "red", size = 1.5, alpha = 0.6) +
    facet_wrap(~ facet_label, scales = "free") +
    geom_text(
      data = metrics_labels, aes(x = -Inf, y = Inf, label = label),
      hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4, fontface = "bold"
    ) +
    theme_bw(base_size = 13) +
    labs(
      title = "True abundance vs. Σ calibrated detection probabilities per image", 
      subtitle = bquote("Dataset: " ~ .(dataset_label) * ", Detection Model: " ~ bold(.(model_name))),
      x = "Σ p(detection)", y = "Manual count"
    ) +
    theme(strip.text = element_text(face = "bold", size = 12))
  
  # B. Residuals
  p_resid <- img_df %>%
    mutate(residual = predicted_number - n_annotations) %>%
    ggplot(aes(x = n_annotations, y = residual)) +
    geom_point(alpha = 0.5, color = "#2C7FB8") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = "Residuals vs True Count", subtitle = model_name, x = "True count", y = "Residual")
  
  # C. False Negatives
  p_fn <- ggplot(img_df, aes(x = n_annotations, y = false_negative)) +
    geom_point(alpha = 0.5, color = "#D95F0E") +
    geom_smooth(method = "gam", color = "black") +
    theme_minimal() +
    labs(title = "False Negatives vs True Abundance", subtitle = model_name, x = "True abundance", y = "False negatives")
  
  # D. Total Abundance Summary
  sum_df <- data.frame(
    metric = factor(c("True count", "Raw detector", "GAM Calibrated", "F1 Cutoff"), 
                    levels = c("True count", "Raw detector", "F1 Cutoff", "GAM Calibrated")),
    value  = c(
      sum(img_df$n_annotations, na.rm = TRUE),
      nrow(det_df),
      sum(det_df$pred_p, na.rm = TRUE),
      sum(img_df$predicted_f1_number, na.rm = TRUE)
    )
  )
  
  p_sum <- ggplot(sum_df, aes(x = metric, y = value, fill = metric)) +
    geom_col(width = 0.6, color = "black") +
    geom_hline(yintercept = sum_df$value[1], linetype = "dashed", color = "black", linewidth = 1) +
    geom_text(aes(label = round(value, 0)), vjust = -0.5, fontface = "bold") +
    scale_fill_manual(values = c("grey40", "salmon", "#D95F0E", "#2C7FB8")) +
    theme_minimal() +
    labs(title = "Total Scallop Abundance Estimate", subtitle = model_name, x = "", y = "Total sum count") +
    theme(legend.position = "none")
  
  return(list(p_zoomed = p_zoomed, p_resid = p_resid, p_fn = p_fn, p_sum = p_sum))
}
# ======================================================================
# 4. Build Synergistic Image-Level Dataset
# ======================================================================
build_synergy_dataset <- function(yolo_eval, cas_eval, meta_df) {
  
  # Find images common to both test sets
  common_ids <- intersect(yolo_eval$img_eval$image_id, cas_eval$img_eval$image_id)
  
  img_combined <- yolo_eval$img_eval %>%
    filter(image_id %in% common_ids) %>%
    select(image_id, region, dataset, any_of("fold"), n_annotations, pred_yolo = predicted_number,
           pred_yolo_insample = any_of("predicted_number_insample")) %>%
    left_join(
      cas_eval$img_eval %>%
        filter(image_id %in% common_ids) %>%
        select(image_id, pred_cascade = predicted_number,
               pred_cascade_insample = any_of("predicted_number_insample")),
      by = "image_id"
    ) %>%
    add_synergy_features()
  
  # Attach environmental metadata
  img_combined <- img_combined %>%
    left_join(
      meta_df %>% select(-any_of(c("region", "dataset", "n_annotations", "imagename", "fold"))), 
      by = "image_id"
    )
  
  return(img_combined)
}

# Derived synergy covariates from pred_yolo / pred_cascade.
# Factored out so the leakage diagnostic can rebuild them from in-sample sums.
add_synergy_features <- function(df) {
  df %>%
    mutate(
      pred_mean = (pred_yolo + pred_cascade) / 2,
      pred_diff = pred_cascade - pred_yolo,
      pred_sum  = pred_yolo + pred_cascade,
      log_yolo_pred = log1p(pred_yolo),
      log_cascade_pred = log1p(pred_cascade),
      log_pred_mean = log1p(pred_mean),
      log_pred_sum = log1p(pred_sum),
      # log_pred_diff = log1p(abs(pred_diff)),
      log_true_number = log1p(n_annotations + 1e-6),
      log_pred_diff = log1p(abs(pred_cascade - pred_yolo)),
      log_ratio = log((pred_yolo + 1) / (pred_cascade + 1))
    )
}

# Swap in the in-sample detection-GAM sums and rebuild features: reproduces the
# OLD (leaky) synergy training inputs for the before/after comparison.
insample_synergy_dataset <- function(img_combined) {
  if (!all(c("pred_yolo_insample", "pred_cascade_insample") %in% names(img_combined))) {
    stop("In-sample columns missing; re-run gamstrat.R with the cross-fitting version of gamfunc.R.")
  }
  img_combined %>%
    mutate(pred_yolo = pred_yolo_insample, pred_cascade = pred_cascade_insample) %>%
    add_synergy_features()
}

# ======================================================================
# 5. Compare Candidate Synergistic Image-Level GAMs (Stratified)
# ======================================================================
# Fit candidate models using Tweedie distribution to match final training
compare_image_gams <- function(img_combined, candidate_formulas, region_focus = "GB") {
  train_data <- img_combined %>% filter(dataset == "test_GAM_train", region == region_focus)
  
  results <- purrr::map_dfr(names(candidate_formulas), function(mod_name) {
    form <- candidate_formulas[[mod_name]]
    fit <- gam(form, family = tw(link = "log"), data = train_data, method = "REML")
    preds <- predict(fit, type = "response")
    
    tibble(
      Model_ID = mod_name,
      Formula  = paste(deparse(form), collapse = " "),
      AIC      = AIC(fit),
      Deviance_Explained = summary(fit)$dev.expl,
      RMSE_Train = sqrt(mean((train_data$n_annotations - preds)^2)),
      R2_Train   = summary(lm(train_data$n_annotations ~ preds))$adj.r.squared
    )
  })
  
  return(results %>% arrange(AIC))
}

# ======================================================================
# 6. Fit and Test Stratified Synergistic GAMs
# ======================================================================
fit_synergy_gams <- function(img_combined, formula_gb, formula_mab) {
  
  train_data <- img_combined %>% filter(dataset == "test_GAM_train")
  
  cat("  Fitting Georges Bank Synergy GAM...\n")
  m_gb <- gam(formula_gb, family = tw(link = "log"), data = train_data %>% filter(region == "GB"), method = "REML")
  
  cat("  Fitting Mid-Atlantic Bight Synergy GAM...\n")
  m_mab <- gam(formula_mab, family = tw(link = "log"),  data = train_data %>% filter(region == "MAB"), method = "REML")
  
  return(list(GB = m_gb, MAB = m_mab))
}

test_synergy_gams <- function(gams, img_combined) {
  
  test_data <- img_combined %>% filter(dataset == "test_GAM_test")
  test_data$pred_synergy <- NA_real_
  
  # Predict safely by region
  gb_idx <- which(test_data$region == "GB")
  if (length(gb_idx) > 0) test_data$pred_synergy[gb_idx] <- predict(gams$GB, newdata = test_data[gb_idx, ], type = "response")
  
  mab_idx <- which(test_data$region == "MAB")
  if (length(mab_idx) > 0) test_data$pred_synergy[mab_idx] <- predict(gams$MAB, newdata = test_data[mab_idx, ], type = "response")
  
  # Calculate holdout metrics
  metrics <- test_data %>%
    group_by(region) %>%
    summarise(
      r2_synergy      = summary(lm(n_annotations ~ pred_synergy))$adj.r.squared,
      r2_1to1_synergy = r2_1to1(n_annotations, pred_synergy),
      rmse_synergy    = sqrt(mean((n_annotations - pred_synergy)^2, na.rm = TRUE)),
      mae_synergy     = mean(abs(n_annotations - pred_synergy), na.rm = TRUE),
      pct_bias_total_synergy = 100 * (sum(pred_synergy) - sum(n_annotations)) / sum(n_annotations),
      .groups = "drop"
    )
  
  return(list(img_eval = test_data, metrics = metrics))
}

# ======================================================================
# 7. Cross-validated synergy formula selection
# ======================================================================
# AIC on the training images favours flexible surfaces (e.g. 3-D tensors) that
# can extrapolate badly for the few high-density images. Using the same block
# folds, each candidate is fit on folds != k and scored on fold k; the
# formula with the lowest CV error is selected. No test_GAM_test data used.
cv_select_synergy <- function(img_combined, candidate_formulas, region_focus,
                              metric = c("cv_rmse", "cv_mae")) {
  metric <- match.arg(metric)
  d <- img_combined %>% filter(dataset == "test_GAM_train", region == region_focus, !is.na(fold))
  if (nrow(d) == 0) stop("No test_GAM_train rows with a fold for ", region_focus, " (fold column missing?).")
  folds <- sort(unique(d$fold))
  
  res <- purrr::map_dfr(names(candidate_formulas), function(id) {
    pred <- rep(NA_real_, nrow(d))
    for (k in folds) {
      tr <- d$fold != k; ho <- d$fold == k
      fit <- tryCatch(gam(candidate_formulas[[id]], family = tw(link = "log"),
                          data = d[tr, ], method = "REML"),
                      error = function(e) NULL)
      if (is.null(fit)) return(tibble(Model_ID = id, cv_rmse = NA_real_, cv_mae = NA_real_,
                                      cv_r2_1to1 = NA_real_, cv_pct_bias_total = NA_real_))
      pred[ho] <- predict(fit, newdata = d[ho, ], type = "response")
    }
    tibble(Model_ID = id,
           cv_rmse = sqrt(mean((d$n_annotations - pred)^2, na.rm = TRUE)),
           cv_mae  = mean(abs(d$n_annotations - pred), na.rm = TRUE),
           cv_r2_1to1 = r2_1to1(d$n_annotations, pred),
           cv_pct_bias_total = 100 * (sum(pred, na.rm = TRUE) - sum(d$n_annotations)) / sum(d$n_annotations))
  })
  res %>% arrange(.data[[metric]])
}

# ======================================================================
# 8. Block-bootstrap uncertainty for holdout metrics
# ======================================================================
# Resamples whole transect blocks (within region) of the holdout set, so
# spatially correlated neighbours move together. Reports 95% intervals for
# each method's RMSE / 1:1 R2 / total % bias and for the RMSE difference
# relative to `ref`.
#   pred_cols: named character vector, label -> column (e.g. c("Synergy" = "pred_synergy"))
block_bootstrap_metrics <- function(df, pred_cols, truth = "n_annotations",
                                    block = "transect_block_id", ref = names(pred_cols)[1],
                                    B = 1000, seed = 2226) {
  set.seed(seed)
  one <- function(d) {
    y <- d[[truth]]
    sapply(pred_cols, function(col) {
      yhat <- d[[col]]
      c(rmse = sqrt(mean((y - yhat)^2, na.rm = TRUE)),
        r2_1to1 = r2_1to1(y, yhat),
        pct_bias_total = 100 * (sum(yhat, na.rm = TRUE) - sum(y)) / sum(y))
    })
  }
  purrr::map_dfr(split(df, df$region), function(dr) {
    idx_by_block <- split(seq_len(nrow(dr)), dr[[block]])
    nb <- length(idx_by_block)
    point <- one(dr)
    boots <- replicate(B, {
      rows <- unlist(idx_by_block[sample.int(nb, nb, replace = TRUE)], use.names = FALSE)
      one(dr[rows, , drop = FALSE])
    }, simplify = "array")                          # 3 x methods x B
    purrr::map_dfr(names(pred_cols), function(m) {
      q <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
      d_rmse <- boots["rmse", m, ] - boots["rmse", ref, ]
      tibble(region = as.character(dr$region[1]), method = m,
             rmse = point["rmse", m], rmse_lo = q(boots["rmse", m, ])[1], rmse_hi = q(boots["rmse", m, ])[2],
             r2_1to1 = point["r2_1to1", m], r2_lo = q(boots["r2_1to1", m, ])[1], r2_hi = q(boots["r2_1to1", m, ])[2],
             pct_bias_total = point["pct_bias_total", m],
             bias_lo = q(boots["pct_bias_total", m, ])[1], bias_hi = q(boots["pct_bias_total", m, ])[2],
             d_rmse_vs_ref = point["rmse", m] - point["rmse", ref],
             d_rmse_lo = q(d_rmse)[1], d_rmse_hi = q(d_rmse)[2],
             n_images = nrow(dr), n_blocks = nb)
    })
  })
}

# ======================================================================
# 9. Image-level detection-structure ("crowding") features
# ======================================================================
# Computed from DETECTOR OUTPUT ONLY (boxes + confidences), so they are
# available for new, unannotated images and carry no label information.
#   n_det_all    detections at any confidence (>= model floor, 0.1)
#   nn_ratio     median nearest-neighbour centre distance / box size (capped);
#                small = scallops packed together
#   overlap_frac share of boxes overlapping another box (IoU > iou_thr);
#                high = possible NMS merging (evaluation used NMS 0.65, match IoU 0.1)
#   lowconf_frac share of detections below the regional F1 threshold
#   med_box_mm   median box size (sqrt area) in mm; small = juvenile-sized
#                animals, which the 2026-10-06 prototype linked to missed scallops
# Prototype (2026-10-06, 202404 images with >= 20 scallops): recall vs nn_ratio
# rho = 0.00, overlap rho = 0.04, box size rho = +0.34, lowconf rho = -0.34.
image_detection_features <- function(det_df, iou_thr = 0.05, nn_cap = 5) {
  one <- function(g) {
    x1 <- g$t_lx; y1 <- g$t_ly; x2 <- g$b_rx; y2 <- g$b_ry
    n  <- length(x1)
    s  <- sqrt(pmax((x2 - x1) * (y2 - y1), 1))
    mmpp <- suppressWarnings(as.numeric(g$millimeter_per_pixel[1]))
    nn <- nn_cap; ov <- 0
    if (n >= 2) {
      cx <- (x1 + x2) / 2; cy <- (y1 + y2) / 2
      D <- sqrt(outer(cx, cx, "-")^2 + outer(cy, cy, "-")^2); diag(D) <- Inf
      nn <- min(median(apply(D, 1, min) / s), nn_cap)
      iw <- outer(x2, x2, pmin) - outer(x1, x1, pmax); iw[iw < 0] <- 0   # keep matrix dims
      ih <- outer(y2, y2, pmin) - outer(y1, y1, pmax); ih[ih < 0] <- 0
      inter <- iw * ih; a <- (x2 - x1) * (y2 - y1)
      iou <- inter / (outer(a, a, "+") - inter + 1e-9); diag(iou) <- 0
      ov <- mean(apply(iou, 1, max) > iou_thr)
    }
    tibble(n_det_all = n, nn_ratio = nn, overlap_frac = ov,
           lowconf_frac = mean(g$conf < g$regional_f1_thresh, na.rm = TRUE),
           med_box_mm = median(s) * mmpp)
  }
  det_df %>%
    filter(!is.na(t_lx), !is.na(b_rx)) %>%
    group_by(image_id) %>%
    group_modify(~ one(.x)) %>%
    ungroup()
}

CROWD_VARS <- c("n_det_all", "nn_ratio", "overlap_frac", "lowconf_frac", "med_box_mm")

# Join features to the synergy table; images without detections get neutral
# values (0 detections, isolated, no overlap, no low-confidence share, regional
# median box size) so they stay in the model rather than dropping as NA.
add_crowding_features <- function(img_df, feats, nn_cap = 5) {
  img_df %>%
    select(-any_of(CROWD_VARS)) %>%
    left_join(feats, by = "image_id") %>%
    group_by(region) %>%
    mutate(n_det_all = coalesce(n_det_all, 0L),
           nn_ratio = coalesce(nn_ratio, nn_cap),
           overlap_frac = coalesce(overlap_frac, 0),
           lowconf_frac = coalesce(lowconf_frac, 0),
           med_box_mm = coalesce(med_box_mm, median(med_box_mm, na.rm = TRUE))) %>%
    ungroup()
}
