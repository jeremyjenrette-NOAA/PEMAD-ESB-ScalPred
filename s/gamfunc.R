library(dplyr)
library(tidyr)
library(mgcv)
library(ggplot2)
library(patchwork)

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
      Formula = deparse(form),
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
fit_calibration_gams <- function(det_df, formula_gb, formula_mab) {
  
  train_data <- det_df %>%
    filter(dataset == "test_GAM_train") %>%
    drop_na(truedetect, conf, bottom_depth, latitude, longitude, altitude, 
            field_of_view_sq_meter, fluorometer_backscatter_ntu, ctd_temperature_celsius,
            boxsize, iu)
  
  cat("  Fitting Georges Bank GAM...\n")
  m_gb <- gam(
    formula_gb,
    family = binomial(),
    data = train_data %>% filter(region == "GB"),
    method = "REML"
  )
  
  cat("  Fitting Mid-Atlantic Bight GAM...\n")
  m_mab <- gam(
    formula_mab,
    family = binomial(),
    data = train_data %>% filter(region == "MAB"),
    method = "REML"
  )
  
  return(list(GB = m_gb, MAB = m_mab))
}
# ======================================================================
# 2. Test GAMs on Holdout Data (Dynamic F1 Thresholds)
# ======================================================================
# ======================================================================
# UPDATED: Safely Predict (Processes all data, evaluates only holdout)
# ======================================================================
test_calibration_gams <- function(gams, det_df, img_df, f1_thresholds) {
  
  # We NO LONGER filter out the training data here! We want predictions for EVERYTHING.
  det_eval <- det_df 
  img_eval_base <- img_df 
  
  # 1. Safely assign regional F1 thresholds
  det_eval$regional_f1_thresh <- ifelse(det_eval$region == "GB", f1_thresholds[["GB"]], f1_thresholds[["MAB"]])
  
  # 2. Safely predict by explicitly subsetting data (prevents silent vector recycling)
  det_eval$pred_p <- NA_real_
  
  gb_idx <- which(det_eval$region == "GB")
  if (length(gb_idx) > 0) {
    det_eval$pred_p[gb_idx] <- predict(gams$GB, newdata = det_eval[gb_idx, ], type = "response", na.action = na.pass)
  }
  
  mab_idx <- which(det_eval$region == "MAB")
  if (length(mab_idx) > 0) {
    det_eval$pred_p[mab_idx] <- predict(gams$MAB, newdata = det_eval[mab_idx, ], type = "response", na.action = na.pass)
  }
  
  # 3. Aggregate predictions to the image level
  img_predictions <- det_eval %>%
    group_by(image_id) %>%
    summarise(
      raw_detection_number = n(),
      predicted_f1_number  = sum(conf >= regional_f1_thresh, na.rm = TRUE),
      predicted_number     = sum(pred_p, na.rm = TRUE),
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
      true_positive_sum    = replace_na(true_positive_sum, 0),
      false_negative       = n_annotations - true_positive_sum
    )
  
  # 5. Calculate metrics STRICTLY on the holdout test set to avoid data leakage
  metrics <- img_eval %>%
    filter(dataset == "test_GAM_test") %>%  # <-- Crucial line!
    group_by(region) %>%
    summarise(
      r2_calibrated   = summary(lm(n_annotations ~ predicted_number))$adj.r.squared,
      rmse_calibrated = sqrt(mean((n_annotations - predicted_number)^2, na.rm = TRUE)),
      r2_f1           = summary(lm(n_annotations ~ predicted_f1_number))$adj.r.squared,
      rmse_f1         = sqrt(mean((n_annotations - predicted_f1_number)^2, na.rm = TRUE)),
      n_images        = n(),
      .groups = "drop"
    )
  
  return(list(img_eval = img_eval, det_eval = det_eval, metrics = metrics))
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
      label = sprintf("R² = %.2f\nRMSE = %.2f", r2_calibrated, rmse_calibrated)
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
    select(image_id, region, dataset, n_annotations, pred_yolo = predicted_number) %>%
    left_join(
      cas_eval$img_eval %>%
        filter(image_id %in% common_ids) %>%
        select(image_id, pred_cascade = predicted_number),
      by = "image_id"
    ) %>%
    mutate(
      pred_mean = (pred_yolo + pred_cascade) / 2,
      pred_diff = pred_cascade - pred_yolo,
      pred_sum  = pred_yolo + pred_cascade,
      log_yolo_pred = log1p(pred_yolo),
      log_cascade_pred = log1p(pred_cascade),
      log_pred_mean = log1p(pred_mean),
      log_pred_sum = log1p(pred_sum),
      log_pred_diff = log1p(abs(pred_diff)),
      log_true_number = log1p(n_annotations + 1e-6)
    )
  
  # Attach environmental metadata
  img_combined <- img_combined %>%
    left_join(
      meta_df %>% select(-any_of(c("region", "dataset", "n_annotations", "imagename"))), 
      by = "image_id"
    )
  
  return(img_combined)
}

# ======================================================================
# 5. Compare Candidate Synergistic Image-Level GAMs (Stratified)
# ======================================================================
compare_image_gams <- function(img_combined, candidate_formulas, region_focus = "GB") {
  
  train_data <- img_combined %>% filter(dataset == "test_GAM_train", region == region_focus)
  
  cat("\nEvaluating", length(candidate_formulas), "synergistic models for", region_focus, "(n =", nrow(train_data), ")...\n")
  
  results <- purrr::map_dfr(names(candidate_formulas), function(mod_name) {
    form <- candidate_formulas[[mod_name]]
    fit <- gam(form, family = nb(), data = train_data, method = "REML")
    
    # Generate predictions once to use for both RMSE and R-squared
    preds <- predict(fit, type = "response")
    
    tibble(
      Model_ID = mod_name,
      Formula = deparse(form),
      AIC = AIC(fit),
      Deviance_Explained = summary(fit)$dev.expl,
      RMSE_Train = sqrt(mean((train_data$n_annotations - preds)^2)),
      R2_Train = summary(lm(train_data$n_annotations ~ preds))$adj.r.squared
    )
  })
  
  # Sort by lowest AIC
  return(results %>% arrange(AIC))
}

# ======================================================================
# 6. Fit and Test Stratified Synergistic GAMs
# ======================================================================
fit_synergy_gams <- function(img_combined, formula_gb, formula_mab) {
  
  train_data <- img_combined %>% filter(dataset == "test_GAM_train")
  
  cat("  Fitting Georges Bank Synergy GAM...\n")
  m_gb <- gam(formula_gb, family = nb(), data = train_data %>% filter(region == "GB"), method = "REML")
  
  cat("  Fitting Mid-Atlantic Bight Synergy GAM...\n")
  m_mab <- gam(formula_mab, family = nb(), data = train_data %>% filter(region == "MAB"), method = "REML")
  
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
      r2_synergy   = summary(lm(n_annotations ~ pred_synergy))$adj.r.squared,
      rmse_synergy = sqrt(mean((n_annotations - pred_synergy)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  return(list(img_eval = test_data, metrics = metrics))
}
