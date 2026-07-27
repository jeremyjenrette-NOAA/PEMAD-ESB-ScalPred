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
  results <- results %>% arrange(Deviance_Explained) # changed from AIC
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
# UPDATED: Multiclass Calibration Testing & Species-Level Evaluation
# ======================================================================
test_calibration_gams <- function(gams, det_df, img_df, f1_thresh_df, 
                                  target_classes = c("jonah_crab", "rock_crab", "cancer_sp")) {
  
  det_eval <- det_df 
  img_eval_base <- img_df 
  
  # 1. Join species- and region-specific F1 thresholds
  det_eval <- det_eval %>%
    left_join(
      f1_thresh_df %>% select(region, species, regional_f1_thresh = conf),
      by = c("region", "spname" = "species")
    ) %>%
    mutate(regional_f1_thresh = replace_na(regional_f1_thresh, 0.5))
  
  # 2. Predict probability of correct detection using regional GAMs
  det_eval$pred_p <- NA_real_
  
  gb_idx <- which(det_eval$region == "GB")
  if (length(gb_idx) > 0) {
    det_eval$pred_p[gb_idx] <- predict(gams$GB, newdata = det_eval[gb_idx, ], type = "response", na.action = na.pass)
  }
  
  mab_idx <- which(det_eval$region == "MAB")
  if (length(mab_idx) > 0) {
    det_eval$pred_p[mab_idx] <- predict(gams$MAB, newdata = det_eval[mab_idx, ], type = "response", na.action = na.pass)
  }
  
  # 3. Aggregate predictions to image level BY SPECIES
  sp_predictions <- det_eval %>%
    group_by(image_id, spname) %>%
    summarise(
      pred_f1  = sum(conf >= regional_f1_thresh, na.rm = TRUE),
      pred_gam = sum(pred_p, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    pivot_wider(
      id_cols = image_id,
      names_from = spname,
      values_from = c(pred_f1, pred_gam),
      names_glue = "{.value}_{spname}",
      values_fill = 0
    )
  
  # 4. Aggregate TOTAL predictions across all species
  tot_predictions <- det_eval %>%
    group_by(image_id) %>%
    summarise(
      raw_detection_number = n(),
      predicted_f1_number  = sum(conf >= regional_f1_thresh, na.rm = TRUE),
      predicted_number     = sum(pred_p, na.rm = TRUE),
      true_positive_sum    = sum(truedetect, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Ensure total annotations exist in base image df
  if (!"n_annotations" %in% names(img_eval_base)) {
    gt_cols <- intersect(paste0("n_", target_classes), names(img_eval_base))
    img_eval_base$n_annotations <- rowSums(img_eval_base[, gt_cols, drop = FALSE], na.rm = TRUE)
  }
  
  # 5. Join predictions back to master image dataset
  img_eval <- img_eval_base %>%
    left_join(tot_predictions, by = "image_id") %>%
    left_join(sp_predictions, by = "image_id") %>%
    mutate(
      across(starts_with("pred_f1_"), ~ replace_na(.x, 0)),
      across(starts_with("pred_gam_"), ~ replace_na(.x, 0)),
      raw_detection_number = replace_na(raw_detection_number, 0),
      predicted_f1_number  = replace_na(predicted_f1_number, 0),
      predicted_number     = replace_na(predicted_number, 0),
      true_positive_sum    = replace_na(true_positive_sum, 0),
      false_negative       = n_annotations - true_positive_sum
    )
  
  # 6. Calculate holdout evaluation metrics (STRICTLY on test set)
  holdout <- img_eval %>% filter(dataset == "test_GAM_test")
  
  # A. Species-Specific Metrics
  metrics_sp <- purrr::map_dfr(target_classes, function(sp) {
    gt_col   <- paste0("n_", sp)
    pred_col <- paste0("pred_gam_", sp)
    f1_col   <- paste0("pred_f1_", sp)
    
    if (!gt_col %in% names(holdout) || !pred_col %in% names(holdout)) return(NULL)
    
    holdout %>%
      group_by(region) %>%
      summarise(
        species         = sp,
        r2_calibrated   = summary(lm(get(gt_col) ~ get(pred_col)))$adj.r.squared,
        rmse_calibrated = sqrt(mean((get(gt_col) - get(pred_col))^2, na.rm = TRUE)),
        r2_f1           = summary(lm(get(gt_col) ~ get(f1_col)))$adj.r.squared,
        rmse_f1         = sqrt(mean((get(gt_col) - get(f1_col))^2, na.rm = TRUE)),
        n_images        = n(),
        .groups         = "drop"
      )
  })
  
  # B. Overall Total Crab Metrics
  metrics_tot <- holdout %>%
    group_by(region) %>%
    summarise(
      species         = "TOTAL_CRABS",
      r2_calibrated   = summary(lm(n_annotations ~ predicted_number))$adj.r.squared,
      rmse_calibrated = sqrt(mean((n_annotations - predicted_number)^2, na.rm = TRUE)),
      r2_f1           = summary(lm(n_annotations ~ predicted_f1_number))$adj.r.squared,
      rmse_f1         = sqrt(mean((n_annotations - predicted_f1_number)^2, na.rm = TRUE)),
      n_images        = n(),
      .groups         = "drop"
    )
  
  metrics <- bind_rows(metrics_sp, metrics_tot)
  
  return(list(img_eval = img_eval, det_eval = det_eval, metrics = metrics))
}

# ======================================================================
# UPDATED: Multiclass Diagnostic Plots with Metrics Overlay & F1 Comparison
# ======================================================================
generate_gam_plots <- function(
    eval_res, 
    model_name, 
    dataset_label = "2024, 2026", 
    target_classes = c("jonah_crab", "rock_crab", "cancer_sp")
) {
  
  # 1. Isolate Holdout Test Data & Metrics
  img_df  <- eval_res$img_eval %>% filter(dataset == "test_GAM_test")
  det_df  <- eval_res$det_eval %>% filter(dataset == "test_GAM_test")
  metrics <- eval_res$metrics
  
  # Enhance region labels with image count (n)
  region_counts <- img_df %>%
    group_by(region) %>%
    summarise(n_img = n(), .groups = "drop") %>%
    mutate(
      region_full = case_when(
        region == "GB" ~ "Georges Bank",
        region == "MAB" ~ "Mid-Atlantic Bight",
        TRUE ~ as.character(region)
      ),
      facet_region_label = paste0(region_full, " (n = ", n_img, " images)")
    )
  
  img_df <- img_df %>% left_join(region_counts, by = "region")
  
  # --------------------------------------------------------------------
  # A. Species-Level Calibration Plot with GAM vs F1 Metrics
  # --------------------------------------------------------------------
  sp_long <- img_df %>%
    pivot_longer(
      cols = all_of(paste0("n_", target_classes)),
      names_to = "gt_sp", names_prefix = "n_", values_to = "true_count"
    ) %>%
    pivot_longer(
      cols = all_of(paste0("pred_gam_", target_classes)),
      names_to = "pred_sp", names_prefix = "pred_gam_", values_to = "pred_count"
    ) %>%
    filter(gt_sp == pred_sp) %>%
    mutate(species_clean = gsub("_", " ", tools::toTitleCase(gt_sp)))
  
  # Format metrics annotation labels per species & region
  sp_metrics_labels <- metrics %>%
    filter(species != "TOTAL_CRABS") %>%
    left_join(region_counts, by = "region") %>%
    mutate(
      species_clean = gsub("_", " ", tools::toTitleCase(species)),
      gt_sp = species,
      label_text = sprintf(
        "GAM: R²=%.2f, RMSE=%.2f\nF1:   R²=%.2f, RMSE=%.2f",
        r2_calibrated, rmse_calibrated, r2_f1, rmse_f1
      )
    )
  
  p_species_fit <- ggplot(sp_long, aes(x = pred_count, y = true_count)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = species_clean), alpha = 0.6, size = 1.8) +
    facet_grid(facet_region_label ~ species_clean, scales = "free") +
    geom_text(
      data = sp_metrics_labels,
      aes(x = -Inf, y = Inf, label = label_text),
      hjust = -0.05, vjust = 1.15, inherit.aes = FALSE,
      size = 3.1, fontface = "bold", family = "mono", color = "grey20"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Species-Level Abundance Calibration: ", model_name),
      subtitle = paste0("Evaluated on Holdout Test Set (Dataset: ", dataset_label, ")"),
      x = "GAM Calibrated Count (Σ p)",
      y = "Manual Ground Truth Count",
      color = "Species"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "none"
    )
  
  # --------------------------------------------------------------------
  # B. Total Crab Abundance Fit with Metrics Overlay
  # --------------------------------------------------------------------
  tot_metrics_labels <- metrics %>%
    filter(species == "TOTAL_CRABS") %>%
    left_join(region_counts, by = "region") %>%
    mutate(
      label_text = sprintf(
        "GAM Calibrated:\n  R² = %.2f | RMSE = %.2f\nF1 Cutoff:\n  R² = %.2f | RMSE = %.2f",
        r2_calibrated, rmse_calibrated, r2_f1, rmse_f1
      )
    )
  
  p_zoomed <- ggplot(img_df, aes(x = predicted_number, y = n_annotations)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(color = "#2C7FB8", size = 2, alpha = 0.6) +
    facet_wrap(~ facet_region_label, scales = "free") +
    geom_text(
      data = tot_metrics_labels,
      aes(x = -Inf, y = Inf, label = label_text),
      hjust = -0.08, vjust = 1.18, inherit.aes = FALSE,
      size = 3.6, fontface = "bold", color = "grey15"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Total Crab Abundance Calibration: ", model_name),
      subtitle = paste0("Evaluated on Holdout Test Set (Dataset: ", dataset_label, ")"),
      x = "Total GAM Calibrated Count (Σ p)",
      y = "Total Manual Ground Truth Count"
    ) +
    theme(strip.text = element_text(face = "bold", size = 11))
  
  # --------------------------------------------------------------------
  # C. Residual Analysis
  # --------------------------------------------------------------------
  p_resid <- img_df %>%
    mutate(residual = predicted_number - n_annotations) %>%
    ggplot(aes(x = n_annotations, y = residual)) +
    geom_point(alpha = 0.5, color = "#2C7FB8") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ facet_region_label, scales = "free_x") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Residuals vs True Abundance: ", model_name),
      subtitle = paste0("Holdout Test Set (Dataset: ", dataset_label, ")"),
      x = "Manual Ground Truth Count",
      y = "Residual (GAM Predicted - True)"
    ) +
    theme(strip.text = element_text(face = "bold"))
  
  # --------------------------------------------------------------------
  # D. Species Macro-Abundance Comparison Bar Chart
  # --------------------------------------------------------------------
  sum_by_sp <- purrr::map_dfr(target_classes, function(sp) {
    gt_col   <- paste0("n_", sp)
    pred_gam <- paste0("pred_gam_", sp)
    pred_f1  <- paste0("pred_f1_", sp)
    
    raw_cnt  <- sum(det_df$spname == sp, na.rm = TRUE)
    gt_cnt   <- sum(img_df[[gt_col]], na.rm = TRUE)
    f1_cnt   <- sum(img_df[[pred_f1]], na.rm = TRUE)
    gam_cnt  <- sum(img_df[[pred_gam]], na.rm = TRUE)
    
    tibble(
      species = gsub("_", " ", tools::toTitleCase(sp)),
      `True Count`     = gt_cnt,
      `Raw Detector`   = raw_cnt,
      `F1 Cutoff`      = f1_cnt,
      `GAM Calibrated` = gam_cnt
    )
  }) %>%
    pivot_longer(-species, names_to = "Approach", values_to = "Count") %>%
    mutate(Approach = factor(Approach, levels = c("True Count", "Raw Detector", "F1 Cutoff", "GAM Calibrated")))
  
  p_sum <- ggplot(sum_by_sp, aes(x = species, y = Count, fill = Approach)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    geom_text(
      aes(label = round(Count, 0)),
      position = position_dodge(width = 0.8),
      vjust = -0.4, size = 3.2, fontface = "bold"
    ) +
    scale_fill_manual(values = c("grey40", "salmon", "#D95F0E", "#2C7FB8")) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("Macro Abundance Estimation by Species: ", model_name),
      subtitle = paste0("Total Count Comparison across Holdout Test Set (n = ", nrow(img_df), " images)"),
      x = "Species",
      y = "Total Abundance Count",
      fill = "Approach"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold")
    )
  
  return(list(
    p_species_fit = p_species_fit,
    p_zoomed      = p_zoomed,
    p_resid       = p_resid,
    p_sum         = p_sum
  ))
}

# ======================================================================
# 4. Build Synergistic Image-Level Dataset
# ======================================================================
build_synergy_dataset_multi <- function(yolo_eval, cas_eval, meta_df, 
                                        target_classes = c("jonah_crab", "rock_crab", "cancer_sp")) {
  
  # 1. Identify common image IDs across both holdout sets
  common_ids <- intersect(yolo_eval$img_eval$image_id, cas_eval$img_eval$image_id)
  
  gt_cols <- paste0("n_", target_classes)
  
  # 2. Extract YOLO predictions & Ground Truth counts
  yolo_sp <- yolo_eval$img_eval %>%
    filter(image_id %in% common_ids) %>%
    select(
      image_id, region, dataset, n_annotations, 
      any_of(gt_cols), 
      starts_with("pred_gam_")
    ) %>%
    rename_with(~ paste0("yolo_", .x), starts_with("pred_gam_"))
  
  # 3. Extract Cascade predictions ONLY (avoids repeating GT cols)
  cas_sp <- cas_eval$img_eval %>%
    filter(image_id %in% common_ids) %>%
    select(image_id, starts_with("pred_gam_")) %>%
    rename_with(~ paste0("cascade_", .x), starts_with("pred_gam_"))
  
  # 4. Join model outputs
  img_combined <- yolo_sp %>%
    left_join(cas_sp, by = "image_id")
  
  # 5. Calculate per-species diffs and log-transforms
  for (sp in target_classes) {
    y_col <- paste0("yolo_pred_gam_", sp)
    c_col <- paste0("cascade_pred_gam_", sp)
    
    if (y_col %in% names(img_combined) && c_col %in% names(img_combined)) {
      img_combined[[paste0("log_yolo_", sp)]]    <- log1p(img_combined[[y_col]])
      img_combined[[paste0("log_cascade_", sp)]] <- log1p(img_combined[[c_col]])
      img_combined[[paste0("diff_", sp)]]        <- img_combined[[c_col]] - img_combined[[y_col]]
    }
  }
  
  # 6. Safely attach environmental metadata (Drop ALL existing columns to prevent .x/.y suffixes)
  existing_cols <- setdiff(names(img_combined), "image_id")
  
  img_combined <- img_combined %>%
    left_join(
      meta_df %>% select(-any_of(existing_cols), -any_of("imagename")), 
      by = "image_id"
    )
  
  return(img_combined)
}

# Fit per-species synergistic GAMs
fit_synergy_gams_multi <- function(img_combined, formula_list_gb, formula_list_mab, 
                                   target_classes = c("jonah_crab", "rock_crab", "cancer_sp")) {
  train_data <- img_combined %>% filter(dataset == "test_GAM_train")
  
  gams <- list()
  
  for (sp in target_classes) {
    cat("Fitting Synergy GAMs for:", sp, "\n")
    form_gb  <- formula_list_gb[[sp]]
    form_mab <- formula_list_mab[[sp]]
    
    m_gb  <- gam(form_gb,  family = nb(), data = train_data %>% filter(region == "GB"),  method = "REML")
    m_mab <- gam(form_mab, family = nb(), data = train_data %>% filter(region == "MAB"), method = "REML")
    
    gams[[sp]] <- list(GB = m_gb, MAB = m_mab)
  }
  
  return(gams)
}

# Test per-species synergistic GAMs on holdout set
test_synergy_gams_multi <- function(gams_multi, img_combined, 
                                    target_classes = c("jonah_crab", "rock_crab", "cancer_sp")) {
  
  test_data <- img_combined %>% filter(dataset == "test_GAM_test")
  
  for (sp in target_classes) {
    pred_col <- paste0("pred_synergy_", sp)
    test_data[[pred_col]] <- NA_real_
    
    gb_idx  <- which(test_data$region == "GB")
    mab_idx <- which(test_data$region == "MAB")
    
    if (length(gb_idx) > 0) {
      test_data[[pred_col]][gb_idx] <- predict(gams_multi[[sp]]$GB, newdata = test_data[gb_idx, ], type = "response")
    }
    if (length(mab_idx) > 0) {
      test_data[[pred_col]][mab_idx] <- predict(gams_multi[[sp]]$MAB, newdata = test_data[mab_idx, ], type = "response")
    }
  }
  
  # Total synergy count = sum of individual species synergy predictions
  synergy_pred_cols <- paste0("pred_synergy_", target_classes)
  test_data$pred_synergy_total <- rowSums(test_data[, synergy_pred_cols, drop = FALSE], na.rm = TRUE)
  
  # Calculate metrics per species & total
  metrics <- purrr::map_dfr(target_classes, function(sp) {
    gt_col   <- paste0("n_", sp)
    pred_col <- paste0("pred_synergy_", sp)
    
    test_data %>%
      group_by(region) %>%
      summarise(
        species      = sp,
        r2_synergy   = summary(lm(get(gt_col) ~ get(pred_col)))$adj.r.squared,
        rmse_synergy = sqrt(mean((get(gt_col) - get(pred_col))^2, na.rm = TRUE)),
        .groups      = "drop"
      )
  })
  
  tot_metrics <- test_data %>%
    group_by(region) %>%
    summarise(
      species      = "TOTAL_CRABS",
      r2_synergy   = summary(lm(n_annotations ~ pred_synergy_total))$adj.r.squared,
      rmse_synergy = sqrt(mean((n_annotations - pred_synergy_total)^2, na.rm = TRUE)),
      .groups      = "drop"
    )
  
  return(list(img_eval = test_data, metrics = bind_rows(metrics, tot_metrics)))
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

# ======================================================================
# FIXED: Multiclass Synergistic Model Comparison & Extraction
# ======================================================================
compare_synergy_gams_multi <- function(
    img_combined, 
    candidate_templates, 
    region_focus = "GB", 
    target_classes = c("jonah_crab", "rock_crab", "cancer_sp")
) {
  
  train_data <- img_combined %>% filter(dataset == "test_GAM_train", region == region_focus)
  
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("Synergy Model Selection: %s (n = %d images)\n", region_focus, nrow(train_data)))
  cat(sprintf("==================================================\n"))
  
  results <- purrr::map_dfr(target_classes, function(sp) {
    gt_var <- paste0("n_", sp)
    
    purrr::map_dfr(names(candidate_templates), function(mod_name) {
      # 1. Substitute species placeholder into template
      template_str <- candidate_templates[[mod_name]]
      form_str     <- gsub("\\{sp\\}", sp, template_str)
      form         <- as.formula(form_str)
      
      # 2. Fit model safely with error handling
      fit <- tryCatch({
        gam(form, family = nb(), data = train_data, method = "REML")
      }, error = function(e) {
        cat(sprintf("   [Warning] Model %s failed for %s: %s\n", mod_name, sp, e$message))
        return(NULL)
      })
      
      if (is.null(fit)) return(NULL)
      
      # 3. Calculate evaluation metrics on training set
      preds   <- predict(fit, type = "response")
      gt_vals <- train_data[[gt_var]]
      
      # FIX: Use deparse1() or collapse deparse() output into a single string
      formula_single_str <- paste(deparse(form), collapse = " ")
      
      tibble(
        region             = region_focus,
        species            = sp,
        Model_ID           = mod_name,
        Formula            = formula_single_str, # Single clean string
        AIC                = AIC(fit),
        Deviance_Explained = summary(fit)$dev.expl * 100,
        RMSE_Train         = sqrt(mean((gt_vals - preds)^2, na.rm = TRUE)),
        R2_Train           = summary(lm(gt_vals ~ preds))$adj.r.squared
      )
    })
  })
  
  # Group by species and rank by lowest AIC
  results <- results %>%
    group_by(species) %>%
    arrange(AIC, .by_group = TRUE) %>%
    ungroup()
  
  return(results)
}

extract_best_formulas <- function(selection_results, metric = "AIC") {
  best_df <- selection_results %>%
    group_by(species) %>%
    {
      if (metric == "AIC") slice_min(., order_by = AIC, n = 1)
      else if (metric == "Deviance") slice_max(., order_by = Deviance_Explained, n = 1)
      else slice_min(., order_by = RMSE_Train, n = 1)
    } %>%
    ungroup()
  
  out_list <- list()
  for (i in seq_len(nrow(best_df))) {
    sp <- best_df$species[i]
    # FIX: Ensure clean single string conversion to formula
    form_str <- paste(best_df$Formula[[i]], collapse = " ")
    out_list[[sp]] <- as.formula(form_str)
  }
  return(out_list)
}