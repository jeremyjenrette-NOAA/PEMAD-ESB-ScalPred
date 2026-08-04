library(dplyr)
library(tidyr)
library(mgcv)
library(ggplot2)
library(patchwork)
library(purrr)

# ======================================================================
# Helper: Autolocate Classes from Datasets
# ======================================================================
get_target_classes <- function(img_data = NULL, det_data = NULL) {
  classes <- character(0)
  
  if (!is.null(img_data)) {
    img_classes <- sub("^n_", "", grep("^n_", names(img_data), value = TRUE))
    classes <- c(classes, img_classes)
  }
  
  if (!is.null(det_data)) {
    det_classes <- c(det_data$gt_label, det_data$spname)
    classes <- c(classes, det_classes)
  }
  
  reserved <- c("Background", "background", "Missed", "missed", "annotations", "auto", "total")
  valid_classes <- unique(classes[!is.na(classes) & classes != "" & !classes %in% reserved])
  valid_classes <- valid_classes[!grepl("^auto_", valid_classes)]
  
  return(valid_classes)
}

# ======================================================================
# 1. Compare Candidate Detection-Level GAM Structures
# ======================================================================
compare_detection_gams <- function(det_df, candidate_formulas, region_focus = "GB") {
  train_data <- det_df
  
  if ("dataset" %in% names(train_data)) {
    train_data <- train_data %>% filter(dataset == "test_GAM_train")
  }
  if ("region" %in% names(train_data) && !is.null(region_focus)) {
    train_data <- train_data %>% filter(region == region_focus)
  }
  
  train_data <- train_data %>% drop_na(truedetect, conf)
  
  cat("\nEvaluating", length(candidate_formulas), "candidate models (n =", nrow(train_data), ")...\n")
  
  results <- purrr::map_dfr(names(candidate_formulas), function(mod_name) {
    form <- candidate_formulas[[mod_name]]
    cat("  Fitting:", mod_name, "...\n")
    
    fit <- tryCatch({
      gam(form, family = binomial(), data = train_data, method = "REML", gamma = 1.4)
    }, error = function(e) {
      cat("   [Warning] Model failed:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(fit)) return(NULL)
    
    tibble(
      Model_ID           = mod_name,
      Formula            = paste(deparse(form), collapse = " "),
      AIC                = AIC(fit),
      Deviance_Explained = summary(fit)$dev.expl,
      UBRE_Score         = fit$gcv.ubre
    )
  })
  
  return(results %>% arrange(AIC))
}

# ======================================================================
# 2. Build Synergistic Image-Level Dataset
# ======================================================================
build_synergy_dataset_multi <- function(
    yolo_eval, 
    cas_eval, 
    meta, 
    target_classes = NULL,
    conf_thresh = 0.10
) {
  if (is.null(target_classes)) {
    target_classes <- get_target_classes(img_data = meta, det_data = yolo_eval$det)
  }
  
  yolo_counts <- yolo_eval$det %>%
    filter(conf >= conf_thresh) %>%
    count(image_id, spname, name = "count") %>%
    pivot_wider(
      id_cols = image_id, 
      names_from = spname, 
      values_from = count, 
      values_fill = 0,
      names_prefix = "yolo_"
    )
  
  cas_counts <- cas_eval$det %>%
    filter(conf >= conf_thresh) %>%
    count(image_id, spname, name = "count") %>%
    pivot_wider(
      id_cols = image_id, 
      names_from = spname, 
      values_from = count, 
      values_fill = 0,
      names_prefix = "cascade_"
    )
  
  img_combined <- meta %>%
    left_join(yolo_counts, by = "image_id") %>%
    left_join(cas_counts, by = "image_id")
  
  for (sp in target_classes) {
    yolo_col <- paste0("yolo_", sp)
    cas_col  <- paste0("cascade_", sp)
    
    if (!yolo_col %in% names(img_combined)) img_combined[[yolo_col]] <- 0
    if (!cas_col %in% names(img_combined))  img_combined[[cas_col]]  <- 0
    
    img_combined[[yolo_col]] <- replace_na(img_combined[[yolo_col]], 0)
    img_combined[[cas_col]]  <- replace_na(img_combined[[cas_col]], 0)
    
    img_combined[[paste0("log_yolo_", sp)]]    <- log1p(img_combined[[yolo_col]])
    img_combined[[paste0("log_cascade_", sp)]] <- log1p(img_combined[[cas_col]])
    img_combined[[paste0("diff_", sp)]]        <- abs(img_combined[[yolo_col]] - img_combined[[cas_col]])
  }
  
  return(img_combined)
}

compare_synergy_gams_multi <- function(
    img_combined, 
    candidate_templates, 
    region_focus = "GB", 
    target_classes = NULL
) {
  if (is.null(target_classes)) {
    target_classes <- get_target_classes(img_data = img_combined)
  }
  
  train_data <- img_combined
  
  if ("dataset" %in% names(train_data)) {
    train_data <- train_data %>% filter(dataset == "test_GAM_train")
  }
  
  if (!is.null(region_focus) && "region" %in% names(train_data)) {
    train_data <- train_data %>% filter(region == region_focus)
  } else {
    region_focus <- if ("region" %in% names(train_data)) "ALL" else "FULL_DATASET"
  }
  
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("Synergy Model Selection: %s (n = %d images)\n", region_focus, nrow(train_data)))
  cat(sprintf("==================================================\n"))
  
  results <- purrr::map_dfr(target_classes, function(sp) {
    gt_var <- paste0("n_", sp)
    
    if (!gt_var %in% names(train_data)) {
      cat(sprintf("   [Warning] Column '%s' not found. Skipping '%s'.\n", gt_var, sp))
      return(NULL)
    }
    
    # Calculate available degrees of freedom in current subset
    u_yolo <- max(3, min(5, dplyr::n_distinct(train_data[[paste0("log_yolo_", sp)]]) - 1))
    u_cas  <- max(3, min(5, dplyr::n_distinct(train_data[[paste0("log_cascade_", sp)]]) - 1))
    u_diff <- max(3, min(5, dplyr::n_distinct(train_data[[paste0("diff_", sp)]]) - 1))
    
    purrr::map_dfr(names(candidate_templates), function(mod_name) {
      template_str <- candidate_templates[[mod_name]]
      
      # Inject species names
      form_str <- gsub("\\{sp\\}", sp, template_str)
      
      # Dynamically inject basis dimension limits if k is not explicitly hardcoded
      if (!grepl("k\\s*=", form_str)) {
        form_str <- gsub(paste0("s\\(log_yolo_", sp, "\\)"), sprintf("s(log_yolo_%s, k = %d)", sp, u_yolo), form_str)
        form_str <- gsub(paste0("s\\(log_cascade_", sp, "\\)"), sprintf("s(log_cascade_%s, k = %d)", sp, u_cas), form_str)
        form_str <- gsub(paste0("s\\(diff_", sp, "\\)"), sprintf("s(diff_%s, k = %d)", sp, u_diff), form_str)
        form_str <- gsub(paste0("te\\(log_yolo_", sp, ", log_cascade_", sp, "\\)"), sprintf("te(log_yolo_%s, log_cascade_%s, k = c(%d, %d))", sp, sp, u_yolo, u_cas), form_str)
        form_str <- gsub(paste0("ti\\(log_yolo_", sp, ", log_cascade_", sp, "\\)"), sprintf("ti(log_yolo_%s, log_cascade_%s, k = c(%d, %d))", sp, sp, u_yolo, u_cas), form_str)
      }
      
      form <- as.formula(form_str)
      
      fit <- tryCatch({
        gam(form, family = nb(), data = train_data, method = "REML")
      }, error = function(e) {
        cat(sprintf("   [Warning] Model %s failed for %s: %s\n", mod_name, sp, e$message))
        return(NULL)
      })
      
      if (is.null(fit)) return(NULL)
      
      preds   <- predict(fit, type = "response")
      gt_vals <- train_data[[gt_var]]
      
      tibble(
        region             = region_focus,
        species            = sp,
        Model_ID           = mod_name,
        Formula            = paste(deparse(form), collapse = " "),
        AIC                = AIC(fit),
        Deviance_Explained = summary(fit)$dev.expl * 100,
        RMSE_Train         = sqrt(mean((gt_vals - preds)^2, na.rm = TRUE)),
        R2_Train           = summary(lm(gt_vals ~ preds))$adj.r.squared
      )
    })
  })
  
  if (nrow(results) == 0 || !"species" %in% names(results)) {
    cat("   [Warning] No GAM models were successfully fit. Returning empty output.\n")
    return(tibble(
      region = character(), species = character(), Model_ID = character(),
      Formula = character(), AIC = numeric(), Deviance_Explained = numeric(),
      RMSE_Train = numeric(), R2_Train = numeric()
    ))
  }
  
  return(results %>% group_by(species) %>% arrange(AIC, .by_group = TRUE) %>% ungroup())
}

# ======================================================================
# 4. Extract Best Formula per Species
# ======================================================================
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
    form_str <- paste(best_df$Formula[[i]], collapse = " ")
    out_list[[sp]] <- as.formula(form_str)
  }
  return(out_list)
}

# ======================================================================
# 5. Fit & Test Multiclass Synergy GAMs
# ======================================================================
fit_synergy_gams_multi <- function(
    img_combined, 
    formula_list_gb, 
    formula_list_mab = NULL, 
    target_classes = NULL
) {
  if (is.null(target_classes)) {
    target_classes <- get_target_classes(img_data = img_combined)
  }
  
  train_data <- img_combined
  if ("dataset" %in% names(train_data)) {
    train_data <- train_data %>% filter(dataset == "test_GAM_train")
  }
  
  has_region <- "region" %in% names(train_data)
  gams <- list()
  
  for (sp in target_classes) {
    cat("Fitting Synergy GAMs for:", sp, "\n")
    form_gb  <- if (is.list(formula_list_gb)) formula_list_gb[[sp]] else formula_list_gb
    
    if (has_region && !is.null(formula_list_mab)) {
      form_mab <- if (is.list(formula_list_mab)) formula_list_mab[[sp]] else formula_list_mab
      m_gb  <- gam(form_gb,  family = nb(), data = train_data %>% filter(region == "GB"),  method = "REML")
      m_mab <- gam(form_mab, family = nb(), data = train_data %>% filter(region == "MAB"), method = "REML")
      gams[[sp]] <- list(GB = m_gb, MAB = m_mab)
    } else {
      m_all <- gam(form_gb, family = nb(), data = train_data, method = "REML")
      gams[[sp]] <- list(ALL = m_all)
    }
  }
  
  return(gams)
}

test_synergy_gams_multi <- function(
    gams_multi, 
    img_combined, 
    target_classes = NULL,
    max_multiplier = 2.5
) {
  if (is.null(target_classes)) {
    target_classes <- get_target_classes(img_data = img_combined)
  }
  
  train_data <- img_combined %>% filter(dataset == "test_GAM_train")
  test_data  <- img_combined %>% filter(dataset == "test_GAM_test")
  
  has_region <- "region" %in% names(test_data)
  
  for (sp in target_classes) {
    pred_col  <- paste0("pred_synergy_", sp)
    yolo_col  <- paste0("yolo_", sp)
    cas_col   <- paste0("cascade_", sp)
    diff_col  <- paste0("diff_", sp)
    lyolo_col <- paste0("log_yolo_", sp)
    lcas_col  <- paste0("log_cascade_", sp)
    
    # 1. Clamp test predictors to training range limits to prevent extrapolation explosion
    feature_cols <- intersect(c(yolo_col, cas_col, diff_col, lyolo_col, lcas_col), names(train_data))
    for (fc in feature_cols) {
      tr_max <- max(train_data[[fc]], na.rm = TRUE)
      tr_min <- min(train_data[[fc]], na.rm = TRUE)
      test_data[[fc]] <- pmax(pmin(test_data[[fc]], tr_max), tr_min)
    }
    
    test_data[[pred_col]] <- NA_real_
    
    # 2. Predict counts
    if (has_region && !is.null(gams_multi[[sp]]$GB)) {
      gb_idx  <- which(test_data$region == "GB")
      mab_idx <- which(test_data$region == "MAB")
      
      if (length(gb_idx) > 0)  test_data[[pred_col]][gb_idx]  <- predict(gams_multi[[sp]]$GB, newdata = test_data[gb_idx, ], type = "response")
      if (length(mab_idx) > 0) test_data[[pred_col]][mab_idx] <- predict(gams_multi[[sp]]$MAB, newdata = test_data[mab_idx, ], type = "response")
    } else if (!is.null(gams_multi[[sp]]$ALL)) {
      test_data[[pred_col]] <- predict(gams_multi[[sp]]$ALL, newdata = test_data, type = "response")
    }
    
    # 3. Post-Processing Cap: Clamp predicted count to a logical ceiling relative to raw model detectors
    if (yolo_col %in% names(test_data) && cas_col %in% names(test_data)) {
      cap_limit <- pmax(test_data[[yolo_col]], test_data[[cas_col]]) * max_multiplier + 2
      test_data[[pred_col]] <- pmin(test_data[[pred_col]], cap_limit)
    }
  }
  
  # Calculate totals and metrics
  synergy_pred_cols <- intersect(paste0("pred_synergy_", target_classes), names(test_data))
  test_data$pred_synergy_total <- rowSums(test_data[, synergy_pred_cols, drop = FALSE], na.rm = TRUE)
  
  metrics <- purrr::map_dfr(target_classes, function(sp) {
    gt_col   <- paste0("n_", sp)
    pred_col <- paste0("pred_synergy_", sp)
    
    if (!gt_col %in% names(test_data) || !pred_col %in% names(test_data)) return(NULL)
    
    df <- if (has_region) test_data %>% group_by(region) else test_data %>% mutate(region = "ALL") %>% group_by(region)
    df %>%
      summarise(
        species      = sp,
        r2_synergy   = summary(lm(get(gt_col) ~ get(pred_col)))$adj.r.squared,
        rmse_synergy = sqrt(mean((get(gt_col) - get(pred_col))^2, na.rm = TRUE)),
        .groups      = "drop"
      )
  })
  
  tot_df <- if (has_region) test_data %>% group_by(region) else test_data %>% mutate(region = "ALL") %>% group_by(region)
  tot_metrics <- tot_df %>%
    summarise(
      species      = "TOTAL_ORGANISMS",
      r2_synergy   = summary(lm(n_annotations ~ pred_synergy_total))$adj.r.squared,
      rmse_synergy = sqrt(mean((n_annotations - pred_synergy_total)^2, na.rm = TRUE)),
      .groups      = "drop"
    )
  
  return(list(img_eval = test_data, metrics = bind_rows(metrics, tot_metrics)))
}