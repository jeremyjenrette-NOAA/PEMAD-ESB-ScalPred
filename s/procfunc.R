# ======================================================================
# procfunc.R
# ======================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# ----------------------------------------------------------------------
# 1. Class-Specific Precision-Recall Curve Evaluator
# ----------------------------------------------------------------------
evaluate_pr_curve <- function(
    calib_df,
    img_df,
    conf_grid = seq(0, 1, by = 0.005), 
    model_id  = "model",
    target_class = "jonah_crab"
) {
  # Calculate total ground truth for this target class across all images
  gt_col <- paste0("n_", target_class)
  total_gt <- sum(img_df[[gt_col]], na.rm = TRUE)
  
  map_dfr(conf_grid, function(th) {
    
    # Isolate predictions above confidence threshold for the target class
    df_th <- calib_df %>% filter(conf >= th, spname == target_class)
    
    # TP: Predicted as target_class, matched to GT, and actual label matches
    TP <- sum(df_th$truedetect == TRUE & df_th$gt_label == target_class, na.rm = TRUE)
    
    # FP: Predicted as target_class but matched to background OR matched to another species
    FP <- sum(df_th$truedetect == FALSE | (df_th$truedetect == TRUE & df_th$gt_label != target_class), na.rm = TRUE)
    
    # FN: Ground truth targets of this species missed by predictions of this species
    FN <- max(0L, total_gt - TP)
    
    precision <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
    recall    <- if (total_gt == 0) NA_real_ else TP / total_gt
    f1        <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_
    else 2 * precision * recall / (precision + recall)
    
    tibble(
      model     = model_id,
      species   = target_class,
      conf      = th,
      TP        = TP,
      FP        = FP,
      FN        = FN,
      precision = precision,
      recall    = recall,
      f1        = f1
    )
  })
}

# ----------------------------------------------------------------------
# 2. Average Precision (mAP) Calculator via Trapezoidal AUC
# ----------------------------------------------------------------------
calculate_ap <- function(recall, precision) {
  valid <- !is.na(recall) & !is.na(precision)
  r <- recall[valid]
  p <- precision[valid]
  
  if (length(r) < 2) return(NA_real_)
  
  # Order by recall ascending to calculate area correctly
  ord <- order(r)
  r <- r[ord]
  p <- p[ord]
  
  sum(diff(r) * (p[-1] + p[-length(p)]) / 2)
}

# ----------------------------------------------------------------------
# 3. Model Wrapper (Stratification and Plot Generation)
# ----------------------------------------------------------------------
evaluate_pr_models <- function(
    models_list, 
    conf_grid = seq(0, 1, by = 0.005),
    stratify_region = FALSE
) {
  target_classes <- c("jonah_crab", "rock_crab", "cancer_sp")
  
  pr_all <- map_dfr(names(models_list), function(mod_name) {
    img_data <- models_list[[mod_name]]$img
    det_data <- models_list[[mod_name]]$det
    
    map_dfr(target_classes, function(sp) {
      if (stratify_region) {
        # Georges Bank Subset
        gb_img <- img_data %>% filter(region == "GB")
        gb_det <- det_data %>% filter(region == "GB")
        pr_gb <- evaluate_pr_curve(gb_det, gb_img, conf_grid, paste0(mod_name, " GB"), sp)
        
        # Mid-Atlantic Bight Subset
        mab_img <- img_data %>% filter(region == "MAB")
        mab_det <- det_data %>% filter(region == "MAB")
        pr_mab <- evaluate_pr_curve(mab_det, mab_img, conf_grid, paste0(mod_name, " MAB"), sp)
        
        return(bind_rows(pr_gb, pr_mab))
      } else {
        return(evaluate_pr_curve(det_data, img_data, conf_grid, mod_name, sp))
      }
    })
  })
  
  # Multi-class PR Plot
  legend_title <- if (stratify_region) "Model Region" else "Model"
  p_pr <- ggplot(pr_all, aes(x = recall, y = precision, color = species, linetype = model)) +
    geom_path(linewidth = 1.1, na.rm = TRUE) +
    theme_minimal() +
    xlim(0.05, 1) +
    labs(
      title = "Precision-Recall",
      x = "Recall",
      y = "Precision",
      color = "Species",
      linetype = legend_title
    )
  
  # Calculate mAP per model and species
  map_summary <- pr_all %>%
    group_by(model, species) %>%
    summarize(Average_Precision = calculate_ap(recall, precision), .groups = "drop")
  
  list(pr_all = pr_all, p_pr = p_pr, map_summary = map_summary)
}

# ----------------------------------------------------------------------
# 4. Multi-Class Confusion Matrix Generator
# ----------------------------------------------------------------------
generate_confusion_matrix <- function(
    det_data, 
    img_data, 
    threshold, 
    target_classes = c("jonah_crab", "rock_crab", "cancer_sp")
) {
  # Filter predictions above the selected threshold
  df_th <- det_data %>% filter(conf >= threshold)
  
  # 1. Matches & Misclassifications (True Positives and Cross-Class Errors)
  matches <- df_th %>%
    filter(truedetect == TRUE) %>%
    mutate(
      predicted = spname,
      actual = gt_label
    ) %>%
    count(actual, predicted)
  
  # 2. Background False Positives (Predicted as a species, but actually Background)
  fps <- df_th %>%
    filter(truedetect == FALSE) %>%
    mutate(
      predicted = spname,
      actual = "Background"
    ) %>%
    count(actual, predicted)
  
  # 3. Missed Ground Truths (False Negatives)
  total_manual <- img_data %>%
    summarise(across(all_of(paste0("n_", target_classes)), sum, na.rm = TRUE)) %>%
    pivot_longer(everything(), names_to = "species", values_to = "total_gt") %>%
    mutate(actual = str_remove(species, "^n_"))
  
  # Count unique GT boxes successfully matched above threshold
  detected_counts <- df_th %>%
    filter(truedetect == TRUE) %>%
    count(gt_label, name = "detected") %>%
    mutate(actual = gt_label)
  
  fns <- total_manual %>%
    left_join(detected_counts, by = "actual") %>%
    mutate(
      detected = replace_na(detected, 0L),
      n = pmax(0L, total_gt - detected),
      predicted = "Missed"
    ) %>%
    select(actual, predicted, n)
  
  # Combine into a structured grid layout
  all_counts <- bind_rows(matches, fps, fns)
  
  grid <- expand_grid(
    actual = c(target_classes, "Background"),
    predicted = c(target_classes, "Missed")
  ) %>%
    filter(!(actual == "Background" & predicted == "Missed")) # Background cannot be missed
  
  grid %>%
    left_join(all_counts, by = c("actual", "predicted")) %>%
    mutate(n = replace_na(n, 0L))
}