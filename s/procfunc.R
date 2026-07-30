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
    target_class = "adult"
) {
  target_class <- tolower(target_class)
  gt_col <- paste0("n_", target_class)
  
  # Determine total ground truth count
  if (gt_col %in% names(img_df)) {
    total_gt <- sum(img_df[[gt_col]], na.rm = TRUE)
  } else {
    # Fallback: Count total ground truths present in detection file
    total_gt <- sum(tolower(calib_df$gt_label) == target_class, na.rm = TRUE)
  }
  
  map_dfr(conf_grid, function(th) {
    
    # Isolate predictions above confidence threshold for the target class
    df_th <- calib_df %>% 
      filter(conf >= th, tolower(spname) == target_class)
    
    # TP: Predicted as target_class, true detection, actual label matches
    TP <- sum(df_th$truedetect == TRUE & tolower(df_th$gt_label) == target_class, na.rm = TRUE)
    
    # FP: Predicted as target_class but background OR matched to wrong species
    FP <- sum(df_th$truedetect == FALSE | (df_th$truedetect == TRUE & tolower(df_th$gt_label) != target_class), na.rm = TRUE)
    
    # FN: Ground truth targets missed
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
# 2. Average Precision (AP) Calculator via Trapezoidal AUC
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
# 3. Model Wrapper (Plot Generation & mAP Calculation)
# ----------------------------------------------------------------------
evaluate_pr_models <- function(
    models_list, 
    conf_grid = seq(0, 1, by = 0.005),
    target_classes = c("adult", "pup"),
    stratify_region = FALSE
) {
  pr_all <- map_dfr(names(models_list), function(mod_name) {
    img_data <- models_list[[mod_name]]$img
    det_data <- models_list[[mod_name]]$det
    
    map_dfr(target_classes, function(sp) {
      if (stratify_region && "region" %in% names(img_data)) {
        gb_img <- img_data %>% filter(region == "GB")
        gb_det <- det_data %>% filter(region == "GB")
        pr_gb <- evaluate_pr_curve(gb_det, gb_img, conf_grid, paste0(mod_name, " GB"), sp)
        
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
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      title = "Seal Detection: Precision-Recall Curve",
      x = "Recall",
      y = "Precision",
      color = "Class",
      linetype = legend_title
    )
  
  # Calculate AP per model and species
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
    target_classes = c("adult", "pup")
) {
  target_classes <- tolower(target_classes)
  df_th <- det_data %>% filter(conf >= threshold)
  
  # 1. True Positives and Cross-Class Errors
  matches <- df_th %>%
    filter(truedetect == TRUE) %>%
    mutate(
      predicted = tolower(spname),
      actual    = tolower(gt_label)
    ) %>%
    count(actual, predicted)
  
  # 2. Background False Positives
  fps <- df_th %>%
    filter(truedetect == FALSE) %>%
    mutate(
      predicted = tolower(spname),
      actual    = "background"
    ) %>%
    count(actual, predicted)
  
  # 3. Missed Ground Truths (False Negatives)
  gt_cols <- paste0("n_", target_classes)
  
  if (all(gt_cols %in% names(img_data))) {
    total_manual <- img_data %>%
      summarise(across(all_of(gt_cols), \(x) sum(x, na.rm = TRUE))) %>%
      pivot_longer(everything(), names_to = "species", values_to = "total_gt") %>%
      mutate(actual = sub("^n_", "", species))
  } else {
    total_manual <- det_data %>%
      filter(tolower(gt_label) %in% target_classes) %>%
      count(gt_label, name = "total_gt") %>%
      mutate(actual = tolower(gt_label))
  }
  
  detected_counts <- df_th %>%
    filter(truedetect == TRUE) %>%
    count(gt_label, name = "detected") %>%
    mutate(actual = tolower(gt_label))
  
  fns <- total_manual %>%
    left_join(detected_counts, by = "actual") %>%
    mutate(
      detected = replace_na(detected, 0L),
      n = pmax(0L, total_gt - detected),
      predicted = "missed"
    ) %>%
    select(actual, predicted, n)
  
  all_counts <- bind_rows(matches, fps, fns)
  
  grid <- expand_grid(
    actual = c(target_classes, "background"),
    predicted = c(target_classes, "missed")
  ) %>%
    filter(!(actual == "background" & predicted == "missed"))
  
  grid %>%
    left_join(all_counts, by = c("actual", "predicted")) %>%
    mutate(n = replace_na(n, 0L))
}