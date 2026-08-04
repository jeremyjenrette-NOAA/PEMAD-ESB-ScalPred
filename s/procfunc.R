# ======================================================================
# procfunc.R
# ======================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# ----------------------------------------------------------------------
# Helper: Autolocate Classes from Image and Detection Datasets
# ----------------------------------------------------------------------
get_target_classes <- function(img_data = NULL, det_data = NULL) {
  img_classes <- character(0)
  det_classes <- character(0)
  
  # 1. Extract potential class names from n_* count columns in img_data
  if (!is.null(img_data)) {
    img_classes <- sub("^n_", "", grep("^n_", names(img_data), value = TRUE))
  }
  
  # 2. Extract species labels from detection ground truth and predictions
  if (!is.null(det_data)) {
    det_classes <- c(det_data$gt_label, det_data$spname)
    det_classes <- unique(det_classes[!is.na(det_classes) & det_classes != ""])
  }
  
  # 3. Exclude non-species metadata labels
  reserved_words <- c("Background", "background", "Missed", "missed", "annotations", "auto", "total")
  img_classes <- setdiff(img_classes, reserved_words)
  det_classes <- setdiff(det_classes, reserved_words)
  
  # 4. Determine final species: Intersect if both datasets exist, otherwise filter
  if (length(img_classes) > 0 && length(det_classes) > 0) {
    # Species MUST exist in detection labels AND have a corresponding n_<class> count column
    target_classes <- intersect(det_classes, img_classes)
  } else if (length(det_classes) > 0) {
    target_classes <- det_classes
  } else {
    target_classes <- img_classes
  }
  
  # Exclude any remaining automated summary prefix columns (e.g., auto_asterias)
  target_classes <- target_classes[!grepl("^auto_", target_classes)]
  
  return(target_classes)
}

# ----------------------------------------------------------------------
# 1. Class-Specific Precision-Recall Curve Evaluator
# ----------------------------------------------------------------------
evaluate_pr_curve <- function(
    calib_df,
    img_df,
    conf_grid = seq(0, 1, by = 0.005), 
    model_id  = "model",
    target_class = NULL
) {
  # If target_class is missing, default to the first autodetected class
  if (is.null(target_class)) {
    target_class <- get_target_classes(img_df, calib_df)[1]
  }
  
  # Calculate total ground truth for this target class across all images
  gt_col <- paste0("n_", target_class)
  total_gt <- if (gt_col %in% names(img_df)) sum(img_df[[gt_col]], na.rm = TRUE) else 0
  
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
    target_classes = NULL,
    conf_grid = seq(0, 1, by = 0.005),
    stratify_region = FALSE
) {
  # Auto-detect target classes across all models if not explicitly passed
  if (is.null(target_classes)) {
    target_classes <- map(models_list, function(m) {
      get_target_classes(m$img, m$det)
    }) %>% unlist() %>% unique()
  }
  
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
    target_classes = NULL
) {
  # Auto-detect target classes if not specified
  if (is.null(target_classes)) {
    target_classes <- get_target_classes(img_data, det_data)
  }
  
  # Filter predictions above the selected threshold
  df_th <- det_data %>% filter(conf >= threshold)
  
  # 1. Matches & Misclassifications
  matches <- df_th %>%
    filter(truedetect == TRUE) %>%
    mutate(
      predicted = spname,
      actual = gt_label
    ) %>%
    count(actual, predicted)
  
  # 2. Background False Positives
  fps <- df_th %>%
    filter(truedetect == FALSE) %>%
    mutate(
      predicted = spname,
      actual = "Background"
    ) %>%
    count(actual, predicted)
  
  # 3. Missed Ground Truths (False Negatives)
  gt_cols <- paste0("n_", target_classes)
  valid_gt_cols <- intersect(gt_cols, names(img_data))
  
  total_manual <- img_data %>%
    summarise(across(all_of(valid_gt_cols), \(x) sum(x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "species", values_to = "total_gt") %>%
    mutate(actual = sub("^n_", "", species))
  
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
  
  # Combine into structured grid layout
  all_counts <- bind_rows(matches, fps, fns)
  
  grid <- expand_grid(
    actual = c(target_classes, "Background"),
    predicted = c(target_classes, "Missed")
  ) %>%
    filter(!(actual == "Background" & predicted == "Missed"))
  
  grid %>%
    left_join(all_counts, by = c("actual", "predicted")) %>%
    mutate(n = replace_na(n, 0L))
}