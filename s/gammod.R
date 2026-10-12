library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(mgcv)
library(forcats) 
library(latex2exp)
source("./gamfunc.R")

# Synergy formula selection: "cv" = block-fold CV error on test_GAM_train
# (robust to over-flexible surfaces); "aic" = previous behaviour.
SYN_SELECT     <- "cv"
SYN_CV_METRIC  <- "cv_rmse"
N_BOOT         <- 1000   # block-bootstrap replicates for holdout intervals

# 1. Load metadata and past evaluation results
meta <- read.csv("../data/raw/dataset_split_2226.csv") %>% 
  janitor::clean_names() %>%
  mutate(image_id = stringr::str_remove(imagename, "\\.[A-Za-z0-9]+$")) %>%
  distinct(image_id, .keep_all = TRUE)

yolo_eval <- readRDS("../data/processed/YOLOv12_deteval_2226.rds")
cas_eval  <- readRDS("../data/processed/CascadeR-CNN_deteval_2226.rds")

# 2. Build Synergistic Dataset & ensure log_ratio is defined for M16_PropDiff
img_combined <- build_synergy_dataset(yolo_eval, cas_eval, meta) %>%
  mutate(log_ratio = log((pred_yolo + 1) / (pred_cascade + 1)))

# Detection-structure features from YOLO's boxes (gamfunc.R section 9). YOLO is
# used as the image descriptor because it is the stronger GB detector and the
# CV-selected GB synergy model is YOLO-based.
crowd_feats  <- image_detection_features(yolo_eval$det_eval)
img_combined <- add_crowding_features(img_combined, crowd_feats)
cat("Crowding features (holdout + GAM-train images):\n")
print(img_combined %>% group_by(region) %>%
        summarise(across(all_of(CROWD_VARS), list(med = ~ median(.x), p90 = ~ quantile(.x, 0.9)))),
      width = Inf)

# ======================================================================
# 2. Define Systematic Synergistic Candidate Formulas
# ======================================================================
image_candidate_forms <- list(
  # STAGE 1: Single Model Raw Baselines
  M01_YoloRaw     = n_annotations ~ s(pred_yolo),
  M02_CascadeRaw  = n_annotations ~ s(pred_cascade),
  M03_MeanRaw     = n_annotations ~ s(pred_mean),
  M035_SumRaw     = n_annotations ~ s(pred_sum),
  
  # STAGE 2: Log-Transformed Single Predictors
  M04_YoloLog     = n_annotations ~ s(log_yolo_pred),
  M05_CascadeLog  = n_annotations ~ s(log_cascade_pred),
  M06_MeanLog     = n_annotations ~ s(log_pred_mean),
  
  # STAGE 3: Disagreement/Difference Covariate
  M07_YoloDiff    = n_annotations ~ s(log_yolo_pred, pred_diff),
  M08_CascadeDiff = n_annotations ~ s(log_cascade_pred, pred_diff),
  M09_MeanDiff    = n_annotations ~ s(log_pred_mean, pred_diff),
  M09_MeanCasDiff = n_annotations ~ s(log_pred_mean, pred_diff) + s(log_cascade_pred, pred_diff),
  M10_YoloLogDiff = n_annotations ~ s(log_yolo_pred, log_pred_diff),
  
  # STAGE 4: Multi-Model Synergy
  M11_Additive    = n_annotations ~ s(log_yolo_pred) + s(log_cascade_pred) + s(log_pred_diff),
  M12_Tensor      = n_annotations ~ te(log_yolo_pred, log_cascade_pred),
  M13_TensorDiff  = n_annotations ~ ti(log_yolo_pred, log_cascade_pred) + log_pred_diff,
  M14_Champion    = n_annotations ~ s(log_yolo_pred, log_pred_diff) + s(log_cascade_pred, log_pred_diff),
  # M15_TensorDiff2 = n_annotations ~ ti(log_yolo_pred, log_cascade_pred) + 
  #   ti(log_cascade_pred, log_pred_diff, k = c(6,6)) + 
  #   ti(log_yolo_pred, log_pred_diff, k = c(6,6)),
  
  # STAGE 5: Advanced ANOVA & Ratio Models
  M15b_FullDecomp = n_annotations ~ s(log_yolo_pred, k = 8) + 
    s(log_cascade_pred, k = 8) + 
    s(log_pred_diff, k = 8) + 
    ti(log_yolo_pred, log_cascade_pred, k = c(6, 6)) + 
    ti(log_cascade_pred, log_pred_diff, k = c(6, 6)) + 
    ti(log_yolo_pred, log_pred_diff, k = c(6, 6)),
  
  M16_PropDiff    = n_annotations ~ s(log_yolo_pred, k = 8) + 
    s(log_cascade_pred, k = 8) + 
    s(log_ratio, k = 8) + 
    ti(log_yolo_pred, log_cascade_pred, k = c(6, 6)) + 
    ti(log_yolo_pred, log_ratio, k = c(6, 6)) + 
    ti(log_cascade_pred, log_ratio, k = c(6, 6)),
  
  M17_Tensor3D    = n_annotations ~ te(log_yolo_pred, log_cascade_pred, log_pred_diff, k = c(5, 5, 5)),
  
  # STAGE 6: Detection-structure ("crowding") terms, added to the two CV
  # winners (M04 GB, M12 MAB). k kept small: few holdout images exceed 100
  # scallops, nearly all from cruise 202404.
  C01_Overlap     = n_annotations ~ s(log_yolo_pred) + s(overlap_frac, k = 5),
  C02_Spacing     = n_annotations ~ s(log_yolo_pred) + s(nn_ratio, k = 5),
  C03_LowConf     = n_annotations ~ s(log_yolo_pred) + s(lowconf_frac, k = 5),
  C04_BoxSize     = n_annotations ~ s(log_yolo_pred) + s(med_box_mm, k = 5),
  C05_teOverlap   = n_annotations ~ te(log_yolo_pred, overlap_frac, k = c(5, 4)),
  C06_teBoxSize   = n_annotations ~ te(log_yolo_pred, med_box_mm, k = c(5, 4)),
  C07_teLowConf   = n_annotations ~ te(log_yolo_pred, lowconf_frac, k = c(5, 4)),
  C08_Structure   = n_annotations ~ s(log_yolo_pred) + s(med_box_mm, k = 5) + s(lowconf_frac, k = 5) +
                      s(overlap_frac, k = 5) + s(nn_ratio, k = 5),
  C09_TensorBox   = n_annotations ~ te(log_yolo_pred, log_cascade_pred) + s(med_box_mm, k = 5) +
                      s(lowconf_frac, k = 5),
  C10_TensorCrowd = n_annotations ~ te(log_yolo_pred, log_cascade_pred) + s(overlap_frac, k = 5) +
                      s(nn_ratio, k = 5)
)

# ======================================================================
# Helper Function: Automated Formula Selection
# ======================================================================
select_best_synergy_model <- function(eval_table, candidate_list, metric = "AIC") {
  # Sort table by designated metric
  sorted_table <- eval_table %>% arrange(.data[[metric]])
  
  best_id   <- sorted_table$Model_ID[1]
  best_score <- round(sorted_table[[metric]][1], 2)
  best_r2    <- round(sorted_table$R2_Train[1], 3)
  
  cat(sprintf("  [Selected] %s (by %s = %g | R² = %g)\n", best_id, metric, best_score, best_r2))
  
  return(candidate_list[[best_id]])
}

# ======================================================================
# Run Evaluation and Dynamically Assign Best Models
# ======================================================================
cat("--- Synergy Selection: Georges Bank ---\n")
selection_gb <- compare_image_gams(img_combined, image_candidate_forms, "GB")
print(selection_gb, n = Inf)

best_syn_gb <- select_best_synergy_model(selection_gb, image_candidate_forms, metric = "AIC")
if (SYN_SELECT == "cv") {
  cat("--- CV synergy selection: Georges Bank ---\n")
  cv_gb <- cv_select_synergy(img_combined, image_candidate_forms, "GB", metric = SYN_CV_METRIC)
  print(cv_gb, n = Inf)
  cat("  [Selected by CV]", cv_gb$Model_ID[1], "\n")
  best_syn_gb <- image_candidate_forms[[cv_gb$Model_ID[1]]]
}

cat("\n--- Synergy Selection: Mid-Atlantic Bight ---\n")
selection_mab <- compare_image_gams(img_combined, image_candidate_forms, "MAB")
print(selection_mab, n = Inf)

best_syn_mab <- select_best_synergy_model(selection_mab, image_candidate_forms, metric = "AIC")
if (SYN_SELECT == "cv") {
  cat("--- CV synergy selection: Mid-Atlantic Bight ---\n")
  cv_mab <- cv_select_synergy(img_combined, image_candidate_forms, "MAB", metric = SYN_CV_METRIC)
  print(cv_mab, n = Inf)
  cat("  [Selected by CV]", cv_mab$Model_ID[1], "\n")
  best_syn_mab <- image_candidate_forms[[cv_mab$Model_ID[1]]]
  write.csv(bind_rows(cv_gb %>% mutate(region = "GB"), cv_mab %>% mutate(region = "MAB")),
            "../data/processed/synergy_cv_selection_2226.csv", row.names = FALSE)
}

# best_syn_gb <- image_candidate_forms$M17_Tensor3D
# best_syn_mab <- image_candidate_forms$M07_YoloDiff
# ======================================================================
# 3. Train and Test Final Synergistic Model
# ======================================================================
syn_gams <- fit_synergy_gams(img_combined, best_syn_gb, best_syn_mab)
syn_eval <- test_synergy_gams(syn_gams, img_combined)

saveRDS(syn_gams, file = paste0("../data/processed/Syngams_2226.rds"))

# ----------------------------------------------------------------------
# Leakage diagnostic: refit the SAME synergy formulas on the old in-sample
# detection-GAM sums and score both on test_GAM_test. The test inputs are
# identical in both versions; only the stage-2 training inputs differ.
# ----------------------------------------------------------------------
img_insample     <- insample_synergy_dataset(img_combined)
syn_gams_insamp  <- fit_synergy_gams(img_insample, best_syn_gb, best_syn_mab)
syn_eval_insamp  <- test_synergy_gams(syn_gams_insamp, img_insample)

leakage_synergy <- bind_rows(
  syn_eval$metrics        %>% mutate(stage1_inputs = "out-of-fold (cross-fit)"),
  syn_eval_insamp$metrics %>% mutate(stage1_inputs = "in-sample (previous)")
) %>% arrange(region, stage1_inputs)
cat("\n--- Synergy GAM holdout metrics: cross-fit vs in-sample stage-1 inputs ---\n")
print(leakage_synergy)
write.csv(leakage_synergy, "../data/processed/synergy_leakage_check_2226.csv", row.names = FALSE)

# ----------------------------------------------------------------------
# Crowding-term report (for the paper; NOT used for selection, which used CV
# on test_GAM_train only). Best CV formula with vs without detection-structure
# terms, both scored on the holdout overall and on dense images.
# ----------------------------------------------------------------------
DENSE_MIN <- 50   # annotations per image counted as "dense" in this report
if (SYN_SELECT == "cv") {
  best_id <- function(cv, crowd) {
    ids <- cv$Model_ID[!is.na(cv$cv_rmse)]           # cv tables are sorted by CV error
    ids <- if (crowd) ids[grepl("^C[0-9]", ids)] else ids[!grepl("^C[0-9]", ids)]
    ids[1]
  }
  crowd_report <- purrr::map_dfr(c(FALSE, TRUE), function(cr) {
    id_gb <- best_id(cv_gb, cr); id_mab <- best_id(cv_mab, cr)
    ev <- test_synergy_gams(fit_synergy_gams(img_combined, image_candidate_forms[[id_gb]],
                                             image_candidate_forms[[id_mab]]), img_combined)$img_eval
    ev %>%
      mutate(err = pred_synergy - n_annotations, dense = n_annotations >= DENSE_MIN) %>%
      group_by(region) %>%
      summarise(variant = if (cr) "best with crowding terms" else "best without crowding terms",
                formula = if (first(region) == "GB") id_gb else id_mab,
                cv_rmse = if (first(region) == "GB") cv_gb$cv_rmse[cv_gb$Model_ID == id_gb] else cv_mab$cv_rmse[cv_mab$Model_ID == id_mab],
                holdout_rmse = sqrt(mean(err^2)),
                holdout_pct_bias = 100 * sum(err) / sum(n_annotations),
                dense_n = sum(dense),
                dense_rmse = sqrt(mean(err[dense]^2)),
                dense_pct_bias = 100 * sum(err[dense]) / sum(n_annotations[dense]),
                .groups = "drop")
  }) %>% arrange(region, variant)
  cat("\n--- Synergy with vs without crowding terms (holdout; dense = >=", DENSE_MIN, "annotations) ---\n")
  print(crowd_report, width = Inf)
  write.csv(crowd_report, "../data/processed/synergy_crowding_report_2226.csv", row.names = FALSE)
}

# Where does any difference come from? RMSE by density bin (holdout only).
leakage_by_bin <- bind_rows(
  syn_eval$img_eval        %>% mutate(stage1_inputs = "out-of-fold (cross-fit)"),
  syn_eval_insamp$img_eval %>% mutate(stage1_inputs = "in-sample (previous)")
) %>%
  group_by(region, density_bin, stage1_inputs) %>%
  summarise(rmse = sqrt(mean((n_annotations - pred_synergy)^2, na.rm = TRUE)),
            mean_bias = mean(pred_synergy - n_annotations, na.rm = TRUE),
            n_images = n(), .groups = "drop") %>%
  arrange(region, density_bin, stage1_inputs)
print(leakage_by_bin, n = Inf)
write.csv(leakage_by_bin, "../data/processed/synergy_leakage_by_bin_2226.csv", row.names = FALSE)

# ----------------------------------------------------------------------
# Block-bootstrap 95% intervals on holdout metrics for every method
# (resampling whole transect blocks). d_rmse_* = RMSE minus synergy RMSE;
# an interval excluding 0 means the difference is distinguishable.
# ----------------------------------------------------------------------
boot_df <- syn_eval$img_eval %>%
  left_join(yolo_eval$img_eval %>% select(image_id, yolo_f1 = predicted_f1_number), by = "image_id") %>%
  left_join(cas_eval$img_eval  %>% select(image_id, cascade_f1 = predicted_f1_number), by = "image_id")
boot_metrics <- block_bootstrap_metrics(
  boot_df,
  pred_cols = c("Synergistic GAM" = "pred_synergy",
                "YOLOv12 det-GAM" = "pred_yolo", "Cascade det-GAM" = "pred_cascade",
                "YOLOv12 F1" = "yolo_f1", "Cascade F1" = "cascade_f1"),
  ref = "Synergistic GAM", B = N_BOOT
)
print(boot_metrics, n = Inf, width = Inf)
write.csv(boot_metrics, "../data/processed/holdout_bootstrap_metrics_2226.csv", row.names = FALSE)

# Holdout table with every method's prediction, for error.R
saveRDS(boot_df, "../data/processed/synergy_holdout_eval_2226.rds")
# ======================================================================
# 4. The 12-Panel Plot Generation
# ======================================================================

# A. Helper function with axis controls and zoom capabilities
make_1to1_plot <- function(data, x_col, region_name, plot_title, r2, rmse, pt_color, 
                           show_x = FALSE, show_y = FALSE, zoom_limit = NULL) {
  
  val_r2 <- round(r2, 3)
  val_rmse <- round(rmse, 3)
  label_text <- sprintf("R²(1:1) = %.3f\nRMSE = %.3f", val_r2, val_rmse)
  
  p <- ggplot(data %>% filter(region == region_name), aes(x = .data[[x_col]], y = n_annotations)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(color = pt_color, size = 1.2, alpha = 0.5) +
    geom_text(
      x = -Inf, y = Inf, label = label_text, 
      hjust = -0.2, vjust = 1.5, size = 3.5, fontface = "plain"
    ) +
    theme_bw(base_size = 12) +
    labs(title = plot_title, x = "Predicted Count", y = "True Count") +
    theme(plot.title = element_text(face = "plain", size = 11, hjust = 0.5))
  
  # Conditionally remove axis titles to reduce clutter
  if (!show_x) p <- p + theme(axis.title.x = element_blank())
  if (!show_y) p <- p + theme(axis.title.y = element_blank())
  
  # Apply zoom if a limit is provided
  if (!is.null(zoom_limit)) {
    p <- p + coord_cartesian(xlim = c(0, zoom_limit), ylim = c(0, zoom_limit))
  }
  
  return(p)
}

# B. Prepare data and calculate dynamics limits/metadata
y_test <- yolo_eval$img_eval %>% filter(dataset == "test_GAM_test")
ym <- yolo_eval$metrics

c_test <- cas_eval$img_eval %>% filter(dataset == "test_GAM_test")
cm <- cas_eval$metrics

s_test <- syn_eval$img_eval
sm <- syn_eval$metrics

# Dynamically calculate the 95th percentile limit for the zoomed plots (ensure minimum limit of 5)
limit_gb <- max(5, quantile(y_test$n_annotations[y_test$region == "GB"], 0.95, na.rm = TRUE))
limit_mab <- max(5, quantile(y_test$n_annotations[y_test$region == "MAB"], 0.95, na.rm = TRUE))

# Calculate n for subcaption
n_gb <- nrow(y_test %>% filter(region == "GB"))
n_mab <- nrow(y_test %>% filter(region == "MAB"))

n_gb_abundance <- sum((y_test %>% filter(region == "GB"))$n_annotations)
n_mab_abundance <- sum((y_test %>% filter(region == "MAB"))$n_annotations)

# C. Build the 12 plots sequentially
p <- list()

# -- Row 1: YOLOv12 (Red Points) --
p[[1]] <- make_1to1_plot(y_test, "predicted_number", "GB", expression("YOLOv12" ~ Sigma ~ "P(detection)"), ym$r2_1to1_calibrated[ym$region=="GB"], ym$rmse_calibrated[ym$region=="GB"], "red", show_x = FALSE, show_y = TRUE)
p[[2]] <- make_1to1_plot(y_test, "predicted_f1_number", "GB", expression("YOLOv12 F"[1] ~ "Threshold"), ym$r2_1to1_f1[ym$region=="GB"], ym$rmse_f1[ym$region=="GB"], "red", show_x = FALSE, show_y = FALSE)
p[[3]] <- make_1to1_plot(y_test, "predicted_number", "MAB", expression("YOLOv12" ~ Sigma ~ "P(detection)"), ym$r2_1to1_calibrated[ym$region=="MAB"], ym$rmse_calibrated[ym$region=="MAB"], "red", show_x = FALSE, show_y = FALSE)
p[[4]] <- make_1to1_plot(y_test, "predicted_f1_number", "MAB", expression("YOLOv12 F"[1] ~ "Threshold"), ym$r2_1to1_f1[ym$region=="MAB"], ym$rmse_f1[ym$region=="MAB"], "red", show_x = FALSE, show_y = FALSE)

# -- Row 2: Cascade R-CNN (Blue Points) --
p[[5]] <- make_1to1_plot(c_test, "predicted_number", "GB", expression("Cascade R-CNN" ~ Sigma ~ "P(detection)"), cm$r2_1to1_calibrated[cm$region=="GB"], cm$rmse_calibrated[cm$region=="GB"], "#2C7FB8", show_x = FALSE, show_y = TRUE)
p[[6]] <- make_1to1_plot(c_test, "predicted_f1_number", "GB", expression("Cascade R-CNN F"[1] ~ "Threshold"), cm$r2_1to1_f1[cm$region=="GB"], cm$rmse_f1[cm$region=="GB"], "#2C7FB8", show_x = FALSE, show_y = FALSE)
p[[7]] <- make_1to1_plot(c_test, "predicted_number", "MAB", expression("Cascade R-CNN" ~ Sigma ~ "P(detection)"), cm$r2_1to1_calibrated[cm$region=="MAB"], cm$rmse_calibrated[cm$region=="MAB"], "#2C7FB8", show_x = FALSE, show_y = FALSE)
p[[8]] <- make_1to1_plot(c_test, "predicted_f1_number", "MAB", expression("Cascade R-CNN F"[1] ~ "Threshold"), cm$r2_1to1_f1[cm$region=="MAB"], cm$rmse_f1[cm$region=="MAB"], "#2C7FB8", show_x = FALSE, show_y = FALSE)

# -- Row 3: Synergy Image-GAM (Purple Points, Full and Zoomed) --
p[[9]]  <- make_1to1_plot(s_test, "pred_synergy", "GB", "Synergistic GAM", sm$r2_1to1_synergy[sm$region=="GB"], sm$rmse_synergy[sm$region=="GB"], "#7570B3", show_x = TRUE, show_y = TRUE)
p[[10]] <- make_1to1_plot(s_test, "pred_synergy", "GB", "Synergistic GAM (95% Data)", sm$r2_1to1_synergy[sm$region=="GB"], sm$rmse_synergy[sm$region=="GB"], "#7570B3", show_x = TRUE, show_y = FALSE, zoom_limit = limit_gb)
p[[11]] <- make_1to1_plot(s_test, "pred_synergy", "MAB", "Synergistic GAM", sm$r2_1to1_synergy[sm$region=="MAB"], sm$rmse_synergy[sm$region=="MAB"], "#7570B3", show_x = TRUE, show_y = FALSE)
p[[12]] <- make_1to1_plot(s_test, "pred_synergy", "MAB", "Synergistic GAM (95% Data)", sm$r2_1to1_synergy[sm$region=="MAB"], sm$rmse_synergy[sm$region=="MAB"], "#7570B3", show_x = TRUE, show_y = FALSE, zoom_limit = limit_mab)


# ======================================================================
# 5. Assemble and Export the 12-Panel Masterpiece
# ======================================================================

# Create Overarching Headers using Grid graphics wrapped for patchwork
gb_header  <- wrap_elements(textGrob("Georges Bank", gp = gpar(fontsize = 16, fontface = "bold")))
mab_header <- wrap_elements(textGrob("Mid-Atlantic Bight", gp = gpar(fontsize = 16, fontface = "bold")))
header_row <- gb_header | mab_header

# Assemble the rows (4 columns wide)
row1 <- p[[1]] | p[[2]] | p[[3]] | p[[4]]
row2 <- p[[5]] | p[[6]] | p[[7]] | p[[8]]
row3 <- p[[9]] | p[[10]]| p[[11]]| p[[12]]

# Create dynamic subcaption string
caption_string <- sprintf("Datasets: 2022-2024, 2026 | Models Evaluated: YOLOv12, Cascade R-CNN, Synergistic GAM\nHoldout Test Images: Georges Bank (n = %d), Mid-Atlantic Bight (n = %d)\nScallop Abundance: Georges Bank (n = %d), Mid-Atlantic Bight (n = %d)", n_gb, n_mab, n_gb_abundance, n_mab_abundance)

# Stitch it all together: Header row on top, followed by the 3 data rows
final_plot <- (header_row / row1 / row2 / row3) + 
  plot_layout(heights = c(0.1, 1, 1, 1)) + # Gives the header a slim vertical footprint
  plot_annotation(
    title = bquote(bold("Comparative Assessment: Optimal F"[1] ~ "vs. Detection GAMs vs. Synergistic Calibration")),
    caption = caption_string,
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.caption = element_text(size = 12, face = "italic", color = "grey30", hjust = 0.5, margin = margin(t = 15))
    )
  )

print(final_plot)

# Save the high-resolution figure
ggsave(
  filename = "../figures/ms_figures/2226_12panel_calibration_summary.png",
  plot = final_plot,
  width = 13,
  height = 12,
  device = "png",
  bg = "white"
)
ggsave(
  filename = "~/saltnoaa/presentations/figures/2226_12panel_calibration_summary.png",
  plot = final_plot,
  width = 17,
  height = 12,
  device = "png",
  bg = "white"
)


###########################################################################

# ======================================================================
# 6. Final Comprehensive Abundance Summary Plot
# ======================================================================

# Ensure we are strictly using the holdout test set
true_val <- sum(s_test$n_annotations, na.rm = TRUE)

# Build the comparison dataframe
sum_df <- data.frame(
  Method = factor(
    c("True Abundance", 
      "YOLOv12 F1", "YOLOv12 Det-GAM", 
      "Cascade F1", "Cascade Det-GAM", 
      "Synergistic GAM"),
    levels = c("True Abundance", 
               "YOLOv12 F1", "YOLOv12 Det-GAM", 
               "Cascade F1", "Cascade Det-GAM", 
               "Synergistic GAM")
  ),
  Value = c(
    true_val,
    sum(y_test$predicted_f1_number, na.rm = TRUE),
    sum(y_test$predicted_number, na.rm = TRUE),
    sum(c_test$predicted_f1_number, na.rm = TRUE),
    sum(c_test$predicted_number, na.rm = TRUE),
    sum(s_test$pred_synergy, na.rm = TRUE)
  )
)

# 1. Define the Labels in the EXACT same order as the factor levels
label_map <- c(
  "True Abundance"  = "True~Abundance",
  "YOLOv12 F1"      = "YOLO~F[1]",
  "YOLOv12 Det-GAM" = "YOLO~Sigma*p[det]", # Moved to 3rd to match factor levels
  "Cascade F1"      = "Cascade~F[1]",       # Moved to 4th to match factor levels
  "Cascade Det-GAM" = "Cascade~Sigma*p[det]",
  "Synergistic GAM" = "Synergistic~GAM"
)

sum_df <- sum_df %>%
  mutate(
    Error_Pct = ((Value - true_val) / true_val) * 100,
    Label_Text = case_when(
      Method == "True Abundance" ~ as.character(round(Value, 0)),
      TRUE ~ sprintf("%d\n(%+.1f%%)", round(Value, 0), Error_Pct)
    )
  )

# Update global text geom defaults once before plotting
ggplot2::update_geom_defaults("text", list(family = "sans"))

# 2. Generate the plot with parsed expressions directly
p_sum_final <- ggplot(sum_df, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.7, color = "black") +
  geom_hline(yintercept = true_val, linetype = "dashed", color = "black", linewidth = 1.0) +
  geom_text(aes(label = Label_Text), y = 10000, 
            vjust = 0, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c(
    "True Abundance"  = "grey40", 
    "YOLOv12 F1"      = "#FB6A4A", 
    "Cascade F1"      = "#6BAED6", 
    "YOLOv12 Det-GAM" = "#CB181D", 
    "Cascade Det-GAM" = "#2171B5", 
    "Synergistic GAM" = "#7570B3"
  )) +
  
  # Parse the properly ordered label map directly here
  scale_x_discrete(labels = parse(text = label_map)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) + 
  theme_minimal(base_size = 14) +
  labs(title = "Validated Abundance Estimates", x = "", y = "Total Sum Count") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11, face = "plain"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
p_sum_final
# Save the final render
ggsave("../figures/ms_figures/2226_final_abundance_summary.png",
       plot = p_sum_final, width = 8, height = 6, dpi = 300)