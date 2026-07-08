library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(mgcv)
source("./gamfunc.R")

# 1. Load your master metadata and previously evaluated model results
meta <- read.csv("../data/raw/dataset_split_2226.csv") %>% 
  janitor::clean_names() %>%
  mutate(image_id = stringr::str_remove(imagename, "\\.[A-Za-z0-9]+$")) %>%
  distinct(image_id, .keep_all = TRUE)

# (Assuming yolo_eval and cas_eval were saved from your gamstrat.R output)
yolo_eval <- readRDS("../data/processed/YOLOv12_deteval_2226.rds")
cas_eval  <- readRDS("../data/processed/CascadeR-CNN_deteval_2226.rds")

# 2. Build Synergistic Dataset
img_combined <- build_synergy_dataset(yolo_eval, cas_eval, meta) 
table(img_combined$dataset)

# ======================================================================
# 2. Define Systematic Synergistic Candidate Formulas (Regionally Stratified)
# ======================================================================
image_candidate_forms <- list(
  # STAGE 1: Single Model Raw Baselines
  M01_YoloRaw     = n_annotations ~ s(pred_yolo),
  M02_CascadeRaw  = n_annotations ~ s(pred_cascade),
  M03_MeanRaw     = n_annotations ~ s(pred_mean),
  
  # STAGE 2: Log-Transformed Single Predictors
  M04_YoloLog     = n_annotations ~ s(log_yolo_pred),
  M05_CascadeLog  = n_annotations ~ s(log_cascade_pred),
  M06_MeanLog     = n_annotations ~ s(log_pred_mean),
  
  # STAGE 3: Adding the Disagreement/Difference Covariate
  M07_YoloDiff    = n_annotations ~ s(log_yolo_pred, pred_diff),
  M08_CascadeDiff = n_annotations ~ s(log_cascade_pred, pred_diff),
  M09_MeanDiff    = n_annotations ~ s(log_pred_mean, pred_diff),
  M09_MeanCasDiff    = n_annotations ~ s(log_pred_mean, pred_diff) + s(log_cascade_pred, pred_diff),
  M10_YoloLogDiff = n_annotations ~ s(log_yolo_pred, log_pred_diff),
  
  
  # STAGE 4: Multi-Model Synergy
  M11_Additive    = n_annotations ~ s(log_yolo_pred) + s(log_cascade_pred) + s(log_pred_diff),
  M12_Tensor      = n_annotations ~ te(log_yolo_pred, log_cascade_pred),
  M13_TensorDiff  = n_annotations ~ ti(log_yolo_pred, log_cascade_pred) + log_pred_diff,
  M14_Champion    = n_annotations ~ s(log_yolo_pred, log_pred_diff) + s(log_cascade_pred, log_pred_diff)
)

cat("--- Synergy Selection: Georges Bank ---\n")
print(compare_image_gams(img_combined, image_candidate_forms, "GB"))

cat("--- Synergy Selection: Mid-Atlantic Bight ---\n")
print(compare_image_gams(img_combined, image_candidate_forms, "MAB"))

# Assign the winners dynamically based on your tests!
# M09_MeanDiff best for total count error
best_syn_gb  <- image_candidate_forms$M09_MeanCasDiff
best_syn_mab <- image_candidate_forms$M09_MeanCasDiff
# ======================================================================
# 3. Train and Test Final Synergistic Model
# ======================================================================
syn_gams <- fit_synergy_gams(img_combined, best_syn_gb, best_syn_mab)
syn_eval <- test_synergy_gams(syn_gams, img_combined)
# ======================================================================
# 4. The 12-Panel Plot Generation
# ======================================================================

# A. Helper function with axis controls and zoom capabilities
make_1to1_plot <- function(data, x_col, region_name, plot_title, r2, rmse, pt_color, 
                           show_x = FALSE, show_y = FALSE, zoom_limit = NULL) {
  
  val_r2 <- round(r2, 3)
  val_rmse <- round(rmse, 3)
  label_text <- sprintf("R² = %.3f\nRMSE = %.3f", val_r2, val_rmse)
  
  p <- ggplot(data %>% filter(region == region_name), aes_string(x = x_col, y = "n_annotations")) +
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
p[[1]] <- make_1to1_plot(y_test, "predicted_number", "GB", "YOLOv12 \u03A3 P(detection)", ym$r2_calibrated[ym$region=="GB"], ym$rmse_calibrated[ym$region=="GB"], "red", show_x = FALSE, show_y = TRUE)
p[[2]] <- make_1to1_plot(y_test, "predicted_f1_number", "GB", expression("YOLOv12 F"[1] ~ "Threshold"), ym$r2_f1[ym$region=="GB"], ym$rmse_f1[ym$region=="GB"], "red", show_x = FALSE, show_y = FALSE)
p[[3]] <- make_1to1_plot(y_test, "predicted_number", "MAB", "YOLOv12 \u03A3 P(detection)", ym$r2_calibrated[ym$region=="MAB"], ym$rmse_calibrated[ym$region=="MAB"], "red", show_x = FALSE, show_y = FALSE)
p[[4]] <- make_1to1_plot(y_test, "predicted_f1_number", "MAB", expression("YOLOv12 F"[1] ~ "Threshold"), ym$r2_f1[ym$region=="MAB"], ym$rmse_f1[ym$region=="MAB"], "red", show_x = FALSE, show_y = FALSE)

# -- Row 2: Cascade R-CNN (Blue Points) --
p[[5]] <- make_1to1_plot(c_test, "predicted_number", "GB", "Cascade R-CNN \u03A3 P(detection)", cm$r2_calibrated[cm$region=="GB"], cm$rmse_calibrated[cm$region=="GB"], "#2C7FB8", show_x = FALSE, show_y = TRUE)
p[[6]] <- make_1to1_plot(c_test, "predicted_f1_number", "GB", expression("Cascade R-CNN F"[1] ~ "Threshold"), cm$r2_f1[cm$region=="GB"], cm$rmse_f1[cm$region=="GB"], "#2C7FB8", show_x = FALSE, show_y = FALSE)
p[[7]] <- make_1to1_plot(c_test, "predicted_number", "MAB", "Cascade R-CNN \u03A3 P(detection)", cm$r2_calibrated[cm$region=="MAB"], cm$rmse_calibrated[cm$region=="MAB"], "#2C7FB8", show_x = FALSE, show_y = FALSE)
p[[8]] <- make_1to1_plot(c_test, "predicted_f1_number", "MAB", expression("Cascade R-CNN F"[1] ~ "Threshold"), cm$r2_f1[cm$region=="MAB"], cm$rmse_f1[cm$region=="MAB"], "#2C7FB8", show_x = FALSE, show_y = FALSE)

# -- Row 3: Synergy Image-GAM (Purple Points, Full and Zoomed) --
p[[9]]  <- make_1to1_plot(s_test, "pred_synergy", "GB", "Synergistic GAM", sm$r2_synergy[sm$region=="GB"], sm$rmse_synergy[sm$region=="GB"], "#7570B3", show_x = TRUE, show_y = TRUE)
p[[10]] <- make_1to1_plot(s_test, "pred_synergy", "GB", "Synergistic GAM (95% Data)", sm$r2_synergy[sm$region=="GB"], sm$rmse_synergy[sm$region=="GB"], "#7570B3", show_x = TRUE, show_y = FALSE, zoom_limit = limit_gb)
p[[11]] <- make_1to1_plot(s_test, "pred_synergy", "MAB", "Synergistic GAM", sm$r2_synergy[sm$region=="MAB"], sm$rmse_synergy[sm$region=="MAB"], "#7570B3", show_x = TRUE, show_y = FALSE)
p[[12]] <- make_1to1_plot(s_test, "pred_synergy", "MAB", "Synergistic GAM (95% Data)", sm$r2_synergy[sm$region=="MAB"], sm$rmse_synergy[sm$region=="MAB"], "#7570B3", show_x = TRUE, show_y = FALSE, zoom_limit = limit_mab)


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

library(forcats) 
library(latex2exp)
library(ggplot2)
library(dplyr)

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

# Save the final render
ggsave("../figures/ms_figures/2226_final_abundance_summary.png", 
       plot = p_sum_final, width = 8, height = 6, dpi = 300)