# ==============================================================================
# Script: measuring.R
# Purpose: Grid-Stratified Empirical Size-Frequency Polygons & Subclass Metrics
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

source("lengthfunc2.R")

# 1. Load Data -----------------------------------------------------------------
gt_2226   <- read.csv("../data/raw/groundtruth2226.csv", stringsAsFactors = FALSE)
dat_split <- read.csv("../data/raw/dataset_split_2226.csv", stringsAsFactors = FALSE)

yolo_eval <- readRDS("../data/processed/YOLOv12_deteval_2226.rds")$det_eval
cas_eval  <- readRDS("../data/processed/CascadeR-CNN_deteval_2226.rds")$det_eval

test_frames <- dat_split %>% filter(dataset == "test_GAM_test") %>% pull(imagename)

# 2. Process Ground Truth Data -------------------------------------------------
gt_clean <- gt_2226 %>%
  mutate(class_name = if_else(class_name == "Scallop", "probable live sea scallop", class_name)) %>%
  filter(imagename %in% test_frames) %>%
  filter(class_name != "Scallop") %>%
  mutate(subclass = assign_scallop_subclass(class_name)) %>%
  group_by(imagename) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  mutate(
    true_pixel_length = parse_gt_length(geometry_text),
    true_length_mm = true_pixel_length * millimeter_per_pixel
  ) %>%
  mutate(
    region = if_else(longitude >= -71, "GB", "MAB") %>% factor(levels = c("MAB", "GB"))
  ) %>%
  filter(!is.na(true_length_mm))

# 3. Process Detections and Compute Error Statistics per Subclass -------------
calc_method <- "average"

matched_yolo <- yolo_eval %>%
  filter(imagename %in% test_frames & truedetect == TRUE) %>%
  inner_join(gt_clean %>% select(imagename, index, subclass, true_length_mm), by = c("imagename", "index")) %>%
  mutate(
    pred_pixel_length = back_calc_length(t_lx, t_ly, b_rx, b_ry, method = calc_method),
    pred_length_mm = pred_pixel_length * millimeter_per_pixel,
    raw_error = pred_length_mm - true_length_mm,
    abs_error = abs(raw_error)
  )

matched_cas <- cas_eval %>%
  filter(imagename %in% test_frames & truedetect == TRUE) %>%
  inner_join(gt_clean %>% select(imagename, index, subclass, true_length_mm), by = c("imagename", "index")) %>%
  mutate(
    pred_pixel_length = back_calc_length(t_lx, t_ly, b_rx, b_ry, method = calc_method),
    pred_length_mm = pred_pixel_length * millimeter_per_pixel,
    raw_error = pred_length_mm - true_length_mm,
    abs_error = abs(raw_error)
  )

# --- REFINED EMPIRICAL LIKENESS METRIC ENGINE (Discrete OVL) ---
# Computes the true intersection of empirical discrete frequency distributions
compute_discrete_ovl <- function(p_manual, p_pred) {
  if (length(p_manual) == 0 || length(p_pred) == 0) return(NA)
  return(sum(pmin(p_manual, p_pred)))
}

# ==============================================================================
# 4. Generate Pre-Calculated Empirical Bins & True Baselines (5mm Resolution)
# ==============================================================================
bin_width <- 5
bin_breaks <- seq(0, 220, by = bin_width)
bin_mids <- bin_breaks[-length(bin_breaks)] + (bin_width / 2)

# Helper function to compute binned percentages uniformly
bin_frequencies <- function(data_vec, weight_vec = NULL) {
  binned <- cut(data_vec, breaks = bin_breaks, include.lowest = TRUE)
  if (is.null(weight_vec)) {
    counts <- tabulate(binned, nbins = length(bin_mids))
  } else {
    counts <- sapply(levels(binned), function(lvl) {
      sum(weight_vec[binned == lvl], na.rm = TRUE)
    })
  }
  total <- sum(counts)
  if (total == 0) return(rep(0, length(bin_mids)))
  return(counts / total)
}

# --- CRITICAL FIX: Extract True Global Manual Population Distributions Upfront ---
global_manual_bins <- gt_clean %>%
  group_by(subclass) %>%
  summarise(p_man = list(bin_frequencies(true_length_mm)), .groups = "drop")

# Compile Subclass Statistics & Panel Labels evaluated against Global Truth
generate_panel_stats <- function(df_matched, model_name) {
  df_matched %>%
    group_by(subclass) %>%
    summarise(
      MAE = mean(abs_error, na.rm = TRUE),
      Bias = mean(raw_error, na.rm = TRUE),
      p_f1  = list(bin_frequencies(pred_length_mm[pred_p >= unique(regional_f1_thresh)])),
      p_wt  = list(bin_frequencies(pred_length_mm, pred_p)),
      .groups = "drop"
    ) %>%
    inner_join(global_manual_bins, by = "subclass") %>% 
    group_by(subclass) %>%
    mutate(
      OVL_F1 = compute_discrete_ovl(p_man[[1]], p_f1[[1]]),
      OVL_GAM = compute_discrete_ovl(p_man[[1]], p_wt[[1]])
    ) %>%
    ungroup() %>%
    mutate(
      Model = model_name,
      Label = paste0("MAE: ", round(MAE, 1), " mm\nBias: ", round(Bias, 1), 
                     " mm\nOVL F1: ", round(OVL_F1, 3), "\nOVL GAM: ", round(OVL_GAM, 3))
    )
}

stats_yolo <- generate_panel_stats(matched_yolo, "YOLOv12")
stats_cas  <- generate_panel_stats(matched_cas, "Cascade R-CNN")
panel_annotations <- bind_rows(stats_yolo, stats_cas)

# Construct Model Predictions Only
build_predictions_df <- function(df_matched, model_name) {
  df_matched %>%
    group_by(subclass) %>%
    do({
      sub_data <- .
      p_f1  <- bin_frequencies(sub_data$pred_length_mm[sub_data$pred_p >= unique(sub_data$regional_f1_thresh)])
      p_wt  <- bin_frequencies(sub_data$pred_length_mm, sub_data$pred_p)
      
      data.frame(
        Bin_Mid = rep(bin_mids, 2),
        Percentage = c(p_f1, p_wt),
        Workflow = rep(c("F1 Threshold", "Weighted P(detection)"), each = length(bin_mids))
      )
    }) %>%
    ungroup() %>%
    mutate(Model = model_name)
}

# Construct Absolute Manual Baseline Dataset directly from gt_clean
manual_baseline_df <- gt_clean %>%
  group_by(subclass) %>%
  do({
    p_man <- bin_frequencies(.$true_length_mm)
    data.frame(
      Bin_Mid = bin_mids,
      Percentage = p_man,
      Workflow = "Manual"
    )
  }) %>%
  ungroup()

# Duplicate the manual track across both model facets so ggplot mirrors them perfectly
manual_faceted_baseline <- bind_rows(
  manual_baseline_df %>% mutate(Model = "YOLOv12"),
  manual_baseline_df %>% mutate(Model = "Cascade R-CNN")
)

# Combine everything into integrated dataset
grid_inventory_df <- bind_rows(
  manual_faceted_baseline,
  build_predictions_df(matched_yolo, "YOLOv12"),
  build_predictions_df(matched_cas, "Cascade R-CNN")
) %>%
  mutate(Workflow = factor(Workflow, levels = c("Manual", "F1 Threshold", "Weighted P(detection)")))

# Extract empirical peaks from global baseline
gt_peaks <- manual_faceted_baseline %>%
  group_by(subclass, Model) %>%
  slice(which.max(Percentage)) %>%
  select(subclass, Model, Peak_Size = Bin_Mid) %>%
  ungroup()

# ==============================================================================
# 5. Build Grid Matrix Size-Frequency Visualization (Empirical Polygons)
# ==============================================================================

plot_grid_ms <- ggplot(grid_inventory_df, aes(x = Bin_Mid, y = Percentage)) +
  # Layer 1: Manual Baseline Area (Shaded area under direct frequency point steps)
  geom_area(data = filter(grid_inventory_df, Workflow == "Manual"),
            fill = "#E6E6E6", color = "black", alpha = 0.6, linewidth = 0.6) +
  
  # Layer 2: Dynamic Empirical Peak Line Markers
  # geom_vline(data = gt_peaks, aes(xintercept = Peak_Size),
  #            color = "red", linetype = "dotted", linewidth = 0.75) +
  
  # Layer 3: Annotate Peak Value Text
  # geom_text(data = gt_peaks, 
  #           aes(x = Peak_Size, y = Inf, label = paste0(round(Peak_Size, 1), " mm")),
  #           color = "red", size = 2.8, fontface = "bold", 
  #           hjust = -0.45, vjust = 1.25, inherit.aes = FALSE) +
  
  # Layer 4: Automated Empirical Frequency Lines (Connected points, no estimation)
  geom_line(data = filter(grid_inventory_df, Workflow != "Manual"),
            aes(color = Workflow), linewidth = 0.65, linetype = "solid") +
  geom_point(data = filter(grid_inventory_df, Workflow != "Manual"),
             aes(color = Workflow), size = 1.0) +
  
  # Layer 5: Text Metadata Summary Blocks
  geom_text(data = panel_annotations,
            aes(x = 210, y = Inf, label = Label),
            hjust = 1, vjust = 1.25, size = 2.8, fontface = "italic", inherit.aes = FALSE) +
  
  facet_grid(subclass ~ Model, scales = "free_y") +
  coord_cartesian(xlim = c(-5, 220)) +
  
  scale_color_manual(
    values = c(
      "F1 Threshold"         = "#0072B2",  
      "Weighted P(detection)"= "#D9A00B"   
    ),
    labels = c(
      "F1 Threshold"         = expression(F[1]~"Threshold"),
      "Weighted P(detection)"= "Weighted P(detection)"
    )
  ) +
  
  scale_y_continuous(
    labels = function(x) paste0(round(x * 100, 1), "%")
  ) +
  
  labs(
    title = "Estimating Size Classes",
    x = "Size Classes (mm)",
    y = "Percentage",
    color = "Detection Method"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "#f2f2f2", color = "#d9d9d9"),
    strip.text = element_text(face = "bold", size = 9.5),
    panel.spacing = unit(1, "lines"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13, hjust = 0)
  )

print(plot_grid_ms)

ggsave(
  filename = "../figures/ms_figures/2226_sizeclasses.pdf",
  plot = plot_grid_ms,
  width = 9.5,
  height = 10,
  units = "in",
  device = "pdf"
)
ggsave(
  filename = "~/saltnoaa/presentations/figures/2226_sizeclasses.png",
  plot = plot_grid_ms,
  width = 11,
  height = 8
)

# ==============================================================================
# 7. Global Distribution Likeness Metrics (Pooled Empirical OVL)
# ==============================================================================

global_ovl_summary <- bind_rows(
  matched_yolo %>% mutate(Model = "YOLOv12"),
  matched_cas %>% mutate(Model = "Cascade R-CNN")
) %>%
  group_by(Model) %>%
  summarise(
    Global_OVL_F1 = compute_discrete_ovl(
      bin_frequencies(true_length_mm),
      bin_frequencies(pred_length_mm[pred_p >= unique(regional_f1_thresh)])
    ),
    Global_OVL_Weighted = compute_discrete_ovl(
      bin_frequencies(true_length_mm),
      bin_frequencies(pred_length_mm, pred_p)
    ),
    .groups = "drop"
  )

print("--- GLOBAL WORKFLOW DISTRIBUTION LIKENESS (EMPIRICAL DISCRETE OVL) ---")
print(global_ovl_summary)
