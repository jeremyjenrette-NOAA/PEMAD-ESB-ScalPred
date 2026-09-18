library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(readr)

# ======================================================================
# 1. Load Master Calibrated Inference Object
# ======================================================================
cat("Loading calibrated inference dataset...\n")
master_inference <- readRDS("../data/processed/eval_inference_calibrated_2226.rds")

# Create output directories if they don't exist
# dir.create("../tables", showWarnings = FALSE, recursive = TRUE)
# dir.create("../figures/ms_figures", showWarnings = FALSE, recursive = TRUE)

syn_df <- master_inference$synergistic

# ======================================================================
# 2. Manuscript Summary Table: Stratified Abundance Estimates
# ======================================================================
cat("Generating stratified abundance summary statistics...\n")

abundance_summary <- syn_df %>%
  group_by(year, region) %>%
  summarise(
    n_images              = n(),
    
    # Raw & F1 Cutoff Workflow Estimates
    YOLO_F1_Total         = sum(yolo_f1_pred, na.rm = TRUE),
    Cascade_F1_Total      = sum(cascade_f1_pred, na.rm = TRUE),
    
    # Detection-Level GAM Calibrated Workflow Estimates
    YOLO_DetGAM_Total     = sum(yolo_pred, na.rm = TRUE),
    Cascade_DetGAM_Total  = sum(cascade_pred, na.rm = TRUE),
    
    # Synergistic GAM Fusion Workflow Estimate
    Synergistic_GAM_Total = sum(final_abundance, na.rm = TRUE),
    
    # Density Metrics (if field_of_view_sq_meter exists)
    total_area_sq_m       = if("field_of_view_sq_meter" %in% colnames(.)) sum(field_of_view_sq_meter, na.rm = TRUE) else NA_real_,
    synergy_density_m2    = if(!is.na(total_area_sq_m[1])) Synergistic_GAM_Total / total_area_sq_m else NA_real_,
    
    .groups = "drop"
  ) %>%
  mutate(
    # Quantify relative divergence between single models and synergistic fusion
    yolo_vs_synergy_pct_diff = ((YOLO_DetGAM_Total - Synergistic_GAM_Total) / Synergistic_GAM_Total) * 100,
    cas_vs_synergy_pct_diff  = ((Cascade_DetGAM_Total - Synergistic_GAM_Total) / Synergistic_GAM_Total) * 100
  )

# Export Summary Table to CSV
# write_csv(abundance_summary, "../tables/supp_abundance_summary_by_year_region.csv")
# cat("Summary table exported to '../tables/supp_abundance_summary_by_year_region.csv'\n")

# ======================================================================
# 3. Figure S1: Detection-Level Calibration Curves (Conf vs. P_det)
# ======================================================================
cat("Generating Figure S1: Detection-Level Calibration Curves...\n")

yolo_det <- master_inference$YOLOv12$det %>% 
  mutate(detector = "YOLOv12") %>% 
  select(conf, pred_p, region, detector)

cas_det <- master_inference$`Cascade R-CNN`$det %>% 
  mutate(detector = "Cascade R-CNN") %>% 
  select(conf, pred_p, region, detector)

combined_det <- bind_rows(yolo_det, cas_det)

fig_s1 <- ggplot(combined_det, aes(x = conf, y = pred_p, color = detector)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.2) +
  facet_wrap(~ region, labeller = as_labeller(c("GB" = "Georges Bank", "MAB" = "Mid-Atlantic Bight"))) +
  scale_color_manual(values = c("YOLOv12" = "#CB181D", "Cascade R-CNN" = "#2171B5")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 13) +
  labs(
    title = "Detection-Level Calibration Curves",
    subtitle = "Mapping raw detector confidence to GAM calibrated detection probability p(det)",
    x = "Detector Raw Confidence (conf)",
    y = "Calibrated Detection Probability p(det)",
    color = "Detector"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
fig_s1
# ggsave("../figures/ms_figures/supp_fig_S1_detection_calibration.png", plot = fig_s1, width = 10, height = 6, dpi = 300)

# ======================================================================
# 4. Figure S2: Spatial Scallop Abundance Density Map
# ======================================================================
cat("Generating Figure S2: Spatial Abundance Map...\n")

fig_s2 <- ggplot(syn_df, aes(x = longitude, y = latitude, z = final_abundance)) +
  stat_summary_hex(fun = "sum", bins = 80) +
  scale_fill_viridis_c(option = "magma", trans = "log1p", name = "Est. Count\n(log scale)") +
  facet_wrap(~ region, scales = "free", labeller = as_labeller(c("GB" = "Georges Bank", "MAB" = "Mid-Atlantic Bight"))) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Spatial Scallop Abundance Density",
    subtitle = "Aggregated Synergistic GAM estimates across survey domain",
    x = "Longitude (°W)",
    y = "Latitude (°N)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
fig_s2
# ggsave("../figures/ms_figures/supp_fig_S2_spatial_abundance_map.png", plot = fig_s2, width = 11, height = 6, dpi = 300)

# ======================================================================
# 5. Figure S3: Abundance Relationships with Survey Covariates
# ======================================================================
cat("Generating Figure S3: Covariate Trends (Depth & Altitude)...\n")

p_depth <- ggplot(syn_df, aes(x = bottom_depth, y = final_abundance)) +
  geom_point(alpha = 0.1, color = "grey30", size = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "#7570B3", linewidth = 1.2) +
  facet_wrap(~ region, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(title = "Estimated Abundance vs. Bottom Depth", x = "Bottom Depth (m)", y = "Synergistic Abundance per Image")

p_alt <- ggplot(syn_df, aes(x = altitude, y = final_abundance)) +
  geom_point(alpha = 0.1, color = "grey30", size = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "#1B9E77", linewidth = 1.2) +
  facet_wrap(~ region, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(title = "Estimated Abundance vs. Vehicle Altitude", x = "Altitude (m)", y = "Synergistic Abundance per Image")

fig_s3 <- p_depth / p_alt + 
  plot_annotation(
    title = "Scallop Abundance Distribution across Survey Covariates",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )
fig_s3
# ggsave("../figures/ms_figures/supp_fig_S3_covariate_trends.png", plot = fig_s3, width = 10, height = 8, dpi = 300)

# ======================================================================
# 6. Figure S4: Yearly & Regional Workflow Comparison
# ======================================================================
cat("Generating Figure S4: Workflow Comparison Bar Plots...\n")

plot_df <- abundance_summary %>%
  select(year, region, YOLO_F1_Total, Cascade_F1_Total, YOLO_DetGAM_Total, Cascade_DetGAM_Total, Synergistic_GAM_Total) %>%
  pivot_longer(
    cols = c(-year, -region),
    names_to = "Workflow",
    values_to = "Abundance"
  ) %>%
  mutate(
    Workflow = case_when(
      Workflow == "YOLO_F1_Total"         ~ "YOLOv12 F1",
      Workflow == "Cascade_F1_Total"      ~ "Cascade F1",
      Workflow == "YOLO_DetGAM_Total"     ~ "YOLOv12 Det-GAM",
      Workflow == "Cascade_DetGAM_Total"  ~ "Cascade Det-GAM",
      Workflow == "Synergistic_GAM_Total" ~ "Synergistic GAM"
    ),
    Workflow = factor(Workflow, levels = c("YOLOv12 F1", "Cascade F1", "YOLOv12 Det-GAM", "Cascade Det-GAM", "Synergistic GAM"))
  )

fig_s4 <- ggplot(plot_df, aes(x = factor(year), y = Abundance, fill = Workflow)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  facet_wrap(~ region, scales = "free_y") +
  scale_fill_manual(values = c(
    "YOLOv12 F1"      = "#FB6A4A", 
    "Cascade F1"      = "#6BAED6", 
    "YOLOv12 Det-GAM" = "#CB181D", 
    "Cascade Det-GAM" = "#2171B5", 
    "Synergistic GAM" = "#7570B3"
  )) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Annual Scallop Abundance Estimates Across Workflows",
    subtitle = "Comparing F1 Cutoffs, Detection GAMs, and Synergistic Fusion",
    x = "Survey Year",
    y = "Total Abundance Estimate",
    fill = "Workflow"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
fig_s4
# ggsave("../figures/ms_figures/supp_fig_S4_annual_workflow_comparison.png", plot = fig_s4, width = 11, height = 6, dpi = 300)

cat("All manuscript supplementary figures and tables successfully generated and saved!\n")