# ======================================================================
# 7. Environmental and Spatial Residual Analysis
# ======================================================================
library(scales) # Required for the squish function

# A. Prepare the Data
# Merge the F1 predictions into the synergistic test dataframe
cov_error_data <- s_test %>%
  # Bring in YOLO F1 predictions
  left_join(y_test %>% select(image_id, yolo_f1 = predicted_f1_number), 
            by = "image_id") %>%
  # Bring in Cascade F1 predictions
  left_join(c_test %>% select(image_id, cascade_f1 = predicted_f1_number), 
            by = "image_id") %>%
  # Calculate residuals: Predicted - True
  mutate(
    error_YOLO_F1 = yolo_f1 - n_annotations,
    error_Cascade_F1 = cascade_f1 - n_annotations,
    error_Synergy = pred_synergy - n_annotations
  ) %>%
  # Pivot to long format for ggplot faceting/coloring
  tidyr::pivot_longer(
    cols = c(error_YOLO_F1, error_Cascade_F1, error_Synergy),
    names_to = "Method",
    values_to = "Count_Error"
  ) %>%
  mutate(
    # Clean up names and set factor levels for consistent legend ordering
    Method = factor(Method, 
                    levels = c("error_YOLO_F1", "error_Cascade_F1", "error_Synergy"),
                    labels = c("YOLOv12 F1", "Cascade F1", "Synergistic GAM"))
  ) %>%
  filter(ctd_temperature_celsius > 0)

# Define the color palette matching your abundance summary plot
method_colors <- c(
  "YOLOv12 F1"      = "#FB6A4A", 
  "Cascade F1"      = "#6BAED6", 
  "Synergistic GAM" = "#7570B3"
)

# B. Helper function to plot 1D continuous covariates
# Added y_limits argument to control visual zoom without dropping data for geom_smooth
# B. Helper function to plot 1D continuous covariates with binned pointranges
plot_covariate_error <- function(df, x_var, x_label, y_limits = c(-3, 3), n_bins = 15) {
  
  # 1. Dynamically pre-calculate binned means and standard errors for the points
  binned_data <- df %>%
    mutate(x_binned = ggplot2::cut_interval(!!sym(x_var), n = n_bins)) %>%
    group_by(Method, x_binned) %>%
    summarize(
      # Find the center of the bin for plotting on the x-axis
      x_mid = mean(!!sym(x_var), na.rm = TRUE),
      mean_err = mean(Count_Error, na.rm = TRUE),
      # Calculate Standard Error (SE) = Standard Deviation / sqrt(n)
      se_err = sd(Count_Error, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    filter(!is.na(x_mid))
  
  # 2. Calculate a dynamic dodge width based on the range of the specific x-axis
  # This ensures the three error bars sit nicely side-by-side regardless of the variable's scale
  x_range <- max(df[[x_var]], na.rm = TRUE) - min(df[[x_var]], na.rm = TRUE)
  dodge_width <- x_range / (n_bins * 2.5) 
  
  # 3. Build the plot
  ggplot(df, aes_string(x = x_var, y = "Count_Error", color = "Method", fill = "Method")) +
    # Add a zero-error reference line
    # geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.5) +
    
    # Smooth trend lines using the RAW continuous dataframe (df)
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, alpha = 0.2, linewidth = 1.2) +
    
    # Binned mean and SE points using the SUMMARIZED dataframe (binned_data)
    geom_pointrange(
      data = binned_data,
      aes(x = x_mid, y = mean_err, ymin = mean_err - se_err, ymax = mean_err + se_err),
      position = position_dodge(width = dodge_width),
      size = 0.4, 
      alpha = 0.9, 
      shape = 21, 
      fill = "white" # Gives the points a clean, hollow look to distinguish them from the lines
    ) +
    
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    # coord_cartesian strictly zooms the view; it does NOT alter the geom_smooth math
    coord_cartesian(ylim = y_limits) + 
    theme_bw(base_size = 12) +
    # Use parsed expressions for standard math formatting
    labs(x = x_label, y = expression(N[pred] - N[true])) +
    theme(
      legend.position = "none", 
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

# C. Generate the 1D Covariate Plots
# You can adjust the c(-3, 3) limits based on the actual distribution of your count errors
p_alt   <- plot_covariate_error(cov_error_data, "altitude", "Altitude (m)", y_limits = c(-3, 3)) + ggtitle("HabCam Vehicle Altitude")
p_depth <- plot_covariate_error(cov_error_data, "bottom_depth", "Bottom Depth (m)", y_limits = c(-3, 3)) + ggtitle("Bottom Depth")
p_temp  <- plot_covariate_error(cov_error_data, "ctd_temperature_celsius", "Temperature (°C)", y_limits = c(-3, 3)) + ggtitle("CTD Temperature")
p_fov   <- plot_covariate_error(cov_error_data, "field_of_view_sq_meter", "Field of View (m²)", y_limits = c(-3, 3)) + ggtitle("Image Field of View")

# p_alt
# p_depth
# p_temp
# p_fov
# D. Generate the 2D Spatial Error Plot
# Faceted by method to show if errors clump geographically
p_space <- ggplot(cov_error_data, aes(x = longitude, y = latitude, color = Count_Error)) +
  geom_point(alpha = 0.85, size = 1.05) +
  facet_wrap(~ Method, ncol = 3) +
  # Diverging palette: Over-prediction (Red), Under-prediction (Blue), Accurate (White/Grey)
  # scales::squish takes values outside c(-2, 2) and caps them at the max/min color, preventing NA blank spots
  scale_color_gradient2(
    low = "#2C7FB8", mid = "grey90", high = "#E34A33", midpoint = 0, 
    limits = c(-2, 2), 
    oob = scales::squish
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Spatial Distribution of Residual Errors",
    x = "Longitude", 
    y = "Latitude",
    color = "Count Error (Pred - True)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )
p_space
# E. Assemble with Patchwork
# Top row: 1D environmental covariates (now a 4-panel row)
# Bottom row: Spatial maps
# Make sure to install and load ggtext
# install.packages("ggtext")
library(ggtext)

# Assemble with Patchwork
cov_assembly <- (p_alt | p_depth | p_temp | p_fov) / p_space + 
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(
    title = "Bias Resolution: F1 Thresholds vs. Synergistic GAM Calibration",
    # Add <br> for a line break, and <span> tags for specific colors
    subtitle = "Top Row Trend Lines:</b> <span style='color:#FB6A4A;'>YOLOv12 F1</span> &bull; <span style='color:#6BAED6;'>Cascade F1</span> &bull; <span style='color:#7570B3;'>Synergistic GAM</span>",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      # IMPORTANT: Swap element_text() for element_markdown() to render the HTML
      plot.subtitle = element_markdown(size = 12, hjust = 0.5, color = "grey30", lineheight = 1.3)
    )
  )

print(cov_assembly)

# F. Save the figure
ggsave(
  filename = "../figures/diag_crab1/2426_crab_covariate_error_comparison.pdf",
  plot = cov_assembly,
  width = 11.5, # Increased width slightly to accommodate 4 columns in the top row
  height = 8
)


# ======================================================================
# Diagnostic: Quantifying Error Volatility and Covariate Stability
# ======================================================================

# 1. Calculate Global Frame-Level Variance and Dispersion
global_variance_metrics <- cov_error_data %>%
  group_by(Method) %>%
  summarize(
    Global_SD   = sd(Count_Error, na.rm = TRUE),
    Global_MAE  = mean(abs(Count_Error), na.rm = TRUE),
    Max_Frame_Error = max(abs(Count_Error), na.rm = TRUE),
    .groups = "drop"
  )

print("--- GLOBAL FRAME-LEVEL ERROR METRICS ---")
print(global_variance_metrics)


# 2. Calculate Environmental Stability (How much error swings across bins)
quantify_environmental_swings <- function(df, x_var, var_label, n_bins = 15) {
  df %>%
    mutate(x_binned = ggplot2::cut_interval(!!sym(x_var), n = n_bins)) %>%
    group_by(Method, x_binned) %>%
    summarize(mean_err = mean(Count_Error, na.rm = TRUE), .groups = "drop_last") %>%
    filter(!is.na(x_binned)) %>%
    summarize(
      Covariate = var_label,
      SD_of_Bin_Means = sd(mean_err, na.rm = TRUE),
      Peak_Bin_Bias   = max(abs(mean_err), na.rm = TRUE)
    )
}

# Bind metrics for all continuous environmental features
environmental_stability <- bind_rows(
  quantify_environmental_swings(cov_error_data, "altitude", "Altitude"),
  quantify_environmental_swings(cov_error_data, "bottom_depth", "Bottom Depth"),
  quantify_environmental_swings(cov_error_data, "ctd_temperature_celsius", "Temperature"),
  quantify_environmental_swings(cov_error_data, "field_of_view_sq_meter", "Field of View")
)

print("--- COVARIATE RESIDUAL STABILITY METRICS ---")
print(environmental_stability %>% arrange(Covariate, Method))


# ======================================================================
# Corrected Diagnostic: True 2D Global Spatial Error Metrics
# ======================================================================

quantify_true_spatial_grid_errors <- function(df, grid_size = 0.1) {
  df %>%
    # 1. Group coordinates into discrete spatial grid cells
    mutate(
      lat_grid = round(latitude / grid_size) * grid_size,
      lon_grid = round(longitude / grid_size) * grid_size
    ) %>%
    # 2. Calculate the mean counting error inside each cell
    group_by(Method, lat_grid, lon_grid) %>%
    summarize(mean_cell_err = mean(Count_Error, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(mean_cell_err)) %>%
    # 3. UNGROUP COMPLETELY to calculate global spatial metrics across the whole grid
    group_by(Method) %>%
    summarize(
      Spatial_Resolution  = paste0("2D Grid (", grid_size, "° cell)"),
      Spatial_Grid_MAE    = mean(abs(mean_cell_err), na.rm = TRUE),
      Spatial_Grid_SD     = sd(mean_cell_err, na.rm = TRUE),
      Max_Localized_Bias  = max(abs(mean_cell_err), na.rm = TRUE),
      .groups = "drop"
    )
}

# Run the true spatial metric
true_spatial_stability <- quantify_true_spatial_grid_errors(cov_error_data, grid_size = 0.1)
print(true_spatial_stability)
