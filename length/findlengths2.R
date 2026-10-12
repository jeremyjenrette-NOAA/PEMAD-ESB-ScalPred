# ==============================================================================
# Script: sampling_simulation.R
# Purpose: Parallel bootstrap simulation of sampling effort for scallop size-classes
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(foreach)
library(doParallel)
source("lengthfunc2.R")

# ------------------------------------------------------------------------------
# 1. Setup & Helper Functions
# ------------------------------------------------------------------------------

gt_2226 = read.csv("../data/raw/groundtruth2226.csv")
yolo_eval <- readRDS("../data/processed/YOLOv12_deteval_2226.rds")$det_eval
cas_eval  <- readRDS("../data/processed/CascadeR-CNN_deteval_2226.rds")$det_eval
dat_split <- read.csv("../data/raw/dataset_split_2226.csv", stringsAsFactors = FALSE)

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

# Re-encapsulating your binning function to be self-contained for parallel workers
# Uniformly bin distribution frequencies
bin_frequencies_sim <- function(data_vec, weight_vec = NULL, bin_breaks) {
  if (length(data_vec) == 0) return(rep(0, length(bin_breaks) - 1))
  
  binned <- cut(data_vec, breaks = bin_breaks, include.lowest = TRUE)
  
  if (is.null(weight_vec)) {
    counts <- tabulate(binned, nbins = length(bin_breaks) - 1)
  } else {
    counts <- sapply(levels(binned), function(lvl) {
      sum(weight_vec[binned == lvl], na.rm = TRUE)
    })
  }
  
  total <- sum(counts)
  if (total == 0) return(rep(0, length(bin_breaks) - 1))
  return(counts / total)
}

# Quantify intersection between manual and AI curves
compute_ovl_sim <- function(p_manual, p_pred) {
  if (length(p_manual) == 0 || length(p_pred) == 0) return(0)
  return(sum(pmin(p_manual, p_pred)))
}

# ------------------------------------------------------------------------------
# 2. Main Simulation Engine
# ------------------------------------------------------------------------------

# ==============================================================================
# Script: sampling_simulation.R (Optimized)
# Purpose: Parallel bootstrap simulation with global pre-calculation
# ==============================================================================

run_sampling_simulation <- function(gt_data, 
                                    det_data, 
                                    sample_sizes = seq(20, 300, by = 20), 
                                    iterations = 100, 
                                    bin_width = 5,
                                    max_size = 220,
                                    calc_method = "average",
                                    num_cores = parallel::detectCores() - 1) {
  
  # ----------------------------------------------------------------------------
  # Step 1: Pre-Process Globally (ON THE MAIN THREAD)
  # This runs exactly once, bypassing parallel environmental scoping issues completely.
  # ----------------------------------------------------------------------------
  message("Pre-processing spatial data and calculating lengths on main thread...")
  
  gt_processed <- gt_data %>%
    mutate(
      true_pixel_length = parse_gt_length(geometry_text),
      true_length_mm = true_pixel_length * millimeter_per_pixel
    ) %>%
    filter(!is.na(true_length_mm)) %>%
    select(year, region, imagename, true_length_mm)
  
  det_processed <- det_data %>%
    mutate(
      pred_pixel_length = back_calc_length(t_lx, t_ly, b_rx, b_ry, method = calc_method),
      pred_length_mm = pred_pixel_length * millimeter_per_pixel
    ) %>%
    select(year, region, imagename, pred_length_mm, pred_p, regional_f1_thresh)
  
  # Define bin breaks
  bin_breaks <- seq(0, max_size, by = bin_width)
  
  # Identify valid strata using the processed data
  strata <- det_processed %>% 
    select(year, region) %>% 
    distinct() %>% 
    drop_na()
  
  # ----------------------------------------------------------------------------
  # Step 2: Initialize Parallel Backend
  # ----------------------------------------------------------------------------
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  message(paste("Running bootstrap simulation across", nrow(strata), "strata using", num_cores, "cores..."))
  
  # ----------------------------------------------------------------------------
  # Step 3: Run the Parallelized Bootstrap
  # ----------------------------------------------------------------------------
  simulation_results <- foreach(
    s = 1:nrow(strata), 
    .combine = rbind, 
    .packages = c("dplyr", "tidyr"),
    .export = c("bin_frequencies_sim", "compute_ovl_sim") # Export simpler helpers
  ) %dopar% {
    
    current_year   <- strata$year[s]
    current_region <- strata$region[s]
    
    # Filter pre-calculated data frames
    stratum_gt  <- gt_processed %>% filter(year == current_year, region == current_region)
    stratum_det <- det_processed %>% filter(year == current_year, region == current_region)
    
    # Define active image pool
    image_pool <- unique(c(stratum_det$imagename, stratum_gt$imagename))
    
    if (length(image_pool) < 10) return(NULL) # Skip if data is too sparse
    
    stratum_results <- list()
    idx <- 1
    
    for (n_images in sample_sizes) {
      for (iter in 1:iterations) {
        
        # 1. Sample images with replacement
        sampled_images <- data.frame(
          bootstrap_id = 1:n_images,
          imagename = sample(image_pool, size = n_images, replace = TRUE)
        )
        
        # 2. Match corresponding pre-calculated scallops
        sim_gt <- sampled_images %>% 
          inner_join(stratum_gt, by = "imagename", relationship = "many-to-many")
        
        sim_det <- sampled_images %>% 
          inner_join(stratum_det, by = "imagename", relationship = "many-to-many")
        
        # 3. Bin distributions
        p_true <- bin_frequencies_sim(sim_gt$true_length_mm, bin_breaks = bin_breaks)
        
        # F1 Thresholded predictions
        f1_thresh <- if(nrow(sim_det) > 0) unique(sim_det$regional_f1_thresh)[1] else 0.5
        sim_det_f1 <- sim_det %>% filter(pred_p >= f1_thresh)
        p_f1 <- bin_frequencies_sim(sim_det_f1$pred_length_mm, bin_breaks = bin_breaks)
        
        # Weighted P predictions
        p_weighted <- bin_frequencies_sim(
          sim_det$pred_length_mm, 
          weight_vec = sim_det$pred_p, 
          bin_breaks = bin_breaks
        )
        
        # 4. Compute Discrete Overlap Metrics (OVL)
        ovl_f1       <- compute_ovl_sim(p_true, p_f1)
        ovl_weighted <- compute_ovl_sim(p_true, p_weighted)
        
        # Append result row
        stratum_results[[idx]] <- data.frame(
          year            = current_year,
          region          = current_region,
          sample_size     = n_images,
          iteration       = iter,
          ovl_f1          = ovl_f1,
          ovl_weighted    = ovl_weighted,
          n_gt_scallops   = nrow(sim_gt),
          n_ai_detections = nrow(sim_det_f1)
        )
        idx <- idx + 1
      }
    }
    
    do.call(rbind, stratum_results)
  }
  
  return(simulation_results)
}


# --- Run the Simulation ---
# (Assumes gt_clean and yolo$det_eval are prepared and loaded in your workspace)
sim_data <- run_sampling_simulation(
  gt_data      = gt_clean,       # Ground truth dataframe with geometry_text
  det_data     = yolo$det_eval,  # Detections dataframe with pred_p
  sample_sizes = seq(20, 5000, by = 20),
  iterations   = 100,            # High numbers yield cleaner confidence bands
  calc_method  = "average"       # Bounding box-to-line calculation method
)

# --- Process Results for Plotting ---
sim_summary <- sim_data %>%
  pivot_longer(
    cols = c(ovl_f1, ovl_weighted), 
    names_to = "Method", 
    values_to = "OVL"
  ) %>%
  mutate(
    Method = recode(Method, "ovl_f1" = "F1 Threshold", "ovl_weighted" = "Weighted P(detection)")
  ) %>%
  group_by(year, region, sample_size, Method) %>%
  summarise(
    mean_ovl   = mean(OVL, na.rm = TRUE),
    sd_ovl     = sd(OVL, na.rm = TRUE),
    lower_95   = quantile(OVL, 0.025, na.rm = TRUE),
    upper_95   = quantile(OVL, 0.975, na.rm = TRUE),
    .groups    = "drop"
  )

# --- Generate the Power Curve Visualization ---
plot_power_curve <- ggplot(sim_summary, aes(x = sample_size, y = mean_ovl, color = Method)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = Method), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  facet_grid(region ~ year) +
  labs(
    title = "Sampling Effort Simulation: Size-Class Replication Power",
    subtitle = "Quantifying minimum image count required to accurately resolve population structure",
    x = "Sample Size (Number of Images)",
    y = "Distribution Likeness (Discrete OVL vs. Ground Truth)",
    fill = "AI Integration Strategy",
    color = "AI Integration Strategy"
  ) +
  scale_color_manual(values = c("F1 Threshold" = "#0072B2", "Weighted P(detection)" = "#D9A00B")) +
  scale_fill_manual(values = c("F1 Threshold" = "#0072B2", "Weighted P(detection)" = "#D9A00B")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "#f2f2f2", color = "#d9d9d9"),
    strip.text = element_text(face = "bold")
  )

print(plot_power_curve)
