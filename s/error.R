# ======================================================================
# error.R - Residual / error diagnostics on the HOLDOUT set (test_GAM_test)
#
# Standalone: reads the table saved by gammod.R
#   ../data/processed/holdout_eval_2226_v2.rds  (n3/n3_export_holdout.R)
# Run order: makedat.R -> gamstrat.R -> gammod.R -> n3/n3_ablation.R -> n3/n3_export_holdout.R -> error.R
#
# What it reports (all holdout only):
#   1. Annotation audit: images where BOTH detectors find many scallops but the
#      annotation says (almost) none - likely missing labels. These dominate
#      RMSE (2026-10-06: one image, 0 annotated / ~186 detected, was about half
#      of the Georges Bank MSE for every method). Listed for manual review and
#      used for a clearly labelled sensitivity analysis - never silently dropped.
#      NB: a flag is a prompt for review, not a verdict. Images of dead shells /
#      clappers are correctly annotated as empty (background) yet attract
#      detections; of 5 flagged images reviewed on 2026-10-06 only 1 was a real
#      labelling error. Confirmed errors go in data/raw/annotation_exclusions_2226.csv.
#   2. Metrics by region x density bin and region x cruise (block bootstrap CIs).
#   3. Error vs covariates: quantile bins (equal n per bin), block-clustered SEs,
#      count and density (per m^2) error scales.
#   4. Spatial: gridded mean error maps + Moran's I of image residuals.
#   5. Stability summaries (SD of bin means, peak bin bias, spatial grid MAE).
# ======================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
source("./gamfunc.R")   # r2_1to1(), block_bootstrap_metrics()

# ----------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------
METHODS <- c("Two-stage YOLOv12"           = "pred_2s_yolo",
             "Two-stage YOLOv12 + Cascade" = "pred_2s_yc",
             "YOLOv12 stage 1"  = "pred_yolo", "Cascade stage 1" = "pred_cascade",
             "YOLOv12 F1"       = "yolo_f1",   "Cascade F1"      = "cascade_f1")
REF_METHOD <- "Two-stage YOLOv12"
METHOD_COLORS <- c("Two-stage YOLOv12" = "#7570B3", "Two-stage YOLOv12 + Cascade" = "#3F007D",
                   "YOLOv12 stage 1"   = "#CB181D", "Cascade stage 1" = "#2171B5",
                   "YOLOv12 F1"        = "#FB6A4A", "Cascade F1"      = "#6BAED6")
PLOT_METHODS <- c("Two-stage YOLOv12", "YOLOv12 F1", "Cascade F1")   # covariate/spatial panels

COVARIATES <- c("altitude" = "Altitude (m)", "bottom_depth" = "Bottom depth (m)",
                "ctd_temperature_celsius" = "Temperature (°C)",
                "field_of_view_sq_meter" = "Field of view (m²)",
                "fluorometer_backscatter_ntu" = "Backscatter (NTU)")
# Crowding panel for Fig. 5 (YOLO boxes; images with >= 2 detections, where spacing is defined)
CROWD_PANEL <- c("nn_ratio" = "Nearest-neighbour spacing / box size")
# Detection-structure features (computed in gammod.R from YOLO boxes)
STRUCTURE <- c("med_box_mm" = "Median box size (mm)", "lowconf_frac" = "Share below F1 threshold",
               "overlap_frac" = "Share of overlapping boxes", "nn_ratio" = "NN spacing / box size")
N_BINS     <- 10        # quantile bins per covariate
N_CAL_BINS <- 15        # bins of PREDICTED count for the calibration check
ERROR_SCALE <- "count"  # "count" (pred - true per image) or "density" (per m²)

AUDIT_MAX_ANNOT <- 1    # annotation count at or below which an image can be flagged
AUDIT_MIN_DET   <- 10   # ...if BOTH detectors' F1-threshold counts are at least this
AUDIT_TOP_N     <- 30   # also list the N largest |two-stage YOLOv12 error| images

GRID_DEG   <- 0.1       # spatial grid cell (degrees)
GRID_MIN_N <- 3         # minimum images per cell to map
MORAN_KM   <- 1         # neighbour distance for Moran's I (km)
MORAN_PERM <- 199
N_BOOT     <- 1000

OUT_DIR <- "../data/processed"
FIG_DIR <- "../figures/ms_figures"
PRES_DIR <- "~/saltnoaa/presentations/figures"

# ----------------------------------------------------------------------
# 1. Load and clean
# ----------------------------------------------------------------------
ev <- readRDS(file.path(OUT_DIR, "holdout_eval_2226_v2.rds")) %>%
  filter(dataset == "test_GAM_test") %>%
  mutate(region = as.character(region),
         cruise = as.character(name_cruise_token),
         # sentinel covariates (-99.99) become NA per variable, so an image is
         # dropped only from the panel where its value is missing
         across(any_of(names(COVARIATES)), ~ if_else(.x <= -90, NA_real_, .x)))
stopifnot(all(METHODS %in% names(ev)))
cat("Holdout images:", nrow(ev), "\n")

# ----------------------------------------------------------------------
# 2. Annotation audit
# ----------------------------------------------------------------------
ev <- ev %>%
  mutate(min_det = pmin(yolo_f1, cascade_f1),
         suspect_missing_annotation = n_annotations <= AUDIT_MAX_ANNOT & min_det >= AUDIT_MIN_DET,
         main_error = pred_2s_yolo - n_annotations)

audit <- ev %>%
  filter(suspect_missing_annotation | rank(-abs(main_error)) <= AUDIT_TOP_N) %>%
  arrange(desc(abs(main_error))) %>%
  select(image_id, region, cruise, transect_block_id, density_bin, n_annotations,
         yolo_f1, cascade_f1, pred_2s_yolo, main_error, suspect_missing_annotation,
         any_of(c("class_name", "is_empty")))
cat("\n--- Annotation audit:", sum(ev$suspect_missing_annotation),
    "suspected missing-annotation image(s); top", AUDIT_TOP_N, "|errors| listed ---\n")
print(audit, n = Inf, width = Inf)
write.csv(audit, file.path(OUT_DIR, "holdout_annotation_audit_2226.csv"), row.names = FALSE)

# Share of each method's squared error coming from flagged images
mse_share <- ev %>%
  pivot_longer(all_of(unname(METHODS)), names_to = "col", values_to = "pred") %>%
  mutate(method = names(METHODS)[match(col, METHODS)], se = (pred - n_annotations)^2) %>%
  group_by(region, method) %>%
  summarise(share_sq_error_flagged = sum(se[suspect_missing_annotation]) / sum(se),
            .groups = "drop")
print(mse_share, n = Inf)

# ----------------------------------------------------------------------
# 3. Metrics: overall (with / without flagged images), by density bin, by cruise
# ----------------------------------------------------------------------
boot_all  <- block_bootstrap_metrics(ev, METHODS, ref = REF_METHOD, B = N_BOOT) %>%
  mutate(subset = "all holdout images")
boot_excl <- if (any(ev$suspect_missing_annotation)) {
  block_bootstrap_metrics(ev %>% filter(!suspect_missing_annotation), METHODS,
                          ref = REF_METHOD, B = N_BOOT) %>%
    mutate(subset = "excluding suspected missing annotations (sensitivity)")
} else NULL
metrics_overall <- bind_rows(boot_all, boot_excl)
cat("\n--- Holdout metrics (95% block-bootstrap intervals) ---\n")
print(metrics_overall %>% select(subset, region, method, rmse, rmse_lo, rmse_hi, r2_1to1,
                                 pct_bias_total, bias_lo, bias_hi, d_rmse_vs_ref, d_rmse_lo, d_rmse_hi),
      n = Inf, width = Inf)
write.csv(metrics_overall, file.path(OUT_DIR, "holdout_metrics_with_sensitivity_2226.csv"), row.names = FALSE)

strata_metrics <- function(df, by) {
  df %>%
    pivot_longer(all_of(unname(METHODS)), names_to = "col", values_to = "pred") %>%
    mutate(method = names(METHODS)[match(col, METHODS)], err = pred - n_annotations) %>%
    group_by(region, .data[[by]], method) %>%
    summarise(n_images = n(), true_total = sum(n_annotations),
              rmse = sqrt(mean(err^2)), mae = mean(abs(err)), mean_bias = mean(err),
              pct_bias_total = 100 * sum(err) / pmax(sum(n_annotations), 1),
              .groups = "drop")
}
by_density <- strata_metrics(ev, "density_bin")
by_cruise  <- strata_metrics(ev, "cruise")
write.csv(by_density, file.path(OUT_DIR, "holdout_metrics_by_density_2226.csv"), row.names = FALSE)
write.csv(by_cruise,  file.path(OUT_DIR, "holdout_metrics_by_cruise_2226.csv"),  row.names = FALSE)
cat("\n--- By density bin ---\n"); print(by_density, n = Inf, width = Inf)
cat("\n--- By cruise ---\n");      print(by_cruise,  n = Inf, width = Inf)

# ----------------------------------------------------------------------
# 4. Long error table
# ----------------------------------------------------------------------
err_long <- ev %>%
  pivot_longer(all_of(unname(METHODS)), names_to = "col", values_to = "pred") %>%
  mutate(Method = factor(names(METHODS)[match(col, METHODS)], levels = names(METHODS)),
         count_error   = pred - n_annotations,
         density_error = count_error / field_of_view_sq_meter,
         error = if (ERROR_SCALE == "density") density_error else count_error)
y_lab <- if (ERROR_SCALE == "density") expression(D[pred] - D[true] ~ (m^-2)) else expression(N[pred] - N[true])

# Block-clustered SE of a mean (images in a transect block are not independent)
cl_se <- function(e, cl) {
  ok <- !is.na(e); e <- e[ok]; cl <- cl[ok]
  n <- length(e); if (n < 2) return(NA_real_)
  u <- tapply(e - mean(e), cl, sum); G <- length(u)
  if (G < 2) return(NA_real_)
  sqrt(G / (G - 1) * sum(u^2)) / n
}

# ----------------------------------------------------------------------
# 5. Error vs covariates
# ----------------------------------------------------------------------
bin_covariate <- function(df, x_var, n_bins = N_BINS) {
  df %>%
    filter(!is.na(.data[[x_var]])) %>%
    group_by(region) %>%
    mutate(x_bin = dplyr::ntile(.data[[x_var]], n_bins)) %>%   # quantile bins within region
    group_by(region, Method, x_bin) %>%
    summarise(x_mid = median(.data[[x_var]]), mean_err = mean(error),
              se_err = cl_se(error, transect_block_id), n = n(), .groups = "drop")
}

# Plotting bins: fewer, and built from unique quantile breaks so that tied values
# (many zeros, or the spacing cap of 5) form one bin instead of stacked duplicates.
# Table S2 keeps the 10 equal-count bins of bin_covariate().
N_PLOT_BINS <- 7
bin_for_plot <- function(df, x_var, n_bins = N_PLOT_BINS) {
  df %>% filter(!is.na(.data[[x_var]])) %>%
    group_by(region) %>%
    mutate(x_bin = {
      x <- .data[[x_var]]
      br <- unique(quantile(x, seq(0, 1, length.out = n_bins + 1), na.rm = TRUE, type = 1))
      if (length(br) < 2) 1L else as.integer(cut(x, br, include.lowest = TRUE))
    }) %>%
    group_by(region, Method, x_bin) %>%
    summarise(x_mid = median(.data[[x_var]]), mean_err = mean(error),
              se_err = cl_se(error, transect_block_id), n = n(), .groups = "drop") %>%
    filter(n >= 10)
}

plot_covariate_error <- function(df, x_var, x_label) {
  d  <- df %>% filter(Method %in% PLOT_METHODS)
  bd <- bin_for_plot(d, x_var)
  # dodge proportional to the typical gap between bin centres, per region
  gap <- bd %>% distinct(region, x_bin, x_mid) %>% group_by(region) %>% arrange(x_mid, .by_group = TRUE) %>%
    summarise(g = median(diff(x_mid)), .groups = "drop")
  # fallback when a region has a single bin: a tenth of the covariate's 5-95% range
  rng <- d %>% group_by(region) %>%
    summarise(r = diff(quantile(.data[[x_var]], c(0.05, 0.95), na.rm = TRUE)), .groups = "drop") %>%
    mutate(r = ifelse(is.finite(r) & r > 0, r / 10, 0.05))
  bd <- bd %>% left_join(gap, by = "region") %>% left_join(rng, by = "region") %>%
    mutate(g = ifelse(is.finite(g) & g > 0, g, r),
           # offset: the smaller of a fifth of the bin gap and 1.5% of the covariate's 5-95% range,
           # so points never look as if they sit at different covariate values
           off = pmin(0.2 * g, 0.15 * r),
           x_plot = x_mid + (as.integer(factor(Method, levels = PLOT_METHODS)) - 2) * off)
  ggplot(bd, aes(x = x_plot, y = mean_err, colour = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(ymin = mean_err - 1.96 * se_err, ymax = mean_err + 1.96 * se_err), linewidth = 0.4) +
    geom_point(size = 1.4) +
    facet_wrap(~ region, scales = "free") +
    scale_colour_manual(values = METHOD_COLORS) +
    theme_bw(base_size = 11) +
    labs(x = x_label, y = y_lab, colour = NULL) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

cov_plots <- lapply(names(COVARIATES), function(v) plot_covariate_error(err_long, v, COVARIATES[[v]]))
cov_plots <- c(cov_plots, lapply(names(CROWD_PANEL), function(v)
  plot_covariate_error(err_long %>% filter(n_det_all >= 2), v, CROWD_PANEL[[v]])))

# Image-level error summary for Table S1
image_error_summary <- err_long %>%
  group_by(region, Method) %>%
  summarise(SD = sd(error), MAE = mean(abs(error)), RMSE = sqrt(mean(error^2)),
            Max_Error = error[which.max(abs(error))], .groups = "drop")
write.csv(image_error_summary, file.path(OUT_DIR, "holdout_image_error_summary_2226.csv"), row.names = FALSE)

# Detection-structure panels: images WITH detections only (images without
# detections carry fill-in values that would pile up in one bin)
struct_vars <- intersect(names(STRUCTURE), names(err_long))
# Box size needs one detection; shares and spacing are only informative with >= 2
# (with one detection overlap = 0, spacing = cap and the low-confidence share is 0 or 1).
struct_plots <- lapply(struct_vars, function(v)
  plot_covariate_error(err_long %>% filter(n_det_all >= if (v == "med_box_mm") 1 else 2), v, STRUCTURE[[v]]))

# ----------------------------------------------------------------------
# 5b. Calibration by PREDICTED count
# ----------------------------------------------------------------------
# Binning by the TRUE count makes any shrinking predictor look biased (over at
# the low end, under at the high end) even when it is well calibrated. The fair
# check: within bins of the prediction, does the mean true count match? The
# calibration slope (true ~ pred) is 1 for a calibrated predictor; > 1 means
# predictions are too compressed (high counts under-predicted).
calib <- err_long %>%
  group_by(region, Method) %>%
  mutate(p_bin = dplyr::ntile(pred, N_CAL_BINS)) %>%
  group_by(region, Method, p_bin) %>%
  summarise(pred_mean = mean(pred), true_mean = mean(n_annotations),
            true_se = cl_se(n_annotations, transect_block_id), n = n(), .groups = "drop")
calib_slope <- err_long %>%
  group_by(region, Method) %>%
  summarise(calibration_slope = unname(coef(lm(n_annotations ~ pred))[2]),
            calibration_intercept = unname(coef(lm(n_annotations ~ pred))[1]),
            slope_dense = {i <- pred >= 20
                           if (sum(i) > 5) unname(coef(lm(n_annotations[i] ~ pred[i]))[2]) else NA_real_},
            .groups = "drop")
cat("\n--- Calibration (true ~ predicted); slope_dense uses predictions >= 20 ---\n")
print(calib_slope, n = Inf)
write.csv(calib,       file.path(OUT_DIR, "holdout_calibration_bins_2226.csv"),  row.names = FALSE)
write.csv(calib_slope, file.path(OUT_DIR, "holdout_calibration_slope_2226.csv"), row.names = FALSE)

p_calib <- ggplot(calib %>% filter(Method %in% PLOT_METHODS),
                  aes(pred_mean, true_mean, colour = Method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(alpha = 0.5, show.legend = FALSE) +
  geom_linerange(aes(ymin = pmax(true_mean - 1.96 * true_se, 0), ymax = true_mean + 1.96 * true_se),
                 linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 1.4, show.legend = FALSE) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 25, 50, 100, 200)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 25, 50, 100, 200)) +
  facet_wrap(~ region, scales = "free") +
  scale_colour_manual(values = METHOD_COLORS) +
  theme_bw(base_size = 11) +
  labs(x = "Mean predicted count (bin)", y = "Mean annotated count",
       colour = NULL) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# ----------------------------------------------------------------------
# 6. Spatial: gridded mean error + Moran's I
# ----------------------------------------------------------------------
grid_df <- err_long %>%
  filter(Method %in% PLOT_METHODS) %>%
  mutate(lat_grid = round(latitude / GRID_DEG) * GRID_DEG,
         lon_grid = round(longitude / GRID_DEG) * GRID_DEG) %>%
  group_by(region, Method, lat_grid, lon_grid) %>%
  summarise(mean_cell_err = mean(error), n = n(), .groups = "drop") %>%
  filter(n >= GRID_MIN_N)

lim_sp <- quantile(abs(grid_df$mean_cell_err), 0.95, na.rm = TRUE)
p_space <- ggplot(grid_df, aes(lon_grid, lat_grid, fill = mean_cell_err)) +
  geom_tile(width = GRID_DEG, height = GRID_DEG) +
  facet_wrap(~ Method, ncol = length(PLOT_METHODS)) +
  scale_fill_gradient2(low = "#2C7FB8", mid = "grey92", high = "#E34A33", midpoint = 0,
                       limits = c(-lim_sp, lim_sp), oob = scales::squish,
                       name = "Mean error\n(pred - true)") +
  coord_quickmap() +
  theme_bw(base_size = 11) +
  labs(x = "Longitude", y = "Latitude",
       subtitle = sprintf("%.1f° cells with ≥ %d holdout images", GRID_DEG, GRID_MIN_N)) +
  theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())

haversine_km <- function(lat, lon) {
  r <- pi / 180; la <- lat * r; lo <- lon * r
  dlat <- outer(la, la, "-"); dlon <- outer(lo, lo, "-")
  a <- sin(dlat / 2)^2 + outer(cos(la), cos(la)) * sin(dlon / 2)^2
  a[] <- pmin(1, a)            # a[] keeps the matrix dimensions (pmin alone drops them)
  2 * 6371 * asin(sqrt(a))
}
moran_i <- function(x, W) {
  z <- x - mean(x); S0 <- sum(W)
  if (S0 == 0) return(NA_real_)
  (length(x) / S0) * sum(z * (W %*% z)) / sum(z^2)
}
moran_tab <- purrr::map_dfr(split(ev, ev$region), function(dr) {
  W <- (haversine_km(dr$latitude, dr$longitude) <= MORAN_KM) * 1; diag(W) <- 0
  purrr::map_dfr(names(METHODS), function(m) {
    e <- dr[[METHODS[[m]]]] - dr$n_annotations
    I <- moran_i(e, W)
    perm <- replicate(MORAN_PERM, moran_i(sample(e), W))
    tibble(region = dr$region[1], method = m, moran_I = I,
           p_perm = (1 + sum(perm >= I)) / (MORAN_PERM + 1),
           mean_neighbours = mean(rowSums(W)), radius_km = MORAN_KM)
  })
})
cat("\n--- Moran's I of image residuals (neighbours within", MORAN_KM, "km) ---\n")
print(moran_tab, n = Inf)
write.csv(moran_tab, file.path(OUT_DIR, "holdout_moran_2226.csv"), row.names = FALSE)

# ----------------------------------------------------------------------
# 7. Stability summaries
# ----------------------------------------------------------------------
environmental_stability <- purrr::map_dfr(names(COVARIATES), function(v) {
  bin_covariate(err_long, v) %>%
    group_by(region, Method) %>%
    summarise(Covariate = COVARIATES[[v]], SD_of_Bin_Means = sd(mean_err),
              Peak_Bin_Bias = mean_err[which.max(abs(mean_err))], .groups = "drop")
})
spatial_stability <- err_long %>%
  mutate(lat_grid = round(latitude / GRID_DEG) * GRID_DEG,
         lon_grid = round(longitude / GRID_DEG) * GRID_DEG) %>%
  group_by(region, Method, lat_grid, lon_grid) %>%
  summarise(cell_err = mean(error), n = n(), .groups = "drop") %>%
  filter(n >= GRID_MIN_N) %>%
  group_by(region, Method) %>%
  summarise(Grid_MAE = mean(abs(cell_err)), Grid_SD = sd(cell_err),
            Max_Cell_Bias = cell_err[which.max(abs(cell_err))], n_cells = n(), .groups = "drop")
cat("\n--- Covariate stability ---\n");  print(environmental_stability %>% arrange(region, Covariate, Method), n = Inf)
cat("\n--- Spatial stability ---\n");    print(spatial_stability, n = Inf)
write.csv(environmental_stability, file.path(OUT_DIR, "holdout_covariate_stability_2226.csv"), row.names = FALSE)
write.csv(spatial_stability,       file.path(OUT_DIR, "holdout_spatial_stability_2226.csv"),   row.names = FALSE)

# ----------------------------------------------------------------------
# 8. Assemble and save
# ----------------------------------------------------------------------
cov_assembly <- (wrap_plots(cov_plots, ncol = 2) / p_space) +
  plot_layout(heights = c(2.2, 1.3), guides = "collect") &
  theme(legend.position = "bottom")
cov_assembly <- cov_assembly +
  plot_annotation(
    title = "Holdout error vs. environmental covariates, crowding and space",
    subtitle = sprintf("Quantile-bin means ± 95%% CI (block-clustered SE); error scale: %s",
                       ERROR_SCALE),
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey30")))

# cairo_pdf embeds the Unicode labels (°, ², ≥) that the default pdf device warns about
ggsave(file.path(FIG_DIR, "2226_covariate_error_comparison.pdf"), cov_assembly, width = 11, height = 14,
       device = if (capabilities("cairo")) cairo_pdf else "pdf")
if (dir.exists(path.expand(PRES_DIR))) {
  ggsave(file.path(PRES_DIR, "2226_covariate_error_comparison.png"), cov_assembly,
         width = 11, height = 14, bg = "white")
}
struct_assembly <- (wrap_plots(struct_plots, ncol = 2) / p_calib) +
  plot_layout(heights = c(2, 1.2), guides = "collect") &
  theme(legend.position = "bottom")
struct_assembly <- struct_assembly + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG_DIR, "2226_structure_error_calibration.pdf"), struct_assembly, width = 10, height = 11,
       device = if (capabilities("cairo")) cairo_pdf else "pdf")
ggsave(file.path(FIG_DIR, "2226_structure_error_calibration.png"), struct_assembly, width = 10, height = 11,
       dpi = 300, bg = "white")
ovl_fig <- path.expand("~/Library/CloudStorage/Dropbox/Apps/Overleaf/NEFSC-scallop-detection/ms_figures")
if (dir.exists(ovl_fig)) file.copy(file.path(FIG_DIR, c("2226_structure_error_calibration.pdf",
                                                        "2226_covariate_error_comparison.pdf")),
                                   ovl_fig, overwrite = TRUE)

cat("\nerror.R complete. Tables in", OUT_DIR, "| figure in", FIG_DIR, "\n")
