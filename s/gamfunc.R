library(tidyverse)
#============================================================#
# Functions: Model-specific Regional GAMs
#============================================================#
get_data <- function(model, model_name, region, imglvl = FALSE) {
  if (imglvl) nm <- paste0(model_name, "_", region, "_img") else {
  nm <- paste0(model_name, "_", region, "_calib") }
  
  if (!nm %in% names(model)) {
    stop("Calibration dataset not found: ", nm)
  }
  
  model[[nm]]
}

fit_calibration_gams <- function(model, model_name) {
  
  gb_calib  <- get_data(model, model_name, "GB")
  mab_calib <- get_data(model, model_name, "MAB")
  
  m_gb <- mgcv::gam(
    y ~ s(conf, bottom_depth, k = 7) +
      s(conf, field_of_view_sq_meter, k = 7) +
      s(latitude, longitude, k = 5),
    family = binomial(),
    data = gb_calib,
    method = "REML"
  )
  
  m_mab <- mgcv::gam(
    y ~ s(conf, bottom_depth, k = 7) +
      s(conf, field_of_view_sq_meter, k = 5) +
      s(conf, latitude, k = 7) +
      s(conf, longitude, k = 7),
    family = binomial(),
    data = mab_calib,
    method = "REML"
  )
  
  list(
    model_name = model_name,
    GB  = m_gb,
    MAB = m_mab
  )
}

fitfn_calibration_gams <- function(model, model_name) {
  
  gb_img  <- get_data(model, model_name, "GB", imglvl = TRUE)
  mab_img <- get_data(model, model_name, "MAB", imglvl = TRUE)
  
  # Replace NA auto counts
  gb_img$n_auto[is.na(gb_img$n_auto)] <- 0
  mab_img$n_auto[is.na(mab_img$n_auto)] <- 0
  
  # Define FN presence
  gb_img$fn_pres  <- as.integer(gb_img$n_auto < gb_img$n_manual)
  mab_img$fn_pres <- as.integer(mab_img$n_auto < mab_img$n_manual)
  
  # Density term
  gb_img$auto_density_log  <- log(gb_img$auto_density + 1e-6)
  mab_img$auto_density_log <- log(mab_img$auto_density + 1e-6)
  
  #------------------------------#
  # GB model
  #------------------------------#
  
  # m_gb_fn <- mgcv::gam(
  #   fn_pres ~ 
  #     s(auto_density_log, k = 7) +
  #     s(field_of_view_sq_meter, k = 7),
  #   family = binomial(),
  #   data = gb_img,
  #   method = "REML"
  # )
  
  #------------------------------#
  # MAB model
  #------------------------------#
  
  # m_mab_fn <- mgcv::gam(
  #   fn_pres ~ 
  #     s(auto_density_log, k = 7) +
  #     s(field_of_view_sq_meter, k = 7),
  #   family = binomial(),
  #   data = mab_img,
  #   method = "REML"
  # )
  
  #------------------------------#
  # Combined region model
  #------------------------------#
  
  m_comb <- mgcv::gam(
    fn_pres ~ 
      s(auto_density_log, k = 7) +
      s(field_of_view_sq_meter, k = 7),
    family = binomial(),
    data = rbind(mab_img, gb_img),
    method = "REML"
  )
  
  list(
    model_name = model_name,
    #GB  = m_gb_fn,
    #MAB = m_mab_fn,
    comb = m_comb
  )
}

make_depth_bins <- function(df, depth_col = "bottom_depth", breaks, digits = 0) {
  stopifnot(length(breaks) >= 2)
  
  # auto labels like "36–56"
  lbl <- paste0(
    round(breaks[-length(breaks)], digits), "–",
    round(breaks[-1], digits)
  )
  
  # named midpoints, names match labels
  mids <- setNames((breaks[-length(breaks)] + breaks[-1]) / 2, lbl)
  
  df <- df |>
    mutate(
      depth_bin = cut(
        .data[[depth_col]],
        breaks = breaks,
        include.lowest = TRUE,
        right = TRUE,
        labels = lbl
      )
    )
  
  list(df = df, labels = lbl, mids = mids, breaks = breaks)
}

predict_calibration_by_depth <- function(gam,
                                         calib_df,
                                         depth_breaks,
                                         region,
                                         model_name,
                                         conf_grid = seq(0, 1, by = 0.01)) {
  
  bins <- make_depth_bins(calib_df, breaks = depth_breaks)
  
  lat  <- median(calib_df$latitude,  na.rm = TRUE)
  lon  <- median(calib_df$longitude, na.rm = TRUE)
  fov  <- median(calib_df$field_of_view_sq_meter, na.rm = TRUE)
  
  pred <- tidyr::crossing(
    conf = conf_grid,
    depth_bin = factor(names(bins$mids), levels = bins$labels)
  ) |>
    mutate(
      bottom_depth = unname(bins$mids[as.character(depth_bin)]),
      latitude = lat,
      longitude = lon,
      field_of_view_sq_meter = fov,
      region = region,
      model = model_name
    )
  
  pr <- predict(gam, pred, type = "link", se.fit = TRUE)
  
  pred |>
    mutate(
      fit   = pr$fit,
      se    = pr$se.fit,
      p_hat = plogis(fit),
      p_lo  = plogis(fit - 1.96 * se),
      p_hi  = plogis(fit + 1.96 * se)
    )
}

plot_calibration_by_depth <- function(pred_df, model_name) {
  ggplot(pred_df,
         aes(conf, p_hat, color = depth_bin, fill = depth_bin)) +
    geom_line(linewidth = 1.1) +
    geom_ribbon(aes(ymin = p_lo, ymax = p_hi),
                alpha = 0.15, color = NA) +
    facet_grid(~ region) +
    theme_minimal(base_size = 13) +
    labs(
      title = model_name,
      x = "Detection confidence",
      y = "P(True detection)",
      color = "Depth (m)",
      fill  = "Depth (m)"
    )
}

compute_image_level_counts <- function(calib_df,
                                       gam,
                                       region,
                                       model_name) {
  
  calib_df <- calib_df |>
    mutate(pred_p = predict(gam, type = "response"))
  
  img <- calib_df |>
    group_by(image_id) |>
    summarise(
      predicted_number = sum(pred_p, na.rm = TRUE),
      true_number      = sum(y, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(region = region, model = model_name)
  
  metrics <- img |>
    summarise(
      r2      = summary(lm(true_number ~ predicted_number))$adj.r.squared,
      rmse    = sqrt(mean((true_number - predicted_number)^2)),
      .groups = "drop"
    )
  
  list(img = img, metrics = metrics)
}

plot_image_level_fit <- function(img_df, metrics_df) {
  
  metrics_df <- metrics_df |>
    mutate(
      label = paste0(
        "R² = ", round(r2, 2), "\n",
        "RMSE = ", round(rmse, 2)
      )
    )
  
  ggplot(img_df, aes(predicted_number, true_number)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey40") +
    geom_point(color = "red", size = 1) +
    facet_grid(~ region, scales = "free") +
    geom_text(
      data = metrics_df,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.5, vjust = 1,
      inherit.aes = FALSE
    ) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste(model_name, " - True vs Σ P(detection)",sep=""),
      x = "Predicted (Σ calibrated p)",
      y = "Manual count"
    )
}

compute_zoom_limits <- function(df, x, y, q = 0.9) {
  tibble(
    x_max = quantile(df[[x]], q, na.rm = TRUE),
    y_max = quantile(df[[y]], q, na.rm = TRUE)
  )
}

plot_image_level_fit_zoom <- function(img_df,
                                      metrics_df,
                                      model_name,
                                      zoom_q = 0.9) {
  
  metrics_df <- metrics_df |>
    mutate(
      label = paste0(
        "R² = ", round(r2, 2), "\n",
        "RMSE = ", round(rmse, 2)
      )
    )
  
  zoom_limits <- img_df |>
    group_by(region) |>
    summarise(
      x_max = quantile(predicted_number, zoom_q, na.rm = TRUE),
      y_max = quantile(true_number, zoom_q, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(img_df, aes(predicted_number, true_number)) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "grey40"
    ) +
    geom_point(color = "red", size = 1, alpha = 0.8) +
    facet_wrap(~ region, scales = "free") +
    geom_text(
      data = metrics_df,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.4, vjust = 1.1,
      inherit.aes = FALSE
    ) +
    coord_cartesian(
      xlim = c(0, max(zoom_limits$x_max)),
      ylim = c(0, max(zoom_limits$y_max))
    ) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste(model_name, "– True vs Σ P(detection)", sep = " "),
      subtitle = paste(zoom_q*100,"%"," of images", sep = ""),
      x = "Predicted (Σ calibrated p)",
      y = "Manual count"
    )
}

save_cal <- function(p, id, outdir = "~/Downloads",
                       width = 6, height = 6) {
  
  ggsave(
    plot = p,
    filename = file.path(outdir, paste0(id, ".jpg")),
    width = width,
    height = height
  )
}
