library(tidyverse)
#============================================================#
# Functions: Model-specific Regional GAMs
#============================================================#
get_data <- function(model,
                     model_name,
                     region = NULL,
                     imglvl = FALSE,
                     pooled_regions = c("GB","MAB")) {
  
  # choose suffix
  suffix <- if (imglvl) "img" else "calib"
  
  # -------------------------------------------------
  # Case 1: return pooled dataset across regions
  # -------------------------------------------------
  if (is.null(region) || region == "ALL") {
    
    nm <- paste0(model_name, "_", pooled_regions, "_", suffix)
    
    missing <- nm[!nm %in% names(model)]
    if (length(missing) > 0) {
      stop("Missing datasets: ", paste(missing, collapse = ", "))
    }
    
    df <- dplyr::bind_rows(
      lapply(nm, function(n) model[[n]])
    )
    
    return(df)
  }
  
  # -------------------------------------------------
  # Case 2: return region-specific dataset (original behavior)
  # -------------------------------------------------
  nm <- paste0(model_name, "_", region, "_", suffix)
  
  if (!nm %in% names(model)) {
    stop("Dataset not found: ", nm)
  }
  
  model[[nm]]
}

fit_calibration_gams <- function(model, model_name) {
  
  gb_calib  <- get_data(model, model_name, "GB")
  mab_calib <- get_data(model, model_name, "MAB")
  comb <- rbind(gb_calib, mab_calib)
  
  m_gb <- mgcv::gam(
    y ~ s(conf, bottom_depth, k = 3),
      # s(conf, field_of_view_sq_meter, k = 7),
      # s(latitude, longitude, k = 5),
    family = binomial(),
    data = gb_calib,
    method = "REML"
  )
  
  m_mab <- mgcv::gam(
    y ~ s(conf, bottom_depth, k = 3),
      # s(conf, field_of_view_sq_meter, k = 7),
      # s(conf, latitude, k = 7) +
      # s(conf, longitude, k = 7),
    family = binomial(),
    data = mab_calib,
    method = "REML"
  )
  
  m_comb <- gam(
    y ~ 
      s(conf, bottom_depth, k = 3),
    family = binomial(),
    data = comb,
    method = "REML",
    select = FALSE
  )
  
  list(
    model_name = model_name,
    GB  = m_gb,
    MAB = m_mab,
    comb = m_comb
  )
}

fitfn_calibration_mod <- function(model, model_name, p_detect) {
  
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
  
  reg_comb = rbind(mab_img, gb_img)
  
  reg_comb <- reg_comb %>%
    inner_join(p_detect, by = "image_id")
  
  #------------------------------#
  # GB model
  #------------------------------#
  
  m_gb_fn <- mgcv::gam(
    fn_pres ~
      s(millimeter_per_pixel, k = 7),
    family = binomial(),
    data = gb_img,
    method = "REML"
  )
  
  #------------------------------#
  # MAB model
  #------------------------------#
  
  m_mab_fn <- mgcv::gam(
    fn_pres ~
      s(millimeter_per_pixel, k = 7),
    family = binomial(),
    data = mab_img,
    method = "REML"
  )
  
  #------------------------------#
  # Combined region model
  #------------------------------#
  
  # m_comb <- mgcv::gam(
  #   fn_pres ~ s(predicted_number),
  #   family = binomial(),
  #   data = reg_comb,
  #   method = "REML"
  # )
  
  m_comb <- gam(
    fn_pres ~ 
      s(predicted_number, k=7),
    family = binomial(),
    data = reg_comb,
    method = "REML"
  )
  
  list(
    model_name = model_name,
    GB  = m_gb_fn,
    MAB = m_mab_fn,
    comb = m_comb,
    data = reg_comb
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
    dplyr::mutate(pred_p = stats::predict(gam, newdata = calib_df, type = "response"))
  
  img <- calib_df |>
    dplyr::group_by(image_id) |>
    dplyr::summarise(
      predicted_number = sum(pred_p, na.rm = TRUE),
      true_number      = sum(y, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(region = region, model = model_name)
  
  metrics <- img |>
    dplyr::summarise(
      r2   = summary(stats::lm(true_number ~ predicted_number))$adj.r.squared,
      rmse = sqrt(mean((true_number - predicted_number)^2)),
      .groups = "drop"
    )
  
  list(img = img, metrics = metrics)
}

# ------------------------------------------------------------
# NEW: combine regional GAM predictions and compute pooled metrics
# ------------------------------------------------------------
compute_image_level_counts_pooled <- function(
    calib_df,
    gam_by_region,                 # named list: list(GB = gam_gb, MAB = gam_mab)
    region_col   = "region",
    model_name   = "MyModel",
    pooled_label = "ALL"
) {
  if (!region_col %in% names(calib_df)) stop("calib_df must contain column: ", region_col)
  
  img_list <- lapply(names(gam_by_region), function(r) {
    df_r <- calib_df |> dplyr::filter(.data[[region_col]] == r)
    if (nrow(df_r) == 0) return(NULL)
    
    # predict with the region GAM
    df_r <- df_r |>
      dplyr::mutate(pred_p = stats::predict(gam_by_region[[r]], newdata = df_r, type = "response"))
    
    # image-level sums
    df_r |>
      dplyr::group_by(image_id) |>
      dplyr::summarise(
        predicted_number = sum(pred_p, na.rm = TRUE),
        true_number      = sum(y, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(region = r, model = model_name)
  })
  
  img_all <- dplyr::bind_rows(img_list)
  if (nrow(img_all) == 0) stop("No rows produced. Check region names vs calib_df[[region_col]].")
  
  metrics_pooled <- img_all |>
    dplyr::summarise(
      r2   = summary(stats::lm(true_number ~ predicted_number))$adj.r.squared,
      rmse = sqrt(mean((true_number - predicted_number)^2)),
      n_images = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(region = pooled_label, model = model_name)
  
  metrics_by_region <- img_all |>
    dplyr::group_by(region, model) |>
    dplyr::summarise(
      r2   = summary(stats::lm(true_number ~ predicted_number))$adj.r.squared,
      rmse = sqrt(mean((true_number - predicted_number)^2)),
      n_images = dplyr::n(),
      .groups = "drop"
    )
  
  list(img_all = img_all, metrics_pooled = metrics_pooled, metrics_by_region = metrics_by_region)
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