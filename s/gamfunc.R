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
  box  <- median(calib_df$boxsize,  na.rm = TRUE)
  lat  <- median(calib_df$latitude,  na.rm = TRUE)
  lon  <- median(calib_df$longitude, na.rm = TRUE)
  fov  <- median(calib_df$field_of_view_sq_meter, na.rm = TRUE)
  altitude <- median(calib_df$altitude, na.rm = TRUE)
  backscatter <- median(calib_df$backscatter, na.rm = TRUE)
  
  pred <- tidyr::crossing(
    conf = conf_grid,
    depth_bin = factor(names(bins$mids), levels = bins$labels)
  ) |>
    mutate(
      bottom_depth = unname(bins$mids[as.character(depth_bin)]),
      latitude = lat,
      longitude = lon,
      boxsize = box,
      field_of_view_sq_meter = fov,
      backscatter = backscatter,
      altitude = altitude,
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
                                       model_name,
                                       f1_thresh,
                                       use_strat = FALSE,
                                       dat_split,
                                       conf_col = "conf") {
  
  # --- apply stratified GAM training filter ---
  if (use_strat) {
    
    # ensure consistent naming
    dat_split$imagename <- basename(dat_split$imagename)
    
    # join + filter helper
    filter_strat <- function(df) {
      df$imagename <- basename(df$imagename)
      
      df2 <- merge(
        df,
        dat_split[, c("imagename", "stratum", "is_test_gam_test")],
        by = "imagename",
        all.x = TRUE
      ) %>%
        mutate(is_test_gam_test = if_else(is_test_gam_test == "True", TRUE, FALSE))
      
      # keep ONLY GAM test set
      df2 <- df2[df2$is_test_gam_test == TRUE, ]
      
      df2$stratum = factor(df2$stratum)
      
      return(df2)
    }
    
    calib_df  <- filter_strat(calib_df)
  }
  
  calib_df <- calib_df |>
    dplyr::mutate(
      pred_p  = stats::predict(gam, newdata = calib_df, type = "response"),
      conf_val = .data[[conf_col]]
    )
  
  img <- calib_df |>
    dplyr::group_by(imagename) |>
    dplyr::summarise(
      raw_detection_number = dplyr::n(),
      predicted_f1_number  = sum(conf_val >= f1_thresh, na.rm = TRUE),
      predicted_number     = sum(pred_p, na.rm = TRUE),
      true_positive        = sum(y, na.rm = TRUE),  # keep this for diagnostics
      .groups = "drop"
    ) |>
    dplyr::left_join(
      dat_split |>
        dplyr::select(imagename, n_annotations),
      by = "imagename"
    ) |>
    dplyr::mutate(
      true_number = n_annotations,
      false_negative = true_number - true_positive
    ) |>
    dplyr::mutate(region = region, model = model_name)
  
  metrics <- img |>
    dplyr::summarise(
      # calibrated GAM counts
      r2_calibrated   = summary(stats::lm(true_number ~ predicted_number))$adj.r.squared,
      rmse_calibrated = sqrt(mean((true_number - predicted_number)^2)),
      
      # F1-threshold counts
      r2_f1   = summary(stats::lm(true_number ~ predicted_f1_number))$adj.r.squared,
      rmse_f1 = sqrt(mean((true_number - predicted_f1_number)^2))
    )
  
  list(img = img, metrics = metrics)
}

compute_combined_metrics <- function(img_df, region_name = "GB_MAB") {
  
  img_df |>
    dplyr::summarise(
      r2_calibrated   = summary(stats::lm(true_number ~ predicted_number))$adj.r.squared,
      rmse_calibrated = sqrt(mean((true_number - predicted_number)^2)),
      r2_f1           = summary(stats::lm(true_number ~ predicted_f1_number))$adj.r.squared,
      rmse_f1         = sqrt(mean((true_number - predicted_f1_number)^2))
    ) |>
    dplyr::mutate(region = region_name)
}


plot_image_level_fit <- function(img_df,
                                 metrics_df,
                                 model_name = "",
                                 include_f1 = FALSE,
                                 free_scales = TRUE,
                                 plot_title) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Build plotting dataframe
  if (include_f1) {
    plot_df <- img_df |>
      dplyr::select(region, true_number, predicted_number, predicted_f1_number) |>
      tidyr::pivot_longer(
        cols = c(predicted_number, predicted_f1_number),
        names_to = "method",
        values_to = "predicted_value"
      ) |>
      dplyr::mutate(
        method = dplyr::recode(
          method,
          predicted_number    = "GAM Calibrated",
          predicted_f1_number = "F1 cutoff"
        )
      )
    
    metrics_plot <- metrics_df |>
      dplyr::transmute(
        region,
        method = "GAM Calibrated",
        r2 = r2_calibrated,
        rmse = rmse_calibrated
      ) |>
      dplyr::bind_rows(
        metrics_df |>
          dplyr::transmute(
            region,
            method = "F1 cutoff",
            r2 = r2_f1,
            rmse = rmse_f1
          )
      ) |>
      dplyr::mutate(
        label = paste0(
          "R² = ", round(r2, 2), "\n",
          "RMSE = ", round(rmse, 2)
        )
      )
    
    p <- ggplot(plot_df, aes(predicted_value, true_number)) +
      geom_abline(
        slope = 1, intercept = 0,
        linetype = "dashed", color = "grey40"
      ) +
      geom_point(color = "red", size = 1) +
      geom_text(
        data = metrics_plot,
        aes(x = -Inf, y = Inf, label = label),
        hjust = -0.5, vjust = 1.5,
        inherit.aes = FALSE
      ) +
      facet_grid(
        rows = vars(method),
        cols = vars(region),
        scales = if (free_scales) "free" else "fixed"
      ) +
      theme_minimal(base_size = 13) +
      labs(
        title = plot_title,
        x = "Predicted count",
        y = "Manual count"
      )
    
  } else {
    plot_df <- img_df |>
      dplyr::mutate(method = "GAM Calibrated",
                    predicted_value = predicted_number)
    
    metrics_plot <- metrics_df |>
      dplyr::transmute(
        region,
        method = "GAM Calibrated",
        r2 = r2_calibrated,
        rmse = rmse_calibrated
      ) |>
      dplyr::mutate(
        label = paste0(
          "R² = ", round(r2, 2), "\n",
          "RMSE = ", round(rmse, 2)
        )
      )
    
    p <- ggplot(plot_df, aes(predicted_value, true_number)) +
      geom_abline(
        slope = 1, intercept = 0,
        linetype = "dashed", color = "grey40"
      ) +
      geom_point(color = "red", size = 1) +
      geom_text(
        data = metrics_plot,
        aes(x = -Inf, y = Inf, label = label),
        hjust = -0.5, vjust = 1,
        inherit.aes = FALSE
      ) +
      facet_grid(
        ~ region,
        scales = if (free_scales) "free" else "fixed"
      ) +
      theme_minimal(base_size = 13) +
      labs(
        title = plot_title,
        x = "Predicted (Σ calibrated p)",
        y = "Manual count"
      )
  }
  
  p
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
                                      zoom_q = 0.9,
                                      plot_title,
                                      model_label) {
  
  metrics_df <- metrics_df |>
    mutate(
      label = paste0(
        "R² = ", round(r2_calibrated, 2), "\n",
        "RMSE = ", round(rmse_calibrated, 2)
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
      hjust = -0.4, vjust = 1.5,
      inherit.aes = FALSE
    ) +
    geom_blank(
      data = zoom_limits,
      aes(x = x_max, y = y_max)
    ) +
    theme_bw(base_size = 13) +
    labs(
      title = plot_title, 
      subtitle = model_label,
      # subtitle = paste(zoom_q*100,"%"," of images", sep = ""),
      x = "Σ p(detection)",
      y = "Manual count"
    ) + 
    theme(plot.title = element_text(size = 17),
          strip.text = element_text(size = 12, face = "bold", color = "black"))
}


save_cal <- function(p,
                     id,
                     outdir = "~/Downloads",
                     format = "jpg",
                     width = 6,
                     height = 6,
                     dpi = 300) {
  
  format <- tolower(format)
  allowed <- c("jpg", "jpeg", "png", "pdf")
  
  if (!format %in% allowed) {
    stop("format must be one of: jpg, jpeg, png, pdf")
  }
  
  filename <- file.path(outdir, paste0(id, ".", format))
  
  ggsave(
    plot = p,
    filename = filename,
    width = width,
    height = height,
    dpi = dpi
  )
}
