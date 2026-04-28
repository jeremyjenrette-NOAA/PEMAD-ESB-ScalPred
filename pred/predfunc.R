compile_calibrated_image_df <- function(
    gams,
    pred_detections,
    pred_allimgs,
    meta,
    f1conf_gb = NULL,
    f1conf_mab = NULL,
    lat_split = 40,
    positive_lat_only = TRUE
) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # --------------------------------------------------------------
  # 1. Standardize metadata names to match GAM training variables
  # --------------------------------------------------------------
  meta_std <- meta %>%
    rename(
      Imagename              = IMAGE_NAME,
      altitude               = ALTIMETER_ALTITUDE_METER,
      altitude2              = ALTIMETER_ALTITUDE2_METER,
      stereo_altitude        = STEREO_ALTITUDE_METER,
      backscatter            = FLUOROMETER_BACKSCATTER_NTU,
      bottom_depth           = FATHOMETER_OCEAN_DEPTH_METER,
      cdom                   = FLUOROMETER_CDOM,
      chlorophyll            = FLUOROMETER_CHLOROPHYLL,
      field_of_view_sq_meter = FIELD_OF_VIEW_SQ_METER,
      field_of_view_source   = FIELD_OF_VIEW_SOURCE,
      heading                = VEHICLE_MAGNETIC_HEADING,
      latitude               = SHIP_LATITUDE,
      longitude              = SHIP_LONGITUDE,
      millimeter_per_pixel   = MILLIMETER_PER_PIXEL,
      o2                     = CTD2_DISSOLVED_OXYGEN,
      pitch                  = VEHICLE_PITCH_ANGLE,
      roll                   = VEHICLE_ROLL_ANGLE,
      s                      = CTD_SALINITY,
      t                      = CTD_TEMPERATURE_CELSIUS,
      fluorometer_signal     = FLUOROMETER_SIGNAL,
      datetime               = IMAGE_TIMESTAMP,
      v_depth                = CTD_VEHICLE_DEPTH_METER,
      gear                   = GEAR,
      cruise_id              = CRUISE_ID,
      habcam_pk              = HABCAM_PK
    ) %>%
    mutate(
      bottom_depth = altitude + v_depth
    ) %>%
    mutate(
      Imagename = basename(Imagename)
    ) %>%
    distinct(Imagename, .keep_all = TRUE)
  
  # --------------------------------------------------------------
  # 2. Standardize prediction/detection names
  # --------------------------------------------------------------
  pred_det_std <- pred_detections %>%
    rename(
      detectid  = Detectid,
      imagename = Imagename,
      tlx       = TLx,
      tly       = TLy,
      brx       = BRx,
      bry       = BRy,
      conf      = Conf,
      spname    = Spname
    ) %>%
    mutate(
      imagename = basename(imagename)
    )
  
  # --------------------------------------------------------------
  # 3. Join detections to metadata
  # --------------------------------------------------------------
  pred_master <- pred_det_std %>%
    left_join(
      meta_std %>% rename(imagename = Imagename),
      by = "imagename"
    )
  
  if (positive_lat_only) {
    pred_master <- pred_master %>%
      filter(!is.na(latitude), latitude > 0)
  }
  
  # --------------------------------------------------------------
  # 4. Keep rows usable by the GAMs
  # --------------------------------------------------------------
  pred_master_gam <- pred_master %>%
    filter(
      !is.na(conf),
      !is.na(latitude),
      !is.na(longitude),
      !is.na(bottom_depth),
      !is.na(altitude),
      !is.na(backscatter)
    )
  
  # --------------------------------------------------------------
  # 5. Split by latitude and predict
  # --------------------------------------------------------------
  pred_gb <- pred_master_gam %>%
    filter(latitude >= lat_split) %>%
    mutate(
      p_detection = predict(gams$GB, newdata = ., type = "response"),
      region = "GB"
    )
  
  pred_mab <- pred_master_gam %>%
    filter(latitude < lat_split) %>%
    mutate(
      p_detection = predict(gams$MAB, newdata = ., type = "response"),
      region = "MAB"
    )
  
  pred_calibrated <- bind_rows(pred_gb, pred_mab)
  
  # --------------------------------------------------------------
  # 6. Summarize calibrated detections to image level
  # --------------------------------------------------------------
  img_pred_sum <- pred_calibrated %>%
    group_by(imagename) %>%
    summarise(
      predicted_number = sum(p_detection, na.rm = TRUE),
      raw_detection_number = n(),
      mean_conf = mean(conf, na.rm = TRUE),
      latitude = first(latitude),
      longitude = first(longitude),
      bottom_depth = first(bottom_depth),
      altitude = first(altitude),
      backscatter = first(backscatter),
      field_of_view_sq_meter = first(field_of_view_sq_meter),
      region = first(region),
      .groups = "drop"
    )
  
  # --------------------------------------------------------------
  # 7. F1-thresholded detections and image-level summaries
  # --------------------------------------------------------------
  pred_f1_filtered <- NULL
  img_pred_sum_f1 <- NULL
  img_level_final_f1 <- NULL
  
  if (!is.null(f1conf_gb) && !is.null(f1conf_mab)) {
    
    pred_f1_filtered <- pred_calibrated %>%
      filter(
        (region == "GB"  & conf >= f1conf_gb) |
          (region == "MAB" & conf >= f1conf_mab)
      )
    
    img_pred_sum_f1 <- pred_f1_filtered %>%
      group_by(imagename) %>%
      summarise(
        predicted_number_f1 = n(),
        raw_detection_number_f1 = n(),
        mean_conf_f1 = mean(conf, na.rm = TRUE),
        latitude = first(latitude),
        longitude = first(longitude),
        bottom_depth = first(bottom_depth),
        altitude = first(altitude),
        backscatter = first(backscatter),
        field_of_view_sq_meter = first(field_of_view_sq_meter),
        region = first(region),
        .groups = "drop"
      )
  }
  
  # --------------------------------------------------------------
  # 8. Build image-level metadata for all processed images
  # --------------------------------------------------------------
  all_imgs_meta <- pred_allimgs %>%
    rename(imagename = Imagename) %>%
    mutate(
      imagename = basename(imagename)
    ) %>%
    left_join(
      meta_std %>% rename(imagename = Imagename),
      by = "imagename"
    )
  
  # --------------------------------------------------------------
  # 9. Join image-level GAM predictions back to all processed images
  # --------------------------------------------------------------
  img_level_final <- all_imgs_meta %>%
    left_join(img_pred_sum, by = "imagename", suffix = c("", "_pred")) %>%
    mutate(
      predicted_number = replace_na(predicted_number, 0),
      raw_detection_number = replace_na(raw_detection_number, 0),
      mean_conf = replace_na(mean_conf, 0),
      region = case_when(
        !is.na(region) ~ region,
        !is.na(latitude) & latitude >= lat_split ~ "GB",
        !is.na(latitude) & latitude < lat_split ~ "MAB",
        TRUE ~ NA_character_
      )
    )
  
  if (positive_lat_only) {
    img_level_final <- img_level_final %>%
      filter(!is.na(latitude), latitude > 0)
  } else {
    img_level_final <- img_level_final %>%
      filter(!is.na(latitude))
  }
  
  img_level_final <- img_level_final %>%
    mutate(
      density_raw = raw_detection_number / field_of_view_sq_meter,
      density_cal = predicted_number / field_of_view_sq_meter
    )
  
  # --------------------------------------------------------------
  # 10. Join F1 image-level predictions back to all processed images
  # --------------------------------------------------------------
  if (!is.null(img_pred_sum_f1)) {
    img_level_final_f1 <- all_imgs_meta %>%
      left_join(img_pred_sum_f1, by = "imagename", suffix = c("", "_pred")) %>%
      mutate(
        predicted_number_f1 = replace_na(predicted_number_f1, 0),
        raw_detection_number_f1 = replace_na(raw_detection_number_f1, 0),
        mean_conf_f1 = replace_na(mean_conf_f1, 0),
        region = case_when(
          !is.na(region) ~ region,
          !is.na(latitude) & latitude >= lat_split ~ "GB",
          !is.na(latitude) & latitude < lat_split ~ "MAB",
          TRUE ~ NA_character_
        )
      )
    
    if (positive_lat_only) {
      img_level_final_f1 <- img_level_final_f1 %>%
        filter(!is.na(latitude), latitude > 0)
    } else {
      img_level_final_f1 <- img_level_final_f1 %>%
        filter(!is.na(latitude))
    }
    
    img_level_final_f1 <- img_level_final_f1 %>%
      mutate(
        density_raw = predicted_number_f1 / field_of_view_sq_meter,
        density_f1 = predicted_number_f1 / field_of_view_sq_meter
      )
  }
  
  # --------------------------------------------------------------
  # 11. Return useful outputs
  # --------------------------------------------------------------
  list(
    meta_std = meta_std,
    pred_master = pred_master,
    pred_master_gam = pred_master_gam,
    pred_calibrated = pred_calibrated,
    pred_f1_filtered = pred_f1_filtered,
    img_pred_sum = img_pred_sum,
    img_pred_sum_f1 = img_pred_sum_f1,
    img_level_final = img_level_final,
    img_level_final_f1 = img_level_final_f1
  )
}

plot_calibration_by_depth_envelopes <- function(
    pred_calibrated,
    depth_envelopes,
    conf_col = "conf",
    p_col = "p_detection",
    depth_col = "bottom_depth",
    point_alpha = 0.08,
    point_size = 1.1,
    line_size = 1.2,
    smooth_se = FALSE,
    point_color = "grey60",
    show_bin_n = TRUE,
    model_name
) {
  library(dplyr)
  library(ggplot2)
  library(purrr)
  
  req_cols <- c(conf_col, p_col, depth_col)
  missing_cols <- req_cols[!req_cols %in% names(pred_calibrated)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  env_tbl <- purrr::map_dfr(depth_envelopes, function(x) {
    if (length(x) != 2) stop("Each envelope must have exactly 2 values.")
    tibble(
      min_depth = x[1],
      max_depth = x[2]
    )
  }) %>%
    arrange(min_depth) %>%
    mutate(
      depth_bin = paste0(min_depth, "-", max_depth)
    )
  
  plot_dat <- purrr::map_dfr(seq_len(nrow(env_tbl)), function(i) {
    pred_calibrated %>%
      filter(
        !is.na(.data[[conf_col]]),
        !is.na(.data[[p_col]]),
        !is.na(.data[[depth_col]]),
        .data[[depth_col]] >= env_tbl$min_depth[i],
        .data[[depth_col]] <  env_tbl$max_depth[i]
      ) %>%
      mutate(depth_bin = env_tbl$depth_bin[i])
  })
  
  # enforce ordered legend
  plot_dat <- plot_dat %>%
    mutate(depth_bin = factor(depth_bin, levels = env_tbl$depth_bin))
  
  # counts per envelope for subtitle and/or legend labels
  bin_counts <- plot_dat %>%
    count(depth_bin, name = "n") %>%
    mutate(depth_bin = as.character(depth_bin))
  
  if (show_bin_n) {
    env_tbl <- env_tbl %>%
      left_join(bin_counts, by = "depth_bin") %>%
      mutate(
        n = ifelse(is.na(n), 0, n),
        legend_label = paste0(depth_bin)
      )
    
    plot_dat <- plot_dat %>%
      left_join(env_tbl %>% select(depth_bin, legend_label), by = "depth_bin") %>%
      mutate(
        legend_label = factor(legend_label, levels = env_tbl$legend_label)
      )
    
    color_mapping <- setNames(env_tbl$depth_bin, env_tbl$legend_label)
    color_aes <- "legend_label"
    
    subtitle_text <- paste(
      paste0(env_tbl$depth_bin, ": n images=", scales::comma(env_tbl$n)),
      collapse = "   |   "
    )
  } else {
    color_mapping <- setNames(env_tbl$depth_bin, env_tbl$depth_bin)
    color_aes <- "depth_bin"
    subtitle_text <- NULL
  }
  
  ggplot(plot_dat, aes(x = .data[[conf_col]], y = .data[[p_col]])) +
    geom_point(
      color = point_color,
      alpha = point_alpha,
      size = point_size
    ) +
    geom_smooth(
      aes(color = .data[[color_aes]]),
      method = "gam",
      formula = y ~ s(x),
      se = smooth_se,
      linewidth = line_size
    ) +
    labs(
      x = paste0(model_name," confidence"),
      y = "P(detection)",
      color = "Depth envelope (m)",
      title = "Calibration function by depth envelope",
      subtitle = subtitle_text
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
}

bb2length <- function(tlx, tly, brx, bry, mm_per_pixel = NA) {
  width_px  <- brx - tlx
  height_px <- bry - tly
  
  # if boxes came from lineseg2bb(), width and height should be the same
  length_px <- (width_px + height_px) / 2
  
  length_mm <- length_px * mm_per_pixel
  
  data.frame(
    width_px = width_px,
    height_px = height_px,
    length_px = length_px,
    length_mm = length_mm
  )
}

################################################################################
# From calpred.R 
################################################################################

load_detections <- function(config) {
  
  file <- paste0(config$paths$detections, "detections_", 
                 config$model, "_", config$year, ".csv")
  
  df <- read.csv(file)
  
  if (config$model == "cas") {
    df <- df %>%
      
      # 1. Drop unwanted columns
      select(
        -X3..Unique.Frame.Identifier,
        -Confidence.Pairs.or.Attributes
      ) %>%
      
      # 2. Rename columns to match YOLO format
      rename(
        Detectid  = X..1..Detection.or.Track.id,
        Imagename = X2..Video.or.Image.Identifier,
        TLx       = X4.7..Img.bbox.TL_x,
        TLy       = TL_y,
        BRx       = BR_x,
        BRy       = BR_y.,
        Conf      = X8..Detection.or.Length.Confidence,
        Spname    = X10.11...Repeated.Species
      ) %>%
      
      # 3. Add missing columns (to match YOLO structure)
      mutate(
        img_path      = NA_character_,
        pred_datetime = NA_character_,
        model         = "cascade_rcnn",  # or whatever label you prefer
        split         = NULL,             # ensure it's not present
        year = config$year
      ) %>%
      
      # 4. Reorder columns to match YOLO exactly (excluding split)
      dplyr::select(
        Detectid, Imagename, TLx, TLy, BRx, BRy,
        Conf, Spname, img_path, pred_datetime,
        model, year
      ) %>%
      slice(-1)
  } else if (config$model == "yolo") {
    df <- df
  }
  
  return(df)
}

load_all_data <- function(config) {
  
  detections <- load_detections(config)
  
  if (config$model == "cas") {
  completed <- read.table(paste0(config$paths$completed, "completed_", config$model, "_", config$year, ".txt")) %>%
       rename(Imagename = V1) %>%
    mutate(
      Imagename = basename(gsub("\\\\", "/", Imagename))
    )
  } else {
  completed <- read.table(paste0(config$paths$completed, "completed_", config$model, "_", config$year, ".txt")) %>%
      rename(Imagename = V1)
  }
  
  # metadata <- read.csv(
  #   paste0(config$paths$metadata)
  # )
  # 
  # inventory <- readRDS(
  #   paste0(config$paths$inventory)
  # )
  
  return(list(
    detections = detections,
    completed  = completed
    # metadata   = metadata,
    # inventory  = inventory
  ))
}

process_predictions <- function(gams, data, config, f1conf_gb, f1conf_mab) {
  
  df <- compile_calibrated_image_df(
    gams = gams, 
    pred_detections = data$detections, 
    pred_allimgs = data$completed, 
    meta   = data$metadata,
    f1conf_gb = f1conf_gb,
    f1conf_mab = f1conf_mab
  )
  
  # ---- Add region ----
  df <- df %>%
    mutate(region = ifelse(latitude > config$region_split_lat, "GB", "MAB"))
  
  # ---- Length estimation ----
  if (config$estimate_length) {
    df <- df %>%
      mutate(
        length_px = ((BRx - TLx) + (BRy - TLy)) / 2,
        length_mm = length_px * millimeter_per_pixel
      )
  }
  
  return(df)
}

aggregate_image_level <- function(df) {
  
  df %>%
    group_by(Imagename, region) %>%
    summarise(
      n_raw = n(),
      n_prob = sum(p_detection, na.rm = TRUE),
      density = sum(p_detection) / mean(field_of_view_sq_meter),
      .groups = "drop"
    )
}

plot_length_distribution <- function(df, config) {
  
  p <- ggplot(df, aes(x = length_mm, weight = p_detection)) +
    geom_histogram(bins = 50, fill = "darkgreen") +
    labs(
      x = "Length (mm)",
      y = "Expected individuals"
    ) +
    theme_minimal()
  
  if (config$facet_by_year) {
    p <- p + facet_wrap(~year)
  }
  
  return(p)
}

plot_depth_effects <- function(df, config) {
  
  p <- ggplot(df, aes(x = bottom_depth, y = p_detection)) +
    geom_point(alpha = 0.2) +
    geom_smooth() +
    theme_minimal()
  
  if (config$facet_by_region) {
    p <- p + facet_wrap(~region)
  }
  
  return(p)
}