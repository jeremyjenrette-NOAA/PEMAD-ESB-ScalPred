# ----------------------------
# helper functions
# ----------------------------

parse_geometry_text <- function(txt) {
  vals <- str_extract_all(txt, "-?\\d+\\.?\\d*")[[1]]
  vals <- as.numeric(vals)
  if (length(vals) != 4) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  names(vals) <- c("x1", "y1", "x2", "y2")
  vals
}

lineseg2bb <- function(x1, y1, x2, y2) {
  radius  <- sqrt((y2 - y1)^2 + (x2 - x1)^2) / 2
  centerx <- mean(c(x1, x2))
  centery <- mean(c(y1, y2))
  tlx <- centerx - radius
  tly <- centery - radius
  brx <- centerx + radius
  bry <- centery + radius
  c(tlx = tlx, tly = tly, brx = brx, bry = bry)
}

get_year <- function(imagename) {
  suppressWarnings(as.integer(substr(imagename, 1, 4)))
}

extract_img_date <- function(imagename) {
  parts <- strsplit(imagename, "\\.")[[1]]
  if (length(parts) < 2) return(as.Date(NA))
  out <- suppressWarnings(as.Date(parts[2], format = "%Y%m%d"))
  out
}

get_shift_px <- function(year, img_date,
                         cutoff_date = as.Date("2024-05-11"),
                         shift_px = 55) {
  case_when(
    year == 2024 & !is.na(img_date) & img_date > cutoff_date ~ shift_px,
    TRUE ~ 0
  )
}

get_annotated_side <- function(year) {
  case_when(
    year == 2022 ~ "left",
    year %in% c(2023, 2024) ~ "right",
    TRUE ~ NA_character_
  )
}

plot_stereo_annotation <- function(imagename, gt_df = yr_gt_src) {
  dat <- gt_df %>% filter(imagename == !!imagename)
  if (nrow(dat) == 0) stop("No sourced annotations found for this imagename.")
  
  img_path <- dat$image_path[1]
  g <- read_raster_any(img_path)
  w <- dat$img_width[1]
  h <- dat$img_height[1]
  
  ggplot() +
    annotation_custom(g, xmin = 0, xmax = w, ymin = h, ymax = 0) +
    geom_segment(
      data = dat,
      aes(x = x1, y = y1, xend = x2, yend = y2, color = annotated_side),
      linewidth = 0.8
    ) +
    coord_fixed(xlim = c(0, w), ylim = c(h, 0), expand = FALSE) +
    labs(
      title = imagename,
      subtitle = "Original stereo image with raw line annotation"
    ) +
    theme_void()
}

read_raster_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "png") {
    img <- png::readPNG(path)
  } else if (ext %in% c("jpg", "jpeg")) {
    img <- jpeg::readJPEG(path)
  } else {
    stop("Unsupported image extension: ", ext)
  }
  grid::rasterGrob(img, interpolate = FALSE)
}

plot_split_annotation <- function(imagename, gt_df = yr_gt_src) {
  dat <- gt_df %>% filter(imagename == !!imagename)
  if (nrow(dat) == 0) stop("No sourced annotations found for this imagename.")
  
  img_path <- dat$image_path[1]
  side <- dat$split_side[1]
  
  img_magick <- magick::image_read(img_path)
  info <- magick::image_info(img_magick)
  w <- info$width[[1]]
  h <- info$height[[1]]
  hw <- w / 2
  
  # crop the relevant stereo half
  crop_geom <- if (side == "left") {
    magick::geometry_area(width = hw, height = h, x_off = 0, y_off = 0)
  } else {
    magick::geometry_area(width = hw, height = h, x_off = hw, y_off = 0)
  }
  
  cropped <- magick::image_crop(img_magick, crop_geom)
  cropped_grob <- grid::rasterGrob(as.raster(cropped), interpolate = FALSE)
  
  ggplot() +
    annotation_custom(cropped_grob, xmin = 0, xmax = hw, ymin = h, ymax = 0) +
    geom_segment(
      data = dat,
      aes(x = x1_split, y = y1_split, xend = x2_split, yend = y2_split),
      color = "cyan", linewidth = 0.8
    ) +
    geom_rect(
      data = dat,
      aes(xmin = gt_TLx, ymin = gt_TLy, xmax = gt_BRx, ymax = gt_BRy),
      color = "yellow", fill = NA, linewidth = 0.4
    ) +
    coord_fixed(xlim = c(0, hw), ylim = c(h, 0), expand = FALSE) +
    labs(
      title = imagename,
      subtitle = paste("Split image:", side, "| shifted line (cyan) and GT box (yellow)")
    ) +
    theme_void()
}

# find original stereo image path
find_image_path <- function(imagename, year,
                            dir_2022 = "/Volumes/PortableSSD/saltnoaa/images/2022tr",
                            dir_2024 = "/Volumes/PortableSSD/saltnoaa/images/2024tr",
                            dir_2023 = NA_character_) {
  if (is.na(year)) return(NA_character_)
  
  candidate <- case_when(
    year == 2022 ~ file.path(dir_2022, imagename),
    year == 2023 ~ if (!is.na(dir_2023)) file.path(dir_2023, imagename) else NA_character_,
    year == 2024 ~ file.path(dir_2024, imagename),
    TRUE ~ NA_character_
  )
  
  if (is.na(candidate)) return(NA_character_)
  if (file.exists(candidate)) candidate else NA_character_
}

# read image dimensions without assuming width
read_img_dims <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble(img_width = NA_integer_, img_height = NA_integer_))
  }
  info <- magick::image_info(magick::image_read(path))
  tibble(
    img_width = info$width[[1]],
    img_height = info$height[[1]]
  )
}

read_img_dims_safe <- function(path) {
  out <- tryCatch({
    info <- magick::image_info(magick::image_read(path))
    tibble(
      img_width = as.integer(info$width[[1]]),
      img_height = as.integer(info$height[[1]]),
      img_readable = TRUE
    )
  }, error = function(e) {
    tibble(
      img_width = NA_integer_,
      img_height = NA_integer_,
      img_readable = FALSE
    )
  })
  
  out
}

box_iou <- function(tlx1, tly1, brx1, bry1, tlx2, tly2, brx2, bry2) {
  inter_tlx <- pmax(tlx1, tlx2)
  inter_tly <- pmax(tly1, tly2)
  inter_brx <- pmin(brx1, brx2)
  inter_bry <- pmin(bry1, bry2)
  
  inter_w <- pmax(0, inter_brx - inter_tlx)
  inter_h <- pmax(0, inter_bry - inter_tly)
  inter_area <- inter_w * inter_h
  
  area1 <- pmax(0, brx1 - tlx1) * pmax(0, bry1 - tly1)
  area2 <- pmax(0, brx2 - tlx2) * pmax(0, bry2 - tly2)
  
  union_area <- area1 + area2 - inter_area
  
  ifelse(union_area > 0, inter_area / union_area, 0)
}

match_tp_to_gt_one_image <- function(tp_sub, gt_sub, iou_tol = 1e-6) {
  
  n_tp <- nrow(tp_sub)
  n_gt <- nrow(gt_sub)
  
  if (n_tp == 0 || n_gt == 0) {
    return(tibble())
  }
  
  # IoU matrix: rows = detections, cols = GT
  iou_mat <- matrix(NA_real_, nrow = n_tp, ncol = n_gt)
  
  for (i in seq_len(n_tp)) {
    for (j in seq_len(n_gt)) {
      iou_mat[i, j] <- box_iou(
        tp_sub$TLx[i], tp_sub$TLy[i], tp_sub$BRx[i], tp_sub$BRy[i],
        gt_sub$gt_TLx[j], gt_sub$gt_TLy[j], gt_sub$gt_BRx[j], gt_sub$gt_BRy[j]
      )
    }
  }
  
  # cost = how far recomputed IoU is from recorded iu
  cost_mat <- abs(iou_mat - tp_sub$iu)
  
  # Hungarian needs nrow <= ncol; if not, transpose logic
  if (n_tp <= n_gt) {
    assignment <- solve_LSAP(cost_mat)
    matched_gt_idx <- as.integer(assignment)
    out <- tibble(
      tp_row = seq_len(n_tp),
      gt_row = matched_gt_idx
    )
  } else {
    assignment <- solve_LSAP(t(cost_mat))
    matched_tp_idx <- as.integer(assignment)
    out <- tibble(
      tp_row = matched_tp_idx,
      gt_row = seq_len(n_gt)
    )
  }
  
  out <- out %>%
    mutate(
      Imagename = tp_sub$Imagename[1],
      Detectid = tp_sub$Detectid[tp_row],
      tp_index_local = tp_row,
      gt_index_local = gt_row,
      iu_recorded = tp_sub$iu[tp_row],
      iu_recomputed = map2_dbl(tp_row, gt_row, ~ iou_mat[.x, .y]),
      iu_abs_diff = abs(iu_recorded - iu_recomputed),
      man = gt_sub$man[gt_row]
    )
  
  # usually these should be extremely close if reconstruction is right
  out %>%
    filter(iu_abs_diff <= iou_tol)
}

read_raster_any_safe <- function(path) {
  out <- tryCatch({
    ext <- tolower(tools::file_ext(path))
    if (ext == "png") {
      img <- png::readPNG(path)
    } else if (ext %in% c("jpg", "jpeg")) {
      img <- jpeg::readJPEG(path)
    } else {
      stop("Unsupported image extension")
    }
    grid::rasterGrob(img, interpolate = FALSE)
  }, error = function(e) {
    NULL
  })
  out
}

# build a line of length L centered at (cx, cy) with angle theta
make_centered_line <- function(cx, cy, L, theta) {
  dx <- (L / 2) * cos(theta)
  dy <- (L / 2) * sin(theta)
  tibble(
    x1 = cx - dx,
    y1 = cy - dy,
    x2 = cx + dx,
    y2 = cy + dy
  )
}

plot_detection_length_debug <- function(imagename,
                                        matched_df,
                                        detectid = NULL,
                                        row_id = NULL,
                                        box_color = "yellow",
                                        gt_color = "cyan",
                                        backcalc_color = "magenta") {
  
  dat <- matched_df %>%
    filter(Imagename == !!imagename)
  
  if (!is.null(detectid)) {
    dat <- dat %>% filter(Detectid == !!detectid)
  }
  
  if (!is.null(row_id)) {
    dat <- dat %>% slice(row_id)
  } else {
    dat <- dat %>% slice(1)
  }
  
  if (nrow(dat) == 0) {
    stop("No matching row found.")
  }
  
  # required fields check
  needed <- c("image_path", "split_side", "img_width", "img_height",
              "TLx", "TLy", "BRx", "BRy",
              "x1_split", "y1_split", "x2_split", "y2_split")
  missing_cols <- setdiff(needed, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing columns in matched_df: ", paste(missing_cols, collapse = ", "))
  }
  
  img_path <- dat$image_path[1]
  side <- dat$split_side[1]
  full_w <- dat$img_width[1]
  full_h <- dat$img_height[1]
  half_w <- full_w / 2
  
  # crop the relevant stereo half
  img_magick <- magick::image_read(img_path)
  crop_geom <- if (side == "left") {
    magick::geometry_area(width = half_w, height = full_h, x_off = 0, y_off = 0)
  } else {
    magick::geometry_area(width = half_w, height = full_h, x_off = half_w, y_off = 0)
  }
  cropped <- magick::image_crop(img_magick, crop_geom)
  cropped_grob <- grid::rasterGrob(as.raster(cropped), interpolate = FALSE)
  
  # GT line geometry
  gt_dx <- dat$x2_split[1] - dat$x1_split[1]
  gt_dy <- dat$y2_split[1] - dat$y1_split[1]
  gt_theta <- atan2(gt_dy, gt_dx)
  
  # predicted box center
  pred_cx <- (dat$TLx[1] + dat$BRx[1]) / 2
  pred_cy <- (dat$TLy[1] + dat$BRy[1]) / 2
  
  # mean-method inferred length
  pred_width_px  <- dat$BRx[1] - dat$TLx[1]
  pred_height_px <- dat$BRy[1] - dat$TLy[1]
  pred_length_mean_px <- (pred_width_px + pred_height_px) / 2
  
  # back-calculated line using GT orientation
  backcalc <- make_centered_line(
    cx = pred_cx,
    cy = pred_cy,
    L = pred_length_mean_px,
    theta = gt_theta
  )
  
  # summary label
  gt_length_px <- sqrt(gt_dx^2 + gt_dy^2)
  
  ggplot() +
    annotation_custom(cropped_grob, xmin = 0, xmax = half_w, ymin = full_h, ymax = 0) +
    
    # predicted box
    geom_rect(
      aes(xmin = dat$TLx[1], ymin = dat$TLy[1], xmax = dat$BRx[1], ymax = dat$BRy[1]),
      color = box_color, fill = NA, linewidth = 0.5
    ) +
    
    # GT line
    geom_segment(
      aes(x = dat$x1_split[1], y = dat$y1_split[1],
          xend = dat$x2_split[1], yend = dat$y2_split[1]),
      color = gt_color, linewidth = 0
    ) +
    
    # back-calculated line
    geom_segment(
      data = backcalc,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      color = backcalc_color, linewidth = 0.8
    ) +
    
    # predicted center
    geom_point(
      aes(x = pred_cx, y = pred_cy),
      color = backcalc_color, size = 1.5
    ) +
    
    coord_fixed(xlim = c(0, half_w), ylim = c(full_h, 0), expand = FALSE) +
    labs(
      title = imagename,
      subtitle = paste0(
        "yellow = predicted box | magenta line = predicted line",
        # "\norientation from GT | GT pixels: ", round(gt_length_px, 2),
        "\nGT pixels: ", round(gt_length_px, 2),
        " | Pred pixels: ", round(pred_length_mean_px, 2)
      )
    ) +
    theme_void()
}
###############################################################
# Length evaluation functions
###############################################################
evaluate_length_method <- function(df, pred_col) {
  pred <- df[[pred_col]]
  true <- df$gt_length_mm
  
  tibble(
    method = pred_col,
    bias   = mean(pred - true, na.rm = TRUE),
    MAE    = mean(abs(pred - true), na.rm = TRUE),
    RMSE   = sqrt(mean((pred - true)^2, na.rm = TRUE)),
    R2     = summary(lm(pred ~ true))$r.squared
  )
}