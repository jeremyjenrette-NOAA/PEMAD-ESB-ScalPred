library(tidyverse)
library(janitor)
library(lubridate)

#------------------------------------------------------------#
# Standardize image identifiers across all data sources
#------------------------------------------------------------#
std_image_id <- function(x) {
  x |>
    as.character() |>
    na_if("") |>
    str_trim() |>
    str_remove("\\.[A-Za-z0-9]+$")  # remove file extension
}

#------------------------------------------------------------#
# Parse timestamp embedded in HabCam image filenames
#------------------------------------------------------------#
parse_timestamp_from_imagename <- function(x, tz = "UTC") {
  
  date_str <- str_extract(x, "\\.\\d{8}\\.") |> str_remove_all("\\.")
  time_str <- str_extract(x, "\\.\\d{9}\\.") |> str_remove_all("\\.")
  
  ymd_hms(
    paste0(
      date_str, " ",
      str_sub(time_str, 1, 2), ":",
      str_sub(time_str, 3, 4), ":",
      str_sub(time_str, 5, 6), ".",
      str_sub(time_str, 7, 9)
    ),
    tz = tz
  )
}

#============================================================#
# Build by-detection and by-image datasets for calibration
#============================================================#
build_detection_tables <- function(
    mt,           # manual annotations (any model version)
    at,           # automatic detections (any model version)
    allimg,       # complete image metadata (constant)
    meta = NULL   # optional raw metadata (unused if allimg complete)
) {
  
  #-------------------------------#
  # 1) Clean and standardize inputs
  #-------------------------------#
  
  mt <- mt |>
    clean_names() |>
    mutate(
      image_id = coalesce(
        std_image_id(image_name),
        std_image_id(image_url_id)
      ),
      spname = tolower(spname)
    )
  
  at <- at |>
    clean_names() |>
    mutate(
      image_id = std_image_id(imagename),
      spname   = tolower(spname)
    ) |>
    # keep only images that appear in manual annotations
    filter(image_id %in% mt$image_id)
  
  #-------------------------------#
  # 2) Prepare image-level metadata
  #-------------------------------#
  
  meta_img <- allimg |>
    clean_names() |>
    mutate(
      image_id        = std_image_id(imagename),
      image_timestamp = parse_timestamp_from_imagename(imagename)
    )
  
  #-------------------------------#
  # 3) Estimate field of view from altitude
  #    (median k approach)
  #-------------------------------#
  
  meta_img <- meta_img |>
    mutate(
      field_of_view_sq_meter = 0.1802864 * altitude^2 # calculated from separate metadata sheet
    )
  
  #-------------------------------#
  # 4) Compute bounding-box geometry
  #-------------------------------#
  
  at <- at |>
    mutate(
      box_width  = abs(b_rx - t_lx),
      box_height = abs(b_ry - t_ly),
      box_area   = box_width * box_height
    )
  
  #-------------------------------#
  # 5) Build by-image summaries
  #-------------------------------#
  
  man_img <- mt |>
    count(image_id, name = "n_manual")
  
  auto_img <- at |>
    count(image_id, name = "n_auto")
  
  valid_image_ids <- union(man_img$image_id, auto_img$image_id)
  
  img_df <- meta_img |>
    filter(image_id %in% valid_image_ids) |>
    left_join(man_img,  by = "image_id") |>
    left_join(auto_img, by = "image_id") |>
    mutate(
      n_manual   = replace_na(n_manual, 0L),
      n_auto_raw = replace_na(n_auto, 0L),
      validated  = image_id %in% man_img$image_id,
      
      # densities (per m^2)
      man_density  = n_manual   / field_of_view_sq_meter,
      auto_density = n_auto_raw / field_of_view_sq_meter
    )
  
  #-------------------------------#
  # 6) Build by-detection calibration table
  #-------------------------------#
  
  calib_df <- at |>
    mutate(
      y = as.integer(truedetect)  # TP = 1, FP = 0
    ) |>
    left_join(
      meta_img |>
        select(
          image_id,
          latitude,
          longitude,
          bottom_depth,
          chlorophyll,
          cdom,
          altitude,
          field_of_view_sq_meter,
          internal_ph,
          t,
          s,
          heading,
          o2,
          v_depth,
          therm,
          pitch,
          roll
        ),
      by = "image_id"
    ) |>
    drop_na(conf, y)
  
  #-------------------------------#
  # 7) Return structured output
  #-------------------------------#
  
  list(
    calib_df = calib_df,
    img_df   = img_df
  )
}

