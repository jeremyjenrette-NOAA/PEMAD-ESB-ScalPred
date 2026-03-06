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
      image_id = std_image_id(image_name),
      spname = tolower(spname)
    ) |>
    filter(spname == 'scallop')
  
  at <- at |>
    clean_names() |>
    mutate(
      image_id = std_image_id(image_name),
      spname   = tolower(spname)
    ) |>
    filter(spname == "scallop") |>
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
  # 3) Estimate field of view from altitude, roll, and pitch
  #-------------------------------#
  if (!("field_of_view_sq_meter" %in% colnames(meta_img))) {
  meta_img$field_of_view_sq_meter <- mapply(
    FOV,
    altitude = meta_img$altitude,
    roll     = meta_img$roll,
    pitch    = meta_img$pitch
  )
  }
  #-------------------------------#
  # 4) Compute bounding-box geometry
  #-------------------------------#
  
  # at <- at |>
  #   mutate(
  #     box_width  = abs(b_rx - t_lx),
  #     box_height = abs(b_ry - t_ly),
  #     box_area   = box_width * box_height
  #   )
  
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
          millimeter_per_pixel,
          backscatter,
          # therm,
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

# calculate FOV
FOV <- function(altitude, roll, pitch) {
  
  # constants
  DTOR <- pi / 180
  focalLength <- 16 * 0.00133
  PIXEL_SIZE <- 0.00000586
  
  # trig
  sP <- sin(pitch * DTOR); cP <- cos(pitch * DTOR)
  sR <- sin(roll  * DTOR); cR <- cos(roll  * DTOR)
  
  # rotation matrix
  m <- matrix(0, nrow = 3, ncol = 3)
  m[1,1] <- cP;        m[1,2] <- 0;   m[1,3] <- -sP
  m[2,1] <- sP*sR;     m[2,2] <- cR;  m[2,3] <- cP*sR
  m[3,1] <- sP*cR;     m[3,2] <- -sR; m[3,3] <- cP*cR
  
  # image corners in sensor coords
  ulX <- -1936/2 * PIXEL_SIZE; ulY <-  1216/2 * PIXEL_SIZE
  urX <- -ulX;                 urY <-  ulY
  llX <-  ulX;                 llY <- -ulY
  lrX <- -ulX;                 lrY <- -ulY
  
  # store corners in order (ul, ll, lr, ur)
  px <- c(ulX, llX, lrX, urX)
  py <- c(ulY, llY, lrY, urY)
  
  # project rays through rotation matrix
  X <- numeric(4); Y <- numeric(4); Z <- numeric(4)
  for (k in 1:4) {
    X[k] <- m[1,1]*px[k] + m[1,2]*py[k] + m[1,3]*(-focalLength)
    Y[k] <- m[2,1]*px[k] + m[2,2]*py[k] + m[2,3]*(-focalLength)
    Z[k] <- m[3,1]*px[k] + m[3,2]*py[k] + m[3,3]*(-focalLength)
  }
  
  # intersect with seabed plane at given altitude
  X <- X * (altitude / Z)
  Y <- Y * (altitude / Z)
  
  # polygon area (shoelace)
  area <- 0
  j <- 4
  for (i in 1:4) {
    area <- area + (X[j] + X[i]) * (Y[j] - Y[i])
    j <- i
  }
  
  -area / 2
}
#============================================================#
# Functions: Subset model and region
#============================================================#
structure_by_region <- function(out,
                                model_name,
                                region_lat_cutoff = 40,
                                save_dir = "../data/processed",
                                save_rdata = TRUE) {
  
  stopifnot(is.list(out), "calib_df" %in% names(out), "img_df" %in% names(out))
  stopifnot(is.character(model_name), length(model_name) == 1)
  
  # ---- 1) add region + region-specific model_name ----
  out2 <- out
  out2$calib_df <- out2$calib_df %>%
    mutate(
      region = if_else(latitude >= region_lat_cutoff, "GB", "MAB") %>% factor(levels = c("MAB", "GB")),
      model_name = paste0(model_name, "_", region)
    )
  
  out2$img_df <- out2$img_df %>%
    mutate(
      region = if_else(latitude >= region_lat_cutoff, "GB", "MAB") %>% factor(levels = c("MAB", "GB")),
      model_name = paste0(model_name, "_", region)
    )
  
  # ---- 2) split ----
  GB_calib  <- out2$calib_df %>% filter(region == "GB")
  GB_img    <- out2$img_df   %>% filter(region == "GB")
  MAB_calib <- out2$calib_df %>% filter(region == "MAB")
  MAB_img   <- out2$img_df   %>% filter(region == "MAB")
  
  # ---- 3) return as requested: c(out, ...) but prefixed and named ----
  res <- c(
    out2,
    list(
      GB_calib  = GB_calib,
      GB_img    = GB_img,
      MAB_calib = MAB_calib,
      MAB_img   = MAB_img
    )
  )
  
  # prefix every element name with model_name (including out2 pieces)
  names(res) <- paste0(model_name, "_", names(res))
  
  # ---- 4) save (optional) ----
  if (isTRUE(save_rdata)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    
    # save_path <- file.path(save_dir, paste0(model_name, ".RData"))
    
    # Save each object in the list as its own named object in the RData
    # list2env(res, envir = environment())
    # save(list = names(res), file = save_path, envir = environment())
    
    file <- file.path(save_dir, paste0(model_name, ".RData"))
    
    tmp <- list()
    tmp[[model_name]] <- res
    
    save(list = model_name, file = file, envir = list2env(tmp, parent = emptyenv()))
    
    invisible(file)
  }
  
  res
}

