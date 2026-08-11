load("../data/raw/allimg_2026.RData") # allimg

library(dplyr)
library(lubridate)

allmeta = readRDS("../data/raw/combined_years_master_img_ann.rds")
allmeta = allmeta$img

allmeta_standardized <- allmeta %>%
  # 1. Convert all column names to lowercase
  rename_with(tolower) %>%
  
  # 2. Rename the specific columns to match dat_split
  rename(
    altitude  = altimeter_altitude_meter,
    latitude  = ship_latitude,
    longitude = ship_longitude,
    imagename = image_name
  ) %>%
  
  # 3. Create bottom_depth and extract the year
  mutate(
    bottom_depth = altitude + ctd_vehicle_depth_meter,
    # Extracts year; parses automatically if image_timestamp is a standard date string or POSIXct
    year         = year(ymd_hms(image_timestamp)) 
  ) %>%
  select(-c(altimeter_correlation_factor, altimeter_correlation_energy, 
            altimeter_internal_temperature_celsius, dvl_altitude_meter,
            vehicle_latitude, vehicle_longitude, ctd_dissolved_oxygen,
            ctd_electrical_conductivity)) 


colnames(allimg)
colnames(allmeta_standardized)


library(dplyr)

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

allimg_standardized <- allimg %>%
  # 1. Drop the therm column
  select(-therm) %>%
  
  # 2. Rename columns to match allmeta_standardized exactly
  rename(
    latitude                    = Latitude,
    longitude                   = Longitude,
    vehicle_magnetic_heading    = Heading,
    vehicle_pitch_angle         = Pitch,
    vehicle_roll_angle          = Roll,
    altitude                    = Altitude,
    altimeter_altitude2_meter   = Alt2,
    ctd_vehicle_depth_meter     = V_Depth,
    ctd_salinity                = s,
    ctd_temperature_celsius     = t,
    ctd2_dissolved_oxygen       = o2,
    fluorometer_cdom            = cdom,
    fluorometer_chlorophyll     = chlorophyll,
    fluorometer_backscatter_ntu = backscatter,
    seafet_internal_ph          = internal_ph,
    seafet_external_ph          = external_ph,
    bottom_depth                = Bottom_Depth
  ) %>%
  
  # 3. Add constants and vectorize-friendly calculations
  mutate(
    year = 2026,
    imagename = sub("\\.tif$", ".png", imagename),
    # Calculate millimeter per pixel using Similar Triangles: (Altitude / FocalLength) * SensorPixelSize * 1000mm/m
    # Using the constants from your FOV function (16mm focal length * 1.33 water refraction index)
    millimeter_per_pixel = (altitude / (16 * 0.00133)) * 0.00000586 * 1000,
    
    # Optional: Fill in missing fields from allmeta_standardized that don't exist in allimg
    field_of_view_source = "calculated"
  ) %>%
  
  # 4. Calculate field of view row-by-row (required because FOV uses matrix math)
  rowwise() %>%
  mutate(
    field_of_view_sq_meter = FOV(altitude, vehicle_roll_angle, vehicle_pitch_angle)
  ) %>%
  ungroup() # Important: remove rowwise grouping when done

# 1. Bind rows (automatically pads missing columns with NA) and filter coordinates
all_combined_standardized <- bind_rows(allimg_standardized, allmeta_standardized) %>%
  filter(
    # Ensure coordinates exist
    !is.na(latitude),
    !is.na(longitude),
    
    # Generous bounding box for the US Eastern Seaboard / NW Atlantic
    # (Approx. Florida up to Maine/Canada border, out to the shelf edge)
    latitude >= 25.0,
    latitude <= 48.0,
    longitude >= -82.0,
    longitude <= -64.0
  )

# Verify the dimensions and new columns
dim(all_combined_standardized)
colnames(all_combined_standardized)


################################################################################
library(dplyr)
library(purrr)

# 1. Define the directory and identify the files for each model
completed_dir <- "../data/raw/completed"

cas_files  <- list.files(completed_dir, pattern = "completed_cas_.*\\.txt$", full.names = TRUE)
yolo_files <- list.files(completed_dir, pattern = "completed_yolo_.*\\.txt$", full.names = TRUE)

# 2. Helper function to read paths, extract just the image name, and find uniques
get_unique_basenames <- function(files) {
  files %>%
    # Read all lines from the list of files
    map(readLines, warn = FALSE) %>%
    unlist() %>%
    # Remove everything up to the last forward or backward slash to isolate the filename
    gsub(".*[/\\\\]", "", .) %>%
    # Keep only unique filenames
    unique()
}

# 3. Extract the unique image names for each model
cas_images  <- get_unique_basenames(cas_files)
yolo_images <- get_unique_basenames(yolo_files)

# 4. Find the master list of images processed by BOTH models
master_completed_images <- intersect(cas_images, yolo_images)

# Optional: Convert to a dataframe if you plan to merge it with allimg_standardized later
master_completed_df <- data.frame(imagename = master_completed_images, stringsAsFactors = FALSE)

library(dplyr)
library(readr)
library(purrr)

# 1. Locate all detection CSV files
det_dir <- "../data/raw/detections"
cas_det_files  <- list.files(det_dir, pattern = "detections_cas_.*\\.csv$", full.names = TRUE)
yolo_det_files <- list.files(det_dir, pattern = "detections_yolo_.*\\.csv$", full.names = TRUE)

# 2. Process YOLO Detections
yolo_dets <- yolo_det_files %>%
  map_df(~read_csv(.x, show_col_types = FALSE)) %>%
  # Select only the needed columns and drop the rest (Detectid, img_path, etc.)
  select(
    imagename = Imagename, # Standardize capitalization for joining
    TLx, TLy, BRx, BRy, Conf, Spname, model
  ) %>%
  # Calculate bounding box area
  mutate(
    boxsize = abs((BRx - TLx) * (BRy - TLy))
  ) %>%
  # Filter to only images present in the master list
  filter(imagename %in% master_completed_images) %>%
  # Attach the physical metadata
  inner_join(all_combined_standardized, by = "imagename")

# 3. Process Cascade Detections
cas_dets <- cas_det_files %>%
  # Read VIAME format: no headers, ignore '#' comment lines
  map_df(~read_csv(.x, col_names = FALSE, comment = "#", show_col_types = FALSE)) %>%
  # Rename the default X columns to match the YOLO structure
  select(
    imagename = X2,
    TLx       = X4,
    TLy       = X5,
    BRx       = X6,
    BRy       = X7,
    Conf      = X8,
    Spname    = X10
  ) %>%
  # Add the model column explicitly and calculate bounding box area
  mutate(
    model = "Cascade",
    boxsize = abs((BRx - TLx) * (BRy - TLy))
  ) %>%
  # Filter to only images present in the master list
  filter(imagename %in% master_completed_images) %>%
  # Attach the physical metadata
  inner_join(all_combined_standardized, by = "imagename")


#-----------

library(dplyr)
library(stringr)

# 1. Create the master image-level dataframe
# This contains ALL processed images (with or without detections) and their metadata
inference_img_df <- master_completed_df %>%
  left_join(all_combined_standardized, by = "imagename") %>%
  mutate(
    # Define region spatially and set as a factor
    region = if_else(longitude >= -71, "GB", "MAB") %>% factor(levels = c("MAB", "GB")),
    
    # Create an image_id (stripping the extension) to match downstream aggregation joins
    image_id = sub("\\.png$", "", imagename),
    
    # Optional: tag the dataset as inference so holdout functions process it correctly
    dataset = "inference" 
  )

# 2. Append the region and image_id to your detection dataframes
yolo_dets_inf <- yolo_dets %>%
  mutate(
    region = if_else(longitude >= -71, "GB", "MAB") %>% factor(levels = c("MAB", "GB")),
    image_id = sub("\\.png$", "", imagename)
  ) %>%
  rename_with(tolower)

cas_dets_inf <- cas_dets %>%
  mutate(
    region = if_else(longitude >= -71, "GB", "MAB") %>% factor(levels = c("MAB", "GB")),
    image_id = sub("\\.png$", "", imagename)
  ) %>%
  rename_with(tolower)

# 3. Structure into the expected pipeline list format
# This perfectly mimics the structure expected by evaluate_pr_models() and test_calibration_gams()
inference_models <- list(
  `YOLOv12` = list(
    img = inference_img_df,
    det = yolo_dets_inf
  ),
  `Cascade R-CNN` = list(
    img = inference_img_df,
    det = cas_dets_inf
  )
)

saveRDS(inference_models, file = "../data/processed/eval_inference_2226.rds")
