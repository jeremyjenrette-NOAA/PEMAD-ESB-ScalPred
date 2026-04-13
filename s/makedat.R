library(dplyr)
library(purrr)
library(tibble)
source("./datfunc.R")
# load("../data/raw/allimg_2024.RData")  # loads `allimg for 2024`, complete metadata
# metagb = read.csv("../data/raw/metagb22.csv")
# metamab = read.csv("../data/raw/metama22.csv")
# missing = read.csv("../data/raw/missingmeta22_1.csv")
# meta22 = rbind(metagb, metamab, missing)
# meta22 <- meta22[!duplicated(meta22$IMAGENAME), ]
meta2224 <- read.csv("../data/processed/metadata2224.csv")
metastar = read.csv("../data/processed/startestmeta.csv") %>%
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
#============================================================#
# read data
#============================================================#
at <- read.csv("../data/raw/2024star_viame_339506/autotest2024_viame_cascade.csv") %>%
  mutate(truedetect = if_else(truedetect == "True", TRUE, FALSE)) %>%
  rename(spname = Spname) %>%
  rename(image_name = Imagename)
mt <- read.csv("../data/raw/2024star_viame_339506/mantest2024_viame_cascade.csv") %>%
  rename(spname = Spname) %>%
  rename(image_name = Imagename)

out <- build_detection_tables(mt, at, metastar)

out$calib_df$spname = "star"

model_name <- "Casv2star"

res <- structure_by_region(out, model_name, 
                           region_lat_cutoff = 40, 
                           save_rdata = TRUE)

names(res)


# same val set
# 2022scallop_viame_257794
# 2022scallop_yolo11n_257793
# full dataset 22-24
# 2224scallop_yolo11n_260064