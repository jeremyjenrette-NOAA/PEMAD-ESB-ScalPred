library(dplyr)
library(purrr)
library(tibble)
source("./datfunc.R")

at_cas <- read.csv("../data/raw/2224scallop_viame_370248/autotest2224_viame_cascade.csv") %>%
  mutate(truedetect = if_else(truedetect == "True", TRUE, FALSE)) %>%
  rename(spname = Spname) %>%
  rename(imagename = Imagename)

# mt_cas <- read.csv("../data/raw/2224scallop_viame_370248/mantest2224_viame_cascade.csv") %>%
#   rename(spname = Spname) %>%
#   rename(imagename = Imagename)

#============================================================#
# read data
#============================================================#
at <- read.csv("../data/raw/2224scallop_yolo12n_370283/eval/autotest2224_yolo12n.csv") %>%
  mutate(truedetect = if_else(truedetect == "True", TRUE, FALSE))
mt <- read.csv("../data/raw/2224scallop_yolo12n_370283/eval/mantest2224_yolo12n.csv")
# meta = read.csv("../data/processed/metadata2224.csv")

all_test_imgs = read.csv("../data/raw/2224scallop_yolo12n_370283/eval/val_images2224_yolo12n.csv") 

out <- build_detection_tables(mt = mt_cas, at = at_cas, meta = meta_std, all_test_imgs = all_test_imgs)

model_name <- "Casv2strat2224"

res <- structure_by_region(out, model_name, 
                           region_lat_cutoff = 40, 
                           save_rdata = TRUE)

names(res)


####
dat_split = read.csv("../data/raw/dataset_split_2224.csv")

meta_std <- dat_split %>%
  rename(
    # Imagename              = imagename,
    altitude               = ALTIMETER_ALTITUDE_METER,
    altitude2              = ALTIMETER_ALTITUDE2_METER,
    stereo_altitude        = STEREO_ALTITUDE_METER,
    backscatter            = FLUOROMETER_BACKSCATTER_NTU,
    bottom_depth           = FATHOMETER_OCEAN_DEPTH_METER,
    cdom                   = FLUOROMETER_CDOM,
    chlorophyll            = FLUOROMETER_CHLOROPHYLL,
    # field_of_view_sq_meter = FIELD_OF_VIEW_SQ_METER,
    field_of_view_source   = FIELD_OF_VIEW_SOURCE,
    heading                = VEHICLE_MAGNETIC_HEADING,
    # latitude               = SHIP_LATITUDE,
    # longitude              = SHIP_LONGITUDE,
    # millimeter_per_pixel   = MILLIMETER_PER_PIXEL,
    o2                     = CTD2_DISSOLVED_OXYGEN,
    pitch                  = VEHICLE_PITCH_ANGLE,
    roll                   = VEHICLE_ROLL_ANGLE,
    s                      = CTD_SALINITY,
    t                      = CTD_TEMPERATURE_CELSIUS,
    fluorometer_signal     = FLUOROMETER_SIGNAL,
    datetime               = image_timestamp,
    v_depth                = CTD_VEHICLE_DEPTH_METER,
    gear                   = GEAR,
    cruise_id              = CRUISE_ID,
    habcam_pk              = HABCAM_PK
  ) %>%
  mutate(
    bottom_depth = altitude + v_depth
  ) %>%
  # mutate(
  #   Imagename = basename(Imagename)
  # ) %>%
  distinct(imagename, .keep_all = TRUE)



