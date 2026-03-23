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
#============================================================#
# read data
#============================================================#
at <- read.csv("../data/raw/2224scallop_viame_260065/autotest2224_viame_cascade.csv") %>%
  mutate(truedetect = if_else(truedetect == "True", TRUE, FALSE)) %>%
  rename(spname = Spname) %>%
  rename(image_name = Imagename)
mt <- read.csv("../data/raw/2224scallop_viame_260065/mantest2224_viame_cascade.csv") %>%
  rename(spname = Spname) %>%
  rename(image_name = Imagename)

out <- build_detection_tables(mt, at, meta2224)

model_name <- "Cas2224"

res <- structure_by_region(out, model_name, 
                           region_lat_cutoff = 40, 
                           save_rdata = TRUE)

names(res)


# same val set
# 2022scallop_viame_257794
# 2022scallop_yolo11n_257793
# full dataset 22-24
# 2224scallop_yolo11n_260064