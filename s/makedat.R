library(dplyr)
library(purrr)
library(tibble)
source("./datfunc.R")
# load("../data/raw/allimg_2024.RData")  # loads `allimg for 2024`, complete metadata
metagb = read.csv("../data/raw/metagb22.csv")
metamab = read.csv("../data/raw/metama22.csv")
missing = read.csv("../data/raw/missingmeta.csv")
meta22 = rbind(metagb, metamab, missing)
meta22 <- meta22[!duplicated(meta22$IMAGENAME), ]
#============================================================#
# read data
#============================================================#
at <- read.csv("../data/raw/2022scallop_yolov10n_242223/eval/autotest2022_yolov10n.csv") %>%
  mutate(truedetect = if_else(truedetect == "True", TRUE, FALSE)) %>%
  rename(spname = Spname) %>%
  rename(image_name = Imagename)
mt <- read.csv("../data/raw/2022scallop_yolov10n_242223/eval/mantest2022_yolov10n.csv") %>%
  rename(spname = Spname) %>%
  rename(image_name = Imagename)

out <- build_detection_tables(mt, at, meta22)

model_name <- "YOLOv102022"

res <- structure_by_region(out, model_name, 
                           region_lat_cutoff = 40, 
                           save_rdata = TRUE)

names(res)
