library(dplyr)
library(purrr)
library(tibble)
source("./datfunc.R")
# load("../data/raw/allimg_2024.RData")  # loads `allimg for 2024`, complete metadata
# metagb = read.csv("../data/raw/metagb22.csv")
# metamab = read.csv("../data/raw/metama22.csv")
# meta22 = rbind(metagb, metamab)
#============================================================#
# read data
#============================================================#
at <- read.csv("../data/raw/autotestyolo22.csv")
mt <- read.csv("../data/raw/mantestyolo22.csv")

out <- build_detection_tables(mt, at, meta22)

model_name <- "YOLO2022"

res <- structure_by_region(out, model_name, 
                           region_lat_cutoff = 40, 
                           save_rdata = TRUE)

names(res)
