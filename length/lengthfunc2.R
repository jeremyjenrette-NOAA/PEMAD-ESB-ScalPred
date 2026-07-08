# ==============================================================================
# Script: lengthfunc2.R
# Purpose: Geometric utility and dynamic taxonomic routing functions
# ==============================================================================

library(stringr)
library(dplyr)

#' Route Raw Class Strings into 4 Operational Subclasses
#' @param class_name Character vector of original annotations
#' @return Character vector mapping to strict study sub-types
assign_scallop_subclass <- function(class_name) {
  class_clean <- tolower(trimws(class_name))
  
  case_when(
    str_detect(class_clean, "swimming") ~ "Swimming Sea Scallop",
    str_detect(class_clean, "inexact")  ~ "Inexact Length", 
    str_detect(class_clean, "width")    ~ "Live/Probable Width",
    TRUE                                ~ "Live/Probable Length"
  )
}

#' Extract True Line Length from Geometry Text
parse_gt_length <- function(geometry_text) {
  coords <- lapply(str_extract_all(geometry_text, "[-+]?[0-9]*\\.?[0-9]+"), as.numeric)
  lengths <- sapply(coords, function(pts) {
    if (length(pts) >= 4) {
      sqrt((pts[3] - pts[1])^2 + (pts[4] - pts[2])^2)
    } else {
      NA
    }
  })
  return(lengths)
}

#' Back-calculate Line Length from Bounding Box Coordinates
back_calc_length <- function(tlx, tly, brx, bry, method = "average") {
  w <- abs(brx - tlx)
  h <- abs(bry - tly)
  
  case_when(
    method == "average"  ~ (w + h) / 2,
    method == "max"      ~ pmax(w, h),
    method == "diagonal" ~ sqrt(w^2 + h^2) / sqrt(2),
    TRUE ~ (w + h) / 2
  )
}