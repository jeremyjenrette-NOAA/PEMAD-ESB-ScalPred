library(dplyr)
library(mgcv)

# 1. Source the function repository
source("inference_functions.R")

# 2. Load the inference dataset and pre-trained GAMs
models <- readRDS("../data/processed/eval_inference_2226.rds")
yolo_gams <- readRDS("../data/processed/YOLOv12_detgams_2226.rds")
cas_gams  <- readRDS("../data/processed/CascadeR-CNN_detgams_2226.rds")

# ======================================================================
# 3. Process YOLOv12
# ======================================================================
cat("Processing YOLOv12 Inference...\n")
yolo_det_calib <- predict_calibration_gams(models$YOLOv12$det, yolo_gams)
yolo_img_calib <- aggregate_image_predictions(yolo_det_calib, models$YOLOv12$img)

# Update the master list
models$YOLOv12$det <- yolo_det_calib
models$YOLOv12$img <- yolo_img_calib

# ======================================================================
# 4. Process Cascade R-CNN
# ======================================================================
cat("Processing Cascade R-CNN Inference...\n")
cas_det_calib <- predict_calibration_gams(models$`Cascade R-CNN`$det, cas_gams)
cas_img_calib <- aggregate_image_predictions(cas_det_calib, models$`Cascade R-CNN`$img)

# Update the master list
models$`Cascade R-CNN`$det <- cas_det_calib
models$`Cascade R-CNN`$img <- cas_img_calib

# ======================================================================
# 5. Export Calibrated Dataset
# ======================================================================
saveRDS(models, file = "../data/processed/eval_inference_calibrated_2226.rds")
cat("Inference calibration complete and saved!\n")