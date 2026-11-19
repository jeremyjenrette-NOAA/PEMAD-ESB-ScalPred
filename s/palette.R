# ======================================================================
# palette.R - one set of workflow and detector colours for every figure.
# source("./palette.R") from s/ (or "../s/palette.R" from length/). palette.py mirrors it.
#   Two-stage YOLOv12 green · YOLOv12 F1 red · YOLOv12 stage 1 / sum of p orange
#   Cascade R-CNN blue (F1 dark, stage 1 light) · RT-DETR-L purple
# Green and red are also told apart by point shape (circle vs triangle) where both appear.
# ======================================================================
PAL <- list(
  two_stage      = "#1B9E77",   # two-stage YOLOv12
  two_stage_yc   = "#0B5D46",   # two-stage YOLOv12 + Cascade R-CNN
  f1_yolo        = "#D7301F",   # YOLOv12 F1 threshold
  stage1_yolo    = "#EB6834",   # YOLOv12 stage 1 only / calibrated sum of p
  f1_cascade     = "#08519C",   # Cascade R-CNN F1 threshold
  stage1_cascade = "#6BAED6",   # Cascade R-CNN stage 1 only
  yolo           = "#1B9E77",   # detector colours (Fig. 3)
  cascade        = "#2171B5",
  rtdetr         = "#FC2808",
  simplest       = "#8C8C8C",   # simplest-form variants
  ink            = "#1A1A1A")
METHOD_SHAPES <- c("Two-stage YOLOv12" = 16, "Two-stage YOLOv12 + Cascade" = 18, "YOLOv12 F1" = 17,
                   "Cascade F1" = 15, "YOLOv12 stage 1" = 1, "Cascade stage 1" = 0)
