library(tidyverse)
source("./procfunc.R")
#============================================================#
load("../data/processed/YOLOv102022.RData")
load("../data/processed/YOLOv112022.RData")
load("../data/processed/YOLOv262022.RData")

load("../data/processed/Cas2022v2.RData")
load("../data/processed/YOLOv72022.RData")
#============================================================#
# Analysis: Precision-Recall
#============================================================#
# out_pr <- evaluate_pr_models(list(YOLOv10_2022, YOLOv11_2022, YOLOv26_2022))
out_pr <- evaluate_pr_models(list(YOLOv102022, YOLOv112022, YOLOv262022), 
                             stratify_region = FALSE)


pr_all <- out_pr$pr_all
p_pr   <- out_pr$p_pr
p_pr
#============================================================#
# Analysis: F1 score
#============================================================#
best_pts <- pr_all %>%
  group_by(model) %>%
  filter(f1 == max(f1, na.rm = TRUE)) %>%
  slice_max(conf, n = 1) %>%   # break ties by confidence
  ungroup()


p_f1 <- ggplot(pr_all, aes(x = conf, y = f1, color = model)) +
  geom_line(linewidth = 1.2) +
  
  # X marks the spot
  geom_point(
    data = best_pts,
    aes(x = conf, y = f1),
    shape = 4,        # X
    size = 4,
    stroke = 1.2
  ) +
  
  # Text annotation
  # geom_text(
  #   data = best_pts,
  #   aes(
  #     label = sprintf("F1 = %.3f, conf = %.2f", f1, conf)
  #   ),
  #   hjust = 1,
  #   vjust = -1,
  #   size = 3.5,
  #   show.legend = FALSE
  # ) +
  
  theme_minimal() +
  labs(
    title = "F1 Score vs Confidence Threshold",
    x = "Confidence threshold",
    y = "F1 score",
    color = "Model"
  )

p_f1

save_pr_f1(p_pr, p_f1, id = "2022yolo_10v11", outdir = "../figures/")
