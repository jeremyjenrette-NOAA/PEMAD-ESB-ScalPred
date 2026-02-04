mt_v1 <- read.csv("../data/raw/mantest24.csv")  |> clean_names()
at_v1 <- read.csv("../data/raw/autotest24.csv") |> clean_names()
mt_v2 <- read.csv("../data/raw/mantest24v2.csv")  |> clean_names()
at_v2 <- read.csv("../data/raw/autotest24v2.csv") |> clean_names()
load("../data/raw/allimg_2024.RData")  # loads `allimg`, complete metadata
allimg <- allimg |> clean_names()

library(dplyr)
library(purrr)
library(tibble)

source("./datfunc.R")

out_v1 <- build_detection_tables(mt_v1, at_v1, allimg)
out_v2 <- build_detection_tables(mt_v2, at_v2, allimg)
#============================================================#
# Precision–Recall + F1 evaluation for one model
#============================================================#
evaluate_pr_curve <- function(
    calib_df,
    img_df,
    conf_grid = seq(0, 1, by = 0.0005),
    model_id  = "model"
) {
  
  map_dfr(conf_grid, function(th) {
    
    # detections above threshold
    df_th <- calib_df |> filter(conf >= th)
    
    # detection-level counts
    TP <- sum(df_th$truedetect)
    FP <- sum(!df_th$truedetect)
    
    # how many manual annotations were successfully 
    # detected by the model at or above confidence th?
    # False positives are explicitly excluded
    # Images with zero detections above threshold are not present,
    # until next step
    auto_true_img <- df_th |>
      filter(truedetect) |>
      count(image_id, name = "n_auto_true")
    
    # excluding FP, then calculate FN per image
    # even from images where TP = 0
    # Total number of scallops that exist but were not detected at confidence ≥ th
    FN <- img_df |>
      left_join(auto_true_img, by = "image_id") |>
      mutate(
        n_auto_true = replace_na(n_auto_true, 0L),
        FN = pmax(n_manual - n_auto_true, 0L)
      ) |>
      summarise(FN = sum(FN)) |>
      pull(FN)
    
    precision <- TP / (TP + FP)
    recall    <- TP / (TP + FN)
    f1        <- 2 * precision * recall / (precision + recall)
    
    tibble(
      model     = model_id,
      conf      = th,
      TP        = TP,
      FP        = FP,
      FN        = FN,
      precision = precision,
      recall    = recall,
      f1        = f1
    )
  })
}

pr_v1 <- evaluate_pr_curve(
  calib_df = out_v1$calib_df,
  img_df   = out_v1$img_df,
  model_id = "autotest v1"
)

pr_v2 <- evaluate_pr_curve(
  calib_df = out_v2$calib_df,
  img_df   = out_v2$img_df,
  model_id = "autotest v2"
)

pr_all <- bind_rows(pr_v1, pr_v2)
glimpse(pr_all)
#============================================================#
# Joint plotting for model comparison
#============================================================#
p_pr <- ggplot(pr_all, aes(x = recall, y = precision, color = model)) +
  geom_path(linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "Precision–Recall Curves", # by detection
    x = "Recall",
    y = "Precision",
    color = "Model"
  )

p_pr

p_f1 <- ggplot(pr_all, aes(x = conf, y = f1, color = model)) +
  geom_line(linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "F1 Score vs Confidence Threshold",
    x = "Confidence threshold",
    y = "F1 score",
    color = "Model"
  )

p_f1

save_pr_f1 <- function(p_pr, p_f1, id, outdir = "~/Downloads",
                       width = 6, height = 6) {
  
  ggsave(
    plot = p_pr,
    filename = file.path(outdir, paste0("pr_", id, ".jpg")),
    width = width,
    height = height
  )
  
  ggsave(
    plot = p_f1,
    filename = file.path(outdir, paste0("f1_", id, ".jpg")),
    width = width,
    height = height
  )
}

save_pr_f1(p_pr, p_f1, id = "total")
#============================================================#
# Regional comparison
#============================================================#
pr_v1 <- evaluate_pr_curve(
  calib_df = out_v1$calib_df %>% filter(latitude < 40),
  img_df   = out_v1$img_df %>% filter(latitude < 40),
  model_id = "autotest v1"
)

pr_v2 <- evaluate_pr_curve(
  calib_df = out_v2$calib_df %>% filter(latitude < 40),
  img_df   = out_v2$img_df %>% filter(latitude < 40),
  model_id = "autotest v2"
)

pr_all <- bind_rows(pr_v1, pr_v2)
glimpse(pr_all)
#============================================================#
# Joint plotting for model comparison
#============================================================#
p_pr <- ggplot(pr_all, aes(x = recall, y = precision, color = model)) +
  geom_path(linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "Precision–Recall Curves", # by detection
    x = "Recall",
    y = "Precision",
    color = "Model"
  )

p_pr

p_f1 <- ggplot(pr_all, aes(x = conf, y = f1, color = model)) +
  geom_line(linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "F1 Score vs Confidence Threshold",
    x = "Confidence threshold",
    y = "F1 score",
    color = "Model"
  )

p_f1


save_pr_f1(p_pr, p_f1, id = "MAB")