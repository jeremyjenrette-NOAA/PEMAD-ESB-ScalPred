library(mgcv)
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
source("./datfunc.R")
source("./fitfunc.R")

# load("../data/processed/Casfixed2022.RData")
# load("../data/processed/YOLOv11fixed2022.RData")
# meta2224 <- read.csv("../data/processed/metadata2224.csv")
dat_split

model1 = YOLOv12strat2224
model_name1 = deparse(substitute(YOLOv12strat2224))

gams_yolo <- fit_calibration_gams(
  model = model1,
  model_name = model_name1,
  dat_split = dat_split,
  use_strat = TRUE
)

model2 = Casv2strat2224
model_name2 = deparse(substitute(Casv2strat2224))

gams_cas <- fit_calibration_gams(
  model = model2,
  model_name = model_name2,
  dat_split = dat_split,
  use_strat = TRUE
)

img_yolo  <- run_model_pipeline(model1, model_name1, gams = gams_yolo, test = FALSE, 
                                dat_split = dat_split, use_strat = TRUE)

img_cas  <- run_model_pipeline(model2, model_name2, gams = gams_cas, test = FALSE,
                               dat_split = dat_split, use_strat = TRUE)

img_yolo_test  <- run_model_pipeline(model1, model_name1, gams = gams_yolo, test = TRUE, 
                                dat_split = dat_split, use_strat = TRUE)

img_cas_test  <- run_model_pipeline(model2, model_name2, gams = gams_cas, test = TRUE,
                               dat_split = dat_split, use_strat = TRUE)

# for some reason, there is a very small number of mismatch images
# this makes sure the gam receives only the matches
common_ids <- intersect(img_yolo$img$imagename, img_cas$img$imagename)

common_ids_test <- intersect(img_yolo_test$img$imagename, img_cas_test$img$imagename)

img_combined <- img_yolo$img %>%
  filter(imagename %in% common_ids) %>%
  select(imagename, region,
         pred_yolo = predicted_number,
         true_number) %>%
  left_join(
    img_cas$img %>%
      filter(imagename %in% common_ids) %>%
      select(imagename,
             pred_cascade = predicted_number),
    by = "imagename"
  ) %>%
  mutate(
    pred_mean = (pred_yolo + pred_cascade)/2,
    pred_diff = pred_cascade - pred_yolo,
    pred_sum  = pred_yolo + pred_cascade,
    log_yolo_pred = log1p(pred_yolo),
    log_cascade_pred = log1p(pred_cascade),
    log_pred_diff = log1p(abs(pred_diff)),
    log_true_number = log1p(true_number + 1e-6)
  )

img_combined_test <- img_yolo_test$img %>%
  filter(imagename %in% common_ids_test) %>%
  select(imagename, region,
         pred_yolo = predicted_number,
         true_number) %>%
  left_join(
    img_cas_test$img %>%
      filter(imagename %in% common_ids_test) %>%
      select(imagename,
             pred_cascade = predicted_number),
    by = "imagename"
  ) %>%
  mutate(
    pred_mean = (pred_yolo + pred_cascade)/2,
    pred_diff = pred_cascade - pred_yolo,
    pred_sum  = pred_yolo + pred_cascade,
    log_yolo_pred = log1p(pred_yolo),
    log_cascade_pred = log1p(pred_cascade),
    log_pred_diff = log1p(abs(pred_diff)),
    log_true_number = log1p(true_number + 1e-6)
  )

img_combined <- attach_metadata_to_images(
  img_df  = img_combined,
  meta_df = dat_split
)

img_combined_test <- attach_metadata_to_images(
  img_df  = img_combined_test,
  meta_df = dat_split
)

models <- list(
  
  m1 = gam(true_number ~ s(pred_yolo, pred_diff), 
           family = nb(), data = img_combined),
  
  m2 = gam(true_number ~ s(pred_cascade, pred_diff), 
           family = nb(), data = img_combined),
  
  m3 = gam(true_number ~ pred_diff, 
           family = nb(), data = img_combined),
  
  m4 = gam(true_number ~ ti(pred_yolo, pred_cascade) + pred_diff, 
           family = nb(), data = img_combined),
  
  # m5 = gam(true_number ~ ti(pred_yolo, pred_cascade) +
  #   s(pred_yolo - pred_cascade) + ti(pred_yolo, pred_cascade, by = pred_diff),
  #   family = nb(), data = img_combined),
  
  m6 = gam(true_number ~ s(pred_yolo, pred_diff, bs = "tp"),
           family = nb(), data = img_combined),
  
  m7 = gam(true_number ~ s(log_yolo_pred, pred_diff),
           family = nb(), data = img_combined),
  
  m75 = gam(true_number ~ region + s(log_yolo_pred, pred_diff),
           family = nb(), data = img_combined),
  
  m8 = gam(true_number ~ s(log_cascade_pred, pred_diff),
           family = nb(), data = img_combined),
  
  m9 = gam(true_number ~ s(log_yolo_pred, pred_diff, bs = "tp"),
            family = nb(), data = img_combined),
  
  m10 = gam(log_true_number ~ s(log_yolo_pred, log_pred_diff),
           family = nb(), data = img_combined),
  
  m11 = gam(log_true_number ~ s(log_yolo_pred, pred_diff),
           family = nb(), data = img_combined)
  )


evaluate_model <- function(model, data) {
  pred <- predict(model, type = "response")
  true <- data$true_number
  
  data.frame(
    Model = deparse(formula(model)),  # <- this is the key line
    AIC = AIC(model),
    RMSE = sqrt(mean((pred - true)^2)),
    MAE = mean(abs(pred - true)),
    Deviance = deviance(model)
  )
}

results <- lapply(models, evaluate_model, data = img_combined)
results_df <- do.call(rbind, results)

rownames(results_df) <- names(models)  # keep m1, m2, etc.
results_df

df_eval_train <- data.frame(
  pred = predict(models$m1, type = "response"),
  true = img_combined$true_number
)

df_eval <- data.frame(
  pred = predict(models$m8, newdata = img_combined_test, type = "response"),
  true = img_combined_test$true_number
)

pred <- df_eval$pred
true <- df_eval$true

rmse <- sqrt(mean((pred - true)^2))
mae  <- mean(abs(pred - true))

# Correlation-based R² (good for visualization)
r2 <- cor(pred, true)^2

# Deviance (from model)
dev <- deviance(models$m8)

# AIC
aic = AIC(models$m8)

plot(residuals(models$m75) ~ img_combined$pred_yolo, main = "YOLO")
abline(h=0, col = "red", lty = 2)
plot(residuals(models$m75) ~ img_combined$pred_cascade, main = "Cascade R-CNN")
abline(h=0, col = "red", lty = 2)
plot(residuals(models$m75) ~ img_combined$pred_diff, main = "Cascade - YOLO")
abline(h=0, col = "red", lty = 2)
########################################################################
## confirm two similar model eval
# library(rsample)
# library(purrr)
# 
# folds <- rsample::vfold_cv(img_combined, v = 5)
# 
# cv_results <- folds %>%
#   mutate(
#     model = map(splits, ~ gam(
#       true_number ~ te(pred_yolo, pred_cascade) + abs(pred_diff),
#       family = nb(),
#       data = rsample::analysis(.x)
#     )),
#     pred = map2(model, splits, ~ predict(.x, newdata = rsample::assessment(.y), type = "response")),
#     truth = map(splits, ~ rsample::assessment(.x)$true_number)
#   )
# cv_results <- cv_results %>%
#   mutate(
#     rmse = map2_dbl(pred, truth, ~ sqrt(mean((.x - .y)^2))),
#     mae  = map2_dbl(pred, truth, ~ mean(abs(.x - .y)))
#   )
# 
# mean(cv_results$rmse)
# mean(cv_results$mae)
########################################################################
library(ggplot2)

disagreement = ggplot(img_combined, aes(x = pred_diff, y = true_number)) +
  geom_point(alpha = 0.25, size = 2, color = "#2C7FB8") +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 7),
    color = "#D95F02",
    linewidth = 1.2,
    se = TRUE
  ) +
  labs(
    title = "2022-2024: Model Disagreement on True Scallop Count",
    # subtitle = "Difference = Cascade prediction - YOLO prediction",
    x = "Prediction Difference (Cascade - YOLO)",
    y = "True Scallop Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40")
  )
disagreement
########################################################################
label_text <- paste0(
  "R² = ", round(r2, 3), "\n",
  "RMSE = ", round(rmse, 3), "\n",
  "MAE = ", round(mae, 3), "\n",
  "AIC = ", round(aic, 1)
)

gammod = ggplot(df_eval, aes(x = pred, y = true)) +
  geom_point(alpha = 0.25, size = 2, color = "black") +
  
  # Perfect prediction line
  # geom_abline(
  #   slope = 1, intercept = 0,
  #   linetype = "dashed",
  #   color = "gray50",
  #   linewidth = 1
  # ) +
  
  # Fitted relationship
  geom_smooth(
    method = "lm",
    color = "black",
    linewidth = 1.2,
    se = TRUE
  ) +
  
  # Metrics annotation
  annotate(
    "label",
    x = Inf, y = -Inf,
    label = label_text,
    hjust = 1.1, vjust = -0.1,
    size = 4.5,
    fill = "white",
    color = "black",
    label.size = 0.3
  ) +
  
  labs(
    title = "2022 - 2024: Predicted vs True Scallop Counts",
    subtitle = paste0("Model: ", "true_number ~ region + s(log_cas_pred, pred_diff)"),
    x = "Predicted Count",
    y = "True Count"
  ) +
  
  coord_equal() +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40")
  )
gammod

true_total <- sum(df_eval$true, na.rm = TRUE)

summary_df <- data.frame(
  metric = c("True", "Cascade R-CNN", "YOLOv12", "s(YOLOv12 - Cascade)", "YOLOv12 F1"),
  value  = c(
    true_total,
    sum(img_cas_test$calib_df$pred_p, na.rm = TRUE),
    sum(img_yolo_test$calib_df$pred_p, na.rm = TRUE),
    sum(df_eval$pred, na.rm = TRUE),
    sum(img_yolo_test$img$predicted_f1_number, na.rm = TRUE)
  )
) %>%
  mutate(
    error = value - true_total,
    pct_error = 100 * error / true_total,
    abs_pct_error = abs(pct_error),
    label = paste0(
      comma(round(value, 0)),
      "\n",
      ifelse(metric == "True", "Groundtruth",
             paste0(ifelse(error > 0, "+", ""), round(pct_error, 1), "%"))
    ),
    metric = fct_reorder(metric, abs_pct_error, .desc = TRUE)
  )

summary_error_df <- summary_df %>%
  filter(metric != "True") %>%
  mutate(
    metric = fct_reorder(metric, pct_error),
    error_label = paste0(
      ifelse(error > 0, "+", ""),
      comma(round(error, 0)),
      " scallops\n",
      ifelse(pct_error > 0, "+", ""),
      round(pct_error, 1),
      "%"
    )
  )

p_error_summary <- ggplot(summary_error_df, aes(x = metric, y = pct_error)) +
  
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_text(
    aes(
      label = error_label,
      hjust = ifelse(pct_error >= 0, -0.05, 1.05)
    ),
    size = 4
  ) +
  
  coord_flip(clip = "off") +
  
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  
  labs(
    title = "Abundance Estimation Error by Method",
    subtitle = paste0("True Abundance = ", comma(round(true_total, 0)), " scallops"),
    x = NULL,
    y = "Percent error"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = margin(10, 45, 10, 10)
  )

p_error_summary

p_resid = df_eval %>%
  mutate(residual = pred - true) %>%
  ggplot(aes(true, residual)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs True Count",
       subtitle = "s(YOLOv12 - Cascade R-CNN)")
p_resid

save_cal(p_sum, id = paste("2224_",model_name, "_sum_strat_all",sep=""), 
         outdir = "../figures/diag2/", format = "png", width = 5, height = 4.5)
save_cal(gammod, id = paste("2224yolov12cas_gammod",sep=""),
         outdir = "../figures/diag2/", width = 6, height = 7, format = "png")
save_cal(p_error_summary, id = paste("2224total_error",sep=""),
         outdir = "../figures/diag2/", width = 9, height = 8, format = "png")
save_cal(p_resid, id = paste("2224total_residual",sep=""),
         outdir = "../figures/diag2/", width = 5, height = 4.5, format = "png")
