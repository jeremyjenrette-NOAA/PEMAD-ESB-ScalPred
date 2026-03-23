library(mgcv)
source("./datfunc.R")
source("./fitfunc.R")

load("../data/processed/Casfixed2022.RData")
load("../data/processed/YOLOv11fixed2022.RData")
meta2224 <- read.csv("../data/processed/metadata2224.csv")

model1 = YOLOv122224
model_name1 = deparse(substitute(YOLOv122224))

model2 = Cas2224
model_name2 = deparse(substitute(Cas2224))

img_yolo  <- run_model_pipeline(model1, model_name1)
img_cascade  <- run_model_pipeline(model2, model_name2)

# for some reason, there is a very small number of mismatch images
# this makes sure the gam receives only the matches
common_ids <- intersect(img_yolo$image_id, img_cascade$image_id)

img_combined <- img_yolo %>%
  filter(image_id %in% common_ids) %>%
  select(image_id, region,
         pred_yolo = predicted_number,
         true_number) %>%
  left_join(
    img_cascade %>%
      filter(image_id %in% common_ids) %>%
      select(image_id,
             pred_cascade = predicted_number),
    by = "image_id"
  ) %>%
  mutate(
    pred_mean = (pred_yolo + pred_cascade)/2,
    pred_diff = pred_cascade - pred_yolo,
    pred_sum  = pred_yolo + pred_cascade
  )

img_combined <- attach_metadata_to_images(
  img_df  = img_combined,
  meta_df = meta2224
)

models <- list(
  
  m1 = gam(true_number ~ s(pred_yolo, pred_diff), 
           family = nb(), data = img_combined),
  
  m2 = gam(true_number ~ s(pred_cascade, pred_diff), 
           family = nb(), data = img_combined),
  
  m3 = gam(true_number ~ pred_diff, 
           family = nb(), data = img_combined),
  
  m4 = gam(true_number ~ ti(pred_yolo, pred_cascade) + pred_diff, 
           family = nb(), data = img_combined)
  
  # m5 = gam(true_number ~ ti(pred_yolo, pred_cascade) +
  #   s(pred_yolo - pred_cascade) + ti(pred_yolo, pred_cascade, by = pred_diff),
  #   family = nb(), data = img_combined)
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

df_eval <- data.frame(
  pred = predict(models$m1, type = "response"),
  true = img_combined$true_number
)

pred <- df_eval$pred
true <- df_eval$true

rmse <- sqrt(mean((pred - true)^2))
mae  <- mean(abs(pred - true))

# Correlation-based R² (good for visualization)
r2 <- cor(pred, true)^2

# Deviance (from model)
dev <- deviance(models$m1)

# AIC
aic = AIC(models$m1)

plot(residuals(models$m1) ~ img_combined$pred_yolo, main = "YOLO")
abline(h=0, col = "red", lty = 2)
plot(residuals(models$m1) ~ img_combined$pred_cascade, main = "Cascade R-CNN")
abline(h=0, col = "red", lty = 2)
########################################################################
## confirm two similar model eval
library(rsample)
library(purrr)

folds <- rsample::vfold_cv(img_combined, v = 5)

cv_results <- folds %>%
  mutate(
    model = map(splits, ~ gam(
      true_number ~ te(pred_yolo, pred_cascade) + abs(pred_diff),
      family = nb(),
      data = rsample::analysis(.x)
    )),
    pred = map2(model, splits, ~ predict(.x, newdata = rsample::assessment(.y), type = "response")),
    truth = map(splits, ~ rsample::assessment(.x)$true_number)
  )
cv_results <- cv_results %>%
  mutate(
    rmse = map2_dbl(pred, truth, ~ sqrt(mean((.x - .y)^2))),
    mae  = map2_dbl(pred, truth, ~ mean(abs(.x - .y)))
  )

mean(cv_results$rmse)
mean(cv_results$mae)
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
  geom_point(alpha = 0.25, size = 2, color = "#1B9E77") +
  
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
    color = "#D95F02",
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
    subtitle = "Model: true_count ~ s(yolo, cascade - yolo)",
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

save_cal(disagreement, id = paste("disagreement_2224yolov12cas",sep=""), 
         outdir = "../figures/", width = 8.5, height = 5, format = "png")

save_cal(gammod, id = paste("gammod_2224yolov12cas",sep=""), 
         outdir = "../figures/", width = 6, height = 7, format = "png")
