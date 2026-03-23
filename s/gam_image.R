library(tidyverse)
source("./gamfunc.R")
source("./fitfunc.R")

load("../data/processed/Cas2024v2.RData")
# load("../data/processed/YOLO2022.RData")

model = Cas2022v2
model_name = deparse(substitute(Cas2022v2))

gams_cal <- fit_calibration_gams(
  model = model,
  model_name = model_name
)

########################################################################
# P(FN) ~ SUM(P(detection))
########################################################################

res_gb <- compute_image_level_counts(
  calib_df  = get_data(model, model_name, "GB"),
  gam       = gams_cal$GB,
  region    = "GB",
  model_name = model_name
)

res_mab <- compute_image_level_counts(
  calib_df  = get_data(model, model_name, "MAB"),
  gam       = gams_cal$MAB,
  region    = "MAB",
  model_name = model_name
)

str(res_gb)

img_all <- bind_rows(
  res_gb$img,
  res_mab$img
)

gams = fitfn_calibration_mod(model = model,
                              model_name = model_name,
                              p_detect = img_all)
summary(gams$comb)

m <- gams$comb
dat = gams$data
########################################################################
# 
########################################################################
ggplot(dat, aes(x = predicted_number, y = fn_pres)) +
  geom_jitter(height = 0.01, alpha = 0.2) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE) +
  labs(
    y = "Observed FN Presence",
    x = "Σ P(detection)",
    title = "Logistic Fit: FN Presence vs Predicted Density"
  ) +
  theme_minimal()

library(pROC)

dat$fn_prob <- predict(gams$comb, type = "response")

roc_obj <- roc(dat$fn_pres, dat$fn_prob)

plot(roc_obj, col = "blue", main = "ROC Curve: FN Model")
auc(roc_obj)

library(dplyr)

calib_df <- dat %>%
  mutate(bin = cut(fn_prob, breaks = seq(0, 1, by = 0.1))) %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(fn_prob),
    obs_rate  = mean(fn_pres),
    n = n()
  )

ggplot(calib_df, aes(mean_pred, obs_rate)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Mean Predicted Probability",
    y = "Observed FN Rate",
    title = "Calibration Curve"
  ) +
  theme_minimal()

# Brier Score (overall probabilistic accuracy)
# 0 = perfect
mean((dat$fn_prob - dat$fn_pres)^2)

ggplot(dat, aes(fn_prob, fn_pres)) +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE) +
  labs(
    x = "Predicted FN Probability",
    y = "Observed FN Frequency",
    title = "Empirical FN Rate vs Predicted Probability"
  ) +
  theme_minimal()

save_cal()