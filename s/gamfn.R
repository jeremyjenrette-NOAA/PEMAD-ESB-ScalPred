library(tidyverse)
source("./gamfunc.R")

load("../data/processed/Cas2024v2.RData")
# load("../data/processed/YOLO2022.RData")

model = Cas2024v2
model_name = deparse(substitute(Cas2024v2))

gams = fitfn_calibration_gams(model = model,
                       model_name = model_name)
summary(gams$comb)

m <- gams$comb

# Extract ranges
d_range  <- range(exp(m$model$auto_density_log) - 1e-6)
fov_range <- range(m$model$field_of_view_sq_meter)

# Build grid
newdat <- expand.grid(
  auto_density = seq(d_range[1], d_range[2], length.out = 100),
  field_of_view_sq_meter = seq(fov_range[1], fov_range[2], length.out = 100)
)

# Convert back to log scale for model
newdat$auto_density_log <- log(newdat$auto_density + 1e-6)

# Predict
newdat$pred <- predict(m, newdat, type = "response")

pred_link <- predict(m, newdat, type = "link", se.fit = TRUE)

newdat$fit  <- pred_link$fit
newdat$se   <- pred_link$se.fit

newdat$upper <- plogis(newdat$fit + 1.96 * newdat$se)
newdat$lower <- plogis(newdat$fit - 1.96 * newdat$se)
newdat$pred  <- plogis(newdat$fit)

library(pROC)

pred_obs <- predict(m, type = "response")
roc_obj <- roc(m$model$fn_pres, pred_obs)

plot(roc_obj)
auc(roc_obj)


