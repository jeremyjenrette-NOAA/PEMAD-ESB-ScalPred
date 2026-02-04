library(mgcv)
set.seed(44)

source("./datfunc.R")
out_v1 <- build_detection_tables(mt_v1, at_v1, allimg)
out_v2 <- build_detection_tables(mt_v2, at_v2, allimg)

m_calib_0 <- gam(
  y ~ s(conf, k = 5) + # y is binary 0,1 of truedetect 
    s(bottom_depth, k = 5) +
    ti(conf, bottom_depth, 
       bs = c("tp", "tp"),
       k = c(5, 5)),
  family = binomial(link = "logit"),
  data = out_v2$calib_df %>% filter(latitude < 40),
  method = "REML"
)
summary(m_calib_0)
gam.check(m_calib_0)


out_v1$calib_df <- out_v1$calib_df |>
  mutate(
    p_hat = predict(m_calib_0, type = "response"),
    conf_bin = cut(conf, breaks = seq(0, 1, by = 0.05))
  )

out_v1$calib_df |>
  group_by(conf_bin) |>
  summarise(
    mean_conf = mean(conf),
    empirical_p = mean(y),
    pred_p = mean(p_hat),
    n = n()
  ) |>
  ggplot(aes(mean_conf, empirical_p)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(
    x = "Mean detector confidence",
    y = "Observed P(true detection)"
  )

at_calibrated <- at |>
  filter(spname == "scallop") |>
  mutate(
    p_true = predict(m_calib_0, type = "response")
  )

image_counts <- at_calibrated |>
  group_by(image_id) |>
  summarise(
    expected_scallops = sum(p_true),
    n_detections = n()
  )

image_counts_m <- at_calibrated |>
  filter(truedetect) |>
  group_by(image_id) |>
  summarise(
    expected_scallops = sum(p_true),
    n_manual = n()
  )

dev.off()
plot(expected_scallops~n_detections, data=image_counts)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
dev.off()
plot(expected_scallops~n_manual, data=image_counts_m)
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

gam.check(m_calib_0)


######################################################################
# P(detection | confidence, depth strata)
######################################################################

depth_bins <- tibble(
  depth_bin = c("70m", "80m", "90m"),
  depth_mid = c(70, 80, 90)
)

conf_grid <- seq(0, 1, by = 0.01)

pred_grid <- expand_grid(
  conf = conf_grid,
  depth_bin = depth_bins$depth_bin
) |>
  left_join(depth_bins, by = "depth_bin") |>
  rename(bottom_depth = depth_mid)

pred_grid <- pred_grid |>
  mutate(
    p_hat = predict(
      m_calib_0,
      newdata = pred_grid,
      type = "response"
    )
  )

pred <- predict(
  m_calib_0,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE
)

pred_grid <- pred_grid |>
  mutate(
    eta = pred$fit,
    se  = pred$se.fit,
    p_lo = plogis(eta - 1.96 * se),
    p_hi = plogis(eta + 1.96 * se)
  )

library(ggplot2)

p_calib <- ggplot(
  pred_grid,
  aes(x = conf, y = p_hat, color = depth_bin)
) +
  geom_line(linewidth = 1.2, lty = 1) +
  theme_minimal() +
  labs(
    title = "Model_v2 - MAB",
    x = "YOLO detection confidence",
    y = "Probability of true detection",
    color = "Depth"
  )

p_calib

ggsave(p_calib, filename = "~/Downloads/p_calib_v2_MAB.jpg", width=6, height=6)

