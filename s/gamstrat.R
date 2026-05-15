source("./fitfunc.R")
source("./gamfunc.R")
source("./procfunc.R")

load("../data/processed/Casv2strat2224.RData")
load("../data/processed/YOLOv12strat2224.RData")

dat_split = read.csv("../data/raw/dataset_split_2224.csv")
table(dat_split$is_test_gam_train)
table(dat_split$is_test_gam_test)
colnames(dat_split)

model = YOLOv12strat2224
model_name = deparse(substitute(YOLOv12strat2224))

# trains GAMs on is_test_gam_train (65% of detection test - stratified space x density)
gams <- fit_calibration_gams(
  model = model,
  model_name = model_name,
  dat_split = dat_split,
  use_strat = TRUE
)

summary(gams$GB); summary(gams$MAB)
AIC(gams$GB); AIC(gams$MAB)

# tests GAMs on is_test_gam_test (35% of detection test - stratified space x density)
pred = run_model_pipeline(model, model_name, gams = gams, test = TRUE,
                               dat_split = dat_split, use_strat = TRUE)

# ---- Step 1: rename regions FIRST ----
img_df <- pred$img %>%
  mutate(
    region = recode(region,
                    "GB" = "Georges Bank",
                    "MAB" = "Mid-Atlantic Bight"
    )
  )

# ---- Step 2: compute counts AFTER renaming ----
region_counts <- img_df %>%
  count(region, name = "n_images")

total_n <- sum(region_counts$n_images)

# ---- Step 3: create labels ----
region_labels <- region_counts %>%
  mutate(region_label = paste0(region, " (n = ", n_images, ")"))

# ---- Step 4: apply labels ----
img_df <- img_df %>%
  left_join(region_labels, by = "region") %>%
  mutate(region = region_label) %>%
  select(-region_label)

metrics <- pred$metrics %>%
  mutate(
    region = recode(region,
                    "GB" = "Georges Bank",
                    "MAB" = "Mid-Atlantic Bight"
    )
  ) %>%
  left_join(region_labels, by = "region") %>%
  mutate(region = region_label) %>%
  select(-region_label)

p_zoomed <- plot_image_level_fit_zoom(
  # img_df = img_all |> dplyr::mutate(region = "All Survey Regions"),
  img_df = img_df, # |> filter(region == "GB") |> mutate(region = "George's Banks"),
  metrics_df = metrics[1:2,], # |> mutate(region = "George's Banks"),
  model_name = model_name,
  zoom_q=1,
  plot_title = "True abundance vs. Σ calibrated detection probabilities per image",
  model_label = expression("Dataset: 2022 - 2024, Detection Model: " * bold("YOLOv12"))
)
p_zoomed

p_resid = pred$img %>%
  mutate(residual = predicted_number - true_number) %>%
  ggplot(aes(true_number, residual)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs True Count",
       subtitle = "YOLOv12")
p_resid
# strong, but how to visualize after calibration?
p_fn = ggplot(pred$img, aes(true_number, false_negative)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  labs(title = "False negatives vs true abundance")
p_fn
# for each break of p(detection), what is the average true positive count?
pred$calib_df %>%
  mutate(bin = cut(pred_p, breaks = seq(0,1,0.1))) %>%
  group_by(bin) %>%
  summarise(
    observed = mean(y),
    predicted = mean(pred_p)
  ) %>%
  mutate(
    r2   = summary(stats::lm(observed ~ predicted))$adj.r.squared,
  ) %>%
  ggplot(aes(predicted, observed)) +
  geom_point() +
  geom_text(
    aes(x = -Inf, y = Inf, label = paste0("R²=", round(r2, 2) ) ),
    hjust = -0.5, vjust = 1.5,
    inherit.aes = FALSE
  ) +
  geom_abline(slope = 1, intercept = 0) +
  labs(title = "Calibration curve",
       subtitle = "for each seq of p(detection), what is the average true positive count?")

summary_df <- data.frame(
  metric = c("True", "Raw detector", "Calibrated", "F1-cutoff"),
  value  = c(
    sum(pred$img$n_annotations, na.rm = TRUE),
    sum(pred$calib_df$conf, na.rm = TRUE),
    sum(pred$calib_df$pred_p, na.rm = TRUE),
    sum(pred$img$predicted_f1_number, na.rm = TRUE)
  )
)

p_sum = ggplot(summary_df, aes(x = metric, y = value, fill = metric)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = summary_df$value[1], linetype = "dashed", color = "black") +
  geom_text(aes(label = round(value, 0)), vjust = -0.5) +
  theme_minimal() +
  labs(
    title = "Total scallop abundance",
    subtitle = "YOLOv12",
    x = "",
    y = "Total count"
  ) +
  theme(legend.position = "none")
p_sum

save_cal(p_zoomed, id = paste("2224_",model_name, "_count_strat",sep=""), 
         outdir = "../figures/diag2/", format = "png", width = 9)
save_cal(p_resid, id = paste("2224_",model_name, "_resid_strat",sep=""), 
         outdir = "../figures/diag2/", format = "png", width = 5, height = 4.5)
save_cal(p_fn, id = paste("2224_",model_name, "_fn_strat",sep=""), 
         outdir = "../figures/diag2/", format = "png", width = 5, height = 4.5)
save_cal(p_sum, id = paste("2224_",model_name, "_sum_strat",sep=""), 
         outdir = "../figures/diag2/", format = "png", width = 5, height = 4.5)


