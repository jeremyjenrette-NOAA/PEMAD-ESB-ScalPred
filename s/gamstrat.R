source("./fitfunc.R")
source("./gamfunc.R")

dat_split = read.csv("../data/raw/dataset_split_2224.csv")
table(dat_split$is_test_gam_train)
table(dat_split$is_test_gam_test)
colnames(dat_split)

model = Casv2strat2224
model_name = deparse(substitute(Casv2strat2224))

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
pred = run_model_pipeline(model, model_name, gams = gams, 
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
  model_label = expression("Dataset: 2022 - 2024, Detection Model: " * bold("Cascade R-CNN"))
)
p_zoomed

save_cal(p_zoomed, id = paste("2224_",model_name, "_count_strat",sep=""), 
         outdir = "../figures/diag2/", format = "png", width = 9)
