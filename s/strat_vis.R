library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(patchwork)
library(maps) # Explicitly loaded for map_data("world")
source("./gamfunc.R")

# Load the newly formatted dataset
dat_split <- read.csv("../data/raw/dataset_split_2226.csv")

# ======================================================================
# 1. Safe Fallbacks for Transect Logic
# Just in case these were excluded from the final CSV export
# ======================================================================
if(!"transect_group" %in% names(dat_split)) {
  dat_split <- dat_split %>% 
    mutate(transect_group = sub("_block_.*", "", transect_block_id))
}

if(!"transect_position" %in% names(dat_split)) {
  dat_split <- dat_split %>%
    group_by(transect_group) %>%
    mutate(transect_position = row_number() - 1) %>%
    ungroup()
}

# ======================================================================
# 2. Prepare Plotting Data
# ======================================================================
dat_plot <- dat_split %>%
  mutate(
    # Directly map our new dataset column to your plot labels
    split_stage = case_when(
      dataset == "train" ~ "Detection train",
      dataset == "test_GAM_train" ~ "Calibration train",
      dataset == "test_GAM_test" ~ "Calibration test",
      TRUE ~ "Other"
    ),
    density_class = case_when(
      n_annotations == 0 ~ "0 scallops",
      n_annotations == 1 ~ "1 scallop",
      n_annotations == 2 ~ "2 scallops",
      n_annotations > 2 ~ ">2 scallops",
      TRUE ~ NA_character_
    ),
    transect_group = as.factor(transect_group),
    transect_block_id = as.factor(transect_block_id)
  ) %>%
  filter(!is.na(transect_group), !is.na(transect_position))

# ======================================================================
# 3. Transect Block Plot (p_blocks)
# ======================================================================
p_blocks <- ggplot(
  dat_plot,
  aes(
    x = transect_position,
    y = fct_reorder(transect_group, as.numeric(transect_group)),
    fill = split_stage
  )
) +
  geom_tile(
    aes(alpha = pmin(n_annotations, 15)),
    height = 0.8,
    width = 0.95
  ) +
  scale_alpha_continuous(
    name = "Manual annotations\n(max = 15+)",
    range = c(0.4, 1)
  ) +
  scale_fill_manual(
    values = c(
      "Detection train" = "black",
      "Calibration train"  = "#2C7FB8",
      "Calibration test"   = "#D95F0E",
      "Other"      = "grey80"
    )
  ) +
  labs(
    title = "Transect-block stratification of annotated HabCam images",
    subtitle = "1 out of every 50 images annotated",
    x = "Image position within transect",
    y = "Transect",
    fill = "Dataset"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

print(p_blocks)

ggsave(p_blocks, filename = "../figures/ms_figures/2226_strat_block.pdf", 
       width = 13, height = 8)
