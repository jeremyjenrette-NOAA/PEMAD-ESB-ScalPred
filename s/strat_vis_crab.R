library(maps)
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(patchwork)
library(marmap)
source("./gamfunc.R")

# ======================================================================
# 1. Load and Prepare Plotting Data
# ======================================================================
dat_split <- read.csv("../data/raw/dataset_split_crab.csv")

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

dat_plot <- dat_split %>%
  mutate(
    split_stage = case_when(
      dataset == "train" ~ "Detection train",
      dataset == "test_GAM_train" ~ "Calibration train",
      dataset == "test_GAM_test" ~ "Calibration test",
      TRUE ~ "Other"
    ),
    region = ifelse(longitude >= -71, "Georges Bank", "Mid-Atlantic Bight"),
    year = as.factor(year)
  ) %>%
  filter(!is.na(transect_group), !is.na(transect_position))

# --- Define a Common Theme ---
my_theme <- theme(
  plot.title = element_text(hjust = 0.5, margin = margin(b = 13), size = 21), 
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
  strip.background = element_rect(color = "black", fill = "grey90", linewidth = 0.5), 
  
  plot.tag = element_text(size = 21, face = "plain"),
  
  # c(x, y) coordinates from 0 to 1 relative to the entire plot area
  # Changing x from 0.02 to 0.08 moves it a bit to the right
  plot.tag.position = c(0.045, 0.97), 
  
  plot.margin = margin(t = 4, r = 5, b = 6, l = 4) 
)

split_colors <- c(
  "Detection train" = "black", 
  "Calibration train" = "#2C7FB8", 
  "Calibration test" = "#D95F0E", 
  "Other" = "grey80"
)

world <- map_data("world")

# ======================================================================
# 2. Regional Geospatial Maps (Panels A & B)
# ======================================================================
create_region_map <- function(target_region, title_text, tag_text, show_legend = TRUE) {
  plot_data <- dat_plot %>% filter(region == target_region)
  
  p <- ggplot() +
    geom_polygon(
      data = world, aes(x = long, y = lat, group = group),
      fill = "grey85", color = "grey60", linewidth = 0.2
    ) +
    geom_point(
      data = plot_data,
      aes(x = longitude, y = latitude, color = split_stage, size = pmin(total_annotations, 15)),
      shape = 15, alpha = 0.75
    ) +
    coord_quickmap(
      xlim = range(plot_data$longitude, na.rm = TRUE) + c(-0.2, 0.2),
      ylim = range(plot_data$latitude, na.rm = TRUE) + c(-0.2, 0.2)
    ) +
    scale_size_continuous(name = "Manual annotations\n(max = 15+)", range = c(0.6, 3.8)) +
    scale_color_manual(values = split_colors) +
    facet_wrap(~ year, nrow = 1) + 
    labs(title = title_text, x = "Longitude", y = "Latitude", color = "Dataset", tag = tag_text) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", margin = margin(t=4, b=4), size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.margin = margin(l = 10, r = 0),
      legend.box.margin = margin(0, 0, 0, 0)
    ) +
    my_theme
  
  if(!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    # Adjust 2: Significantly enlarge the legend text, keys, and symbols
    p <- p + theme(
      legend.title = element_text(size = 20, face = "plain"),
      legend.text = element_text(size = 18),
      legend.key.size = unit(0.9, "cm"),
      legend.spacing.y = unit(0.4, "cm")
    ) +
      guides(
        color = guide_legend(order = 1, override.aes = list(size = 8)),
        size = guide_legend(order = 2)
      )
  }
  
  return(p)
}

p_map_mab <- create_region_map("Mid-Atlantic Bight", "Mid-Atlantic Bight", tag_text = "A", show_legend = TRUE)
p_map_gb  <- create_region_map("Georges Bank", "Georges Bank", tag_text = "B", show_legend = FALSE)

# ======================================================================
# 3. Density Violin Plot (Panel C)
# ======================================================================
p_density <- ggplot(
  dat_plot,
  aes(x = split_stage, y = total_annotations, fill = split_stage)
) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.alpha = 0.15) +
  scale_y_continuous(
    trans = "pseudo_log",
    breaks = c(0, 1, 2, 5, 10, 25, 50, 100, 250)
  ) +
  scale_fill_manual(
    values = c(
      "Detection train" = "#4D4D4D",
      "Calibration train"  = "#2C7FB8",
      "Calibration test"   = "#D95F0E",
      "Other"      = "grey80"
    )
  ) +
  labs(
    title = "Scallop density across datasets",
    x = "",
    y = "Annotations per image",
    tag = "C"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  my_theme

# ======================================================================
# 4. High-Detail Bathymetric Overview Map (Panel D)
# ======================================================================
states <- map_data("state") 

bathy_data <- getNOAA.bathy(lon1 = -78, lon2 = -64, lat1 = 36, lat2 = 44, resolution = 1, keep = TRUE)
bathy_df <- fortify.bathy(bathy_data)

gb_bounds <- dat_plot %>% filter(region == "Georges Bank") %>%
  summarize(xmin = min(longitude, na.rm=T) - 0.2, xmax = max(longitude, na.rm=T) + 0.2,
            ymin = min(latitude, na.rm=T) - 0.2, ymax = max(latitude, na.rm=T) + 0.2)

mab_bounds <- dat_plot %>% filter(region == "Mid-Atlantic Bight") %>%
  summarize(xmin = min(longitude, na.rm=T) - 0.2, xmax = max(longitude, na.rm=T) + 0.2,
            ymin = min(latitude, na.rm=T) - 0.2, ymax = max(latitude, na.rm=T) + 0.2)

p_overview <- ggplot() +
  geom_raster(data = bathy_df %>% filter(z <= 0), aes(x = x, y = y, fill = z)) +
  scale_fill_gradientn(
    colors = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#E0F3F8"),
    values = scales::rescale(c(-4000, -2000, -500, -100, 0)),
    guide = "none" 
  ) +
  geom_contour(data = bathy_df, aes(x = x, y = y, z = z), 
               breaks = c(-50, -200, -500, -1000, -2000), color = "white", alpha = 0.3, linewidth = 0.2) +
  geom_contour(data = bathy_df, aes(x = x, y = y, z = z), 
               breaks = -100, color = "#08519C", linewidth = 0.4, alpha = 0.8) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey85", color = NA) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill = NA, color = "grey50", linewidth = 0.2) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = NA, color = "grey40", linewidth = 0.4) + 
  geom_rect(data = gb_bounds, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "red", fill = NA, linewidth = 0.5) +
  geom_rect(data = mab_bounds, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "red", fill = NA, linewidth = 0.5) +
  annotate("text", x = gb_bounds$xmin + 0.5, y = gb_bounds$ymax + 0.4, label = "Georges Bank", color = "red", size = 6, fontface = "plain", hjust = 0) +
  annotate("text", x = mab_bounds$xmin + 0.5, y = mab_bounds$ymin - 0.4, label = "Mid-Atlantic Bight", color = "red", size = 6, fontface = "plain", hjust = 0) +
  coord_quickmap(xlim = c(-78, -64), ylim = c(36, 44), expand = FALSE) +
  labs(x = "", y = "", tag = "D") +
  theme_minimal(base_size = 12) +
  ggtitle("Study Regions") +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(), 
    axis.title = element_blank()
  ) +
  my_theme

# ======================================================================
# 5. Final Patchwork Assembly and Output Save
# ======================================================================
# Adjust 2: Pair Panel A with an empty spacer to pull it and its legend to the left
row1 <- p_map_mab + plot_spacer() + plot_layout(widths = c(3.2, 0.25))

bottom_row <- p_density + p_overview + plot_layout(widths = c(1.4, 1.6))

# Adjust 3: Optimized row heights (0.5 for row 2 safely eliminates the GB title whitespace gap)
combined_manuscript_plot <- row1 / p_map_gb / bottom_row + 
  plot_layout(heights = c(1.5, 1, 1.9))
combined_manuscript_plot
# Export file as requested
ggsave(combined_manuscript_plot, filename = "../figures/diag_crab1/2426_spatial_strat_final.pdf",
       width = 14, height = 15)
