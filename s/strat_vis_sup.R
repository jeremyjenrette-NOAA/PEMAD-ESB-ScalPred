library(maps)
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(patchwork)
library(marmap)
source("./gamfunc.R")

# ======================================================================
# 1. Prepare Plotting Data
# ======================================================================
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
  plot.title = element_text(hjust = 0.5, margin = margin(b = 5)), 
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
  strip.background = element_rect(color = "black", fill = "grey90", linewidth = 0.5), 
  plot.tag = element_text(size = 16, face = "bold"),
  # This aggressively trims the default white space around every plot:
  plot.margin = margin(t = 2, r = 5, b = 2, l = 2) 
)

split_colors <- c(
  "Detection train" = "black", 
  "Calibration train" = "#2C7FB8", 
  "Calibration test" = "#D95F0E", 
  "Other" = "grey80"
)

# ======================================================================
# 2. Detailed Geospatial Maps (GB and MAB)
# ======================================================================
world <- map_data("world")

create_region_map <- function(target_region, title_text, tag_text, show_legend = TRUE) {
  plot_data <- dat_plot %>% filter(region == target_region)
  
  p <- ggplot() +
    geom_polygon(
      data = world, aes(x = long, y = lat, group = group),
      fill = "grey85", color = "grey60", linewidth = 0.2
    ) +
    geom_path(
      data = plot_data %>% arrange(transect_group, transect_position),
      aes(x = longitude, y = latitude, group = transect_group),
      color = "grey50", linewidth = 0.3, alpha = 0.5
    ) +
    geom_point(
      data = plot_data,
      aes(x = longitude, y = latitude, color = split_stage, size = pmin(n_annotations, 15)),
      shape = 15, alpha = 0.75
    ) +
    coord_quickmap(
      xlim = range(plot_data$longitude, na.rm = TRUE) + c(-0.2, 0.2),
      ylim = range(plot_data$latitude, na.rm = TRUE) + c(-0.2, 0.2)
    ) +
    scale_size_continuous(name = "Manual annotations\n(max = 15+)", range = c(0.4, 3.2)) +
    scale_color_manual(values = split_colors) +
    facet_wrap(~ year, nrow = 1) + 
    labs(title = title_text, x = "Longitude", y = "Latitude", color = "Datasets", tag = tag_text) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", margin = margin(t=4, b=4)),
      axis.text.x = element_text(angle = 45, hjust = 1),
      # Pulls the legend tightly into the plot boundary to reduce right-side white space:
      legend.margin = margin(l = 0, r = 0),
      legend.box.margin = margin(0, 0, 0, -5)
    ) +
    my_theme
  
  if(!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# Tag C for GB (Bottom), Tag B for MAB (Top)
p_map_gb  <- create_region_map("Georges Bank", "Georges Bank", tag_text = "C", show_legend = FALSE)
p_map_mab <- create_region_map("Mid-Atlantic Bight", "Mid-Atlantic Bight", tag_text = "B", show_legend = TRUE)

# ======================================================================
# 3. High-Detail Bathymetric Overview Map
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
  annotate("text", x = gb_bounds$xmin + 2.2, y = gb_bounds$ymax + 0.4, label = "GB", color = "red", size = 4, fontface = "bold") +
  annotate("text", x = mab_bounds$xmin + 2.2, y = mab_bounds$ymin - 0.4, label = "MAB", color = "red", size = 4, fontface = "bold") +
  coord_quickmap(xlim = c(-78, -64), ylim = c(36, 44), expand = FALSE) +
  labs(x = "", y = "", tag = "A") +
  theme_minimal(base_size = 12) +
  ggtitle("Study Regions") +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(), 
    axis.title = element_blank()
  ) +
  my_theme

# ======================================================================
# 4. Final Patchwork Assembly
# ======================================================================
# Increased the overview map ratio (from 2:4 to 1.2:1.8) to make it larger and reduce gap
top_row <- p_overview + p_map_mab + plot_layout(widths = c(1.2, 1.8), guides = "collect")

map_assembly <- top_row / p_map_gb 

print(map_assembly)

save_cal(map_assembly, id = paste0("2226_spatial_strat"), format = "png", 
         outdir = "../figures/diag3/", width = 17, height = 10)
