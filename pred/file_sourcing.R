library(dplyr)
library(stringr)
library(readr)
library(tidyr)

inv_dir <- "img_inventory_out"
out_dir <- "img_inventory_out"

#-----------------------------
# Helper: read one yearly txt
#-----------------------------
read_inventory_txt <- function(year, inv_dir = "img_inventory_out") {
  f <- file.path(inv_dir, paste0("img_inventory_", year, ".txt"))
  
  if (!file.exists(f)) {
    stop("Missing inventory file: ", f)
  }
  
  read_lines(f, progress = FALSE) |>
    tibble(actual_path = _) |>
    mutate(
      year_folder = as.character(year),
      imagename   = basename(actual_path)
    )
}

#-----------------------------
# Parse imagename
# Example:
# 202203.20220530.205541302.129435.png
#-----------------------------
parse_imagename_vec <- function(imagename) {
  x <- str_trim(imagename)
  x <- gsub("\\\\", "/", x)
  x <- basename(x)
  
  parts <- str_split_fixed(x, "\\.", 5)
  
  tibble(
    imagename      = x,
    prefix         = na_if(parts[, 1], ""),
    date_yyyymmdd  = na_if(parts[, 2], ""),
    time_token     = na_if(parts[, 3], ""),
    frame_token    = na_if(parts[, 4], ""),
    extension      = na_if(parts[, 5], "")
  ) |>
    mutate(
      year   = str_sub(date_yyyymmdd, 1, 4),
      month  = str_sub(date_yyyymmdd, 5, 6),
      day    = str_sub(date_yyyymmdd, 7, 8),
      hour   = str_sub(time_token, 1, 2),
      minute = str_sub(time_token, 3, 4),
      second = str_sub(time_token, 5, 6),
      msec   = str_sub(time_token, 7, 9),
      date   = suppressWarnings(as.Date(date_yyyymmdd, format = "%Y%m%d")),
      year_month = if_else(
        !is.na(year) & !is.na(month),
        paste0(year, "-", month),
        NA_character_
      ),
      datetime_str = if_else(
        !is.na(date_yyyymmdd) & !is.na(time_token) & nchar(time_token) >= 6,
        paste0(
          date_yyyymmdd, " ",
          hour, ":", minute, ":", second, ".",
          if_else(is.na(msec), "000", msec)
        ),
        NA_character_
      ),
      frame_num = suppressWarnings(as.numeric(frame_token))
    )
}

#-----------------------------
# Read and combine inventories
#-----------------------------
img_inventory <- bind_rows(
  read_inventory_txt(2022, inv_dir),
  read_inventory_txt(2023, inv_dir),
  read_inventory_txt(2024, inv_dir)
) |>
  distinct(actual_path, .keep_all = TRUE)

parsed <- parse_imagename_vec(img_inventory$imagename)

img_inventory <- img_inventory |>
  left_join(parsed, by = "imagename")

#-----------------------------
# Order images within day and
# mark every 10th image for processing
#-----------------------------
img_inventory <- img_inventory |>
  arrange(date, time_token, frame_num, imagename) |>
  group_by(date) |>
  mutate(
    image_index_in_day = row_number(),
    process_image = image_index_in_day %% 10 == 1
  ) |>
  ungroup()

#-----------------------------
# Summaries
#-----------------------------
summary_year <- img_inventory |>
  count(year, name = "n_images") |>
  arrange(year)

summary_month <- img_inventory |>
  count(year, year_month, name = "n_images") |>
  arrange(year, year_month)

summary_day_sampling <- img_inventory |>
  group_by(date, year, year_month) |>
  summarise(
    n_images = n(),
    n_to_process = sum(process_image, na.rm = TRUE),
    prop_to_process = n_to_process / n_images,
    .groups = "drop"
  ) |>
  arrange(date)

summary_year_sampling <- img_inventory |>
  group_by(year) |>
  summarise(
    n_images = n(),
    n_to_process = sum(process_image, na.rm = TRUE),
    prop_to_process = n_to_process / n_images,
    .groups = "drop"
  ) |>
  arrange(year)

#-----------------------------
# Print summaries
#-----------------------------
cat("\nImages by year:\n")
print(summary_year, n = Inf)

cat("\nImages by month:\n")
print(summary_month, n = Inf)

cat("\nSampling summary by year:\n")
print(summary_year_sampling, n = Inf)

#-----------------------------
# Write outputs
#-----------------------------
write_tsv(
  img_inventory,
  file.path(out_dir, "img_inventory_2022_2024_compiled.tsv")
)

write_lines(
  img_inventory$actual_path,
  file.path(out_dir, "img_inventory_2022_2024_fullpaths.txt")
)

write_lines(
  img_inventory$actual_path[img_inventory$process_image],
  file.path(out_dir, "img_inventory_2022_2024_to_process_fullpaths.txt")
)

write_csv(summary_year, file.path(out_dir, "img_summary_by_year.csv"))
write_csv(summary_month, file.path(out_dir, "img_summary_by_month.csv"))
write_csv(summary_day_sampling, file.path(out_dir, "img_summary_by_day_sampling.csv"))
write_csv(summary_year_sampling, file.path(out_dir, "img_summary_sampling_by_year.csv"))

saveRDS(img_inventory, file.path(out_dir, "img_inventory_2022_2024_compiled.rds"))
