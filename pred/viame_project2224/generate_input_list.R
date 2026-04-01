library(dplyr)
library(readr)

# Load inventory
img_inventory <- readRDS("../img_inventory_out/img_inventory_2022_2024_compiled.rds")

# 1) subset to images flagged for processing
img_inventory_proc <- img_inventory %>%
  filter(process_image == TRUE)

# 2) convert Linux-style /z/... paths to Windows-style Z:\...
to_windows_z_path <- function(x) {
  x %>%
    gsub("^/z/", "Z:/", ., ignore.case = TRUE) %>%
    gsub("/", "\\\\", ., fixed = TRUE)
}

img_inventory_proc <- img_inventory_proc %>%
  mutate(actual_path_windows = to_windows_z_path(actual_path))

# 3) function to generate an input_list text file for a supplied year
write_viame_input_list <- function(year,
                                   img_inventory_df,
                                   out_dir = ".",
                                   only_processable = TRUE,
                                   filename_prefix = "input_list") {
  
  df <- img_inventory_df
  
  if (only_processable && "process_image" %in% names(df)) {
    df <- df %>% filter(process_image == TRUE)
  }
  
  df_year <- df %>%
    filter(year == !!year) %>%
    mutate(actual_path_windows = to_windows_z_path(actual_path))
  
  if (nrow(df_year) == 0) {
    stop(paste("No images found for year", year))
  }
  
  out_file <- file.path(out_dir, paste0(filename_prefix, "_", year, ".txt"))
  
  write_lines(df_year$actual_path_windows, out_file)
  
  message("Wrote ", nrow(df_year), " image paths to: ", out_file)
  
  return(invisible(df_year))
}

# View the subset you care about
# View(subset(img_inventory, process_image == TRUE))

# Write 2024 input list
write_viame_input_list(
  year = 2024,
  img_inventory_df = img_inventory,
  out_dir = "../viame_project2224"
)
