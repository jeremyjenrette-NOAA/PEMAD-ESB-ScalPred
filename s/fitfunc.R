fit_calibration_gams <- function(model, model_name, dat_split, use_strat = FALSE) {
  
  # --- pull raw calibration data ---
  gb_calib  <- get_data(model, model_name, "GB")
  mab_calib <- get_data(model, model_name, "MAB")
  
  # --- apply stratified GAM training filter ---
  if (use_strat) {
    
    # ensure consistent naming
    dat_split$imagename <- basename(dat_split$imagename)
    
    # join + filter helper
    filter_strat <- function(df) {
      df$imagename <- basename(df$imagename)
      
      df2 <- merge(
        df,
        dat_split[, c("imagename", "stratum", "is_test_gam_train")],
        by = "imagename",
        all.x = TRUE
      ) %>%
        mutate(is_test_gam_train = if_else(is_test_gam_train == "True", TRUE, FALSE))
      
      # keep ONLY GAM training set
      df2 <- df2[df2$is_test_gam_train == TRUE, ]
      
      df2$stratum = factor(df2$stratum)
      
      return(df2)
    }
    
    gb_calib  <- filter_strat(gb_calib)
    mab_calib <- filter_strat(mab_calib)
    
    cat("GB n:", nrow(gb_calib), "\n")
    cat("MAB n:", nrow(mab_calib), "\n")
    
    # --- fit STRATIFIED GAMs ---
    m_gb <- mgcv::gam(
      y ~ s(conf, bottom_depth, k = 5) +
        # s(conf, boxsize, k = 7) +
        s(latitude, longitude, k = 7) + 
        s(altitude, millimeter_per_pixel, k = 7),
      family = binomial(),
      data = gb_calib,
      method = "REML"
    )
    
    m_mab <- mgcv::gam(
      y ~ s(conf, bottom_depth, k = 5) +
        # s(conf, boxsize, k = 7) +
        s(latitude, longitude, k = 7) + 
        s(altitude, millimeter_per_pixel, k = 7),
      family = binomial(),
      data = mab_calib,
      method = "REML"
    )
    
    m_comb <- mgcv::gam(
      y ~ s(conf, bottom_depth, k = 5) +
        # s(conf, boxsize, k = 7) +
        s(latitude, longitude, k = 7) + 
        s(altitude, backscatter, k = 7),
      family = binomial(),
      data = bind_rows(gb_calib,mab_calib),
      method = "REML"
    )
    
  } else {
  
  # --- sanity check ---
  cat("GB n:", nrow(gb_calib), "\n")
  cat("MAB n:", nrow(mab_calib), "\n")
  
  # --- fit GAMs ---
  m_gb <- mgcv::gam(
    y ~ s(conf, bottom_depth, k = 5) +
      s(altitude, backscatter, k = 7) +
      s(latitude, longitude, k = 7), # + s(stratum, bs = "re"),
    family = binomial(),
    data = gb_calib,
    method = "REML"
  )
  
  m_mab <- mgcv::gam(
    y ~ s(conf, bottom_depth, k = 5) +
      s(altitude, backscatter, k = 7) +
      s(latitude, longitude, k = 7), # + s(stratum, bs = "re"),
    family = binomial(),
    data = mab_calib,
    method = "REML"
  )
  
  }
  
  list(
    model_name = model_name,
    GB  = m_gb,
    MAB = m_mab
  )
}

# DEPRECATED!!!
fitimage_calibration_mod <- function(model, model_name, p_detect) {
  
  gb_img  <- get_data(model, model_name, "GB", imglvl = TRUE)
  mab_img <- get_data(model, model_name, "MAB", imglvl = TRUE)
  
  # Replace NA auto counts
  gb_img$n_auto[is.na(gb_img$n_auto)] <- 0
  mab_img$n_auto[is.na(mab_img$n_auto)] <- 0
  
  # Define FN presence
  gb_img$fn_pres  <- as.integer(gb_img$n_auto < gb_img$n_manual)
  mab_img$fn_pres <- as.integer(mab_img$n_auto < mab_img$n_manual)
  
  # Density term
  gb_img$auto_density_log  <- log(gb_img$auto_density + 1e-6)
  mab_img$auto_density_log <- log(mab_img$auto_density + 1e-6)
  
  reg_comb = rbind(mab_img, gb_img)
  
  reg_comb <- reg_comb %>%
    inner_join(p_detect, by = "image_id")
  
  #------------------------------#
  # GB model
  #------------------------------#
  
  m_gb_fn <- mgcv::gam(
    fn_pres ~
      s(millimeter_per_pixel, k = 7),
    family = binomial(),
    data = gb_img,
    method = "REML"
  )
  
  #------------------------------#
  # MAB model
  #------------------------------#
  
  m_mab_fn <- mgcv::gam(
    fn_pres ~
      s(millimeter_per_pixel, k = 7),
    family = binomial(),
    data = mab_img,
    method = "REML"
  )
  
  #------------------------------#
  # Combined region model
  #------------------------------#
  
  # m_comb <- mgcv::gam(
  #   fn_pres ~ s(predicted_number),
  #   family = binomial(),
  #   data = reg_comb,
  #   method = "REML"
  # )
  
  m_comb <- gam(
    fn_pres ~ 
      s(predicted_number, k=7),
    family = binomial(),
    data = reg_comb,
    method = "REML"
  )
  
  list(
    model_name = model_name,
    GB  = m_gb_fn,
    MAB = m_mab_fn,
    comb = m_comb,
    data = reg_comb
  )
}


run_model_pipeline <- function(model, model_name, gams, test = TRUE,
                               use_strat = FALSE, dat_split = NULL) {
  
  # gams <- fit_calibration_gams(
  #   model = model,
  #   model_name = model_name
  # )
  
  # pred_gb <- predict_calibration_by_depth(
  #   gam = gams$GB,
  #   calib_df = get_data(model, model_name, "GB"),
  #   depth_breaks = c(40, 80, 120),
  #   region = "GB",
  #   model_name = model_name
  # )
  # 
  # pred_mab <- predict_calibration_by_depth(
  #   gam = gams$MAB,
  #   calib_df = get_data(model, model_name, "MAB"),
  #   depth_breaks = c(40, 50, 60, 70),
  #   region = "MAB",
  #   model_name = model_name
  # )
  
  dat_split_gb <- dat_split %>% filter(latitude > 40)
  dat_split_mab <- dat_split %>% filter(latitude < 40)
  
  out_pr <- evaluate_pr_models(list(model), stratify_region = TRUE)
  best_pts <- out_pr$pr_all %>%
    group_by(model) %>%
    filter(f1 == max(f1, na.rm = TRUE)) %>%
    slice_max(conf, n = 1) %>%
    ungroup()
  
  f1conf_gb  <- best_pts$conf[1]
  f1conf_mab <- best_pts$conf[2]
  
  res_gb <- compute_image_level_counts(
    calib_df  = get_data(model, model_name, "GB", imglvl = FALSE),
    gam       = gams$GB,
    region    = "GB",
    model_name = model_name,
    f1_thresh = f1conf_gb,
    dat_split = dat_split_gb,
    use_is_test_gam_test = test,
    use_strat = use_strat
  )
  
  res_mab <- compute_image_level_counts(
    calib_df  = get_data(model, model_name, "MAB", imglvl = FALSE),
    gam       = gams$MAB,
    region    = "MAB",
    model_name = model_name,
    f1_thresh = f1conf_mab,
    dat_split = dat_split_mab,
    use_is_test_gam_test = test,
    use_strat = use_strat
  )
  
  img = bind_rows(res_gb$img, res_mab$img)
  calib_df = bind_rows(res_gb$calib_df, res_mab$calib_df)
  
  metrics = bind_rows(res_gb$metrics, res_mab$metrics) %>%
    mutate(region = c("GB", "MAB"))
  metrics = bind_rows(metrics, 
                      compute_combined_metrics(img,region_name = "All Survey Regions"))
  
  list(img = img, metrics = metrics, calib_df = calib_df)
}

attach_metadata_to_images <- function(img_df, meta_df) {
  
  library(dplyr)
  library(stringr)
  
  # 1. Standardize column names to lowercase
  # meta_df <- meta_df %>%
  #   rename_with(tolower) %>%
  #   mutate(image_id = std_image_id(imagename))
  
  # 3. Ensure img_df also has lowercase names (safe)
  img_df <- img_df %>%
    rename_with(tolower)
  
  # 4. Select only desired metadata columns (if they exist)
  keep_cols <- c(
    "imagename",
    "latitude",
    "longitude",
    "bottom_depth",
    "chlorophyll",
    "cdom",
    "altitude",
    "field_of_view_sq_meter",
    "t",
    "s",
    "heading",
    "o2",
    "v_depth",
    "millimeter_per_pixel",
    "backscatter",
    "pitch",
    "roll"
  )
  
  meta_selected <- meta_df %>%
    select(any_of(keep_cols))
  
  # 5. Join
  out <- img_df %>%
    left_join(meta_selected, by = "imagename")
  
  return(out)
}
