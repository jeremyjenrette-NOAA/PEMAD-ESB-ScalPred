library(tidyverse)
#============================================================#
# Functions: Precision–Recall, F1 evaluation, save plot 
#============================================================#
evaluate_pr_curve <- function(
    calib_df,
    img_df,
    conf_grid = seq(0, 1, by = 0.0005),
    model_id  = "model"
) {
  
  purrr::map_dfr(conf_grid, function(th) {
    
    df_th <- calib_df |> dplyr::filter(conf >= th)
    
    TP <- sum(df_th$truedetect, na.rm = TRUE)
    FP <- sum(!df_th$truedetect, na.rm = TRUE)
    
    auto_true_img <- df_th |>
      dplyr::filter(truedetect) |>
      dplyr::count(image_id, name = "n_auto_true")
    
    FN <- img_df |>
      dplyr::left_join(auto_true_img, by = "image_id") |>
      dplyr::mutate(
        n_auto_true = tidyr::replace_na(n_auto_true, 0L),
        FN = pmax(n_manual - n_auto_true, 0L)
      ) |>
      dplyr::summarise(FN = sum(FN, na.rm = TRUE)) |>
      dplyr::pull(FN)
    
    precision <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
    recall    <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
    f1        <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_
    else 2 * precision * recall / (precision + recall)
    
    tibble::tibble(
      model     = model_id,
      conf      = th,
      TP        = TP,
      FP        = FP,
      FN        = FN,
      precision = precision,
      recall    = recall,
      f1        = f1
    )
  })
}

evaluate_pr_models <- function(
    res_list,                         # one res OR list of res
    conf_grid = seq(0, 1, by = 0.0005),
    title = "Precision–Recall Curves",
    stratify_region = TRUE,           # NEW: if FALSE, pool GB+MAB into one PR curve per model
    pooled_label = ""             # NEW: suffix for pooled model_id
) {
  # allow passing a single res directly
  if (!is.list(res_list)) stop("res_list must be a list")
  
  is_single_res <- all(c("calib_df","img_df") %in% names(res_list)) ||
    any(grepl("_GB_calib$|_MAB_calib$", names(res_list)))
  
  if (is_single_res) res_list <- list(res_list)
  
  # helper: get a named element by exact name or by suffix match
  get_el <- function(res, key_exact, key_suffix) {
    if (key_exact %in% names(res)) return(res[[key_exact]])
    hit <- names(res)[grepl(paste0(key_suffix, "$"), names(res))]
    if (length(hit) == 1) return(res[[hit]])
    stop("Couldn't uniquely find: ", key_exact, " or suffix: ", key_suffix)
  }
  
  pr_all <- purrr::map_dfr(res_list, function(res) {
    # infer model base name from names like "Cas2024v2_GB_calib"
    model_base <- sub("_.*$", "", names(res)[1])
    
    # If someone passes already-pooled objects (calib_df/img_df), respect that:
    if (all(c("calib_df", "img_df") %in% names(res)) && !any(grepl("_GB_calib$|_MAB_calib$", names(res)))) {
      pr <- evaluate_pr_curve(
        calib_df  = res[["calib_df"]],
        img_df    = res[["img_df"]],
        conf_grid = conf_grid,
        model_id  = paste0(model_base, pooled_label)
      )
      return(pr)
    }
    
    GB_calib  <- get_el(res, "GB_calib",  "_GB_calib")
    GB_img    <- get_el(res, "GB_img",    "_GB_img")
    MAB_calib <- get_el(res, "MAB_calib", "_MAB_calib")
    MAB_img   <- get_el(res, "MAB_img",   "_MAB_img")
    
    if (isTRUE(stratify_region)) {
      pr_gb <- evaluate_pr_curve(
        calib_df  = GB_calib,
        img_df    = GB_img,
        conf_grid = conf_grid,
        model_id  = paste0(model_base, "_GB")
      )
      
      pr_mab <- evaluate_pr_curve(
        calib_df  = MAB_calib,
        img_df    = MAB_img,
        conf_grid = conf_grid,
        model_id  = paste0(model_base, "_MAB")
      )
      
      return(dplyr::bind_rows(pr_gb, pr_mab))
    }
    
    # pooled across regions
    calib_all <- dplyr::bind_rows(GB_calib, MAB_calib)
    img_all   <- dplyr::bind_rows(GB_img,   MAB_img)
    
    pr_all_model <- evaluate_pr_curve(
      calib_df  = calib_all,
      img_df    = img_all,
      conf_grid = conf_grid,
      model_id  = paste0(model_base, pooled_label)
    )
    
    pr_all_model
  })
  
  legend_title <- if (isTRUE(stratify_region)) "Mdl_Region" else "Model"
  
  p_pr <- ggplot2::ggplot(pr_all, ggplot2::aes(x = recall, y = precision, color = model)) +
    ggplot2::geom_path(linewidth = 1.2, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = "Recall",
      y = "Precision",
      color = legend_title
    )
  
  list(pr_all = pr_all, p_pr = p_pr)
}

save_pr_f1 <- function(p_pr, p_f1, id, outdir = "~/Downloads",
                       width = 6, height = 6) {
  
  ggsave(
    plot = p_pr,
    filename = file.path(outdir, paste0("pr_", id, ".jpg")),
    width = width,
    height = height
  )
  
  ggsave(
    plot = p_f1,
    filename = file.path(outdir, paste0("f1_", id, ".jpg")),
    width = width,
    height = height
  )
}