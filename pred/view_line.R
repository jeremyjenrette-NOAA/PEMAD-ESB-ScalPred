library(dplyr)
library(magick)

visualize_backcalculated_lengths_images <- function(
    pred_calibrated,
    image_names,
    img_dir,
    split = c("right", "left"),
    out_dir = "length_backcalc_check",
    draw_boxes = TRUE,
    draw_lines = TRUE,
    box_color = "red",
    line_color = "cyan",
    processed_border_color = "yellow",
    label_lines = TRUE
) {
  
  split <- match.arg(split)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  req <- c("imagename", "tlx", "tly", "brx", "bry")
  miss <- req[!req %in% names(pred_calibrated)]
  if (length(miss) > 0) {
    stop("pred_calibrated is missing required columns: ",
         paste(miss, collapse = ", "))
  }
  
  image_names <- basename(image_names)
  
  det_sub <- pred_calibrated %>%
    mutate(imagename = basename(imagename)) %>%
    filter(imagename %in% image_names) %>%
    mutate(
      tlx = as.numeric(tlx),
      tly = as.numeric(tly),
      brx = as.numeric(brx),
      bry = as.numeric(bry),
      width_px  = brx - tlx,
      height_px = bry - tly,
      length_px = (width_px + height_px) / 2,
      center_x  = (tlx + brx) / 2,
      center_y  = (tly + bry) / 2,
      line_x1   = center_x - length_px / 2,
      line_x2   = center_x + length_px / 2,
      line_y1   = center_y,
      line_y2   = center_y,
      length_mm = if ("millimeter_per_pixel" %in% names(.)) {
        length_px * millimeter_per_pixel
      } else {
        NA_real_
      }
    )
  
  if (nrow(det_sub) == 0) {
    stop("No matching processed detections found for the supplied image_names.")
  }
  
  split_image_exact <- function(img, side = "right") {
    info <- image_info(img)
    w <- info$width[1]
    h <- info$height[1]
    half_w <- w %/% 2
    
    if (side == "left") {
      image_crop(img, paste0(half_w, "x", h, "+0+0"))
    } else {
      image_crop(img, paste0(half_w, "x", h, "+", half_w, "+0"))
    }
  }
  
  draw_backcalc_magick <- function(img, dets, label_text = NULL) {
    info <- image_info(img)
    w <- info$width[1]
    h <- info$height[1]
    
    img_draw <- image_draw(img)
    
    rect(
      xleft = 1, ybottom = 1, xright = w - 1, ytop = h - 1,
      border = processed_border_color, lwd = 8
    )
    
    if (!is.null(dets) && nrow(dets) > 0) {
      for (j in seq_len(nrow(dets))) {
        
        if (draw_boxes) {
          rect(
            xleft = dets$tlx[j],
            ybottom = dets$tly[j],
            xright = dets$brx[j],
            ytop = dets$bry[j],
            border = box_color,
            lwd = 3
          )
        }
        
        if (draw_lines) {
          segments(
            x0 = dets$line_x1[j],
            y0 = dets$line_y1[j],
            x1 = dets$line_x2[j],
            y1 = dets$line_y2[j],
            col = line_color,
            lwd = 4
          )
          
          points(
            x = dets$center_x[j],
            y = dets$center_y[j],
            col = line_color,
            pch = 16,
            cex = 0.7
          )
        }
        
        if (label_lines) {
          conf_txt <- if ("conf" %in% names(dets) && !is.na(dets$conf[j])) {
            sprintf(" conf=%.3f", dets$conf[j])
          } else ""
          
          mm_txt <- if ("length_mm" %in% names(dets) && !is.na(dets$length_mm[j])) {
            sprintf(" | %.1f mm", dets$length_mm[j])
          } else ""
          
          text(
            x = dets$tlx[j],
            y = max(20, dets$tly[j] - 10),
            labels = paste0("L=", round(dets$length_px[j], 1), " px", mm_txt, conf_txt),
            col = line_color,
            cex = 0.9,
            adj = c(0, 1)
          )
        }
      }
    }
    
    if (!is.null(label_text)) {
      text(
        x = 20,
        y = 35,
        labels = label_text,
        col = "white",
        bg = "black",
        cex = 1.0,
        adj = c(0, 0.5)
      )
    }
    
    dev.off()
    img_draw
  }
  
  output_files <- character(0)
  
  for (img_name in unique(image_names)) {
    img_path <- file.path(img_dir, img_name)
    
    if (!file.exists(img_path)) {
      warning(sprintf("Missing image: %s", img_path))
      next
    }
    
    img_full  <- image_read(img_path)
    img_split <- split_image_exact(img_full, split)
    
    dets_i <- det_sub %>%
      filter(imagename == img_name)
    
    split_info <- image_info(img_split)
    split_w <- split_info$width[1]
    split_h <- split_info$height[1]
    
    dets_i <- dets_i %>%
      filter(
        tlx >= 0, tly >= 0,
        brx <= split_w, bry <= split_h,
        brx > tlx, bry > tly,
        line_x1 >= 0, line_x2 <= split_w,
        line_y1 >= 0, line_y2 <= split_h
      )
    
    label_text <- paste0(
      img_name,
      " | split=", split,
      " | n=", nrow(dets_i)
    )
    
    img_annot <- draw_backcalc_magick(
      img = img_split,
      dets = if (nrow(dets_i) > 0) dets_i else NULL,
      label_text = label_text
    )
    
    out_file <- file.path(out_dir, paste0("annotated_", img_name))
    image_write(img_annot, path = out_file, format = tools::file_ext(out_file))
    
    output_files <- c(output_files, out_file)
  }
  
  invisible(list(
    detections_used = det_sub,
    image_files = output_files,
    split = split
  ))
}

img_dir = "/Volumes/PortableSSD/saltnoaa/images/overlap_2023/"

img_files <- list.files(
  img_dir,
  pattern = "\\.png$",   # adjust if needed
  full.names = FALSE
)

img_inventory_subset <- img_inventory %>%
  filter(imagename %in% img_files) %>%
  filter(process_image == TRUE)

res2 <- visualize_backcalculated_lengths_images(
  pred_calibrated = pred_calibrated %>% filter(conf > 0.05),
  image_names = img_inventory_subset$imagename,
  img_dir = img_dir,
  split = "right",
  out_dir = "~/saltnoaa/PEMAD-ESB-ScalPred/data/images/length_backcalc_multi"
)
