library(dplyr)
library(magick)

visualize_overlap_sequence <- function(
    view_overlap,
    pred_detections,
    img_dir,
    start_index = 1,
    n_frames = 20,
    split = c("right", "left"),
    out_dir = "overlap_check",
    montage_file = "montage.png",
    box_color = "red",
    processed_border_color = "yellow",
    label_boxes = TRUE,
    tile = "4x5",
    geometry = "400x300+4+4"
) {
  
  split <- match.arg(split)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  stopifnot("imagename" %in% names(view_overlap))
  stopifnot(all(c("Imagename", "TLx", "TLy", "BRx", "BRy") %in% names(pred_detections)))
  
  # Match ordering used in the sequence
  vo <- view_overlap
  if ("datetime" %in% names(vo)) {
    vo <- vo %>% arrange(datetime, frame_num)
  } else if ("frame_num" %in% names(vo)) {
    vo <- vo %>% arrange(frame_num)
  }
  
  end_index <- min(start_index + n_frames - 1, nrow(vo))
  seq_df <- vo[start_index:end_index, , drop = FALSE]
  
  if (nrow(seq_df) == 0) {
    stop("No rows found for the requested sequence window.")
  }
  
  seq_df <- seq_df %>%
    mutate(processed = imagename %in% unique(pred_detections$Imagename))
  
  det_sub <- pred_detections %>%
    filter(Imagename %in% seq_df$imagename) %>%
    mutate(
      TLx = as.numeric(TLx),
      TLy = as.numeric(TLy),
      BRx = as.numeric(BRx),
      BRy = as.numeric(BRy)
    )
  
  # Reproduce Python split exactly:
  # left  -> img.crop((0, 0, w // 2, h))
  # right -> img.crop((w // 2, 0, w, h))
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
  
  draw_boxes_magick <- function(img, dets, processed = FALSE, label_text = NULL) {
    info <- image_info(img)
    w <- info$width[1]
    h <- info$height[1]
    
    img_draw <- image_draw(img)
    
    if (processed) {
      rect(
        xleft = 1, ybottom = 1, xright = w - 1, ytop = h - 1,
        border = processed_border_color, lwd = 8
      )
    }
    
    if (!is.null(dets) && nrow(dets) > 0) {
      for (j in seq_len(nrow(dets))) {
        rect(
          xleft = dets$TLx[j],
          ybottom = dets$TLy[j],
          xright = dets$BRx[j],
          ytop = dets$BRy[j],
          border = box_color,
          lwd = 4
        )
        
        if (label_boxes) {
          conf_txt <- if ("Conf" %in% names(dets) && !is.na(dets$Conf[j])) {
            sprintf(" %.3f", dets$Conf[j])
          } else {
            ""
          }
          
          text(
            x = dets$TLx[j],
            y = max(20, dets$TLy[j] - 10),
            labels = paste0("scallop", conf_txt),
            col = box_color,
            cex = 1.0,
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
  
  output_files <- vector("character", nrow(seq_df))
  annotated_imgs <- vector("list", nrow(seq_df))
  
  for (i in seq_len(nrow(seq_df))) {
    img_name <- seq_df$imagename[i]
    img_path <- file.path(img_dir, img_name)
    
    if (!file.exists(img_path)) {
      warning(sprintf("Missing image: %s", img_path))
      next
    }
    
    img_full <- image_read(img_path)
    img_split <- split_image_exact(img_full, split)
    
    # IMPORTANT:
    # detections were generated on the split image itself,
    # so DO NOT remap x coordinates and DO NOT rescale for imgsz=1024
    dets_i <- det_sub %>%
      filter(Imagename == img_name)
    
    split_info <- image_info(img_split)
    split_w <- split_info$width[1]
    split_h <- split_info$height[1]
    
    # Keep only boxes that live inside the split image canvas
    dets_i <- dets_i %>%
      filter(
        TLx >= 0, TLy >= 0,
        BRx <= split_w, BRy <= split_h,
        BRx > TLx, BRy > TLy
      )
    
    label_text <- paste0(
      "Frame ", i, "/", nrow(seq_df),
      " | ", img_name,
      " | split=", split,
      if (seq_df$processed[i]) " | PROCESSED" else ""
    )
    
    img_annot <- draw_boxes_magick(
      img = img_split,
      dets = if (nrow(dets_i) > 0) dets_i else NULL,
      processed = seq_df$processed[i],
      label_text = label_text
    )
    
    out_file <- file.path(out_dir, sprintf("%02d_%s", i, img_name))
    image_write(img_annot, path = out_file, format = tools::file_ext(out_file))
    
    output_files[i] <- out_file
    annotated_imgs[[i]] <- img_annot
  }
  
  annotated_imgs <- annotated_imgs[!vapply(annotated_imgs, is.null, logical(1))]
  
  if (length(annotated_imgs) > 0) {
    montage <- image_montage(
      do.call(c, annotated_imgs),
      tile = "6x",
      geometry = "400x300+4+4"
    )
    
    montage_path <- file.path(out_dir, montage_file)
    image_write(montage, path = montage_path, format = tools::file_ext(montage_path))
  } else {
    montage_path <- NA_character_
  }
  
  invisible(list(
    sequence_df = seq_df,
    detections_used = det_sub,
    image_files = output_files,
    montage_file = montage_path,
    split = split
  ))
}


res <- visualize_overlap_sequence(
  view_overlap = view_overlap,
  pred_detections = pred_detections %>% filter(Conf > 0.3),
  img_dir = "/Volumes/PortableSSD/saltnoaa/images/overlap_2023/",
  start_index = 79,
  n_frames = 21,
  split = "right",
  out_dir = "~/Downloads/overlap_check_seq1",
  montage_file = "seq1_montage.png"
)
