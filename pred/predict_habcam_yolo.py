#!/usr/bin/env python3

import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from PIL import Image
import argparse
import sys
import time
from pathlib import Path
from typing import List, Dict, Any

import torch
from ultralytics import YOLO
import pandas as pd

def chunk_list(items, batch_size):
    for i in range(0, len(items), batch_size):
        yield items[i:i + batch_size]

def load_and_split_image(image_path: str, split: str = "none"):
    img = Image.open(image_path).convert("RGB")
    w, h = img.size

    if split == "none":
        return img, 0
    elif split == "left":
        return img.crop((0, 0, w // 2, h)), 0
    elif split == "right":
        return img.crop((w // 2, 0, w, h)), w // 2
    else:
        raise ValueError(f"Invalid split option: {split}")

def normalize_windows_path(path: str) -> str:
    """
    Convert Git Bash/MSYS paths into Windows-native paths for Python.
    """
    if not path:
        return path

    p = str(path).strip()

    # MSYS/Git-Bash mounted drive path: /z/... -> Z:/...
    if p.startswith("/z/"):
        return "Z:/" + p[3:]

    # General /x/... -> X:/...
    if len(p) >= 4 and p[0] == "/" and p[2] == "/":
        drive = p[1].upper()
        if drive.isalpha():
            return f"{drive}:/{p[3:]}"

    return p

def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"{ts} | {msg}", flush=True)


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def load_inventory(inventory_path: Path, year: str | None = None, process_col: str = "process_image") -> pd.DataFrame:
    """
    Load compiled inventory and retain only images flagged for processing.
    Supports .tsv, .csv, or .rds only if already converted upstream.
    """
    suffix = inventory_path.suffix.lower()

    if suffix == ".tsv":
        df = pd.read_csv(inventory_path, sep="\t")
    elif suffix == ".csv":
        df = pd.read_csv(inventory_path)
    else:
        raise ValueError(
            f"Unsupported inventory format: {inventory_path}. "
            "Use TSV or CSV exported from R."
        )

    required = {"actual_path", "imagename", process_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Inventory missing required columns: {sorted(missing)}")

    df = df.copy()
    df[process_col] = df[process_col].astype(bool)

    df["actual_path_raw"] = df["actual_path"].astype(str)
    df["actual_path"] = df["actual_path_raw"].apply(normalize_windows_path)

    if "year" in df.columns:
        df["year"] = df["year"].astype(str)
    else:
        # fallback from path if year column absent
        df["year"] = df["actual_path"].astype(str).str.extract(r"/(20\d{2})/")[0]

    df = df[df[process_col]].copy()

    if year is not None:
        df = df[df["year"] == str(year)].copy()

    df = df.drop_duplicates(subset=["actual_path"]).reset_index(drop=True)
    return df


def load_completed(completed_file: Path) -> set[str]:
    if not completed_file.exists():
        return set()
    with completed_file.open("r", encoding="utf-8") as f:
        return {line.strip() for line in f if line.strip()}


def append_completed(completed_file: Path, image_path: str) -> None:
    ensure_parent(completed_file)
    with completed_file.open("a", encoding="utf-8") as f:
        f.write(image_path + "\n")


def append_detection_rows(csv_file: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return

    ensure_parent(csv_file)
    df = pd.DataFrame(rows)

    write_header = not csv_file.exists()
    df.to_csv(csv_file, mode="a", header=write_header, index=False)


def build_detection_rows(
    result,
    image_path: str,
    imagename: str,
    x_offset: int = 0,
    split_mode: str = "none",
    model_name_des: str = "none",
    default_label: str = "scallop"
) -> List[Dict[str, Any]]:
    """
    Convert one Ultralytics result into row-wise detections.
    Coordinates are xyxy pixel coordinates.
    """
    rows: List[Dict[str, Any]] = []

    boxes = result.boxes
    if boxes is None or len(boxes) == 0:
        return rows

    names = result.names if hasattr(result, "names") else {}

    xyxy = boxes.xyxy.detach().cpu().numpy()
    confs = boxes.conf.detach().cpu().numpy() if boxes.conf is not None else [None] * len(xyxy)
    clss = boxes.cls.detach().cpu().numpy().astype(int) if boxes.cls is not None else [None] * len(xyxy)

    for i, (coords, conf, cls_id) in enumerate(zip(xyxy, confs, clss), start=1):
        tlx, tly, brx, bry = coords.tolist()
        label = names.get(int(cls_id), default_label) if cls_id is not None else default_label

        rows.append({
            "Detectid": f"{imagename}__det{i}",
            "Imagename": imagename,
            "TLx": float(tlx),
            "TLy": float(tly),
            "BRx": float(brx),
            "BRy": float(bry),
            "Conf": float(conf) if conf is not None else None,
            "Spname": label,
            "img_path": image_path,
            "split": split_mode,
            "pred_datetime": time.strftime("%Y-%m-%d %H:%M:%S"),
            "model": model_name_des
        })

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Resumable YOLO inference over HabCam inventory.")
    parser.add_argument("--split", default="none", choices=["none", "left", "right"], help="Use whole image or split stereo pair before inference")
    parser.add_argument("--model", required=True, help="Path to YOLO .pt weights")
    parser.add_argument("--inventory", required=True, help="Compiled inventory TSV/CSV")
    parser.add_argument("--model_name", required=True, help="YOLOv12")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--year", default=None, help="Year to process, e.g. 2022")
    parser.add_argument("--device", default="0", help="CUDA device, e.g. 0 or 1")
    parser.add_argument("--conf", type=float, default=0.25, help="Confidence threshold")
    parser.add_argument("--nms_iou", type=float, default=0.50, help="NMS IoU threshold")
    parser.add_argument("--imgsize", type=int, default=1024, help="Inference image size")
    parser.add_argument("--max_detections", type=int, default=300, help="Maximum detections per image")
    parser.add_argument("--save", action="store_true", help="Save annotated prediction images")
    parser.add_argument("--save_every", type=int, default=100, help="Log progress every N images")
    parser.add_argument("--batch_size", type=int, default=16, help="Number of images per inference batch")
    parser.add_argument(
    "--process_col",
    default="process_image",
    help="Column name indicating which images to process"
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    detections_csv = outdir / f"detections_{args.year if args.year else 'all'}.csv"
    completed_txt = outdir / f"completed_{args.year if args.year else 'all'}.txt"
    errors_txt = outdir / f"errors_{args.year if args.year else 'all'}.txt"
    run_log = outdir / f"run_{args.year if args.year else 'all'}.log"

    def file_log(msg: str) -> None:
        log(msg)
        with run_log.open("a", encoding="utf-8") as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} | {msg}\n")

    file_log("Starting inference job")
    file_log(f"CUDA available: {torch.cuda.is_available()}")
    file_log(f"Requested device: {args.device}")

    inventory = load_inventory(
    Path(args.inventory),
    year=args.year,
    process_col=args.process_col
    )
    file_log(f"Inventory rows eligible for processing: {len(inventory)}")

    completed = load_completed(completed_txt)
    if completed:
        file_log(f"Previously completed images found: {len(completed)}")

    inventory = inventory[~inventory["actual_path"].isin(completed)].copy().reset_index(drop=True)
    file_log(f"Remaining images to process this run: {len(inventory)}")

    model = YOLO(args.model)
    file_log(f"Model loaded: {args.model}")
    
    file_log(f"Split mode: {args.split}")

    processed_this_run = 0
    detections_this_run = 0

    records = inventory.to_dict("records")

    for batch_idx, batch_rows in enumerate(chunk_list(records, args.batch_size), start=1):
        batch_images = []
        batch_meta = []

        for row in batch_rows:
            image_path = str(row["actual_path"])
            image_path_raw = str(row["actual_path_raw"]) if "actual_path_raw" in row else image_path
            imagename = str(row["imagename"])

            if not os.path.exists(image_path):
                msg = f"MISSING\t{image_path_raw}\tWINDOWS_PATH\t{image_path}"
                with errors_txt.open("a", encoding="utf-8") as f:
                    f.write(msg + "\n")
                append_completed(completed_txt, image_path)
                continue

            try:
                img_for_pred, x_offset = load_and_split_image(
                    image_path,
                    split=args.split
                )

                batch_images.append(img_for_pred)
                batch_meta.append({
                    "image_path": image_path,
                    "imagename": imagename,
                    "x_offset": x_offset
                })

            except Exception as e:
                msg = f"ERROR\t{image_path}\t{repr(e)}"
                with errors_txt.open("a", encoding="utf-8") as f:
                    f.write(msg + "\n")
                continue

        if len(batch_images) == 0:
            continue

        try:
            results = model.predict(
                source=batch_images,
                device=args.device,
                conf=args.conf,
                iou=args.nms_iou,
                imgsz=args.imgsize,
                max_det=args.max_detections,
                save=args.save,
                verbose=False,
            )

            for res, meta in zip(results, batch_meta):
                det_rows = build_detection_rows(
                    res,
                    image_path=meta["image_path"],
                    imagename=meta["imagename"],
                    x_offset=meta["x_offset"],
                    split_mode=args.split,
                    model_name_des=args.model_name
                )

                append_detection_rows(detections_csv, det_rows)
                append_completed(completed_txt, meta["image_path"])

                processed_this_run += 1
                detections_this_run += len(det_rows)

            if processed_this_run % args.save_every == 0:
                file_log(
                    f"Processed {processed_this_run} images this run | "
                    f"batch {batch_idx} | total detections written: {detections_this_run}"
                )

        except Exception as e:
            msg = f"BATCH_ERROR\tbatch={batch_idx}\t{repr(e)}"
            with errors_txt.open("a", encoding="utf-8") as f:
                f.write(msg + "\n")
            continue

    file_log(f"Run complete | images processed: {processed_this_run} | detections written: {detections_this_run}")
    file_log(f"Detections CSV: {detections_csv}")
    file_log(f"Completed checkpoint: {completed_txt}")
    file_log(f"Errors file: {errors_txt}")


if __name__ == "__main__":
    main()
