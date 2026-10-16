#!/usr/bin/env python3
"""
Build the tabular files for the Zenodo data record (run from s/).

Outputs (data/zenodo/, gitignored):
  image_metadata.csv     one row per annotated image: partition, cruise, region, covariates (2026 FOV corrected)
  annotations.csv        one row per manual scallop annotation (line segment + derived square box, size in mm)
  detections_<model>.csv raw detector outputs on the Calibration set (conf >= 0.1), matched to annotations
  missed_<model>.csv     annotations not matched by any detection (false negatives), per model
  run_summary_<model>.csv evaluation settings and counts
  checksums_tables.md5

Personal identifiers (annotator user IDs, assignment IDs, annotation timestamps) are dropped.
"""
import hashlib, re
from pathlib import Path
import numpy as np
import pandas as pd

RAW, RM = Path("../data/raw"), Path("../data/raw_model")
OUT = Path("../data/zenodo"); OUT.mkdir(parents=True, exist_ok=True)
MAX_SH_MM = 180
PARTITION = {"train": "Detector-training", "test_GAM_train": "Calibration-fitting", "test_GAM_test": "Holdout test"}
MODELS = {
    "yolov12n": RM / "2226data_yolo12n_slurm_598585",
    "cascade_rcnn": RM / "2226data_cascade_slurm_598695",
    "rtdetr_l": RM / "2226data_rtdetr_l_gcp_vit_20260924_190833",
}

split = pd.read_csv(RAW / "dataset_split_2226.csv", low_memory=False)
excl = pd.read_csv(RAW / "annotation_exclusions_2226.csv")
excl_ids = set(excl.loc[excl["status"].astype(str).str.strip().eq("confirmed"), "imagename"])
split["partition"] = split["dataset"].map(PARTITION)
split["excluded_after_review"] = split["imagename"].isin(excl_ids)
split["region"] = np.where(split["longitude"] >= -71, "GB", "MAB")
split["image_in_record"] = split["partition"].isin(["Calibration-fitting", "Holdout test"]) & ~split["excluded_after_review"]

meta_cols = {
    "imagename": "imagename", "year": "year", "cruise_id": "cruise_id", "region": "region",
    "partition": "partition", "transect_block_id": "transect_block_id", "image_timestamp_parsed": "image_timestamp_utc",
    "latitude": "latitude", "longitude": "longitude", "bottom_depth": "bottom_depth_m", "altitude": "altitude_m",
    "vehicle_roll_angle": "roll_deg", "vehicle_pitch_angle": "pitch_deg",
    "field_of_view_sq_meter": "field_of_view_m2", "millimeter_per_pixel": "mm_per_pixel",
    "ctd_temperature_celsius": "temperature_c", "fluorometer_backscatter_ntu": "backscatter_ntu",
    "ctd_salinity": "salinity_psu", "n_annotations": "n_scallops", "is_empty": "is_empty",
    "excluded_after_review": "excluded_after_review", "image_in_record": "image_in_record",
}
meta = split[list(meta_cols)].rename(columns=meta_cols)
meta.to_csv(OUT / "image_metadata.csv", index=False)

gt = pd.read_csv(RAW / "groundtruth2226.csv", low_memory=False)
def seg(t):
    n = [float(v) for v in re.findall(r"-?\d+(?:\.\d+)?", str(t))]
    return n[:4] if len(n) >= 4 else [np.nan] * 4
xy = np.array([seg(t) for t in gt["geometry_text"]])
ann = pd.DataFrame({
    "annotation_id": np.arange(1, len(gt) + 1), "imagename": gt["imagename"],
    "class_id": gt["class_id"], "class_name": gt["class_name"],
    "x1_px": xy[:, 0], "y1_px": xy[:, 1], "x2_px": xy[:, 2], "y2_px": xy[:, 3],
    "box_tlx_px": gt["tlx"], "box_tly_px": gt["tly"], "box_brx_px": gt["brx"], "box_bry_px": gt["bry"],
})
ann["length_px"] = np.hypot(ann.x2_px - ann.x1_px, ann.y2_px - ann.y1_px)
ann = ann.merge(meta[["imagename", "mm_per_pixel", "partition", "year", "region"]], on="imagename", how="left")
ann["length_mm"] = ann["length_px"] * ann["mm_per_pixel"]
ann["size_validation_excluded"] = ann["length_mm"] > MAX_SH_MM
ann.drop(columns="mm_per_pixel").to_csv(OUT / "annotations.csv", index=False)

keep_det = ["Imagename", "TLx", "TLy", "BRx", "BRy", "Conf", "boxsize", "truedetect", "iu"]
for tag, d in MODELS.items():
    det = pd.read_csv(next(d.glob("autotest*.csv")), low_memory=False)[keep_det]
    det.columns = ["imagename", "tlx_px", "tly_px", "brx_px", "bry_px", "confidence", "box_area_px2", "true_positive", "iou_with_match"]
    det.to_csv(OUT / f"detections_{tag}.csv", index=False)
    fn = pd.read_csv(next(d.glob("fn*.csv")), low_memory=False)[["Imagename", "TLx", "TLy", "BRx", "BRy", "boxsize"]]
    fn.columns = ["imagename", "tlx_px", "tly_px", "brx_px", "bry_px", "box_area_px2"]
    fn.to_csv(OUT / f"missed_{tag}.csv", index=False)
    rs = pd.read_csv(d / "run_summary.csv")
    rs.drop(columns=[c for c in ("run_dir", "pred_dir", "gt_json", "train_json", "weights", "data_root") if c in rs]).to_csv(OUT / f"run_summary_{tag}.csv", index=False)

with open(OUT / "checksums_tables.md5", "w") as fh:
    for f in sorted(OUT.glob("*.csv")):
        fh.write(f"{hashlib.md5(f.read_bytes()).hexdigest()}  {f.name}\n")
print(f"image_metadata.csv: {len(meta):,} images ({int(meta.image_in_record.sum()):,} with image files in the record)")
print(f"annotations.csv: {len(ann):,} annotations; {int(ann.size_validation_excluded.sum())} > {MAX_SH_MM} mm flagged")
print(f"{len(list(OUT.glob('*')))} files, {sum(f.stat().st_size for f in OUT.glob('*'))/1e6:.1f} MB in {OUT}")
