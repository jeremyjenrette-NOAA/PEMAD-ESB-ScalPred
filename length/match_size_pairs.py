#!/usr/bin/env python3
"""
Match each true-positive detection to the annotation it was matched to, for shell-height validation.

Why: the evaluator's `index` column is the detection's rank within the image, not the matched
annotation's id (eval_yolo_detections_hung.py writes "index": i + 1 and leaves "man" empty).
Joining detections to annotations on it pairs them arbitrarily in multi-scallop images.
This script repeats the evaluator's Hungarian assignment (IoU >= 0.1, maximise total IoU) between
each Holdout-test image's true-positive boxes and its annotation boxes, and writes the pairs.

Run from length/ (needs numpy, pandas, scipy):
  python match_size_pairs.py
Output: ../data/processed/size_pairs_2226.csv (one row per matched pair, all detectors)
Check: the number of pairs per detector must equal its true-positive count (printed).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

RAW = Path("../data/raw"); MODELS = Path("../data/raw_model")
DETS = {
    "YOLOv12":       MODELS / "2226data_yolo12n_slurm_598585/autotest2226_yolo12n.csv",
    "Cascade R-CNN": MODELS / "2226data_cascade_slurm_598695/autotest2226_viame_cascade.csv",
    "RT-DETR-L":     MODELS / "2226data_rtdetr_l_gcp_vit_20260924_190833/autotest2226_rtdetr_l.csv",
}
MATCH_IOU = 0.1
OUT = Path("../data/processed/size_pairs_2226.csv")

split = pd.read_csv(RAW / "dataset_split_2226.csv", usecols=["imagename", "dataset"])
hold = set(split.imagename[split.dataset == "test_GAM_test"])

gt = pd.read_csv(RAW / "groundtruth2226.csv",
                 usecols=["imagename", "tlx", "tly", "brx", "bry", "label", "geometry_text",
                          "millimeter_per_pixel", "longitude", "class_name"])
gt = gt[gt.imagename.isin(hold)].copy()
gt = gt[gt.label.astype(str).str.strip().str.lower() == "scallop"].dropna(subset=["tlx", "tly", "brx", "bry"])
gt["gt_row"] = gt.groupby("imagename").cumcount() + 1          # order within image in groundtruth2226.csv
nums = gt.geometry_text.astype(str).str.findall(r"[-+]?[0-9]*\.?[0-9]+")
gt["true_px"] = nums.apply(lambda p: float(np.hypot(float(p[2]) - float(p[0]), float(p[3]) - float(p[1])))
                           if len(p) >= 4 else np.nan)
gt["true_mm"] = gt.true_px * gt.millimeter_per_pixel
gt["region"] = np.where(gt.longitude >= -71, "GB", "MAB")
x1, x2 = np.minimum(gt.tlx, gt.brx), np.maximum(gt.tlx, gt.brx)
y1, y2 = np.minimum(gt.tly, gt.bry), np.maximum(gt.tly, gt.bry)
gt["gx1"], gt["gy1"], gt["gx2"], gt["gy2"] = x1, y1, x2, y2
by_img = {k: v for k, v in gt.groupby("imagename")}


def iou(a, b):
    ix = np.maximum(0, np.minimum(a[:, None, 2], b[None, :, 2]) - np.maximum(a[:, None, 0], b[None, :, 0]))
    iy = np.maximum(0, np.minimum(a[:, None, 3], b[None, :, 3]) - np.maximum(a[:, None, 1], b[None, :, 1]))
    inter = ix * iy
    aa = (a[:, 2] - a[:, 0]) * (a[:, 3] - a[:, 1]); bb = (b[:, 2] - b[:, 0]) * (b[:, 3] - b[:, 1])
    return inter / (aa[:, None] + bb[None, :] - inter + 1e-9)


rows = []
for model, path in DETS.items():
    if not path.exists():
        print(f"skip {model}: {path} not found"); continue
    d = pd.read_csv(path)
    d = d[d.Imagename.isin(hold)]
    tp = d[d.truedetect.astype(str).str.lower() == "true"]
    n_before = len(rows)
    for im, sub in tp.groupby("Imagename"):
        G = by_img.get(im)
        if G is None:
            continue
        a = sub[["TLx", "TLy", "BRx", "BRy"]].to_numpy(float)
        b = G[["gx1", "gy1", "gx2", "gy2"]].to_numpy(float)
        m = iou(a, b); cost = 1 - m; cost[m < MATCH_IOU] = 1e6
        r, c = linear_sum_assignment(cost)
        for i, j in zip(r, c):
            if m[i, j] < MATCH_IOU:
                continue
            det, g = sub.iloc[i], G.iloc[j]
            pred_px = (abs(det.BRx - det.TLx) + abs(det.BRy - det.TLy)) / 2      # 'average' back-calculation
            rows.append(dict(model=model, imagename=im, detectid=int(det.Detectid), det_index=int(det["index"]),
                             gt_row=int(g.gt_row), iou=float(m[i, j]), conf=float(det.Conf),
                             class_name=g.class_name, region=g.region, mm_per_px=g.millimeter_per_pixel,
                             true_px=g.true_px, true_mm=g.true_mm, pred_px=pred_px,
                             pred_mm=pred_px * g.millimeter_per_pixel))
    n = len(rows) - n_before
    print(f"{model}: {len(tp):,} true positives in Holdout test, {n:,} matched pairs"
          + ("  OK" if n == len(tp) else "  MISMATCH - check annotation boxes"))

out = pd.DataFrame(rows)
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False)
same = (out.det_index == out.gt_row).mean()
print(f"Wrote {OUT} ({len(out):,} rows). Pairs where detection index == matched annotation row: {same:.1%}")
