#!/usr/bin/env python3
"""
Stage the pipeline inputs for the Zenodo data record (run from s/).

Writes ../data/zenodo/pipeline_inputs.zip. Unzipped in the repository root, it recreates
data/raw/ and data/raw_model/ exactly as the scripts read them, so a reader can run
run_pipeline.sh and the figure/table scripts without renaming anything. The cleaned
tables (image_metadata.csv, annotations.csv, detections_<model>.csv) remain the
documented, human-readable version of the same data.

Personal identifiers (annotator user IDs, assignment IDs, annotation timestamps) and
local file paths in run summaries are dropped. Also writes checksums_pipeline_inputs.md5.
"""
import hashlib, io, zipfile
from pathlib import Path
import pandas as pd

ROOT = Path("..")
OUT = ROOT / "data/zenodo"; OUT.mkdir(parents=True, exist_ok=True)
DROP_ID = ["assignment_id", "annotation_timestamp", "annotator_user_id"]
DROP_PATH = ["run_dir", "pred_dir", "gt_json", "train_json", "weights", "data_root"]

RAW = ["data/raw/dataset_split_2226.csv", "data/raw/annotation_exclusions_2226.csv",
       "data/raw/groundtruth2226.csv", "data/raw/dataset_split_tierb2026.csv", "data/raw/allimg_2026.RData"]
MODEL_DIRS = ["data/raw_model/2226data_yolo12n_slurm_598585", "data/raw_model/2226data_cascade_slurm_598695",
              "data/raw_model/2226data_rtdetr_l_gcp_vit_20260924_190833",
              "data/raw_model/2224scallop_yolo12n_slurm_598580/tierb2026"]

def payload(rel):
    p = ROOT / rel
    if p.suffix == ".csv":
        d = pd.read_csv(p, low_memory=False)
        drop = [c for c in DROP_ID + (DROP_PATH if p.name.startswith("run_summary") else []) if c in d.columns]
        if drop:
            buf = io.StringIO(); d.drop(columns=drop).to_csv(buf, index=False); return buf.getvalue().encode()
    return p.read_bytes()

files = list(RAW)
for d in MODEL_DIRS:
    files += [str(f.relative_to(ROOT)) for f in sorted((ROOT / d).glob("*")) if f.is_file() and f.suffix != ".json"]
missing = [f for f in files if not (ROOT / f).exists()]
if missing: raise SystemExit(f"Missing inputs: {missing}")

zpath = OUT / "pipeline_inputs.zip"
sums = []
with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
    for rel in files:
        b = payload(rel); z.writestr(rel, b)
        sums.append(f"{hashlib.md5(b).hexdigest()}  {rel}")
        print(f"{len(b)/1e6:8.1f} MB  {rel}")
(OUT / "checksums_pipeline_inputs.md5").write_text("\n".join(sums) + "\n")
print(f"Wrote {zpath} ({zpath.stat().st_size/1e6:.0f} MB, {len(files)} files)")
