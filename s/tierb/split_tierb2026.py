#!/usr/bin/env python3
"""
Tier B (independent 2026 survey test): split the annotated 2026 images into
Calibration-fitting and Holdout test.

The detector (YOLOv12n, job 598580) was trained on 2022-2024 only, so every
2026 image is new to it and there is no Detector-training level. This script
reproduces the second (GAM) level of build_strat_labels_v2.py on the 2026 set:

  * transect_group = cruise token + date token from the image name
    (e.g. 202605_20260523, 202614_Polaris_20260704)
  * images ordered along the transect (time token, frame number)
  * consecutive chunks of CHUNK (default 20) images; within each chunk 60% go to
    Calibration-fitting and 40% to Holdout test, with the holdout picks spread
    across the chunk's density order (systematic sample with a random start;
    see systematic_pick for why this differs from build_strat_labels_v2.py)
  * --mode whole assigns whole chunks instead (spatially stricter sensitivity)

Dataset labels keep the names the R pipeline expects:
  test_GAM_train = Calibration-fitting, test_GAM_test = Holdout test.
transect_block_id = the chunk id (used for CV folds and the block bootstrap).

Inputs : tierb_manifest.csv from build_tierb_eval.py (excluded rows are dropped)
Outputs: dataset_split_tierb2026.csv + a printed summary (and _summary.txt)

  python s/tierb/split_tierb2026.py \
      --manifest data/raw_model/2224scallop_yolo12n_slurm_598580/tierb2026/tierb_manifest.csv \
      --out data/raw/dataset_split_tierb2026.csv
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

import numpy as np
import pandas as pd

# Same density-bin edges as the 2022-2026 split (pandas qcut produced
# (-0.001, 1.0], (1.0, 2.0], (2.0, max]); fixed here because 2026 is ~90% empty
# and qcut would collapse. assign_folds() balances the top bin across folds.
DENSITY_EDGES = [-0.001, 1.0, 2.0]


def parse_tokens(name: str) -> dict:
    parts = Path(str(name)).stem.split(".")
    out = {"name_cruise_token": None, "name_date_token": None,
           "name_time_token": None, "name_frame_number": np.nan}
    if len(parts) >= 4:
        out.update(name_cruise_token=parts[0], name_date_token=parts[1], name_time_token=parts[2])
        try:
            out["name_frame_number"] = int(parts[3])
        except ValueError:
            pass
    return out


def systematic_pick(idx_sorted: list, n_pick: int, rng: random.Random) -> list:
    """Systematic sample of n_pick items from a density-sorted list, random start.

    Every density rank gets the same holdout probability n_pick/n. (The
    linspace rule in build_strat_labels_v2.py always includes the densest
    image of each chunk; with ~90% empty 2026 chunks that tilted scallops
    toward the holdout: 49% of primary scallops in a 40% holdout.)
    """
    n = len(idx_sorted)
    if n_pick <= 0 or n == 0:
        return []
    h, u = n_pick / n, rng.random()
    return [i for j, i in enumerate(idx_sorted)
            if int(np.floor((j + 1) * h + u)) - int(np.floor(j * h + u)) == 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk", type=int, default=20, help="images per chunk (GAM_BLOCK_SIZE)")
    ap.add_argument("--fit_frac", type=float, default=0.6, help="Calibration-fitting share")
    ap.add_argument("--mode", choices=["interleave", "whole"], default="interleave")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    man = pd.read_csv(a.manifest, keep_default_na=False, na_values=[""])
    n0 = len(man)
    excl = man["excluded"].fillna("").astype(str).str.len() > 0 if "excluded" in man else pd.Series(False, index=man.index)
    df = man[~excl].copy()
    df["n_annotations"] = pd.to_numeric(df["n_annotations"], errors="coerce").fillna(0).astype(int)
    df["is_empty"] = df["n_annotations"] == 0
    df = df.drop_duplicates("imagename").reset_index(drop=True)

    tok = pd.DataFrame([parse_tokens(n) for n in df["imagename"]])
    df = pd.concat([df, tok], axis=1)
    bad = df["name_cruise_token"].isna()
    if bad.any():
        print(f"WARNING: {bad.sum()} image names could not be parsed; they are dropped:")
        print(df.loc[bad, "imagename"].head().to_string())
        df = df[~bad].copy()
    if "platform" not in df:
        df["platform"] = np.where(df["imagename"].str.lower().str.contains("polaris"), "polaris",
                          np.where(df["imagename"].str.lower().str.contains("stella"), "stella", "primary"))

    df["image_timestamp_parsed"] = pd.to_datetime(
        df["name_date_token"] + df["name_time_token"].str[:6], format="%Y%m%d%H%M%S", errors="coerce")
    df["transect_group"] = df["name_cruise_token"] + "_" + df["name_date_token"]
    df = df.sort_values(["transect_group", "name_time_token", "name_frame_number", "imagename"]).reset_index(drop=True)
    df["transect_position"] = df.groupby("transect_group").cumcount()
    df["transect_block"] = df["transect_position"] // a.chunk
    df["transect_block_id"] = df["transect_group"] + "_chunk_" + df["transect_block"].astype(str)

    top = max(int(df["n_annotations"].max()), 3)
    df["density_bin"] = pd.cut(df["n_annotations"], DENSITY_EDGES + [float(top)]).astype(str)

    df["dataset"] = None
    if a.mode == "interleave":
        for _, sub in df.groupby("transect_block_id", sort=False):
            idx = list(sub.index)
            n = len(idx)
            n_fit = min(max(int(round(n * a.fit_frac)), 0), n)
            idx_sorted = sorted(idx, key=lambda i: (df.at[i, "n_annotations"], rng.random()))
            hold = systematic_pick(idx_sorted, n - n_fit, rng)
            df.loc[[i for i in idx if i not in hold], "dataset"] = "test_GAM_train"
            df.loc[hold, "dataset"] = "test_GAM_test"
    else:
        # whole chunks: within each transect_group, assign chunks 60/40 with the
        # holdout chunks spread across the group's chunk-density order
        for _, sub in df.groupby("transect_group", sort=False):
            ch = sub.groupby("transect_block_id")["n_annotations"].sum()
            ids = sorted(ch.index, key=lambda c: (ch[c], rng.random()))
            n_hold = len(ids) - min(max(int(round(len(ids) * a.fit_frac)), 0), len(ids))
            if len(ids) == 1:   # a single chunk cannot be split: randomise it
                n_hold = int(rng.random() >= a.fit_frac)
            hold = set(systematic_pick(ids, n_hold, rng))
            df.loc[sub.index, "dataset"] = np.where(sub["transect_block_id"].isin(hold),
                                                    "test_GAM_test", "test_GAM_train")

    df["partition"] = df["dataset"].map({"test_GAM_train": "Calibration-fitting",
                                         "test_GAM_test": "Holdout test"})
    df["is_test"] = True
    df["is_train"] = False
    df["dataset_level1"] = "test"
    df["year"] = 2026
    df["split_mode"] = a.mode

    # Along-transect spacing between consecutive annotated frames (seconds),
    # for the spatial-independence argument (distances are computed in R once
    # latitude/longitude are joined).
    dt = df.groupby("transect_group")["image_timestamp_parsed"].diff().dt.total_seconds()
    df["sec_since_prev_annotated"] = dt

    lines = []
    p = lambda s="": lines.append(str(s))
    p(f"Manifest rows: {n0:,}; excluded: {int(excl.sum()):,}; used: {len(df):,}  (mode = {a.mode}, chunk = {a.chunk})")
    p("\nImages / scallop images / scallops by platform x partition:")
    tab = df.groupby(["platform", "partition"]).agg(images=("imagename", "size"),
                                                   scallop_images=("is_empty", lambda s: int((~s).sum())),
                                                   scallops=("n_annotations", "sum"),
                                                   max_per_image=("n_annotations", "max"))
    tab["share_images"] = (tab["images"] / tab.groupby(level=0)["images"].transform("sum")).round(3)
    tab["share_scallops"] = (tab["scallops"] / tab.groupby(level=0)["scallops"].transform("sum")).round(3)
    p(tab.to_string())
    p("\nEmpty-image share by platform x partition:")
    p(df.groupby(["platform", "partition"])["is_empty"].mean().round(3).unstack().to_string())
    p("\nDensity bin x partition (all platforms):")
    p(pd.crosstab(df["density_bin"], df["partition"]).to_string())
    p("\nChunks: " + str(df["transect_block_id"].nunique()) + "; transect groups: " + str(df["transect_group"].nunique()))
    p("Seconds between consecutive annotated frames (median, by platform):")
    p(df.groupby("platform")["sec_since_prev_annotated"].median().round(1).to_string())

    cols = ["imagename", "platform", "n_annotations", "is_empty", "dataset", "partition", "dataset_level1",
            "is_train", "is_test", "year", "transect_group", "transect_position", "transect_block",
            "transect_block_id", "density_bin", "name_cruise_token", "name_date_token", "name_time_token",
            "name_frame_number", "image_timestamp_parsed", "sec_since_prev_annotated", "split_mode",
            "img_w", "img_h", "n_raw_pos"]
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df[[c for c in cols if c in df.columns]].to_csv(out, index=False)
    p(f"\nWrote {out}")
    txt = "\n".join(lines)
    print(txt)
    out.with_name(out.stem + "_summary.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
