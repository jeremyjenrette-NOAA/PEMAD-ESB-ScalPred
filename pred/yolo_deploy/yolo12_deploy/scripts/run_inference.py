from ultralytics import YOLO
from pathlib import Path
import pandas as pd
from PIL import Image
import os
MODEL_PATH = "model/best.pt"
IMAGE_LIST = "image_list.txt"
OUTPUT_DIR = "output"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load model
model = YOLO(MODEL_PATH)

# Read image paths
with open(IMAGE_LIST, "r") as f:
    image_paths = [line.strip() for line in f.readlines()]

all_dets = []

for img_path in image_paths:

    results = model.predict(
        source=img_path,
        conf=0.01,
        imgsz=1024,
        save=True,
        project=OUTPUT_DIR,
        name="annotated",
        exist_ok=True
    )

    for r in results:

        boxes = r.boxes

        if boxes is None:
            continue

        for box in boxes:

            xyxy = box.xyxy.cpu().numpy()[0]
            conf = float(box.conf.cpu().numpy()[0])
            cls = int(box.cls.cpu().numpy()[0])

            all_dets.append({
                "image": Path(img_path).name,
                "xmin": xyxy[0],
                "ymin": xyxy[1],
                "xmax": xyxy[2],
                "ymax": xyxy[3],
                "confidence": conf,
                "class": cls
            })

# Save detections
pd.DataFrame(all_dets).to_csv(
    f"{OUTPUT_DIR}/detections.csv",
    index=False
)

print("Inference complete.")