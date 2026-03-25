python - <<EOF
from ultralytics import YOLO
import torch

# Paths
model_path = "../models/2224scallop_yolo12n_261070/weights/best.pt"
# model_path = "../models/2224scallop_viame_260065/best_snapshot.pt"
img_path   = "../data/test2.jpg"

# Check GPU
print("CUDA available:", torch.cuda.is_available())
print("Device:", "cuda" if torch.cuda.is_available() else "cpu")

# Load model
model = YOLO(model_path)
print("Model loaded")

# Run inference
results = model.predict(
    source=img_path,
    device=0,          # use GPU 0
    conf=0.25,
    save=True          # saves output image with boxes
)

print("Inference complete")
EOF