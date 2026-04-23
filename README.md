# PEMAD-ESB-ScalPred

## Maintainer

Name: Jeremy J

Email: jeremy.jenrette@noaa.gov

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

## Description
PEMAD-ESB-ScalPred is a framework for automated detection, abundance estimation, and size inference of Atlantic sea scallops (***Placopecten magellanicus***) from large-scale benthic imagery collected by the HabCam system.

The project integrates the latest object detection models (YOLO and Cascade R-CNN with VIAME) with statistical calibration workflows to address well-known sources of bias in image-based ecological surveys, including false positives, false negatives, and environmental variability. The resulting outputs provide spatially and temporally scalable estimates of scallop abundance and size structure, supporting fisheries assessment and ecosystem monitoring.

This repository contains the code necessary to evaluate trained raw model outputs with image metadata. Here, we can train Generative Additive Models (GAMs), perform inference, predict calibrated abundance estimates, and predict length classes.  Included are figures and diagnostic plots. To see the **Detection Training** repository, visit...

### Project
PEMAD-ESB-ScalPred provides:

* Automating scallop detection across multi-year survey datasets (2022–2024)
* Integrating predictions with environmental and survey metadata
* Calibrating detection outputs using Generative Additive Models
* Converting raw detections into biologically meaningful quantities:

  * **Density (n m⁻²)**
  * **Length distributions (mm)**
  * **Spatial abundance patterns**

```
.
├── data
│   ├── images # diagnostic images
│   │   ├── 202203.20220515.135703301.59728.png
│   │   └── 202203.20220515.224949946.9971.png
│   │
│   ├── processed
│   │   ├── Cas2224.RData # VIAME trained + evaluated model paired with metadata
│   │   ├── YOLOv122224.RData # YOLOv12 trained + evaluated paired with metadata
│   │   ├── completed/ # these are all images processed by the detection model
│   │   │   └── completed_yolo_2024.txt
│   │   ├── detections/ # these are all images WITH DETECTIONS
│   │   │   └── detections_yolo_2024.csv
│   │   └── metapred/ # metadata folder
│   │       └── meta_all2224.csv
│   │
│   └── raw # evaluated models and raw weights, also contains other unprocessed data
│       ├── 2022scallop_viame_257794/
│       │   └── run_summary.csv
│       ├── 2022scallop_yolo11n_242222/
│       │   ├── results.csv
│       │   └── weights/
│       │       └── best.pt
│       └── groundtruth2224.csv
│
├── doc
│
├── figures
│   ├── 2224_caldense.jpg
│   ├── 2224_lengthdist.jpg
│   ├── 2224_totalcal.jpg
│   ├── 2224yolov12cas_gammod.png
│   ├── NEFSC-scal-flowchart.png
│   └── diag1/
│       └── 2224_yolo_compare_m7.jpg
│
├── length # for calculating length from boxes and other diagnostics
│   ├── findlengths.R
│   └── lengthfunc.R
│
├── pred # for dataset-wide inference using YOLO and Cascade R-CNN
│   ├── calpred.R
│   ├── predfunc.R
│   ├── predict_habcam_yolo.py
│   ├── imglist.sh
│   ├── view_overlap.R
│   └── viame_project2224/
│       ├── category_models/
│       │   └── detector.pipe
│       ├── pipelines/
│       │   └── detector_project_folder.pipe
│       ├── predict_habcam_viame.py
│       ├── run_trained_model.sh
│       └── README_PROJECT_FOLDER.md
│
├── s # data-processing and visualization scripts
│   ├── datfunc.R
│   ├── fitfunc.R
│   ├── gammod.R
│   ├── gamcal.R
│   ├── performance.R
│   └── procfunc.R
│
├── tex # LaTeX
│   ├── main.tex
│   ├── main.pdf
│   └── compile_main.sh
│
├── README.md
├── requirements.txt
└── LICENSE
```

## Workflow

The pipeline follows a structured, end-to-end workflow:

### 1. Data Ingestion and Preparation

* Raw HabCam imagery and metadata are compiled
* Ground-truth annotations are harmonized (`groundtruth2224.csv`)
* Image inventories and environmental covariates are standardized

---

### 2. Model Training and Evaluation

* Object detection models trained using:

  * **YOLO (Ultralytics)**
  * **Cascade R-CNN (VIAME framework)**
* Model outputs include:

  * Bounding boxes
  * Confidence scores
* Performance evaluated using held-out datasets (F1, precision-recall)

---

### 3. Dataset-wide Inference

* Trained models deployed across full survey datasets
* Outputs stored as:

  * `detections_*.csv` (detections only)
  * `completed_*.txt` (all processed images)

---

### 4. Prediction Calibration

* Detection-level predictions are calibrated using **Generalized Additive Models (GAMs)**
* Calibration accounts for:

  * Environmental conditions (e.g., depth, backscatter, turbidity)
  * Detection uncertainty (confidence scores)
* Produces **probabilistic detection estimates** (p(detection))

---

### 5. Abundance Estimation

* Image-level abundance computed as:

[
\hat{N}_{image} = \sum_i p_i
]

* Density derived using field-of-view metadata:

[
\hat{D} = \frac{\sum_i p_i}{\text{area}_{image}}
]

---

### 6. Length Estimation

* Bounding boxes converted to scallop length using pixel-to-mm scaling
* Produces:

  * Individual size estimates
  * Population-level length distributions

---

### 7. Visualization and Diagnostics

* Model comparisons (YOLO vs VIAME)
* Calibration diagnostics
* Spatial and depth-based trends
* Length-frequency distributions

Outputs are stored in `/figures` and used for manuscript development.

---

## Results

This framework demonstrates that:

* Deep learning models can reliably detect scallops across heterogeneous benthic environments
* **Calibration is essential** to correct systematic detection bias
* Probabilistic predictions outperform threshold-based counting approaches
* Automated workflows can produce:

  * **Accurate abundance estimates**
  * **Robust size distributions**
  * **Scalable survey analytics**

Across the 2022–2024 datasets, the pipeline enables:

* Processing of **millions of images**
* Generation of **survey-wide density estimates**
* Consistent comparison across models and years

---

## Key Contributions

* End-to-end pipeline linking **computer vision → statistical ecology**
* Integration of **multiple detection frameworks** in a unified workflow
* Explicit handling of **false positives and false negatives**
* Scalable methodology for **image-based fisheries assessment**

---

## Contact