\
# viame_project2224

This project folder was scaffolded from the provided Windows runtime files.

## Included
- `pipelines/detector_project_folder.pipe`
- `run_trained_model.bat`
- `setup_viame.bat`
- `process_video.py`
- `category_models/detector.pipe.template` (copied from embedded_default.pipe)
- `category_models/detector.pipe` (placeholder stub)

## Still needed before a real test run
1. Copy your trained model bundle to `category_models/trained_detector.zip`
2. Replace `category_models/detector.pipe` with a valid VIAME detector pipe
   that loads `trained_detector.zip`
3. Populate `input_list.txt` with sample image paths

## Why the placeholder exists
`process_video.py` explicitly checks for `category_models/detector.pipe`
when `pipelines/detector_project_folder.pipe` is used. The exact rendered
pipe produced by native `viame train` was not available in the supplied files,
so this scaffold preserves the expected project structure while keeping the
missing piece isolated to one file.

## Intended Windows location
Copy this folder to:
`C:\Users\jeremy.jenrette\PEMAD-ESB-ScalPred\pred\viame_project2224`
