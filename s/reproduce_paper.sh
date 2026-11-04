#!/bin/bash
# ======================================================================
# Reproduce every calibration, abundance and size result in the paper.
#
# Before running: unzip pipeline_inputs.zip from the Zenodo data record
# (https://doi.org/10.5281/zenodo.23481212) in the repository root, so that
# data/raw/ and data/raw_model/ exist. Then, from s/:
#   ./reproduce_paper.sh            # everything (about 8-10 h on a laptop)
#   ./reproduce_paper.sh n3_ablation  # start at a step (earlier outputs must exist)
# Each step logs to s/logs/<timestamp>_<step>.log; the script stops at the first
# failure and names the step, so you can restart from it.
# Table/figure map: see "Reproducing the manuscript" in README.md.
# ======================================================================
set -o pipefail
cd "$(dirname "$0")" || exit 1
command -v caffeinate >/dev/null && caffeinate -w $$ &

STEPS=(
  "makedat|Rscript makedat.R"                       # eval_2226.rds
  "gamstrat|Rscript gamstrat.R"                     # detection-level GAMs, Table 1 CV (~2 h)
  "gammod|Rscript gammod.R"                         # image-level inputs
  "performance|Rscript performance.R"               # Fig. 3, Table S1
  "table1|python3 tables/make_table2.py"            # Table 1 (detection-level selection)
  "n3_ablation|Rscript n3/n3_ablation.R"            # ablation, image-level CV, concurvity (~1-1.5 h)
  "n3_tables|python3 n3/make_n3_tables.py"          # Table 2, Table S5
  "n3_figure|python3 n3/plot_n3_ablation.py"        # Fig. 4
  "n3_export|Rscript n3/n3_export_holdout.R"        # holdout_eval_2226_v2.rds
  "error|Rscript error.R"                           # Fig. 6, Fig. S9
  "tableS2|python3 tables/make_tableS1.py"          # Table S2
  "fig_holdout|Rscript fig_holdout.R"               # Fig. 5, Fig. S2
  "depth|Rscript n11_depth.R"                       # Fig. S8
  "density|Rscript n7_density.R"                    # Table S7
  "ml|Rscript n3/n6_ml_compare.R"                   # Table S6, Fig. S10 (~30 min)
  "final_gams|Rscript n3/plot_final_gams.R"         # Fig. S11
  "covariates|Rscript tables/make_tableS7_covariates.R"  # Table S8
  "tierb|Rscript tierb/tierb_calibration.R"         # Table S3, Fig. S3 (2026 split shipped in the record)
  "tierb_effort|Rscript tierb/tierb_effort.R"       # Fig. S4a (~1-1.5 h)
  "tierb_effort_cv|Rscript tierb/tierb_effort_cv.R" # Fig. S4b (~2 h)
  "tiera|Rscript tiera/tiera_loyo.R"                # Fig. S5 (~2 h)
  "revision_figs|python3 revision_figs.py"          # Figs. S4, S5
  "size_pairs|cd ../length && python3 match_size_pairs.py"   # matched shell heights
  "size|cd ../length && Rscript size_revision.R"    # Fig. 7, Figs. S6-S7, Table S4
  "design|Rscript fig2_design.R"                    # Fig. 2 (downloads NOAA bathymetry once)
  "strat|Rscript strat_vis.R"                       # Fig. S1
)
START="${1:-makedat}"
mkdir -p logs; STAMP=$(date +%Y%m%d_%H%M%S); started=0; T0=$SECONDS
for s in "${STEPS[@]}"; do
  name="${s%%|*}"; cmd="${s#*|}"
  [[ "$name" == "$START" ]] && started=1
  (( started )) || continue
  echo "=== $(date +%H:%M:%S) $name: $cmd"
  if ! ( eval "$cmd" ) 2>&1 | tee "logs/${STAMP}_${name}.log"; then
    echo "FAILED at step '$name' (log: logs/${STAMP}_${name}.log). Restart with: ./reproduce_paper.sh $name"; exit 1
  fi
done
(( started )) || { echo "Unknown step '$START'"; exit 1; }
Rscript -e 'sessionInfo()' > "logs/${STAMP}_sessionInfo.txt" 2>&1
echo "Done in $(( (SECONDS - T0) / 60 )) min. Figures: figures/ms_figures/; tables: s/tables/."
