#!/bin/bash
# ======================================================================
# Full calibration workflow: makedat -> gamstrat -> gammod -> error
#
# Usage (from s/):
#   ./run_pipeline.sh              # run everything
#   ./run_pipeline.sh gammod       # start at a step (skips earlier ones)
#   ./run_pipeline.sh --check      # run the random-effect self-check first
#
# Each step's output is saved to s/logs/<timestamp>_<step>.log and echoed to
# the terminal. The script stops at the first failing step and says which
# one, so you can fix it and restart from there.
# Runtime: makedat ~1 min, gamstrat ~2 h (formula CV is cached and refit only
# when eval_2226.rds changes), gammod ~5 min, error ~2 min.
# ======================================================================
set -o pipefail
cd "$(dirname "$0")" || exit 1

# keep the Mac awake until this script exits
caffeinate -w $$ &

STEPS=(makedat gamstrat gammod error)
START="${STEPS[0]}"
RUN_CHECK=0
for arg in "$@"; do
  case "$arg" in
    --check) RUN_CHECK=1 ;;
    *) START="$arg" ;;
  esac
done
if [[ ! " ${STEPS[*]} " =~ " ${START} " ]]; then
  echo "Unknown step '$START'. Choose one of: ${STEPS[*]}"; exit 1
fi

mkdir -p logs
STAMP=$(date +%Y%m%d_%H%M%S)
BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
echo "Pipeline started $(date) | branch: $BRANCH | from step: $START"
T0=$SECONDS

run_step () {
  local name="$1" script="$2"
  local log="logs/${STAMP}_${name}.log"
  local t=$SECONDS
  echo
  echo "=================== $name ($(date +%H:%M:%S)) ==================="
  if ! Rscript "$script" 2>&1 | tee "$log"; then
    echo
    echo "!! $name FAILED after $(( (SECONDS - t) / 60 )) min. Log: s/$log"
    echo "   Fix it, then resume with: ./run_pipeline.sh $name"
    exit 1
  fi
  echo "-- $name done in $(( (SECONDS - t) / 60 )) min (log: s/$log)"
}

[[ $RUN_CHECK -eq 1 ]] && run_step "self_check" "test_re_marginal.R"

started=0
for step in "${STEPS[@]}"; do
  [[ "$step" == "$START" ]] && started=1
  [[ $started -eq 1 ]] && run_step "$step" "${step}.R"
done

echo
echo "Pipeline complete in $(( (SECONDS - T0) / 60 )) min. Logs in s/logs/ (prefix $STAMP)."
