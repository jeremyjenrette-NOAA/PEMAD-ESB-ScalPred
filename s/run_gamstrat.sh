#!/bin/bash

# This line ensures the Mac stays awake until this shell script completes
caffeinate -w $$ &

# --- Your script logic goes here ---
echo "Starting gamstrat processing..."
Rscript ./gammod.R
echo "Done!"
