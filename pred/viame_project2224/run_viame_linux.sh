#!/bin/bash
#SBATCH -J viame_inf
#SBATCH --account=sharkpulse
#SBATCH --partition=a30_normal_q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --output=/projects/sharkpulse/archived/PEMAD-ESB-ScalPred/pred/viame_project2224/log/viame_inf_%j.out
#SBATCH --error=/projects/sharkpulse/archived/PEMAD-ESB-ScalPred/pred/viame_project2224/log/viame_inf_%j.err

BASE_DIR="/projects/sharkpulse/archived/PEMAD-ESB-ScalPred/pred/viame_project2224"
cd "${BASE_DIR}"

bash run_trained_model.sh