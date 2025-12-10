#!/bin/bash

#SBATCH --job-name=model_vis
#SBATCH --output=model_vis_log.out
#SBATCH --error=model_vis_log.err
#SBATCH --time=12:00:00
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1536g

eval "$(conda shell.bash hook)"
conda activate python-env

python model_visualization.py