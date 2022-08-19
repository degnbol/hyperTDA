#!/usr/bin/env zsh
mkdir -p CNN
sbatch -p mig -c 8 --mem=64G -t 8-00:00:00 -e CNN/slurm_134.err -o CNN/slurm_134.out ./CNN_134.sh
