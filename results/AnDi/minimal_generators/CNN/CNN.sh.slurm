#!/usr/bin/env zsh
mkdir -p CNN
sbatch -p mig -c 4 --mem=64G -t 8-00:00:00 -e CNN/slurm.err -o CNN/slurm.out ./CNN.sh
