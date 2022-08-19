#!/usr/bin/env zsh
sbatch -p mig -c 16 --mem=64G -t 8-00:00:00 -e slurm_134.err -o slurm_134.out ./CNN_134.sh
