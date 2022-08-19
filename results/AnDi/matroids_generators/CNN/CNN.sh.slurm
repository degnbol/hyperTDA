#!/usr/bin/env zsh
sbatch -p mig -c 16 --mem=64G -t 8-00:00:00 -e slurm.err -o slurm.out ./CNN.sh
