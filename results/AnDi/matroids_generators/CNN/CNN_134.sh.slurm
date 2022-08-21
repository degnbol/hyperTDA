#!/usr/bin/env zsh
sbatch -p mig-gpu -c 1 -G 1 -J 'gpu134' --mem=32G -t 10-00:00:00 -e slurm_134.err -o slurm_134.out ./CNN_134.sh
