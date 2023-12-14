#!/usr/bin/env zsh
sbatch -p mig-gpu -c 1 -G 1 -J 'gpu' --mem=32G -t 10-00:00:00 -e slurm_0V.err -o slurm_0V.out ./CNN_0V.sh
