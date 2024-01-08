#!/usr/bin/env zsh
sbatch -p mig-gpu -J 'gpu134Mini' -G 1 -c 1 --mem=32G -t 8-00:00:00 -e slurm_0H_134.err -o slurm_0H_134.out ./CNN_0H_134.sh
