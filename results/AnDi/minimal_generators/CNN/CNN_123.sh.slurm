#!/usr/bin/env zsh
sbatch -p mig-gpu -J 'gpu123Mini' -G 1 -c 1 --mem=32G -t 8-00:00:00 -e slurm_123.err -o slurm_123.out ./CNN_123.sh
