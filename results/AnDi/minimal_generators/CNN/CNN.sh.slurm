#!/usr/bin/env zsh
sbatch -p mig-gpu -J 'gpuMini' -G 1 -c 1 --mem=32G -t 8-00:00:00 -e slurm.err -o slurm.out ./CNN.sh
