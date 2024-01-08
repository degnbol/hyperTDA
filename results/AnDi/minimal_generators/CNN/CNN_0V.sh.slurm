#!/usr/bin/env zsh
sbatch -p mig-gpu -J 'gpuMini' -G 1 -c 1 --mem=32G -t 8-00:00:00 -e slurm_0V.err -o slurm_0V.out ./CNN_0V.sh
