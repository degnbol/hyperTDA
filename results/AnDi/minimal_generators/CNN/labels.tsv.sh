#!/usr/bin/env zsh
cut -d';' -f2 ../../labels.txt | cut -f1 -d'.' | cat <(echo diffusion) - |  mlr -t cat -n then put '$n = $n-1' > labels.tsv
