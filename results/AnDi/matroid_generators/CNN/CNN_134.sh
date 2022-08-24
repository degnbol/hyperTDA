#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

date

hypergraph_CNN.jl --cv -m 8 -k 2 5 10 15 20 -f 64 -F 32 -e 500 -E 100 -H ../Model_{1,3,4}/H/ -V ../Model_{1,3,4}/nodeCents/ --pred=pred_134.tsv --save-model=model_134.bson

