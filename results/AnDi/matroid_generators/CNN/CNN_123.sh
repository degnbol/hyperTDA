#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

date

hypergraph_CNN.jl --cv -m 8 -k 2 5 10 15 20 -f 64 -F 32 -e 500 -E 100 -H ../Model_{1,2,3}/H/ -V ../Model_{1,2,3}/nodeCents/ --pred=pred_123.tsv --save-model=model_123.bson

