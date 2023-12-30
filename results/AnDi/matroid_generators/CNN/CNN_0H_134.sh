#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

date

mkdir -p models/
# has to be julia v1.8.3
$HOME/packages/julias/julia-1.8/bin/julia $ROOT/src/hypergraph_CNN.jl --cv -m 8 -k 2 5 10 15 20 -f 64 -F 32 -e 500 -E 100 --0H -H ../Model_{1,3,4}/H/ -V ../Model_{1,3,4}/nodeCents/ --pred=pred_0V_134.tsv --save-model models/model_0V_134.bson

