#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

date

mkdir -p models/
# has to be julia v1.8.3
$HOME/packages/julias/julia-1.8/bin/julia $ROOT/src/hypergraph_CNN.jl --cv -m 8 -k 2 5 10 15 20 -f 64 -F 32 -e 500 -E 100 --0H -H ../Model_?/H/ -V ../Model_?/nodeCents/ --pred=pred_0H.tsv --save-model models/model_0H.bson

