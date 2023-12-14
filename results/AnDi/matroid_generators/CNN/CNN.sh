#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

date

mkdir -p models/
hypergraph_CNN.jl --cv -m 8 -k 2 5 10 15 20 -f 64 -F 32 -e 500 -E 100 -H ../Model_?/H/ -V ../Model_?/nodeCents/ --pred=pred.tsv --save-model models/model.bson

