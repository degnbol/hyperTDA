#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

# d1=`sed '1!d' diffusions.txt`
d2=`sed '2!d' diffusions.txt`
d3=`sed '3!d' diffusions.txt`
d4=`sed '4!d' diffusions.txt`
# d5=`sed '5!d' diffusions.txt`

hypergraph_CNN.jl --cv -m 8 -k 2 5 10 15 20 -f 64 -F 32 -e 500 -E 100 -H ../H/{$d2,$d3,$d4}.csv -V ../nodeCents_uninterp/{$d2,$d3,$d4}.tsv --pred=pred_123.tsv --save-model=model_123.bson

