#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

# unzip if not already done
for file in ../Model_?/{H,nodeCents}.zip; do
    if [ ! -d $file:r ]; then
        unzip $file -d $file:h
        rm -r $file:h/__MACOSX 2> /dev/null
    fi
done

hypergraph_CNN.jl --cv -m 8 -k 2 5 10 15 20 -f 64 -F 32 -e 500 -H ../Model_{1,3,4}/H/ -V ../Model_{1,3,4}/nodeCents/ --pred=pred_134.tsv --save-model=model_134.bson

