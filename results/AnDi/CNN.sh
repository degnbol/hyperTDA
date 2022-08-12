#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

# unzip if not already done
for file in Model_?/{H,nodeCents}.zip; do
    if [ ! -d $file:r ]; then
        unzip $file -d $file:h
        rm -r $file:h/__MACOSX 2> /dev/null
    fi
done

mkdir -p CNN

hypergraph_CNN.jl --cv -m 8 -k 2 5 10 15 20 -f 64 -F 32 -H Model_?/H/ -V Model_?/nodeCents/ --pred=CNN/pred.tsv --save-model CNN/model.bson

