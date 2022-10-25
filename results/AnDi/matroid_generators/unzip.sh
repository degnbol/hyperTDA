#!/usr/bin/env zsh
for dir in Model_?; do
    cd $dir
    unzip H.zip
    unzip nodeCents.zip
    cd -
done
