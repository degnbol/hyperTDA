#!/usr/bin/env zsh
for i in {0..4}; do
    diffusion=`mlr -t --ho --from labels.tsv filter '$diffusion == '$i then cut -f n | tr '\n' ,`
    echo 2_\{${diffusion:0:-1}\} # print without terminal ,
done > diffusions.txt
