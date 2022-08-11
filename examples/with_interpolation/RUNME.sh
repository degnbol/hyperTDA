#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"
cd $0:h

# start with a folder pointClouds/ with .npy or .tsv files where each file is a 
# matrix with 3 columns.

echo "interpolate between points"
interpolate.jl pointClouds/ --out interpolated/ --points=500 --nans=0

echo "calculate persistent homology on interpolated points"
xyz2PH.jl interpolated/ PH/

echo "get node centralities from hypergraphs created from persistent homology"
HG_nodeCent.jl PH/ nodeCents/

echo "community detection"
louvain.py PH/ interpolated/ communities.json

echo "uninterpolate node values"
uninterpolate.jl nodeCents/ nodeCents_uninterp/ --interp=interpolated/

echo "uninterpolate hypergraph matrix created from PH"
uninterpolate.jl PH/ H/ --interp=interpolated/ --entry=representatives

echo "then open jupyter notebook and view communities and centrality"
jupyter notebook result_visualisation.ipynb

