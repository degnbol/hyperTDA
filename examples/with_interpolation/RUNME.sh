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
uninterpolate.jl --interp=interpolated/ nodeCents/ nodeCents_uninterp/

echo "uninterpolate hypergraph matrix created from PH"
uninterpolate.jl --interp=interpolated/ PH/ H/ --entry=representatives

echo "uninterpolate community associations"
uninterpolate.jl --interp=interpolated/ communities.json communities_uninterp/ --categorical

echo "then open jupyter notebook and view communities and centrality"
jupyter notebook result_visualisation.ipynb

