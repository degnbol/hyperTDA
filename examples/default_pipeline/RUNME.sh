#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"
cd $0:h

# example without interpolation and uninterpolation.

# start with a folder pointClouds/ with .npy or .tsv files where each file is a 
# matrix with 3 columns.

echo "persistent homology"
xyz2PH.jl pointClouds/ PH/

echo "get hypergraphs"
PH2hypergraph.jl PH/ H/

echo "get node centralities from hypergraphs created from persistent homology"
HG_nodeCent.jl PH/ nodeCents/

echo "community detection"
louvain.py PH/ pointClouds/ communities.json

echo "then open jupyter notebook and view communities and centrality"
jupyter notebook result_visualisation.ipynb

