#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

# example without interpolation and uninterpolation.

# start with a folder pointClouds/ with .npy or .tsv files where each file is a 
# matrix with 3 columns.

# calculate persistent homology
julia --threads=8 xyz2PH.jl pointClouds/ PH/

# get node centralities from hypergraphs created from persistent homology
HG_nodeCent.jl PH/ nodeCents/

# get incidence matrices
PH2hypergraph.jl PH/ H/

# Run CNN on uninterpolated H and V to predict e.g. diffusion model
# TODO set $pred, $V and $H, were we doing any example weighing edges?
julia --threads=8 hypergraph_CNN.jl -m 8 -f 64 -F 32 -k 2 5 10 15 20 -e 500 --pair-files --cv \
    -p $pred -V 'nodeCents_uninterp/*' -H 'H/*'

