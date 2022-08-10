#!/usr/bin/env zsh
ROOT=`git root`
export PATH="$PATH:$ROOT/src"

# start with a folder pointClouds/ with .npy or .tsv files where each file is a 
# matrix with 3 columns.

# interpolate between points
interpolate.jl pointClouds/ interpolated/ --points=100 --nans=0

# calculate persistent homology
julia --threads=8 xyz2PH.jl interpolated/ PH/

# get node centralities from hypergraphs created from persistent homology
HG_nodeCent.jl PH/ nodeCents/

# uninterpolate node values in case of interpolation
uninterpolate.jl nodeCents/ nodeCents_uninterp/ --interp=interpolated/
# uninterpolate incidence matrix created from PH in case of interpolation
uninterpolate.jl PH/ H/ --interp=interpolated/ --entry=representatives

# Run CNN on uninterpolated H and V to predict e.g. diffusion model
# TODO set $pred, $V and $H, were we doing any example weighing edges?
julia --threads=8 hypergraph_CNN.jl -m 8 -f 64 -F 32 -k 2 5 10 15 20 -e 500 --pair-files --cv \
    -p $pred -V 'nodeCents_uninterp/*' -H 'H/*'

