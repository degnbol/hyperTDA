#!/usr/bin/env zsh
# on spartan to use gurobi
# Gurobi.jl also can only be installed after loading module
module load Gurobi

PLPH=~/PLPH
$PLPH/src/xyz2PH_minimal.jl example_3.npy PH_minimal/

