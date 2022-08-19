#!/usr/bin/env zsh
# Load gurobi, e.g. if available as a module on a computer cluster
module load gurobi
# Locate license file for gurobi
export GRB_LICENSE_FILE=/usr/local/easybuild/software/Gurobi/gurobi.lic

julia -t 8 `git root`/src/xyz2PH_minimal.jl 'interpolated/2_*.tsv' PH_minimal/

