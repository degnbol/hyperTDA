#!/usr/bin/env julia
# USAGE: ./install.jl [--global]
if length(ARGS) == 0
    env = true
elseif ARGS[1] == "--global"
    env = false
else
    error("Unknown option(s) $ARGS")
end

using Pkg

if env
    Pkg.activate(".")
    Pkg.instantiate()
else
    Pkg.add([
    "Flux",
    "Eirene",
    "ArgParse",
    "BSON",
    "Chain",
    "Combinatorics",
    "DelimitedFiles",
    "Distances",
    "Glob",
    "JSON",
    "NPZ",
    "Plots",
    "Printf",
    "Revise",
    "SharedArrays",
    "SparseArrays",
    "StatsBase",
    "PaddedViews",
    "Suppressor",
    "MKL"])
    Pkg.add(url="https://github.com/diegozea/ROC.jl")
end
