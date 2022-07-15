#!/usr/bin/env julia
using ROC # https://github.com/diegozea/ROC.jl
# works for both vectors and matrices. 
AUC(ŷ, y) = ROC.roc(ŷ, y) |> ROC.AUC
