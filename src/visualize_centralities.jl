#!/usr/bin/env julia
using LinearAlgebra
using SparseArrays
using StatsBase
using PyPlot

"""
Scatterplot edge centrality vs edge weights
"""
plot_edge_centrality_weight(W::AbstractMatrix, y::Vector) = plot_edge_centrality_weight(diag(W), y)
function plot_edge_centrality_weight(diagW::Vector, y::Vector)
    ftlabels = 15
    ftticks = 12
    fttitle = 18
    
    f = figure()
    suptitle("edge centrality vs edge weights", fontsize=fttitle)
    plot(y./maximum(y), diagW./maximum(diagW), "o")
    ylabel("edge weight", fontsize=ftlabels)
    xlabel("edge centrality", fontsize=ftlabels)
    xticks(fontsize=ftticks)
    yticks(fontsize=ftticks)
end

