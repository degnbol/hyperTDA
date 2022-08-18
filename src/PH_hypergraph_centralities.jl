#!/usr/bin/env julia
using LinearAlgebra
using StatsBase
using JSON
include("PH2hypergraph.jl") # read_PH2hypergraph
ROOT = readchomp(`git root`)
TUDSICO = "$ROOT/Publications/Tudisco_2021/node-edge-hypergraph-centrality"
include("$TUDSICO/centrality_tools.jl") # compute_centrality
using Suppressor # kill stdout and stderr

fgϕψs = Dict("linear"  => (x -> x, x -> x,         x -> x,       x -> x        ), 
             "log-exp" => (x -> x, x -> x.^(1/10), x -> log.(x), x -> exp.(x)  ),
             "max"     => (x -> x, x -> x.^(1/5),  x -> x.^15,   x -> x.^(1/15)))

function compute_centrality(B, fgϕψ::String; kwargs...)
    @suppress compute_centrality(B, fgϕψs[lowercase(fgϕψ)]; kwargs...)
end

"""
Compute node and edge centralities from a hypergraph made from PH.
- path: path to PH .json
- fgϕψ: string "linear", "log-exp", "max" or a 4 tuple of functions (f, g, ϕ, ψ)
return: node centrality x, edge centrality y
"""
function read_PH2centralities(path::AbstractString, fgϕψ="linear"; edgesAsHypernodes::Bool=false)
    B, diagW = read_PH2hypergraph(path; edgesAsHypernodes=edgesAsHypernodes)
    compute_centrality(B, fgϕψ; edge_weights=diagW)
end

"""
Read PH, make hypergraph, calculate node and edge centralities,
write to json with entries "node_centrality" and "hyperedge_centrality".
"""
function readwrite_PH2centralities(infname::AbstractString, outfname::AbstractString, fgϕψ="linear"; edgesAsHypernodes::Bool=false)
    x, y = read_PH2centralities(infname, fgϕψ; edgesAsHypernodes=edgesAsHypernodes)
    open(outfname, "w") do io JSON.print(io, Dict("node_centrality"=>x, "hyperedge_centrality"=>y), 1) end
end

