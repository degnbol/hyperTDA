#!/usr/bin/env julia
using JSON
using SparseArrays
include("hypergraph_utils.jl")

# testpath = readchomp(`git root`) * "/Data/artificialSampleData/PH/generated_curve_0.json"

function representatives2hyperedges(representatives, edgesAsHypernodes::Bool=false)
    if !edgesAsHypernodes
        # Each representative is a generator where we take all nodes to form a hyperedge.
        return [Set(i for e in es for i in e) for es in representatives]
    else
        N = maximum(i for es in representatives for e in es for i in e)
        # each 2-tuple becomes a hypernode 
        # so we will assign new int indexes for each element in a N×N adjacency matrix.
        ij = LinearIndices((N,N))
        return [Set(ij[e...] for e in es) for es in representatives]
    end
end

"""
- path: string path to json file with field "representatives".
- edgesAsHypernodes: should nodes or edges in the curve be made into "hyper" nodes (nodes in the hypegraph)?
    Default: nodes become nodes in the hypergraph. This means a representative is reduced from a set of 2-tuples to a set of nodes.
return: B, diag(W)
Their read functions also returned hyperedges as dict mapping from sets of node ids => diag(W) elements
but that is never used.
"""
function read_PH2hypergraph(path::AbstractString; edgesAsHypernodes::Bool=false)
    PH = JSON.parsefile(path)
    
    barcodes = hcat(PH["barcode"]...)
    @assert size(barcodes, 2) == 2 "Update code to handle transposed barcodes"
    # persistences will be used as hyper edge weights
    persistences = barcodes[:, 2] - barcodes[:, 1]
    @assert all(persistences .> 0)

    # representatives are each 2-tuples (edges). 
    representatives = PH["representatives"]
    hyperedges = representatives2hyperedges(representatives, edgesAsHypernodes)
    hyperedges2B(hyperedges), persistences
end

