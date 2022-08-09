#!/usr/bin/env julia
using JSON
using StatsBase
using Printf
using NPZ
using SparseArrays
using Chain: @chain
using DataFrames, CSV
using DelimitedFiles
ROOT = readchomp(`git root`)
include("hypergraph_utils.jl") # hyperedges2B, get_A_H, get_A_repr


norm2(xyzs::Matrix) = sum(xyzs .^ 2; dims=2) .|> sqrt
distances(xyz::AbstractVector, xyzs::Matrix) = norm2(xyz' .- xyzs) |> vec

"""
Get the index of a xyz point in the list (::Matrix) of interpolated points.
"""
function un2interp(xyz::AbstractVector, xyzs_i::Matrix)::Int
    @chain xyz' .== xyzs_i all(dims=2) vec findall only
end

"""
Default function to get matrix D.
Rolling weighted average sparse conversion matrix.
Each interpolated point gets two weight values, one for each nearest real neighbor point. 
The weights are based on eucledian distance, i.e. w1 = 1 - d1/(d1+d2) = d2/(d1+d2)
The real points get an initial weight=1.
By initial, it is meant that weights will be normalized (otherwise points with more interpolation neighbors would get falsely emphasized).
- xyzs: real xyz points ordered along the curve
- xyzs_i: real mixed with interpolated xyz points, still ordered along the curve
return: matrix that converts from interpolated xyz to uninterpolated xyz. 
E.g matrix is D, we want to convert node values b, we do D*b. 
Similar if we want to convert hypergraph membership with D*B where B has shape (#interp, #hyperedges).
"""
function uninterp_average_weight(xyzs::Matrix, xyzs_i::Matrix)
    @assert !any(isnan.(xyzs)) "NaN in xyzs at $(findall(any(isnan.(xyzs); dims=2)))"
    # build sparse weights from xyz_i to xyz
    xyz2interp = un2interp.(eachrow(xyzs), Ref(xyzs_i))
    uninterp_average_weight(xyz2interp, xyzs_i)
end
"""
-interp_idx: true for interpolation points and false for real points.
"""
function uninterp_average_weight(interp_idx::Union{BitVector,Vector{Bool}}, xyzs_i::Matrix)
    uninterp_average_weight(findall(.!interp_idx), xyzs_i)
end
"""
-xyz2interp: int index for the location of each real point within xyzs_i.
"""
function uninterp_average_weight(xyz2interp::Vector{Int}, xyzs_i::Matrix)
    @assert !any(isnan.(xyzs_i)) "NaN in xyzs_i at $(findall(any(isnan.(xyzs_i); dims=2)))"
    # make sure interp_idx (should be boolean) wasn't accidentally given as ints ∈ {0, 1}.
    if all((xyz2interp .== 0) .| (xyz2interp .== 1))
        @info "Assuming ints ∈ {0,1} given are supposed to be bool."
        return uninterp_average_weight(BitVector(xyz2interp), xyzs_i)
    end
    
    Is, Js, ws = Int[], Int[], Float64[]
    for i in 1:length(xyz2interp)-1
        # slice range for interp xyzs between real xyz i and i+1 (inclusive in both ends).
        _js = xyz2interp[i]:xyz2interp[i+1]
        xyz1to2 = xyzs_i[_js, :]
        xyz1 = xyz1to2[begin, :]
        xyz2 = xyz1to2[end, :]
        # distance to the neighbor real points
        d1s = sum((xyz1to2 .- xyz1') .^ 2; dims=2)
        d2s = sum((xyz1to2 .- xyz2') .^ 2; dims=2)
        # weightings toward the neighbor real points
        w1s = @. d2s / (d1s + d2s)
        w2s = @. d1s / (d1s + d2s)
        # building sparse matrix
        append!(Is, fill(i, length(w1s)))
        append!(Is, fill(i+1, length(w2s)))
        append!(Js, _js)
        append!(Js, _js)
        append!(ws, w1s)
        append!(ws, w2s)
    end
    
    # use max instead of add as combine function since we are adding real points twice 
    # (except for the very first and last)
    D = sparse(Is, Js, ws, length(xyz2interp), size(xyzs_i, 1), max)
    # norm so weights sum to 1 for each real node
    D = D ./ sum(D; dims=2) |> dropzeros
end
"-df: DataFrame with x, y, z columns and bool column indicating interpolation."
function uninterp_average_weight(df::DataFrame; interp::AbstractString="interp",
        x::AbstractString="x", y::AbstractString="y", z::AbstractString="z")
    uninterp_average_weight(BitVector(xyzs_i[:, interp]), Matrix(xyzs_i[:, [x,y,z]]))   
end

