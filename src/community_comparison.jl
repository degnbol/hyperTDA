#!/usr/bin/env julia
using DelimitedFiles
using DataFrames, CSV
using JSON
using SimpleHypergraphs
using HyperModularity
using Statistics
using Conda, PyCall
# INSTALL:
# Conda.pip_interop(true)
# Conda.pip("install", "leidenalg")
@pyinclude "leiden.py"

ROOT = `git root` |> readchomp

mats = Matrix{Float64}[]
nodeCents = Vector{Float64}[]
persistences = Vector{Float64}[]

for (interp, fname, centDir) in zip(["default_pipeline", "with_interpolation"], ["curve", "AnDi"], ["nodeCents", "nodeCents_uninterp"])
    folder = "$ROOT/examples/$interp/"
    for iCurve in 1:2
        # H has dim #nodes × #hyperedges
        push!(mats, readdlm("$folder/H/$(fname)_$iCurve.csv", ','))
        push!(nodeCents, readdlm("$folder/$centDir/$(fname)_$iCurve.tsv") |> vec)
        barcode = hcat(JSON.parsefile("$folder/PH/$(fname)_$iCurve.json")["barcode"]...)
        push!(persistences, barcode[:, 2] .- barcode[:, 1])
    end
end

# first try SimpleHypergraphs

"""
use nothing instead of zero for non-edges
since it is what SimpleHypergraphs expects
"""
function nothingify(m::Matrix{T}) where T<:Real
    m = Matrix{Union{T,Nothing}}(m)
    m[m .== 0] .= nothing
    m
end
function nothingify(m::BitMatrix)
    m = Matrix{Union{Bool,Nothing}}(m)
    m[m .== 0] .= nothing
    m
end

function SimpleHypergraphsSets2Comm(bp)
    N = maximum.(bp) |> maximum
    comm = zeros(Int, N)
    for (i, c) in enumerate(bp)
        if length(c) > 1
            comm[collect(c)] .= i
        end
    end
    comm
end

function CFMCNMLike(H)
    best = nothing
    bestbm = 0.
    # individual runs will vary, so we have multiple attempts and keep the best.
    # The best usually has the same or higher number of non-singleton communities.
    for _ in 1:100
        res = findcommunities(H, CFModularityCNMLike(1_000))
        if res.bm ≥ bestbm
            println(sum(length.(res.bp) .> 1))
            best = res
            bestbm = best.bm
        end
    end
    best.bp |> SimpleHypergraphsSets2Comm
end

function CFMRandom(H)
    best = nothing
    bestbm = 0 # n=1 always gives bm==0
    for n in 1:20
        res = findcommunities(H, CFModularityRandom(n, 1_000))
        if res.bm ≥ bestbm
            println(n)
            best = res
            bestbm = best.bm
        else
            break
        end
    end
    # always n=2 is best, which is obvi, since n=1 always has bm=0 for CFModularityRandom
    best.bp |> SimpleHypergraphsSets2Comm
end


function HyperModularity.hypergraph(mat::Matrix)
    N = 1:size(mat,1)
    D = sum(mat .!= 0.; dims=2) |> vec
    E = Dict{Int,Dict}()
    for he in eachcol(mat)
        n = sum(he .!= 0.)
        if !haskey(E, n) E[n] = Dict() end
        E[n][findall(he .!= 0.)] = mean(he)
    end
    HyperModularity.hypergraph(N, E, D)
end

function HyperModularity.CliqueExpansionModularity(mat::Matrix)
    HyperModularity.CliqueExpansionModularity(HyperModularity.hypergraph(mat))
end

"""
- Z_: initial clustering. Default = trivial. Solution from clique expansion Louvain can be used here.
"""
function HyperModularity.AON_Louvain(mat::Matrix, Z_::Vector{Int}=collect(1:size(mat,1)))
    H = HyperModularity.hypergraph(mat)
    # all or nothing aggregator: p -> [length(p) == 1, sum(p)]
    Ω = estimateΩEmpirically(H, Z_; aggregator = p -> [length(p) == 1, sum(p)])
    AON_Louvain(H, Ω; α=0)
end


function CliqueExpansion(mat)
    N = size(mat, 1)
    ex = zeros(N, N)
    for he in eachcol(mat)
        Is = findall(he .!= 0)
        for i in Is
            for j in Is
                if i != j
                    ex[i,j] = (he[i] + he[j]) / 2
                end
            end
        end
    end
    ex
end

leiden(mat) = py"leiden"(CliqueExpansion(mat))

# for the results
df = DataFrame(dataset=String[], dataseti=Int[], method=String[], E=Bool[], V=Bool[], residue=Int[], community=Int[])

for (i, (dataset, dataseti)) in enumerate([("default_pipeline", 1), ("default_pipeline", 2), ("with_interpolation", 1), ("with_interpolation", 2)])
    println("#", dataset, " ", dataseti)
    mat = mats[i]
    nodeCent = nodeCents[i]
    persistence = persistences[i]

    for E in [false, true]
        for V in [false, true]
            _mat = E ? mat .* persistence' : mat
            _mat = V ? _mat .* nodeCent : _mat

            for (method, fMethod) in [("CFMCNMLike", CFMCNMLike), ("CFMRandom", CFMRandom)]
                println(method)
                H = SimpleHypergraphs.Hypergraph(nothingify(_mat))
                comm = fMethod(H)
                append!(df, DataFrame(dataset=dataset, dataseti=dataseti, method=method, E=E, V=V, residue=1:length(comm), community=comm))
            end

            println("CliqueLouvain")
            comm = CliqueExpansionModularity(_mat)
            append!(df, DataFrame(dataset=dataset, dataseti=dataseti, method="CliqueLouvain", E=E, V=V, residue=1:length(comm), community=comm))

            println("AON_Louvain")
            comm = AON_Louvain(_mat)
            append!(df, DataFrame(dataset=dataset, dataseti=dataseti, method="AON_Louvain", E=E, V=V, residue=1:length(comm), community=comm))

            println("Leiden")
            comm = leiden(_mat)
            append!(df, DataFrame(dataset=dataset, dataseti=dataseti, method="Leiden", E=E, V=V, residue=1:length(comm), community=comm))
        end
    end
end

CSV.write("community_comparison.tsv.gz", df; delim='\t', compress=true)


