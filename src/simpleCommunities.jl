#!/usr/bin/env julia
using SimpleHypergraphs
using DelimitedFiles
using JSON

ROOT = `git root` |> readchomp
folder, fname = "$ROOT/examples/default_pipeline/", "curve"
# folder, fname = "$ROOT/examples/with_interpolation/", "AnDi"
iCurve = 1
folder = expanduser(folder)

# H has dim #nodes × #hyperedges
mat = readdlm("$folder/H/$(fname)_$iCurve.csv", ',');
nodeCents = readdlm("$folder/nodeCents/$(fname)_$iCurve.tsv") |> vec
barcode = hcat(JSON.parsefile("$folder/PH/$(fname)_$iCurve.json")["barcode"]...)
persistences = barcode[:, 2] .- barcode[:, 1]

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


H = SimpleHypergraphs.Hypergraph{Float64,Float64}(nothingify(mat); v_meta=nodeCents, he_meta=persistences)
H = SimpleHypergraphs.Hypergraph(nothingify(mat))
H = SimpleHypergraphs.Hypergraph(nothingify(mat .* nodeCents))
H = SimpleHypergraphs.Hypergraph(nothingify(mat .* persistences'))

best = findcommunities(H, CFModularityRandom(1, 100_000))
bestbm = best.bm
for n in 2:20
    res = findcommunities(H, CFModularityRandom(n, 100_000))
    if res.bm ≥ bestbm
        println(n)
        best = res
        bestbm = best.bm
    else
        break
    end
end
# always n=2 is best, which is obvi, since n=1 always has bm=0 for CFModularityRandom

sum(length.(findcommunities(H, CFModularityCNMLike(100_000)).bp) .> 1)
# mostly trivial. It only gives us 1 or 2 non singleton communities.


## Let's try HyperModularity
using HyperModularity
using Statistics

function HyperModularity.hypergraph(mat::Matrix)
    N = 1:size(mat,1)
    D = sum(mat .!= 0.; dims=2) |> vec
    E = Dict{Int,Dict}()
    for he in eachcol(mat)
        n = sum(he .!= 0.)
        if !haskey(E, n) E[n] = Dict() end
        E[n][findall(he .!= 0)] = mean(he)
    end
    HyperModularity.hypergraph(N, E, D)
end

function HyperModularity.CliqueExpansionModularity(mat::Matrix)
    HyperModularity.CliqueExpansionModularity(HyperModularity.hypergraph(mat))
end

CliqueExpansionModularity(mat)
# looks good, as expected. Should be what we already have done.

n = size(mat,1)
Z_ = collect(1:n) # trivial clustering

# all or nothing aggregator: p -> [length(p) == 1, sum(p)]
# This gives a starter estimate for Ω, from a trivial clustering Z_
Ω = estimateΩEmpirically(H, Z_; aggregator = p -> [length(p) == 1, sum(p)])

Z = AON_Louvain(H,Ω; α=0) |> unique |> length

# Alternatively, one can learn Ω from graph Louvain solution Z_g
Ω = estimateΩEmpirically(H, Z_g; aggregator = p -> [length(p) == 1, sum(p)])
Z = AON_Louvain(H,Ω; α=0) |> unique |> length

# In conclusion, all the results that aren't what we already do (clique expans.) are trivial.

