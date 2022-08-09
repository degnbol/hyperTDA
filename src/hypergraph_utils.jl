#!/usr/bin/env julia
using SparseArrays

"""
B is sparse integer matrix with dim (#nodes, #hyperedges).
- hyperedges: each hyperedge is a set of node indices
- nNodes: optionally specify the total number of nodes
"""
function hyperedges2B(hyperedges::Vector{Set{Int}}, nNodes::Union{Int,Nothing}=nothing)
    Is = [n for     h  in           hyperedges  for n in h] 
    Js = [j for (j, h) in enumerate(hyperedges) for n in h]
    if nNodes === nothing
        sparse(Is, Js, 1)
    else
        sparse(Is, Js, 1, nNodes, length(hyperedges))
    end
end


"""
A_H defined in node and hyperedge centrality paper.
The adjacency matrix created from a hypergraph.
Symmetric.
"""
function get_A_H(
        hyperedges::Vector{Set{Int}}, 
        hyperedge_weights::Vector=ones(length(hyperedges)),
        nNodes::Union{Int,Nothing}=nothing)
    w = hyperedge_weights
    Is = [i    for (k, he) in enumerate(hyperedges) for i in he for j in he]
    Js = [j    for (k, he) in enumerate(hyperedges) for i in he for j in he]
    ws = [w[k] for (k, he) in enumerate(hyperedges) for i in he for j in he]
    ij = Is .!= Js
    # sparse has default combine function sum, 
    # which is what we want to do when the same ij appears in multiple hyperedges.
    if nNodes === nothing
        sparse(Is[ij], Js[ij], ws[ij])
    else
        sparse(Is[ij], Js[ij], ws[ij], nNodes, nNodes)
    end
end

"""
Just like A_H but only summing for the two element edges present in the representatives, 
rather than all pairs in a given hyperedge.
So, if hyperedge_weights is left as default then it is simply an adjancecy 
matrix created from the persistent homology cycle representatives.
"""
function get_A_repr(
        representatives::Vector{Vector{Vector{Int}}}, 
        hyperedge_weights::Vector=ones(length(representatives)),
        nNodes::Union{Int,Nothing}=nothing)
    w = hyperedge_weights
    # each e has two nodes. We run through them forwards and backwards to get symmetry.
    Is = [n    for (k, he) in enumerate(representatives) for e in he for n in e]
    Js = [n    for (k, he) in enumerate(representatives) for e in he for n in e[end:-1:1]]
    ws = [w[k] for (k, he) in enumerate(representatives) for e in he for n in e]
    ij = Is .!= Js
    # sparse has default combine function sum, 
    # which is what we want to do when the same ij appears in multiple hyperedges.
    if nNodes === nothing
        sparse(Is[ij], Js[ij], ws[ij])
    else
        sparse(Is[ij], Js[ij], ws[ij], nNodes, nNodes)
    end
end

"Add a default that is different from nothing"
function Base.findfirst(A::AbstractVector, default::Real)
    out = findfirst(A)
    out === nothing ? default : out
end
"Add a default that is different from nothing"
function Base.findlast(A::AbstractVector, default::Real)
    out = findlast(A)
    out === nothing ? default : out
end

"""
Align first nonzero entry in each hyperedge.
Hyperedges are also truncated to remove trailing zero entries. To avoid this, 
set len, e.g. to length of input hyperedges.
This function is not well defined for hyperedges with all zeros, i.e. no members.
They are currently handled by "aligning" a single zero entry, i.e. the 
resulting aligned hyperedge will be all zeros just like the input.
- H: hypergraph as an incidence matrix with size=(#nodes, #hyperedges).
"""
function nzalign(H::AbstractMatrix, len::Union{Int,Nothing}=nothing)
    nzs = eachcol(H .!= 0)
    
    # deal with all zero hyperedges by having a default value 1 instead of 
    # nothing
    starts = findfirst.(nzs, 1)
    stops = findlast.(nzs, 1)
    
    if len === nothing
        len = maximum(stops .- starts) + 1
    end
    
    out = zeros(len, length(nzs))
    
    for (iE, (start, stop)) in enumerate(zip(starts, stops))
        out[1:stop-start+1, iE] = H[start:stop, iE]
    end
    
    out
end
"Same as nzalign for single hypergraph H, except length truncation is the same for all Hs."
function nzalign(Hs::Vector)
    nzs = [eachcol(H .!= 0) for H ∈ Hs]
    
    # deal with all zero hyperedges by having a default value 1 instead of 
    # nothing
    starts = [findfirst.(nz, 1) for nz ∈ nzs]
    stops = [findlast.(nz, 1) for nz ∈ nzs]
    
    len = maximum(vcat(stops...) .- vcat(starts...)) + 1
    
    nzalign.(Hs, len)
end
nzalign(Hs::Vector, len::Int) = nzalign.(Hs, len)


"Vector wrapper version of the more general function that takes a matrix."
function nzranges(H::AbstractMatrix, W::AbstractVector, len::Union{Int,Nothing}=nothing)
    H, W = nzranges(H, add_dim(W)', len)
    H, vec(W)
end
"""
Generate hyperedges for all ranges (1:end, 1:end-1, ..., end-1:end) which 
contains nonzero elements and doesn't exclude any nonzero elements.
- H: hypergraph incidence matrix, with dim #nodes×#hEdges.
- W: hyperedge weights or other edge value(s). Second dim is #hEdges, but first dim can be any number ∈ [0,1,...]
return: H, W with cycled hyperedges.
"""
function nzranges(H::AbstractMatrix, W::AbstractMatrix, len::Union{Int,Nothing}=nothing)
    @assert size(H,2) == size(W,2)
    nzs = eachcol(H .!= 0)
    
    starts = findfirst.(nzs, 1)
    stops = findlast.(nzs, 1)
    if len === nothing
        len = maximum(stops .- starts) + 1
    end
    
    copies = [len - (b-a) for (a,b) ∈ zip(starts, stops)]
    
    W_out = hcat((hcat(fill(W[:,i], nCopies)...) for (i, nCopies) ∈ enumerate(copies))...)
    
    Is, Js, Vs = Int[], Int[], eltype(H)[]
    # out_col_group: column index in output matrix for the start of each group 
    # where a group contains multiple copies of the same input column with some 
    # kind of shift in the first dimension.
    # in_col: column index in input matrix.
    # nCopies: number of copies to make at different shifts for current 
    # hyperedge.
    out_col_group = 0
    for (in_col, (a, b, nCopies)) ∈ enumerate(zip(starts, stops, copies))
        for iCopy ∈ 1:nCopies
            append!(Is, iCopy:iCopy+(b-a))
            append!(Js, (out_col_group+iCopy for _ ∈ a:b))
            append!(Vs, H[a:b, in_col])
        end
        out_col_group += nCopies
    end
    # fmt = row wise since each row will be a hyperedge, and we need to look at 
    # each hyperedge one at a time. Use it if you implement code for CuVectors 
    # (GPU)
    # H_out = sparse(Is, Js, Vs, len, sum(copies); fmt=:csr)
    H_out = sparse(Is, Js, Vs, len, sum(copies))
    
    H_out, W_out
end

