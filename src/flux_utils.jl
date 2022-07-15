#!/usr/bin/env zsh
using Flux
using Random # shuffle and shuffle!

"""
Onehot encoding. Class (label) along first dim as expected by 
(logit)crossentropy and other flux code.
BitMatrix type.
"""
function onehotBit(labels::Vector{Int}, uLabels=sort(unique(labels)))::BitMatrix
    labels' .== uLabels
end
"""
Onehot encoding. Class (label) along first dim as expected by 
(logit)crossentropy and other flux code.
Flux sparse int matrix type.
"""
function onehot(labels::Vector{Int}, uLabels=sort(unique(labels)))::Flux.OneHotMatrix
    Flux.onehotbatch(labels, uLabels)
end
"Inverse of the onehot functions."
function unonehot(labels::AbstractMatrix, uLabels=1:size(labels,1))::Vector{Int}
    iLabels = zeros(Int, size(labels, 2))
    for i ∈ uLabels
        iLabels[labels[i,:]] .= i
    end
    iLabels
end


get_batches(set::BitVector, batchsize::Int) = get_batches((1:length(set))[set], batchsize)
get_batches(set::Vector, batchsize::Int) = begin
    shuffle!(set)
    (set[i:i+batchsize] for i ∈ 1:batchsize:length(set)-batchsize)
end

"just run the main CV function with integer labels"
CV(labels::Union{BitMatrix,Flux.OneHotMatrix}, parts::Int=5) = CV(unonehot(labels), parts)
"just run the main CV function with integer labels"
CV(labels::BitVector, parts::Int=5) = CV(Int.(labels), parts)
"""
Create "parts" partitions by assigning each datapoint an integer in range 1:"parts".
Partitioning picks balanced amounts from each label.
"""
function CV(labels::Vector{Int}, parts::Int=5)::Vector{Int}
    inds = [shuffle(findall(labels .== i)) for i in unique(labels)]
    # number of observations of each label
    Ns = length.(inds)
    # size of each of the "parts" partitions.
    # The parts will not always add up perfectly to the total number of 
    # datapoints. We use ceil and later min(...) to make sure we use large 
    # enough parts to get all values (@assert !any(CV(labels) .== 0)).
    # The last part may be a few datapoints smaller than the others.
    set_sizes = ceil.(Int, Ns ./ parts)
    # fill in an integer for each datapoint in range 1:"parts"
    sets = zeros(Int, sum(Ns)) # prealloc
    for part in 1:parts
        for (_inds, set_size, N) in zip(inds, set_sizes, Ns)
            ind_range = (part-1)*set_size+1 : min(part*set_size, N)
            sets[_inds[ind_range]] .= part
        end
    end
    sets
end

