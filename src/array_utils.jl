#!/usr/bin/env julia

# e.g. reshape is supposedly faster than cat but doesn't work for sparse 
# arrays. The cat method densifies a sparse array.
add_dim(x::Array) = reshape(x, (size(x)...,1))
add_dim(x::AbstractArray) = cat(x; dims=ndims(x)+1)


