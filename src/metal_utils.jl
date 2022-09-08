#!/usr/bin/env julia
# Make Flux.gpu function use Metal GPU, if found (and Metal is installed)
try using Metal
catch e
    if e isa ArgumentError
        return # Metal not installed, not checking for metal GPU then
    else
        rethrow(e)
    end
end
# Any metal GPUs?
length(Metal.devices()) > 0 || return

using Flux
using Flux: Adapt

metal(x::AbstractArray{Float64}) = MtlArray(Float32.(x))
metal(x::AbstractArray{Float32}) = MtlArray(x)
metal(x::Int) = x # has to stay 64 bit for Flux signatures
metal(x::Float64) = Float32(x)
metal(x) = x # for anything else such as AdaptiveMaxPool and other layers
# inspired by Functor.jl and Adapt.jl
struct FluxMetalAdaptor end
Adapt.adapt_storage(to::FluxMetalAdaptor, x) = metal(x)
# redefine
Flux.gpu(x) = fmap(x -> Adapt.adapt(FluxMetalAdaptor(), x), x)

