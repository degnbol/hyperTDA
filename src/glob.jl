#!/usr/bin/env julia
using Glob
using Glob: glob

"""
For some reason they didn't implement a glob that allows for absolute paths.
This extension makes that possible.
https://githubmemory.com/repo/vtjnash/Glob.jl/issues/2
"""
function Glob.glob(pattern::String)
    first(pattern) == '/' || return glob(Glob.GlobMatch(pattern))
    glob(pattern[2:end], "/")
end

"""
Extend glob from Glob, so far only by expanding {a,b,c,...}.
"""
function eglob(pattern::String)
    patterns = expandCurly(pattern)
    vcat(glob.(patterns)...)
end

"""
Expand zero or more curly brackets as a shell would.
Not tested for nested brackets.
"""
function expandCurly(s::String)::Vector{String}
    ranges = findall(r"{[^{}]*}", s)
    length(ranges) > 0 || return [s]
    pieces = []
    lastStop = 0
    for rang in ranges
        push!(pieces, [s[lastStop+1:rang.start-1]])
        # add each string between commas within a set of curls.
        push!(pieces, split(s[rang.start+1:rang.stop-1], ','))
        lastStop = rang.stop
    end
    push!(pieces, [s[ranges[end].stop+1:end]])
    Iterators.product(pieces...) .|> prod |> vec
end

