#!/usr/bin/env julia
using Glob

"""
For some reason they didn't implement a glob that allows for absolute paths. This extension makes that possible.
https://githubmemory.com/repo/vtjnash/Glob.jl/issues/2
"""
function Glob.glob(pattern::String)
    first(pattern) == '/' || return Glob.glob(Glob.GlobMatch(pattern))
    Glob.glob(pattern[2:end], "/")
end

