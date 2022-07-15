#!/usr/bin/env julia

"""
Get longest prefix in common for a vector of strings.
"""
function commonPrefixLength(ss::Vector)::Int
    len = minimum(length.(ss))
    for i in 1:len
        if !all(startswith.(ss, ss[1][1:i]))
            return i-1
        end
    end
    len
end

"""
Get longest suffix in common for a vector of strings.
"""
function commonSuffixLength(ss::Vector)::Int
    len = minimum(length.(ss))
    for i in 1:len
        if !all(endswith.(ss, ss[1][end-i+1:end]))
            return i-1
        end
    end
    len
end

"""
Given strings
a = ["prefixA...3...suffixA", "prefixA...1...suffixA", "prefixA...2...suffixA"]
b = ["prefixB...2...suffixB", "prefixB...1...suffixB"]
Get the intersecting pairs containing 1 and 2, as well as the indexing.
"""
function prefixSuffixPairs(a::Vector, b::Vector)
    prefixA = commonPrefixLength(a)
    suffixA = commonSuffixLength(a)
    prefixB = commonPrefixLength(b)
    suffixB = commonSuffixLength(b)
    
    middleA = [s[prefixA+1:end-suffixA] for s ∈ a]
    middleB = [s[prefixB+1:end-suffixB] for s ∈ b]
    
    idxA = []
    idxB = []
    for (i, mA) ∈ enumerate(middleA)
        j = findfirst(middleB .== mA)
        if j !== nothing
            push!(idxA, i)
            push!(idxB, j)
        end
    end
    
    a[idxA], b[idxB], idxA, idxB
end

