#!/usr/bin/env julia
using DelimitedFiles
include("PH_hypergraph_centralities.jl")
include("glob.jl")

# write node centrality column vector (type=max) given PH jsons.

INDIR, OUTDIR = ARGS

mkpath(OUTDIR)

for file in glob("$INDIR/*.json")
    fname = splitext(basename(file))[1]
    # print(fname)
    nodeCents, _ = read_PH2centralities("$file", "max") 
    writedlm("$OUTDIR/$fname.tsv", nodeCents)
end

