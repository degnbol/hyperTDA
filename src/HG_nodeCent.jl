#!/usr/bin/env julia
using DelimitedFiles
include("PH_hypergraph_centralities.jl")
include("glob.jl")

# write node centrality column vector (type=max) given PH jsons.

INDIR, OUTDIR = ARGS

mkpath(OUTDIR)

infiles = glob("$INDIR/*.json")

# filter out empty files
emptyIdx = filesize.(infiles) .== 0
if any(emptyIdx)
    @warn "skipping empty files" infiles[emptyIdx]
    infiles = infiles[.!emptyIdx]
end

for file in infiles
    fname = splitext(basename(file))[1]
    println(fname)
    nodeCents, _ = read_PH2centralities("$file", "max")
    writedlm("$OUTDIR/$fname.tsv", nodeCents)
end

