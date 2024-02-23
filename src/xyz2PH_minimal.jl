#!/usr/bin/env julia
# USE: `git root`/src/xyz2PH_minimal.jl 'xyz*.tsv' OUTDIR/
include("compute_PH_minimal.jl")
using .Threads: @threads
include("glob.jl") # glob that handles abs path as well

infile_glob, outdir = ARGS

mkpath(outdir)

@threads for infile in glob(infiles_glob)
    outfile = splitext(basename(infile))[1] * ".json"
    PH_minimal.compute_PH_minimal_generators(infile, outfile)
end

