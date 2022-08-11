#!/usr/bin/env julia
using ArgParse
using JSON
using DelimitedFiles
using DataFrames, CSV
using PaddedViews
include("glob.jl")
include("uninterpolation.jl")
include("PH2hypergraph.jl") # representatives2hyperedges, and from hypergraph_utils.jl hyperedges2B

parser = ArgParseSettings(description="""
Take array and uninterpolate values, or given a json, provide an entry to uninterpolate.
E.g. representatives from PH to create an uninterpolated incidence matrix H.
""")
@add_arg_table parser begin
    "infiles"
    required = true
    arg_type = String
    help = "Directory or glob pattern with tables (.tsv) or PH (.json)."
    "outdir"
    required = true
    arg_type = String
    help = "Directory name for writing outputs. Will be created if nonexistent."
    "--interp", "-i"
    arg_type = String
    default = "interpolated/"
    help = "Directory with interpolated tables (.tsv) matching infiles naming. Columns: x, y, z, interp."
    "--entry", "-e"
    arg_type = String
    default = "representatives"
    help = "If jsons are given, specify which entry to uninterpolate."
end

args = parse_args(ARGS, parser, as_symbols=true) |> NamedTuple

infiles = glob(args.infiles)
if length(infiles) == 1 && isdir(infiles[1])
    indir, = infiles
    infiles = joinpath.(indir, readdir(indir))
end

noext(s::AbstractString) = splitext(s)[1]
getext(s::AbstractString) = splitext(s)[end]
names = noext.(basename.(infiles))
ext = getext.(infiles) |> unique
@assert length(ext) == 1 "All infiles should have the same extension: $ext"
ext = only(ext)

if ext == ".json"
    @assert args.entry !== nothing "For json, -e/--entry should be given."
    arrays = [j[args.entry] for j in JSON.parsefile.(infiles)]
    ext, delim = ".csv", ',' # for output
elseif ext == ".tsv"
    delim = '\t'
    arrays = readdlm.(infiles, delim)
elseif ext == ".csv"
    delim = ','
    arrays = readdlm.(infiles, delim)
else
    error("Unknown file format: $ext")
end

fnames_interp = joinpath.(args.interp, names .* ".tsv")
# create matrix D that converts between interpolated and uninterpolated indexes
dfs_interp = CSV.read.(fnames_interp, DataFrame)
nInterps = nrow.(dfs_interp)
Ds = uninterp_average_weight.(dfs_interp)

if all(eltype.(arrays) .<: Real)
    # add zeros to end of node values. The last entries are missing because 
    # those interp points were never added to a feauture so no hyperedges 
    # contain them.
    arrays = PaddedView.(zero.(eltype.(arrays)), arrays, tuple.(nInterps, size.(arrays, 2)))
elseif all(typeof.(arrays) .<: Vector) # typeof. instead of eltype since JSON will set eltype to Any
    arrays = hyperedges2B.(representatives2hyperedges.(arrays), nInterps)
else
    uTypes = typeof.(arrays) |> unique
    error("Unknown array format(s): $uTypes")
end

mkpath(args.outdir)
outnames = args.outdir .* '/' .* names .* ext
writedlm.(outnames, Ds .* arrays, delim)

