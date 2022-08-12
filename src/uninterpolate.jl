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
Take array in a table or json and uninterpolate values.
If vector of vectors of vectors, then it is assumed to be representatives from PH (use -e/--entry=representatives),
where an uninterpolated incidence matrix H will be written.
""")
@add_arg_table parser begin
    "infiles"
    required = true
    arg_type = String
    help = "Directory or glob pattern with tables (.tsv) or (.json)."
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
    help = "If jsons are given, specify which entry to uninterpolate.
    Default is each entry gets its own outfile, in which case only a single infile is supported."
    "--categorical", "-c"
    action = :store_true
    help = "Set flag to treat values as categorical, which will they will be 
    one-hot encoded before multiplying by D. Only implemented for integer vector infiles.
    The categories will be the sorted unique values found in each infile."
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
    if args.entry === nothing
        @assert length(infiles) == 1 "Using all json entries only supported for a single infile"
        d = JSON.parsefile(only(infiles))
        names = keys(d) # using json entry keys as outfile names instead of infile names
        arrays = values(d)
    else
        arrays = [j[args.entry] for j in JSON.parsefile.(infiles)]
    end
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

"""
eltype for when eltype entry isn't set (it will be Any) but it is assumed to be 
a single type.
"""
elementtype(x) = typeof.(x) |> unique |> only

if all(eltype.(arrays) .<: Real)
    # add zeros to end of node values. The last entries are missing because 
    # those interp points were never added to a feauture so no hyperedges 
    # contain them.
    arrays = PaddedView.(zero.(eltype.(arrays)), arrays, tuple.(nInterps, size.(arrays, 2)))
# For simplicity assume Vector{Vector} is always Vector{Vector{Vector}}
elseif all(elementtype.(arrays) .<: Vector)
    arrays = hyperedges2B.(representatives2hyperedges.(arrays), nInterps)
elseif all(elementtype.(arrays) .<: Real)
    # This scenario occurs for JSONs where all entries are vectors of numbers.
else
    uTypes = typeof.(arrays) |> unique
    error("Unknown array format(s): $uTypes")
end

if args.categorical
    @assert all(ndims.(arrays) .== 1) "Categorical only implemented for vectors."
    categories = unique.(arrays) .|> sort
    onehot(values::Vector, categories::Vector) = values .== categories'
    arrays = onehot.(arrays, categories)
end

mkpath(args.outdir)
outnames = args.outdir .* '/' .* names .* ext
writedlm.(outnames, Ds .* arrays, delim)

