#!/usr/bin/env julia
using ArgParse
parser = ArgParseSettings(autofix_names=true, description="""
Reads .npy or .tsv files with xyz coordinates and writes persistent homology 
barcodes and representatives to .json files in OUTDIR/.
""")
@add_arg_table parser begin
    "paths"
    arg_type = String
    nargs = '+'
    help = "INDIR/ OUTDIR/ or INFILES... OUTDIR/"
    "--dist", "-d"
    arg_type = Float64
    help = "Provide the typical distance between points (interpolation distance 
    choice) in order to improve performance. If provided, filtration time will 
    increment by the given value."
end

# if run as script
if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
else
    # if interactive, do an example TODO
    sArgs = ""
    args = split(sArgs, ' ')
end
args = parse_args(args, parser, as_symbols=true)
# for (k,v) in args println("# $k = $v") end
args = NamedTuple(args)

# load packages after potential -h/--help call.
using Eirene
using Distances
using CSV, DataFrames
using NPZ
using DelimitedFiles
using JSON
using .Threads: @threads

if length(ARGS) == 2 && isdir(ARGS[1])
    indir, outdir = ARGS
    fnames = indir * '/' * readdir(indir)
else
    outdir = ARGS[end]
    @assert isdir(outdir) "Last arg $outdir has to be a dir."
    fnames = ARGS[1:end-1]
end

npys = fnames[endswith.(fnames, ".npy")]
tsvs = fnames[endswith.(fnames, ".tsv")]
println([npys; tsvs])
if isempty(tsvs) && !isempty(npys)
    pointclouds = npzread.(npys)
elseif isempty(npys) && !isempty(tsvs)
    pointclouds = CSV.read.(tsvs, DataFrame; delim='\t')
    # assume there exist 3 columns named x,y,z
    pointclouds = [pc[:, ["x", "y", "z"]] for pc in pointclouds] .|> Matrix
else
    error("Both .npy and .tsv files given.")
end

"Eirene wants points along second dim so make sure #dim1 = 3"
function xyz_dim1(mat::AbstractMatrix)
    if size(mat, 1) == 3 && size(mat, 2) > 3
        return mat
    elseif size(mat, 2) == 3 && size(mat, 1) > 3
        return mat'
    else
        error("Matrix with size=$(size(mat)) not looking like xyz values.")
    end
end
pointclouds = xyz_dim1.(pointclouds)


function get_diameter(pointcloud::AbstractMatrix)
end


@threads for (pointcloud, fname) in collect(zip(pointclouds, fnames))
    fname = joinpath(outdir, splitext(fname)[1] * ".json")
    # same approach as in eirene source code
    d = pairwise(Euclidean(), pointcloud, dims=2)
    maxrad = maximum(d)
    numrad = args.dist === nothing ? Inf : ceil(Int, maxrad / args.dist)
    PH = eirene(d, maxdim=1, minrad=0, maxrad=maxrad, numrad=numrad)
    b = barcode(PH, dim=1)
    representatives = [classrep(PH, class=i, dim=1) for i in 1:size(b, 1)]
    
    open(fname, "w") do io
    	JSON.print(io, Dict("barcode" => b, "representatives" => representatives))
    end
end

