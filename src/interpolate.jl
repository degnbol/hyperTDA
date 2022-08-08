#!/usr/bin/env julia
using ArgParse
using DelimitedFiles
using NPZ
include("glob.jl") # glob that accepts abspath

parser = ArgParseSettings(description="""
Take xyz curve and interpolate linearly between points, given a fixed number of 
required output points.
""")

@add_arg_table parser begin
    "infiles"
    arg_type = String
    nargs = '+'
    help = "Infiles, indir, or glob pattern for infiles. 
    .npy, .csv, or .tsv without headers."
    "--out", "-o"
    arg_type = String
    help = "Output directory. 
    Files are written with the same names as infiles except extension will be .tsv. 
    If it doesn't exist it will be created."
    "--points", "-n"
    arg_type = Int
    help = "Number of points in each output file."
    "--nans", "-N"
    arg_type = Int
    default = 0
    help = "Max number of allowed NaNs in a curve. Curves with more are skipped. 
    NaNs are replaced by points interpolated between the point before and after.
    Curves with NaN in first or last point are skipped."
end

# if run as script
if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
else
    # if interactive, do an example. You may need to unzip pointClouds.zip first.
    ROOT = readchomp(`git root`)
    sArgs = "$ROOT/Data/HCT116_chr21-28-30Mb_untreated/pointClouds -o lars -n 100 -N 1"
    args = split(sArgs, ' ')
end
args = parse_args(args, parser, as_symbols=true)
for (k,v) in args println("# $k = $v") end
args = NamedTuple(args)

# if indir given
infiles = glob(args.infiles)
if isdir(infiles) infiles = glob(infiles[1] .* "/*") end

# a consistent fileformat should be given
extMatches = [endswith.(infiles, ext) for ext in [".npy", ".tsv", ".csv"]] .|> any |> sum
@assert extMatches == 1 "$extMatches extension matches. There should be exactly 1 match (.npy, .tsv or .csv)."

# read
npzs = infiles[endswith.(infiles, ".npy")]
tsvs = infiles[endswith.(infiles, ".tsv")]
csvs = infiles[endswith.(infiles, ".csv")]
names = basename.([npzs; tsvs; csvs])
curves = [npzread.(npzs); readdlm.(tsvs, '\t'); readdlm.(csvs, ',')]

@assert all(size.(curves, 2) .== 3) "All files should have 3 columns."

nIn = maximum(size.(curves, 1))
@assert nIn .<= args.points "-n/--points=$(args.points) should be more than curve lengths=$nIn."

# NaN filters
anynan(v::AbstractArray) = any(isnan.(v))
anynan(v::AbstractArray, dims) = any(isnan.(v); dims=dims)
nnan(curve::AbstractMatrix) = sum(anynan(curve, 2))
hasTerminalNan(curve::AbstractMatrix) = anynan(curve[[1,end], :])

nanFilt = nnan.(curves) .<= args.nans
curves = curves[nanFilt]
names = names[nanFilt]
@info "NaNs > $(args.nans) filter: $(length(nanFilt)) -> $(length(curves))"

terminalNan = hasTerminalNan.(curves)
curves = curves[.!terminalNan]
names = names[.!terminalNan]
@info "terminal NaN filter: $(length(terminalNan)) -> $(length(curves))"

"""
Get point linearly interpolated between point at index "start" and "stop" with 
index "between" where between=start will return point start and between=stop 
will return point stop.
"""
function interp(curve::AbstractMatrix, start::Int, stop::Int, between::Int)
    interp(curve, start, stop, (between - start) / (stop - start))
end
function interp(curve::AbstractMatrix, start::Int, stop::Int, frac::Float64)
    (1-frac) * curve[start,:] + frac * curve[stop,:]
end

"""
Replace points containing NaN with linear interpolation of neighbors.
curve: xyz
"""
function interpNan(curve::Matrix)
    anynan(curve[[1,end], :]) && error("Not implemented.")
    
    n = size(curve, 1)
    
    out = []
    for (i, point) ∈ enumerate(eachrow(curve))
        if !anynan(point)
            push!(out, point)
        else
            # find nearest non-nan neighbors
            start = (1:i-1)[vec(.!anynan(curve[1:i-1, :], 2))] |> last
            stop  = (i+1:n)[vec(.!anynan(curve[i+1:n, :], 2))] |> first
            push!(out, interp(curve, start, stop, i))
        end
    end
    
    hcat(out...)'
end 

if any(anynan.(curves))
    nanRows = anynan.(curves, 2)
    curves = interpNan.(curves)
else
    nanRows = nothing
end

"Euclidean distance between each consecutive point (row) in a curve."
function consecutiveEuclidean(curve::AbstractMatrix)
    sqrt.(sum(diff(curve; dims=1) .^ 2; dims=2))
end

"""
Interpolate between points in a curve of xyzs resulting in n point.
"""
function interp(curve::AbstractMatrix, nOut::Int)
    nIn = size(curve, 1)
    nNew = nOut - nIn
    
    # decide how many points to add between each real point based on the 
    # distances between points
    dists = consecutiveEuclidean(curve)
    nBetween = zeros(Int, nIn - 1)
    for i in 1:nNew
        gaps = dists ./ (nBetween .+ 1) |> vec
        nBetween[argmax(gaps)] += 1
    end
    
    # indicator vector of interpolation
    interpIdx = trues(nOut)
    interpIdx[cumsum(1 .+ [0; nBetween])] .= false
    
    # add the real and interpolated points iteratively
    out = []
    for i in 1:nIn-1
        push!(out, curve[i, :])
        nb = nBetween[i]
        for ii in 1:nb
            push!(out, interp(curve, i, i+1, ii/(nb+1)))
        end
    end
    push!(out, curve[end, :])
    
    hcat(out...)', interpIdx
end

interps = interp.(curves, args.points)
# nanRows are only for uninterped points, update the indices to the 
# interpolated output size.
if nanRows !== nothing
    tmp = nanRows
    nanRows = [falses(args.points) for _ in 1:length(interps)]
    for (i, (_, interpIdx)) in enumerate(interps)
        nanRows[i][.!interpIdx] .= tmp[i]
    end
end

# write
delim, ext = '\t', ".tsv"
header = ["x", "y", "z", "interp"]
nanRows === nothing || push!(header, "nan")
header = join(header, delim) * '\n'

ispath(args.out) || mkdir(args.out)

for i in 1:length(interps)
    interped, interpIdx = interps[i]
    open(args.out * "/" * splitext(names[i])[1] * ext; write=true) do f
        write(f, header)
        mat = nanRows === nothing ? Real[interped interpIdx] : Real[interped interpIdx nanRows[i]]
        writedlm(f, mat)
    end
end

