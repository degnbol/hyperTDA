#!/usr/bin/env julia
ROOT = readchomp(`git root`)
include("path_utils.jl") # glob that accepts abspath
using DelimitedFiles
using NPZ
using Distances
using ArgParse
parser = ArgParseSettings(autofix_names=true, description="""
Train a simple dense net on adjancecy matrices as a comparison to CNNs trained on hypergraphs.
""")
@add_arg_table parser begin
    "adjacencies"
    arg_type = String
    nargs = '+'
    help = "Adjacency matrices or xyz coordinates that will be converted to adjacencies.
    Format can be .csv, .tsv, or .npy 
    where .csv and .tsv can have a header (set -H/--header) with x,y,z[,interp].
    Multiple glob patterns for filenames. 
    Each glob pattern is a separate label. 
    Two glob patterns means boolean classifier, more means categorical."
    "--header", "-H"
    action = :store_true
    help =  "If .tsv or .csv, set flag to indicate files has headers."
    "--nans", "-N"
    arg_type = Int
    default = 0
    help = "Number of NaN entries to be allowed in each input file."
    "--pred", "-p"
    arg_type = String
    help = "Write predictions to given file."
    "--params", "-P"
    arg_type = String
    help = "Write params to given file. Reshaped to match shape of adjacencies, with the bias written in entry (1,1)."
    "--epochs", "-e"
    arg_type = Int
    default = 1000
    help = "Number of epochs of training."
    "--early-stop", "-E"
    arg_type = Int
    default = 200
    help = "If the best testset AUC was this many epochs ago, then stop training. Set to zero to disable."
    "--nodes", "-n"
    arg_type = Int
    default = 1
    help = "Number of nodes in second dense layer. A single layer model is created if -n/--nodes < #labels."
    "--cv"
    action = :store_true
    help = "Perform 5-fold cross-validation instead of a single train-test-eval run for more stable scoring."
end

# if run as script
if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
else
    # if interactive, do an example
    # fnames = ["$ROOT/Data/IMR90_chr21-28-30Mb_cellCycle/interpolated/*$cycle.tsv" for cycle in ["G1", "G2", "S"]]
    # fnames = ["$ROOT/Data/HCT116_chr21-34-37Mb_$treat/interpolated/*" for treat in ["untreated", "6hAuxin"]]
    fnames = ["$ROOT/Data/HCT116_chr21-28-30Mb_$treat/interpolated/*" for treat in ["untreated", "6hAuxin"]]
    args = [fnames; split("-H -N 4 --params lars.tsv --cv", ' ')]
end
args = parse_args(args, parser, as_symbols=true)

for (k,v) in args println("# $k = $v") end
args = NamedTuple(args)

fnames = []
labels = Int[]
for (lab, pattern) in enumerate(args.adjacencies)
    fs = glob.(pattern)
    println("# $(length(fs)) files found for pattern $pattern")
    append!(fnames, fs)
    append!(labels, [lab for _ ∈ 1:length(fs)])
end

# read

if all(endswith.(fnames, ".npy"))
    mats = npzread.(fnames)
else
    if all(endswith.(fnames, ".tsv"))
        mats = readdlm.(fnames, '\t'; header=args.header)
    elseif all(endswith.(fnames, ".csv"))
        mats = readdlm.(fnames, ','; header=args.header)
    else
        error("All filenames should end with either .npy, .tsv, or .csv.")
    end
    if args.header
        # columns should be named "x", "y", "z" and optionally "interp" where 
        # we remove interpolated points.
        dfs = mats
        mats = []
        for (mat, header) in dfs
            header = vec(header)
            interpInd = findfirst(header .== "interp")
            rowIdx = interpInd === nothing ? trues(size(mat,1)) : mat[:, interpInd] .== 0
            xyzInd = findfirst.(header .== l for l in split("xyz", ""))
            mat = mat[:, xyzInd][rowIdx, :]
            push!(mats, mat)
        end
    end
end


nNaNs = [sum(any(isnan.(mat); dims=2)) for mat in mats]
nanKeep = nNaNs .<= args.nans
mats = mats[nanKeep]
labels = labels[nanKeep]
println("# File count reduction due to #NaNs > $(args.nans): $(length(nanKeep)) -> $(sum(nanKeep))")


# interpolate NaNs
if any(any(isnan.(mat)) for mat in mats)
    error("Interpolation not implemented yet, try running on interp/ files instead for now.")
end

# assume adj if square. Could also check with issymmetric from LinearAlgebra 
# but I worry about rounding error.
issquare(mat::Matrix) = size(mat, 1) == size(mat, 2)
if all(issquare.(mats))
    println("# File type detected: adj")
    # adjs = Symmetric(mats)
    adjs = mats
else
    println("# File type detected: xyz")
    adjs = pairwise.(Ref(Euclidean()), mats; dims=1)
end

# wait till after possible -h/--help call and file assertions
using Statistics: mean
include("AUC.jl")
include("flux_utils.jl")

# boolean classifier, or categorical?
uLabels = labels |> unique |> sort
for uLabel in uLabels
    println("# num label $uLabel = ", sum(labels .== uLabel))
end
nLabels = length(uLabels)
if nLabels == 2
    labels = labels .== uLabels[2]
elseif nLabels > 2
    labels = onehot(labels)
else
    error("Multiple patterns (labels) must be provided")
end

# Different number of nodes not implemented yet.
nNodes = unique(size(adj, 1) for adj in adjs)
if length(nNodes) > 1
    error("Different numbers of nodes not implemented yet")
end
nNodes = only(nNodes)

function new_model()
    # input is a symmetric matrix with zeros in the diagonal so we only need 
    # less than half of the entries.
    nIn = (nNodes ^ 2 - nNodes) / 2 |> Int
    # binary classifier?
    if nLabels == 2
        if args.nodes < nLabels
            Chain(Dense(nIn => 1), vec)
        else
            Chain(
                Dense(nIn => args.nodes, relu),
                Dense(args.nodes => 1),
                vec)
        end
    else
        if args.nodes < nLabels
            Dense(nIn => nLabels)
        else
            Chain(
                Dense(nIn => args.nodes, relu),
                Dense(args.nodes => nLabels))
        end
    end
end

"Get entries from lower triangle of matrix, excluding the diagonal."
lowTri(mat::Matrix) = [mat[i,j] for i ∈ 1:size(mat,1) for j ∈ 1:i-1]
"Inverse of function lowTri (except only returns the lower triangle filled)."
fromLowTri(v::Vector) = begin
    # solution to (nNodes ^ 2 - nNodes) / 2 
    n = (√(8*length(v) + 1) + 1) / 2 |> Int
    out = zeros(n, n)
    ii = 1
    for i ∈ 1:nNodes
        for j ∈ 1:i-1
            out[i,j] = v[ii]
            ii += 1
        end
    end
    out
end

features = hcat(lowTri.(adjs)...) .|> Float32


# logit versions used instead of binCE(σ) and CE(softmax) for numerical 
# stability: 
# https://fluxml.ai/Flux.jl/stable/models/losses/#Flux.Losses.logitbinarycrossentropy
if nLabels == 2
    predict(set, model) = σ.(model(features[:, set]))
    CE(set, model) = Flux.logitbinarycrossentropy(model(features[:, set]), labels[set])
    auc(set, model) = AUC(predict(set, model), labels[set])
    acc(set, model) = mean((predict(set, model) .> 0.5) .== labels[set])
else
    predict(set, model) = softmax(model(features[:,set]))
    CE(set, model) = Flux.logitcrossentropy(model(features[:,set]), labels[:, set])
    auc(set, model) = begin
        # just do average AUC for each class vs all other classes
        preds = predict(set, model)
        aucs = [AUC(preds[lab, :], labels[lab, set] .== 1) for lab in 1:nLabels]
        mean(aucs)
    end
    acc(set, model) = begin
        preds = predict(set, model)
        # using argmax means we assume labels are some range 1:n
        mean(view(labels, :, set)[onehot(argmax.(eachcol(preds)), uLabels)])
    end
end

opt = ADAMW()
λ = 1.
function train(set, params, model)
    grads = Flux.gradient(params) do
        CE(set, model) + λ * sum(sum(abs.(p)) for p in params)
    end
    Flux.Optimise.update!(opt, params, grads)
end

bold(v) = "\033[1m$v\033[0m"

evalAUCs    = []
evalAccs    = []
preds       = []
pred_labels = []
best_params = []

# use function for the return statement
function training!(cv_parts=5)
    cv_sets = CV(labels, cv_parts)
    for iTest ∈ 1:cv_parts
        for iEval ∈ 1:cv_parts
            iTest != iEval || continue
            println("# test=$(iTest), eval=$(iEval)")
            testset  = cv_sets .== iTest
            evalset  = cv_sets .== iEval
            trainset = .!(testset .| evalset)
            
            model = new_model()
            params = Flux.params(model)
            
            println("Epoch\tCE\tAcc\tAUC\tTestAUC")
            push!(best_params, deepcopy(collect(params)))
            best_auc = auc(testset, model)
            best_epoch = 0
            println(join((0, CE(trainset, model), acc(trainset, model), auc(trainset, model), best_auc), '\t'))
            for epoch in 1:args.epochs
                for batchset in get_batches(trainset, 20)
                    train(batchset, params, model)
                end
                if epoch % 10 == 0
                    auc_test = auc(testset, model)
                    if auc_test > best_auc
                        best_auc = auc_test
                        best_epoch = epoch
                        best_params[end] = params |> collect |> deepcopy
                        auc_test = bold(auc_test) # bold if best
                        best_auc ≈ 1. && break
                    elseif args.early_stop > 0 && epoch - best_epoch > args.early_stop
                        break # early stop if it's been too long
                    end
                    # use join on tuple instead of list since list will have a 
                    # type, and make epoch int->float sometimes.
                    println(join((epoch, CE(trainset, model), acc(trainset, model), auc(trainset, model), auc_test), '\t'))
                end
            end
            
            Flux.loadparams!(model, best_params[end])
            push!(evalAUCs, auc(evalset, model))
            push!(evalAccs, acc(evalset, model))
            println("eval AUC = ", evalAUCs[end])
            println("eval Acc = ", evalAccs[end])
            
            if args.pred !== nothing
                if nLabels == 2
                    append!(preds, predict(evalset, model))
                    append!(pred_labels, labels[evalset])
                else
                    # using argmax means we assume labels are some range 1:?
                    append!(preds, argmax.(predict(evalset, model)))
                    append!(pred_labels, unonehot(labels[:, evalset]))
                end
            end
            
            args.cv || return
        end
    end
end

training!()

if args.cv
    println("avg eval AUC = ", mean(evalAUCs))
    println("avg eval Acc = ", mean(evalAccs))
    inParamCor = [cor(vec(best_params[i][1]), vec(best_params[j][1]))
         for i in 1:length(best_params)
         for j in 1:i-1]
    println("avg cor among in-layer params = ", mean(inParamCor))
end

if args.pred !== nothing
    delim = '\t'
    open(args.pred; write=true) do f
        write(f, join(["Pred", "Label"], delim) * '\n')
        writedlm(f, [preds pred_labels])
    end
end


function reshape_params(params::Vector)
    weights, (bias,) = params
    out = fromLowTri(vec(weights))
    out[1,1] = bias
    out
end

"writedlm, but zeros are written as empty string."
function writedlm0(f, A::Matrix, delim='\t')
    open(f; write=true) do fh
        for row in eachrow(A)
            row = Vector{Any}(row)
            row[row .== 0] .= ""
            write(fh, join(row, delim) * '\n')
        end
    end
end

if args.params !== nothing
    delim = '\t'
    if length(best_params) == 1
        writedlm0(args.params, reshape_params(only(best_params)))
    else
        outfname, ext = splitext(args.params)
        for (i, params) in enumerate(best_params)
            writedlm0("$outfname-CV$i$ext", reshape_params(params))
        end
    end
end


