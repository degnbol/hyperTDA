#!/usr/bin/env julia
using MKL # potential speedup over openBLAS
using CUDA # use GPU if available
include("glob.jl") # glob that accepts abspath
include("string_utils.jl") # prefixSuffixPairs
using DelimitedFiles
using Random # shuffle
using BSON: @load, @save
using ArgParse
using Dates

parser = ArgParseSettings(autofix_names=true, description="""
Take incidence matrix, e.g. 1s and 0s H or uninterpolated incidence matrix 
with partial memberships. Run conv filters along hyperedges E, with conv 
size=1 since E order is arbitrary. The other dim should be the max number of 
nodes found among the hypergraphs. The hypergraphs with less than that number 
of nodes are either a truncated set of nodes, where the missing nodes exists 
for other curves that have trajectories going further, or a smaller subset of 
interpolated points. We initially implement assuming truncated number of 
points. The missing points are simply zero, i.e. they never appear in any of 
the hyperedges.
""")
@add_arg_table parser begin
"--hypergraphs", "-H"
required = true
arg_type = String
nargs = '+'
help = "Multiple glob patterns for filenames or directories. 
Files should be comma delimited without header (or row index). 
Incidence matrix for hypergraphs, i.e. nodes along first dim, hyperedges along second. 
Each glob pattern is a separate label. 
Two glob patterns means boolean classifier, more means categorical."
"--edge-values", "--weights", "-W"
arg_type = String
nargs = '+'
help = "Multiple glob patterns for filenames or directories. 
Files are each simply a column vector, so a single scalar on each line, no header or row index.
Hyperedge weight columns for hypergraphs. 
Same number of entries as --hypergraphs/-H. 
Alternatively, supply a find and replace string, 
which when used on a -H filename will yield the equivalent -W filename."
"--node-values", "-V"
arg_type = String
nargs = '+'
help = "Same as --edge-values/-W except length should match first dim of --hypergraphs/-H rather than second. 
Multiple node values are also supported, in which case use commas two separate them on each line.
Currently only used by model 6."
"--weigh-edges", "--we"
action = :store_true
help = "Preprocess each H by weighing them (multiplying) by edge values given by --edge-values/-W."
"--weigh-nodes", "--wn"
action = :store_true
help = "Preprocess each H by weighing them (multiplying) by node values given by --node-values/-V."
"--pred", "-p"
arg_type = String
help = "Optionally write predictions to a given filename."
"--load-model"
arg_type = String
help = "Optionally load model(s) to continue training (.bson)."
"--save-model", "-M"
arg_type = String
help = "Optionally write model(s) to a given filename (.bson)."
"--epochs", "-e"
arg_type = Int
default = 500
help = "Number of epochs of training."
"--early-stop", "-E"
arg_type = Int
default = 100
help = "If the best testset AUC was this many epochs ago, then stop training. Set to zero to disable."
"--model", "-m"
arg_type = Int
default = 8
# We should use both H and W. We can feed them as parallel features but the 
# filter only does a linear combination of features, so we would have to give 
# each feature multiple output values, that can then be combined in a 
# subsequent dense layer. Alternatively, we can multiply H and W but then we 
# assume a linear and equal importance of weight for each hyperedge hidden 
# feature.
help = "Model 1 is a single conv layer, which is simply used as a trick to allow variable number of hyperedge input. 
Model 2 has two conv layers. 
Model 3 has two but the edge values are only fed to the second. This should be a meaningful network for finding global features.
Model 4 is a baseline dummy model, an identity layer instead of the conv layers. 
Model 5 uses actual convolutions, scanning along nodes, assumed to be ordered. 
It finds features locally rather than globally, although it will not include all entries of a hyperedge so may only be meaningful in specific cases.
Model 6 is like model 3 except node values are utilized with convolutions then pooling before concatenating with the second layer output from H and W. 
Note that this means there is no connection between individual node values and which hyperedges that specific node is a part of. 
Model 7 runs convolutions on both H and V and combines the pooled conv features before concatenating W. 
Model 8 like model 7 runs convolutions on both H and V but then concatenates them without pooling nodes, but pooling edges instead. This means W is NOT utilized."
"--nFilters", "-f"
arg_type = Int
default = 8
help = "Number of filters in first conv layer."
"--nFilters2", "-F"
arg_type = Int
default = 4
help = "Number of filters in second conv layer, if applicable."
"--kernel-sizes", "-k"
arg_type = Int
nargs = '+'
default = [5, 15]
help = "Sizes of kernels for model 5 and 6."
"--cv"
action = :store_true
help = "Perform 5-fold cross-validation instead of a single train-test-eval run for more stable scoring."
"--pre", "-P"
arg_type = String
help = "Preprocessesing action to perform on each H. \"nzalign\" to align first nonzero entries for each hyperedge. 
\"nzranges\" to generate hyperedges for all ranges (1:end, 1:end-1, ..., end-1:end) which contains nonzero elements and doesn't exclude nonzero elements."
"--sample", "-s"
arg_type = Int
range_tester = x -> x > 0
help = "Randomly sample no more than this number of datapoints (graphs) from each label."
"--pair-files"
action = :store_true
help = "Relax expectation that files are sorted and paired, and find the 
filename pairs between H, W, and V when there may be inconcsistent file count 
and order."
end

# if run as script
if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
else
    # if interactive, do an example
    ROOT = readchomp(`git root`)
    fnames_H = ["$ROOT/results/AnDi/Model_$i/H" for i in 0:4]
    fnames_V = ["$ROOT/results/AnDi/Model_$i/nodeCents" for i in 0:4]
    fnames_H = join(fnames_H, ' ')
    fnames_V = join(fnames_V, ' ')
    args = split("-m 8 -H $fnames_H -V $fnames_V", ' ')
end
args = parse_args(args, parser, as_symbols=true)

if length(args[:edge_values]) == 2 && all(isempty.(glob.(args[:edge_values])))
    # seems a search and replace string has been provided
    find, sub = args[:edge_values]
    args[:edge_values] = [replace(hPattern, find => sub) for hPattern ∈ args[:hypergraphs]]
    @assert length(args[:edge_values]) > 0
end

if length(args[:node_values]) == 2 && all(isempty.(glob.(args[:node_values])))
    # seems a search and replace string has been provided
    find, sub = args[:node_values]
    args[:node_values] = [replace(hPattern, find => sub) for hPattern ∈ args[:hypergraphs]]
    @assert length(args[:node_values]) > 0
end

edgeVals = length(args[:edge_values]) > 0
nodeVals = length(args[:node_values]) > 0

for (k,v) in args println("# $k = $v") end
args = NamedTuple(args)

# assert equal number of labels
if edgeVals
    @assert(   length(args.hypergraphs)  ==   length(args.edge_values),
            "$(length(args.hypergraphs)) != $(length(args.edge_values))")
end
if nodeVals
    @assert(   length(args.hypergraphs)  ==   length(args.node_values),
            "$(length(args.hypergraphs)) != $(length(args.node_values))")
end
@assert args.pre ∈ [nothing, "nzalign", "nzranges"] "$(args.pre) ∉ [nzalign, nzranges]"

"Function to check if a directory has been provided, and if so return all files in it."
function dir2files(paths::AbstractVector)
    length(paths) == 1 || return paths
    path = only(paths)
    isdir(path) || return paths
    # this glob will not include hidden files, as opposed to a readdir call.
    glob(path * "/*")
end

# read
Hs, Ws, Vs, labels = Matrix[], Matrix[], Matrix[], Int[] 
for (label, glob_H) ∈ enumerate(args.hypergraphs)
    fnames_H = eglob(glob_H) |> dir2files
    nLabel = length(fnames_H)
    
    if edgeVals
        glob_W = args.edge_values[label]
        fnames_W = eglob(glob_W) |> dir2files
        if args.pair_files
            nFnames_W = length(fnames_W)
            fnames_H, fnames_W, = prefixSuffixPairs(fnames_H, fnames_W)
            println("# Pairing up files.")
            println("H: $nLabel -> $(length(fnames_H))")
            println("W: $nFnames_W -> $(length(fnames_W))")
            nLabel = length(fnames_H)
        end
        @assert(length(fnames_H) == length(fnames_W),
                "$(length(fnames_H)) $(glob_H) != $(length(fnames_W)) $(glob_W)")
    end
    if nodeVals
        glob_V = args.node_values[label]
        fnames_V = eglob(glob_V) |> dir2files
        if args.pair_files
            nFnames_V = length(fnames_V)
            fnames_H, fnames_V, fnames_H_idx, = prefixSuffixPairs(fnames_H, fnames_V)
            println("# Pairing up files.")
            println("H: $nLabel -> $(length(fnames_H))")
            println("V: $nFnames_V -> $(length(fnames_V))")
            if edgeVals
                fnames_W = fnames_W[fnames_H_idx] 
                println("W: $nLabel -> $(length(fnames_W))")
            end
            nLabel = length(fnames_H)
        end
        @assert(length(fnames_H) == length(fnames_V),
                "$(length(fnames_H)) $(glob_H) != $(length(fnames_V)) $(glob_V)")
    end
    
    println("# Reading label $label from $nLabel hypergraphs matching:")
    println("# $glob_H")
    edgeVals && println("# $glob_W")
    nodeVals && println("# $glob_V")
    
    if args.sample !== nothing
        lSample = min(args.sample, nLabel)
        println("# Sampling $lSample files")
        sample = shuffle(1:nLabel)[1:lSample]
        fnames_H = fnames_H[sample]
        if edgeVals fnames_W = fnames_W[sample] end
        if nodeVals fnames_V = fnames_V[sample] end
        nLabel = lSample
    end
    
    _Hs = readdlm.(fnames_H, ',')
    # If they are not given we can represent them in many ways. To avoid a lot 
    # of if clauses and to have the ML work and to allow for more than a single 
    # feature it is nice to then also allow for 0-dim.
    if edgeVals
        _Ws = readdlm.(fnames_W)
        for (H, W) in zip(_Hs, _Ws)
            @assert(       size(H,2)  ==   size(W,1),
            "#hyperedges $(size(H,2)) != $(size(W,1))")
        end
        # match direction of H
        _Ws = [W' for W ∈ _Ws]
    else
        _Ws = [ones(0, size(H,2)) for H ∈ _Hs]
    end
    if nodeVals
        _Vs = readdlm.(fnames_V)
        for (H, V) in zip(_Hs, _Vs)
            @assert(  size(H,1)  ==   size(V,1),
            "#nodes $(size(H,1)) != $(size(V,1))")
        end
    else
        _Vs = [ones(size(H,1), 0) for H ∈ _Hs]
    end
    append!(Hs, _Hs)
    append!(Ws, _Ws)
    append!(Vs, _Vs)
    append!(labels, fill(label, nLabel))
end

# wait till after possible -h/--help call and file assertions
using Statistics # mean, cor
using Printf
include("flux_utils.jl")
include("array_utils.jl") # add_dim
include("hypergraph_utils.jl") # for args.pre
include("AUC.jl")
# not implemented for GPU
AUC(ŷ::CuArray, y::CuArray) = AUC(cpu(ŷ), cpu(y))

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

# number of edge and node values given
nW = unique(size(W,1) for W ∈ Ws) |> only
nV = unique(size(V,2) for V ∈ Vs) |> only

# enforce same number of nodes by assuming truncation in cases with fewer 
# nodes.
nNodes = [size(H, 1) for H in Hs]
nNodesMax = maximum(nNodes)
nZeroNodes = nNodesMax .- nNodes |> unique |> sort
if nZeroNodes != [0]
    println("# Appending zero nodes: ", nZeroNodes)
    Hs = [vcat(H, zeros(nNodesMax-size(H,1), size(H,2))) for H in Hs]
end

# Preprocessesing
# for column vectors V and W the weightings will be meaningful, but if multiple 
# node values are given for each node then the first node value is the 
# weighting used here.
if args.weigh_nodes
    Hs = [H .* V[:,1] for (H, V) ∈ zip(Hs, Vs)]
end
if args.weigh_edges
    Hs = [H .* W[1,:]' for (H, W) ∈ zip(Hs, Ws)]
end
if args.pre !== nothing
    if args.pre == "nzalign"
        Hs = nzalign.(Hs, nNodesMax)
    elseif args.pre == "nzranges"
        Hs_ranges = []
        Ws_ranges = []
        for (i, (H, W)) ∈ enumerate(zip(Hs, Ws))
            # print progress is nice if interactive session
            # print("# nzranges: $i/$(length(Hs))\r")
            _H, _W = nzranges(H, W, nNodesMax)
            push!(Hs_ranges, _H)
            push!(Ws_ranges, _W)
        end
        Hs = Hs_ranges
        Ws = Ws_ranges
    end
end

# Flux uses Float32
# and we consistently use hyperedges as first dim, so nodes can sometimes be 
# "channels" for the ?Conv.
Hs = [Float32.(H') for H ∈ Hs]
Ws = [Float32.(W') for W ∈ Ws]
Vs = [Float32.(V') for V ∈ Vs]

# add singleton dim since we can't provide multiple batches due to varying size 
# since hyperedges vary in count.
# We would want to concat all features along the last dim which is for batches. 
# However, since they vary in length we cannot do that and instead uses 
# multiple calls with batchsize=1, hence the add_dim.
# The convolution only has size 1, since hyperedge order is meaningless.
# We look for #nFilters hidden features to be detected by separate filters.
convTriv(nIn) = Chain(add_dim, Conv((1,), nIn=>args.nFilters, relu))

convEdges() = begin
    convEdge(kernelSize) = begin
        Chain(
            Conv((1, kernelSize), 1=>args.nFilters, relu, pad=SamePad()),
            # we can't pool hEdge dim yet since we need to combine with W,
            # and we pool nodes since we would otherwise have to complicate the 
            # next conv in order to not conv across W concated to H,
            # and since it should possible to capture similar patterns with a 
            # single real conv layer, if there are enough filters.
            MaxPool((1, nNodesMax)), 
            x -> dropdims(x; dims=2)
        )
    end
    Chain(add_dim, Parallel(hcat, (convEdge(k) for k ∈ args.kernel_sizes)...))
end

convNodes() = begin
    convNode(kernelSize) = begin
        Chain(
            Conv((1, kernelSize), 1=>args.nFilters, relu, pad=SamePad()),
            AdaptiveMaxPool((1, nNodesMax)), 
            x -> dropdims(x; dims=1)
        )
    end
    Chain(add_dim, Parallel(hcat, (convNode(k) for k ∈ args.kernel_sizes)...))
end


if args.model == 1
    features = [hcat(W, H) for (H, W) in zip(Hs, Ws)]
    nIn = nNodesMax+nW
    convLayers() = convTriv(nIn)
    nDenseIn = args.nFilters
elseif args.model == 2
    features = [hcat(W, H) for (H, W) in zip(Hs, Ws)]
    nIn = nNodesMax+nW
    convLayers() = Chain(
        convTriv(nIn),
        # Weight is combined linearly with H before relu in model1. We can add 
        # another Conv layer to combine them in a more nonlinear way.
        Conv((1,), args.nFilters=>args.nFilters2, relu)
    )
    nDenseIn = args.nFilters2
elseif args.model == 3
    # "channel" is nNodesMax
    features = zip(Ws, Hs) |> collect
    # same as above except we are very deliberate about H and W being different 
    # things so instead of concat from start, we first find features in H, then 
    # combine those hidden features with their weight.
    # This means it is fed W,H instead of vcat(W,H')
    convLayers() = Chain(
        Parallel(hcat,
            identity, # pass W through unmodified
            convTriv(nNodesMax), # find linear simple features in H
        ),
        # Combine W and H in a conv layer.
        # using same number of features for simplicity.
        Conv((1,), args.nFilters+nW=>args.nFilters2, relu)
    )
    nDenseIn = args.nFilters2
elseif args.model == 4
    features = [hcat(W, H) for (H, W) in zip(Hs, Ws)]
    nIn = nNodesMax+nW
    convLayers() = add_dim
elseif args.model == 5
    # will become 4 dims so 1 channel and H is a 2dim "image".
    features = [(W, add_dim(H)) for (H, W) in zip(Hs, Ws)]
    convLayers() = begin
        nConvIn2 = nW + length(args.kernel_sizes)*args.nFilters
        Chain(
              Parallel(hcat, identity, convEdges()),
              # we need to do pooling across first dim (hEdges) to deal with the 
              # varying dim size. If we do it now we would just get max weight 
              # and max of each matched pattern. We do the second fake conv layer 
              # to look for simple patterns in each hyperedge where there is 
              # still a connection between individual weights and the hEdge.
              Conv((1,), nConvIn2=>args.nFilters2, relu)
        )
    end
    nDenseIn = args.nFilters2
elseif args.model == 6
    # H,W tuple so they are fed together and separate from V that will become 4 
    # dims so 1 channel, it is a 2dim "image".
    # We run convolutions on node values, then pool, so there is no connnection 
    # between node values and H.
    features = [((H, W), add_dim(V)) for (H, W, V) in zip(Hs, Ws, Vs)]
    convLayers() = begin
        Parallel(hcat,
            Chain(
                Parallel(hcat, convTriv(nNodesMax), identity),
                Conv((1,), args.nFilters+nW=>args.nFilters2, relu),
                GlobalMaxPool(),
            ),
            convEdges()
        )
    end
    nDenseIn = args.nFilters2 + length(args.kernel_sizes)*args.nFilters
elseif args.model == 7
    # append zero(s) to start of W so dim will match when V is concated onto H 
    features = [([zeros(Float32,nV,nW); W], (add_dim(V), add_dim(H))) for (H, W, V) in zip(Hs, Ws, Vs)]
    convLayers() = begin
        Chain(
            Parallel(hcat,
                 identity, # W
                 Parallel(vcat,
                          convEdges(), # V
                          convEdges()  # H
                         )
                ),
            Conv((1,), args.nFilters*2+nW=>args.nFilters2, relu))
    end
    nDenseIn = args.nFilters2
elseif args.model == 8
    # Doesn't use W at all, can be used with Preprocessing weighing of H.
    features = [(add_dim(H), add_dim(V)) for (H, V) in zip(Hs, Vs)]
    # Same as model 7 but run conv on V and H WITHOUT pooling so that we are 
    # concating by node rather than by conv feature.
    convLayers() = Chain(Parallel(hcat, convNodes(), convNodes()),
        # this second conv looks at each node without context of surrounding 
        # nodes, in order to simply combine features coming from H and V.
        Conv((1,), 2*length(args.kernel_sizes)*args.nFilters=>args.nFilters2, relu))
    nDenseIn = args.nFilters2
else
    error("Unknown model option: $(args.model)")
end

function denseLayer(nLabels, nDenseIn)
    if nLabels == 2
        Chain(Dense(nDenseIn => 1), only)
    else
        Chain(Dense(nDenseIn => nLabels), vec)
    end
end

new_model() = Chain(
    convLayers(),
    GlobalMaxPool(),
    x -> dropdims(x; dims=1),
    denseLayer(nLabels, nDenseIn)
) |> gpu
# convert to CUDA if avail.
# Both parameters when generating a model, and the feature data.
features = gpu(features)
labels = gpu(labels)

# logit versions used instead of binCE(σ) and CE(softmax) for numerical 
# stability: 
# https://fluxml.ai/Flux.jl/stable/models/losses/#Flux.Losses.logitbinarycrossentropy
if nLabels == 2
    predict(set, model) = σ.(model.(features[set]))
    CE(set, model) = Flux.logitbinarycrossentropy(model.(features[set]), labels[set])
    auc(set, model) = AUC(predict(set, model), labels[set])
    acc(set, model) = mean((predict(set, model) .> 0.5) .== labels[set])
else
    predict(set, model) = softmax.(model.(features[set]))
    CE(set, model) = Flux.logitcrossentropy(hcat(model.(features[set])...), labels[:, set])
    auc(set, model) = begin
        # just do average AUC for each class vs all other classes
        preds = hcat(predict(set, model)...)
        aucs = [AUC(preds[lab, :], labels[lab, set] .== 1) for lab in 1:nLabels]
        mean(aucs)
    end
    acc(set, model) = begin
        preds = predict(set, model)
        # using argmax means we assume labels are some range 1:n
        # the onehot indexing doesn't seem to work with GPU, so we call cpu
        mean(cpu(labels[:, set])[onehot(argmax.(preds), uLabels)])
    end
end

opt = ADAMW()
function train(set, params, model)
    grads = Flux.gradient(params) do
        CE(set, model)
    end
    Flux.Optimise.update!(opt, params, grads)
end

bold(v) = "\033[1m$v\033[0m"
fmtval(v::Int) = v
fmtval(v::String) = v
fmtval(v::DateTime) = string(v)
# plenty of precision yet less than default float32 print.
fmtval(v::Real) = @sprintf "%.10f" v
printrow(vs...) = begin
    join((fmtval(v) for v ∈ vs), '\t') |> println
    flush(stdout)
end

cv_parts = 5
iTests = args.cv ? (1:cv_parts) : [1]
iEvals = args.cv ? (1:cv_parts) : [2]

if args.load_model !== nothing
    if args.cv
        models = [@load(args.load_model * "-$i") for i ∈ 1:cv_parts^2-cv_parts]
    else
        models = [@load(args.load_model)]
    end
else
    models = [new_model() for iTest ∈ iTests for iEval ∈ iEvals if iTest != iEval]
end

evalAUCs    = []
evalAccs    = []
preds       = []
pred_labels = []
best_params = []

cv_sets = CV(labels, cv_parts)
iModel = 0
for iTest ∈ iTests
    for iEval ∈ iEvals
        iTest != iEval || continue
        args.cv && println("# test=$(iTest), eval=$(iEval)")
        testset  = cv_sets .== iTest
        evalset  = cv_sets .== iEval
        trainset = .!(testset .| evalset)
        global iModel += 1
        model = models[iModel]
        params = Flux.params(model)

        println("Epoch\tCE\tAcc\tAUC\tTestAUC\tTime")
        push!(best_params, deepcopy(collect(params)))
        best_auc = auc(testset, model)
        best_epoch = 0
        printrow(0, CE(trainset, model), acc(trainset, model), auc(trainset, model), best_auc, now())
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
                    auc_test = bold(fmtval(auc_test)) # bold if best
                    best_auc ≈ 1. && break
                elseif args.early_stop > 0 && epoch - best_epoch > args.early_stop
                    break # early stop if it's been too long
                end
                # use join on tuple instead of list since list will have a 
                # type, and make epoch int->float sometimes.
                printrow(epoch, CE(trainset, model), acc(trainset, model), auc(trainset, model), auc_test, now())
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
    end
end

if args.cv
    println("avg eval AUC = ", mean(evalAUCs))
    println("avg eval Acc = ", mean(evalAccs))
    best_params = cpu(best_params)
    inParamCor = [cor(vec.([abs.(best_params[i][1]), abs.(best_params[j][1])])...)
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

if args.save_model !== nothing
    if args.cv
        for (i, model) ∈ enumerate(models)
            @save args.save_model * "-$i" model
        end
    else
        model = only(models) # has to use macro @save with named variable, not just models[1]
        @save args.save_model model
    end
end

