#!/usr/bin/env julia
using DelimitedFiles
using JSON3
using PersistenceDiagrams
using Glob: glob
ROOT = `git root` |> readchomp
include("$ROOT/src/flux_utils.jl")
include("$ROOT/src/AUC.jl")
using Statistics: mean, cor
using Dates # now, DateTime etc
using Printf # @sprintf

labels = readdlm("../labels.txt", ';')[:, 2] .|> Int
phs = JSON3.read.("../matroid_generators/PH/2_" .* string.(0:length(labels)-1) .* ".json")
# explicit type to handle trivial case.
# Otherwise we get an error calling PersistenceDiagram
barcodes::Vector{Vector{Tuple{Float64,Float64}}} = [collect(zip(ph[:barcode]...)) for ph in phs]
diagrams = PersistenceDiagram.(barcodes)
img = PersistenceImage(diagrams)
images = img.(diagrams)
features = vec.(images) .|> Vector{Float32}

# 0-indexed
uLabels = 0:4
# NOTE: COMMENT this line to use all
uLabels = [1, 3, 4]
nLabels = length(uLabels)
if nLabels < 5
    idx = labels .∈ Ref(uLabels)
    labels = labels[idx]
    features = features[idx]
end
labels = onehot(labels)

new_model() = Dense(length(features[1]) => nLabels)

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
    mean(labels[:, set][onehot(argmax.(preds), 1:nLabels)])
end

cv_parts = 5
args = (epochs = 500, early_stop = 100)

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

models = [new_model() for iTest ∈ 1:cv_parts for iEval ∈ 1:cv_parts if iTest != iEval]

evalAUCs    = []
evalAccs    = []
preds       = []
pred_labels = []
best_params = []

cv_sets = CV(labels, cv_parts)
iModel = 0
for iTest ∈ 1:cv_parts
    for iEval ∈ 1:cv_parts
        iTest != iEval || continue
        println("# test=$(iTest), eval=$(iEval)")
        testset  = cv_sets .== iTest
        evalset  = cv_sets .== iEval
        trainset = .!(testset .| evalset)
        global iModel += 1
        model = models[iModel]
        params = Flux.params(model)

        println("Epoch\tCE\tAcc\tAUC\tTestAUC\tTime")
        flush(stdout)
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
    end
end

println("avg eval AUC = ", mean(evalAUCs))
println("avg eval Acc = ", mean(evalAccs))
best_params = cpu(best_params)
inParamCor = [cor(vec.([abs.(best_params[i][1]), abs.(best_params[j][1])])...)
    for i in 1:length(best_params)
    for j in 1:i-1]
println("avg cor among in-layer params = ", mean(inParamCor))

