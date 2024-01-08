#!/usr/bin/env julia
using CSV, DataFrames
using Chain: @chain

# read avg eval accuracies
avgEvalAcc = @chain `./avgEvalAcc.sh` readchomp split(_, '\n') split.(_, '\t')
avgEvalAccGen = [v[1] for v in avgEvalAcc]
avgEvalAccLab = [v[2] for v in avgEvalAcc]
avgEvalAcc = @chain [v[3] for v in avgEvalAcc] parse.(Float64, _)

# parse legend labels
avgEvalModels = fill("all 5", length(avgEvalAcc))
avgEvalModels[match.(r"134", avgEvalAccLab) .!== nothing] .= "1,3,4"
avgEval0 = fill("H+V", length(avgEvalAcc))
avgEval0[match.(r"0H", avgEvalAccLab) .!== nothing] .= "V"
avgEval0[match.(r"0V", avgEvalAccLab) .!== nothing] .= "H"

df = DataFrame(AvgEvalAcc=avgEvalAcc, Generators=avgEvalAccGen, Models=avgEvalModels, Input=avgEval0)
CSV.write("avgEvalAcc.tsv", df; delim='\t')
