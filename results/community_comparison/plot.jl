#!/usr/bin/env julia
using DataFrames, CSV
using DelimitedFiles
using JSON
using Clustering
using Chain
using NPZ
ROOT = `git root` |> readchomp
include("$ROOT/src/plotlyjs_util.jl")

mkpath("figures/")

df = CSV.read("./community_comparison.tsv.gz", DataFrame; delim='\t')

louvainDefault = JSON.parsefile("../../examples/default_pipeline/communities.json")
# louvain was done on the interp H, but we now did all analyses on uninterp H
# louvainInterp = JSON.parsefile("../../examples/with_interpolation/communities.json")

for dataseti in 1:2
    # +1 for python->julia indexing
    comm = louvainDefault["curve_$dataseti"] .+ 1
    append!(df, DataFrame(dataset="default_pipeline", dataseti=dataseti, method="Louvain", E=true, V=false, residue=1:length(comm), community=comm))
end

xyzs = [
    npzread("../../examples/default_pipeline/pointClouds/curve_1.npy"),
    npzread("../../examples/default_pipeline/pointClouds/curve_2.npy"),
    readdlm("../../examples/with_interpolation/pointClouds/AnDi_1.tsv"),
    readdlm("../../examples/with_interpolation/pointClouds/AnDi_2.tsv")
]
dfPCs = DataFrame.(xyzs, Ref([:x, :y, :z]))
for dfPC in dfPCs dfPC[!, "residue"] = 1:nrow(dfPC) end
dfPCs[1][!, "dataset"] .= "default_pipeline"
dfPCs[2][!, "dataset"] .= "default_pipeline"
dfPCs[3][!, "dataset"] .= "with_interpolation"
dfPCs[4][!, "dataset"] .= "with_interpolation"
dfPCs[1][!, "dataseti"] .= 1
dfPCs[2][!, "dataseti"] .= 2
dfPCs[3][!, "dataseti"] .= 1
dfPCs[4][!, "dataseti"] .= 2
dfPC = vcat(dfPCs...)

dfMethods = df[!, [:method, :E, :V]] |> unique
dfCross = crossjoin(rename(x->x*"_i", dfMethods), rename(x->x*"_j", dfMethods))
dfCross = leftjoin(dfCross, rename(df, :community=>:community_i); on=[:method_i=>:method, :E_i=>:E, :V_i=>:V])
dfCross = leftjoin(dfCross, rename(df, :community=>:community_j); on=[:method_j=>:method, :E_j=>:E, :V_j=>:V, :dataset, :dataseti, :residue])
dropmissing!(dfCross)

dfMetric = @chain dfCross groupby(Not([:residue, :community_i, :community_j])) combine([:community_i, :community_j] => vmeasure => "vmeasure")
dfMetric[!, "label_i"] = dfMetric.method_i
dfMetric[!, "label_j"] = dfMetric.method_j
dfMetric[dfMetric.E_i, "label_i"] .*= "E"
dfMetric[dfMetric.V_i, "label_i"] .*= "V"
dfMetric[dfMetric.E_j, "label_j"] .*= "E"
dfMetric[dfMetric.V_j, "label_j"] .*= "V"
dfMetricAvg = @chain dfMetric groupby(Not([:dataset, :dataseti, :vmeasure])) combine(:vmeasure => mean)

p = heatmap(dfMetricAvg; x=:label_i, y=:label_j, z=:vmeasure_mean)
heatmap!(p)
savefig(p, "figures/vmeasureAvg.pdf", width=750, height=700)
# then delete additional page generated.

dfMetrics = @chain dfMetric groupby([:dataset, :dataseti]) collect
titles = [only(unique(_df.dataset))*string(only(unique(_df.dataseti))) for _df in dfMetrics]

p = make_subplots(rows=2, cols=2, shared_xaxes=false, horizontal_spacing=.15, vertical_spacing=.15)
add_trace!(p, heatmap(dfMetrics[1]; x=:label_i, y=:label_j, z=:vmeasure), row=1, col=1)
add_trace!(p, heatmap(dfMetrics[2]; x=:label_i, y=:label_j, z=:vmeasure), row=1, col=2)
add_trace!(p, heatmap(dfMetrics[3]; x=:label_i, y=:label_j, z=:vmeasure), row=2, col=1)
add_trace!(p, heatmap(dfMetrics[4]; x=:label_i, y=:label_j, z=:vmeasure), row=2, col=2)
heatmap!(p; width=750*1.5, height=700*1.5)
savefig(p, "figures/vmeasure.pdf", width=round(Int,750*1.5), height=round(Int,700*1.5))

function singleton0!(df)
    # community 0 is for singletons. Place them here unless there already is a community 0.
    comm = df.community
    uniq = comm |> unique
    singletons = uniq[[sum(comm .== u) == 1 for u in uniq]]
    if !isempty(singletons)
        @assert 0 ∉ comm
        df.community[comm .∈ Ref(singletons)] .= 0
    end
    df
end

function trace3d(df)
    palette = ["gray", "#09ffff", "#e763fa", "#ab63fa", "blue", "red", "maroon", "green"]
    colorscale = setdiff(unique(df.community), [0])
    colorscale = [0; colorscale ./ maximum(colorscale; init=1)] # colorscale [0, 1]
    colorscale = zip(colorscale, palette) |> collect

    scatter3d(df; x=:x, y=:y, z=:z, marker=attr(
        autocolorscale=false, # use manual scale
        color=df.community, # has to be numerical array as per help for autocolorscale
        colorscale=colorscale
    ))
end

# annotate xyz
leftjoin!(df, dfPC, on=[:dataset, :dataseti, :residue])

for dataseti in 1:2
    for method in unique(df.method)
        E = method == "Louvain"
        p = df[
            (df.dataset .== "default_pipeline") .&
            (df.dataseti .== dataseti) .&
            (df.method .== method) .&
            (df.E .== E) .&
            (df.V .== false), :] |> singleton0! |> trace3d |> plot
        savefig(p, "figures/default_pipeline$dataseti-$method$(E ? "E" : "").html")
    end
end

