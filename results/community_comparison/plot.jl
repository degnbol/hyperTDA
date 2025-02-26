#!/usr/bin/env julia
using DataFrames, CSV
using DelimitedFiles
using JSON
using Clustering
using Chain
using NPZ
using Statistics: mean, median
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

"""
Relabel communities for plotting.
0 will refer to all singletons in output and is understood as such if given in the input.
1... will be the rest, where the new label is sorted according to the median along the spatial curve for each cluster.
Assumes input is ordered along the curve.
"""
function order_comms(comms::Vector{Int})
    out = zeros(Int, length(comms))
    uComm = comms |> unique
    singletons = uComm[[sum(comms .== u) == 1 for u in uComm]]
    uComm = setdiff(uComm, [singletons; 0]) # 0 to take each input 0 as a singleton

    medians = [median(findall(comms .== u)) for u in uComm]
    # sort according to medians
    uComm = uComm[sortperm(medians)]
    
    for (i, u) in enumerate(uComm)
        out[comms .== u] .= i
    end
    out
end

"""
- comms: manual override of colorscale
"""
function trace3d(_df; comms=nothing)
    if comms === nothing comms = sort(unique(_df.community)) end
    comms = setdiff(comms, [0])
    palette = ["#09ffff", "#e763fa", "blue", "red", "maroon", "green"][comms]
    # plotly does colorscales weird, so we need colorscale to be [[0, "gray"], [...], ..., [1, COLOUR]]
    colorscale = [[0, "gray"], [[i/length(palette), col] for (i,col) in enumerate(palette)]...]
    scatter3d(_df; x=:x, y=:y, z=:z, marker=attr(
        autocolorscale=false, # use manual scale
        color=_df.community, # has to be numerical array as per help for autocolorscale
        colorscale=colorscale,
    ))
end
function plot_trace3d(_df::DataFrame; comms=nothing)
    plot(trace3d(_df; comms=comms), Layout(
        template="plotly_white",
        paper_bgcolor="rgba(1,1,1,0)", # transparent bg doesn't work in pdf
        scene=attr(
            xaxis_zeroline=false,
            yaxis_zeroline=false,
            zaxis_zeroline=false,
        ),
    ))
end
function plot_trace3d(method::T, dataseti::Int; comms=nothing) where T<:AbstractString
    E = method == "Louvain"
    _df = df[
        (df.dataset .== "default_pipeline") .&
        (df.dataseti .== dataseti) .&
        (df.method .== method) .&
        (df.E .== E) .&
        (df.V .== false), :]
    @assert nrow(_df) > 0 "No data selected"
    _df.community = order_comms(_df.community)
    plot_trace3d(_df; comms=comms)
end
function save_trace3d(method::T, dataseti::Int=2; comms=nothing) where T<:AbstractString
    E = method == "Louvain"
    fig = plot_trace3d(method, dataseti; comms=comms)
    # half of paperwidth minus margins, 200 DPI
    width = round(Int, (21 - 3) / 2 * 0.4 * 200)
    savefig(
        fig, "figures/default_pipeline$dataseti-$method$(E ? "E" : "").pdf";
        width=width,
        height=round(Int, width * .72),
        scale=2,
    )
end

# annotate xyz
leftjoin!(df, dfPC, on=[:dataset, :dataseti, :residue])

for dataseti in 1:2
    for method in unique(df.method)
        fig = plot_trace3d(method, dataseti)
        # savefig(fig, "figures/default_pipeline$dataseti-$method$(E ? "E" : "").html")
    end
end

save_trace3d("Leiden")
save_trace3d("Louvain")
save_trace3d("CliqueLouvain"; comms=[1, 5, 3, 4, 2])
save_trace3d("CFMRandom"; comms=[5])
save_trace3d("CFMCNMLike"; comms=[1, 3, 4])
save_trace3d("AON_Louvain"; comms=[2, 5, 4])

