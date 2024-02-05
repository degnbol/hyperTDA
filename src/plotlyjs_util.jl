#!/usr/bin/env julia
using PlotlyJS
a = … = PlotlyJS.attr

function PlotlyJS.plot(fig::PlotlyJS.SyncPlot, config::PlotConfig)
    plot(fig.plot.data, fig.plot.layout, config=config)
end

axes(indices, attrib) = (
    _axes('x', indices, attrib)...,
    _axes('y', indices, attrib)...
)
axes(indices; kwargs...) = (
    _axes('x', indices; kwargs...)...,
    _axes('y', indices; kwargs...)...
)

"""
Conveniently setting attributes for multiple xaxes.
Example: `Layout(;xaxes(1:4, [attr(title=i) for i in 1:4])...)`
Note: use `Ref` for setting the same vector entry multiple times.
"""
xaxes(indices, attrib) = _axes('x', indices, attrib)
"""
Conveniently setting attributes for multiple yaxes.
Example: `Layout(;yaxes(1:4, [attr(title=i) for i in 1:4])...)`
Note: use `Ref` for setting the same vector entry multiple times.
"""
yaxes(indices, attrib) = _axes('y', indices, attrib)
"""
Conveniently setting attributes for multiple xaxes.
Example: `Layout(;xaxes(1:4, title=1:4)...)`
Note: use `Ref` for setting the same vector entry multiple times.
"""
xaxes(indices; kwargs...) = _axes('x', indices; kwargs...)
"""
Conveniently setting attributes for multiple yaxes.
Example: `Layout(;yaxes(1:4, title=1:4)...)`
Note: use `Ref` for setting the same vector entry multiple times.
"""
yaxes(indices; kwargs...) = _axes('y', indices; kwargs...)
function _axes(xy::Char, indices, attribute::PlotlyBase.PlotlyAttribute)
    _axes(xy, indices, [attribute for _ in indices])
end
function _axes(xy::Char, indices, attributes::Vector)
    indices = _indices(indices)
    (Symbol("$(xy)axis$i") => a for (i,a) in zip(indices, attributes))
end
function _axes(xy::Char, indices; kwargs...)
    indices = _indices(indices)
    ret = Pair{Symbol,Any}[]
    for (k,v) in kwargs
        if !(v isa AbstractVector)
            v = v isa Ref ? v.x : v
            append!(ret, (Symbol("$(xy)axis$(i)_$k") => v for i in indices))
        else
            append!(ret, (Symbol("$(xy)axis$(i)_$k") => _v for (i,_v) in zip(indices,v)))
        end
    end
    ret
end
_indices(index::Int) = range(index, index) |> _indices
function _indices(indices)
    indices = indices |> collect .|> string
    indices[indices .== "1"] .= ""
    indices
end

"""
Heatmap plot with square cells, no background, and convenience for messing with layout.
width should be higher than height to make room for colorbar.
If not, there will be gaps appearing between tick labels and plotted data.
margins (b and l) can be adjusted to make sure tick labels are read in full.
"""
heatmapKwargs(subplot=1, width::Int=750, height::Int=700, b::Int=120, l::Int=120, kwargs...) = (
    autosize=false, # square
    template=:plotly_white,
    width=width,
    height=height,
    margin=attr(b=b, l=l),
    axes(subplot; showgrid=false, ticks="")...,
    # showgrid=false -> hide background
    # ticks="" -> hide ticks that are shown by the template
    # scaleanchor="x" for yaxis makes plot square
    yaxes(subplot; scaleanchor="x" .* _indices(subplot))...,
    kwargs...
)

function heatmap!(p; kwargs...)
    N = p.plot.data |> length
    relayout!(p; heatmapKwargs(range(1,N))..., kwargs...)
end

