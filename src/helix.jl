#!/usr/bin/env julia

"""
Get x,y,z for a right-handed helix.
- ts: vector of time points
- radius: helix radius
- pitch: distance traveled along helix (z-axis) after 1 rotation
"""
function helix(ts::AbstractVector, radius=1, pitch=2π)
    hcat(radius .* cos.(ts), radius .* sin.(ts), pitch / (2π) .* ts)
end
"""
- n: number of points
- nRot: number of rotations
"""
function helix(n::Int, nRot::Int, radius=1, pitch=2π)
    step = nRot / n
    ts = 0:step:nRot-step
    helix(ts .* 2π, radius, pitch)
end

"Take e.g. a right-handed helix and make it left-handed."
mirrorX(h::Matrix) = hcat(-h[:, 1], h[:, 2], h[:, 3])



