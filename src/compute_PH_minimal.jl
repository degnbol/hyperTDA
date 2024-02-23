"""computes PH with minimal generators.
To use Gurobi, do the following before running this script:

1. Download Gurobi and obtain a free academic license from 
https://www.gurobi.com/academia/academic-program-and-licenses/

2. Before running this script, open a Julia REPL, and enter
`julia> ENV["GUROBI_HOME"] = PATH_TO_GUROBI` 

written by Iris Yoon
iris.hr.yoon@gmail.com
"""
module PH_minimal

using Pkg

using SparseArrays
using MultivariateStats, Distances, DelimitedFiles, JuMP, JLD, LinearAlgebra, Eirene, Distributions, GLPK
include("minimal_cycles_rational/computePH.jl")
include("minimal_cycles_rational/utilFunctions.jl")
include("minimal_cycles_rational/edge-loss.jl")
include("minimal_cycles_rational/triangle-loss.jl")
include("minimal_cycles_rational/outputFunctions.jl")
include("minimal_cycles_rational/find_sub_bdr_matrix.jl")
include("Eirene_var.jl")
using Plots
using StatsBase
using NPZ
using Gurobi
using Random
using DataFrames
using JSON
using CSV
using LinearAlgebra

export compute_PH_minimal_generators,
        chain_to_index,
        simplex_to_index,
        get_vertex_perm,
        check_homologous_cycles,
        bounding_chain,
        select_odd_count,
        minimize_cycle_jumps,
        plot_curve,
        get_minimal_cycles_rational,
        get_jump_minimized_cycles,
        plot_PD,
        plot_cycle_3D
        

function get_cyclerep_in_C_epsilon(cyclerep, C_epsilon)
    cyclerep_v = [sort(cyclerep[:,j]) for j = 1:size(cyclerep, 2)]
    cyclerep_Cepsilon = chain_to_index(cyclerep_v, C_epsilon)
    return cyclerep_Cepsilon
end

"""
    plot_curve(P)
Plots curve in 3D. 

### Inputs
- `P::Array` of size (n,3), where `n` is the number of points 

### Outputs
- `p`: plots object
"""
function plot_curve(P; title = "")
   # check intput size
    if size(P,2) != 3
        P = Array(Transpose(P))
    end

    # plot P
    p = Plots.scatter3d(P[:,1], P[:,2], P[:,3], 
                markersize = 0.5, markercolor = "grey", label = "",
                xaxis = nothing, yaxis = nothing, zaxis = nothing,
                title = title)
    
    for i =1:size(P,1)-1
        v1_x, v1_y, v1_z = P[i, 1], P[i,2], P[i,3]
        v2_x, v2_y, v2_z = P[i+1, 1], P[i+1, 2], P[i+1,3]
        plot!(p, [v1_x, v2_x], [v1_y, v2_y], [v1_z, v2_z], label = "", linewidth = 4, linecolor = "grey")
    end
    
    return p
end

"""
    plot_PD(barcode; <keyword arguments>)
Plots persistence diagram.
"""
function plot_PD(barcode; 
        highlight = [], 
        highlight_color = :deeppink2, 
        pd_min = nothing, 
        pd_max = nothing,  
        diagonal_lw = 5,
        kwargs...)
    points = barcode
    

    if size(barcode,1) == 0
        # plot diagonal line
        p = Plots.plot([0, 1], [0, 1], 
        labels ="", 
        linewidth = diagonal_lw,
        framestyle = :box,
        xlims = (0,1),
        ylims = (0,1),
        aspect_ratio = :equal,
        color = "grey"; 
        kwargs...)
        return p
    end
    # find index of points with death parameter == death
    idx = findall(x -> x == Inf, points[:,2])
    
    # plot points with death < Inf
    idx2 = [i for i in 1:size(points,1) if i ∉ idx]
    p = Plots.scatter(points[idx2,1], points[idx2,2]; kwargs..., color = "grey", labels = "", hover = idx2, alpha = 0.5)
    
    # find max death value
    max_death = maximum(points[idx2, 2])
    
    # plot points with death parameter == Inf
    death_plot = ones(size(idx,1)) * max_death
    
    Plots.scatter!(p, points[idx,1], death_plot, aspect_ratio = :equal, legend=:bottomright, hover = idx, labels="", color ="red"; kwargs...)

    # plot diagonal line
    if pd_max == nothing
        
        min_birth = minimum(points[:,1]) * 0.8
        max_death = max_death * 1.1
        Plots.plot!(p, [min_birth, max_death], [min_birth, max_death], 
            labels ="", 
            linewidth = diagonal_lw,
            framestyle = :box,
            xlims = (min_birth, max_death),
            ylims = (min_birth, max_death),
            color = "grey"; 
            kwargs...)
    else
        max_death = pd_max
        min_birth = pd_min
        Plots.plot!(p, [min_birth, max_death], [min_birth, max_death], 
            labels ="", 
            linewidth = diagonal_lw,
            framestyle = :box,
            xlims = (min_birth, max_death),
            ylims = (min_birth, max_death),
            color = "grey"; 
            kwargs...)
    end
    
     # if highlight is provided, color specific points with the given color
    if highlight != []
         Plots.scatter!(points[highlight,1], points[highlight,2]; kwargs..., color = highlight_color, labels = "", hover = highlight)
    end
    
    return p
end

"""
    plot_cycle_3D(P, cycle)
Plots specific cycle in 3-dimensions. P must have size (n, 3), where n is the number of points.
"""
function plot_cycle_3D(P, cycle; title = "")
    # check size of P
     if size(P,2) != 3
        P = Array(Transpose(P))
    end

    # plot P
    p = Plots.scatter3d(P[:,1], P[:,2], P[:,3], 
                markersize = 1, markercolor = "grey", label = "",
                xaxis = nothing, yaxis = nothing, zaxis = nothing,
                title = title)
    
    # plot cycle (just color the points )
    cycle_v = vec(hcat(cycle...))
    Plots.scatter3d!(p, P[cycle_v,1], P[cycle_v,2], P[cycle_v,3], markercolor = "red", markersize = 2, label ="" )  
    
    # convert cycle to vector if needed 
    if (cycle isa Vector) == false
        cycle = [sort(cycle[:,i]) for i=1:size(cycle, 2)]
    end 
    
    # plot 1-simplices
    for i =1:size(cycle,1)
        v1, v2 = cycle[i]
        v1_x, v1_y, v1_z = P[v1, 1], P[v1,2], P[v1,3]
        v2_x, v2_y, v2_z = P[v2, 1], P[v2, 2], P[v2,3]
        plot!(p, [v1_x, v2_x], [v1_y, v2_y], [v1_z, v2_z], label = "", linewidth = 4, linecolor = "red")
    end
    
    return p
end


function chain_to_vertex(chain, C_epsilon; dim = 1)
    chain_v = [Eirene_var.incidentverts(C_epsilon["farfaces"],C_epsilon["firstv"],dim+1,[item]) for item in chain]
    chain_v = [sort(C_epsilon["nvl2ovl"][item]) for item in chain_v];
    return chain_v
end




"""
    find_adjacent_two_simplices(cycle, rv)
Given a cycle, find all two simplices that share a boundary.

### Inputs
- `cycle`: list of indices
- `rv`: rv of boundary matrix from Eirene.

### Outputs
"""
function find_adjacent_two_simplices(cycle, rv)
    
    two_simplices = []
    for simplex in cycle
        # example: find all 2-simplices (in the boundary matrix C_epsilon ) that contain 1-simplex "2"
        rv_idx = findall(x -> x == simplex, rv)
        # find the index of 2-simplex corresponding to these
        two_simplex_idx = [ceil(Int32, x / 3) for x in rv_idx]
        two_simplices = vcat(two_simplices, two_simplex_idx)
    end
    return two_simplices
end


"""
    select_odd_count(orig_list)
Given an array, return the elements that occur odd number of times.

### Inputs
- `orig_list`: N-element array

### Outputs
- `new_list`: M-element array containing only elements that occur odd number of times in `orig_list`

"""
function select_odd_count(
    orig_list::Array)

    count = countmap(orig_list)
    new_list = [item for item in orig_list if count[item] % 2 != 0]
    return unique(new_list)
end


"""
    bounding_chain()
Check if a given chain is a boundary.
"""
function bounding_chain(C;chain=zeros(Int64,0),dim=1)
    ##### CREDITS: This function was written by Greg Henselman-Petrusek. To be included in a future version of Eirene #####
    ##### https://github.com/Eetion/Eirene.jl #####
    
	if isempty(chain)
		return zeros(Int64,0)
	elseif !isempty(Eirene_var.chainboundary(C,chain=chain,dim=dim))
		print("there is no bounding chain"); return nothing
	else
		sd 			= 	dim+1;
		Lrv 		= 	C["Lrv"][sd+1]
		Lcp 		= 	C["Lcp"][sd+1]
		Rrv 		= 	C["Rrv"][sd+1]
		Rcp 		= 	C["Rcp"][sd+1]
		numrows 	= 	Eirene_var.boundarycorank(C,dim=dim)
		translator 	= 	zeros(Int64,Eirene_var.complexrank(C,dim=dim))
		translator[C["tid"][sd+1]] = 1:numrows
		rv 			= 	findall(!iszero,translator[chain])
		rv 			=  	translator[chain[rv]]
		cp 			= 	[1,length(rv)+1]
		rv,cp 		=  	Eirene_var.spmmF2silentLeft(Lrv,Lcp,rv,cp,numrows)
		rv,cp 		=  	Eirene_var.spmmF2silentLeft(Rrv,Rcp,rv,cp,numrows)
		#
		# recall that plo = the (ordered) subvector consisting of
		# the first (rank(boundary operator)) elements of tid
		#
		if maximum(rv) > length(C["plo"][sd+1])
			#print("there is no bounding chain");  # NOTE: Iris removed print
			return nothing
		else
			ocg2rad = 	C["ocg2rad"]
			grain 	= 	C["grain"][sd+1]
			names 	= 	C["phi"][sd+1]
			return 	names[rv]
		end
	end
end


"""
    check_homologous_cycles(v_chain, w_chain, C)
Check if two given cycles [v], [w] are homologous in C

### Intputs
- `v_chain`: (array) of indices of simplices of [v]. Note that simplex indexing must be consistent with C.
- `w_chain`: (array) of indices of simplices of [w]. Note that simplex indexing must be consistent with C. 
- `C``: (dict) output of eirene. When running, eirene, must use the argument `record = "all"`

### Outputs 
- bool: true if [v] = [w] in C
            flase if [v] != [w] in C
"""
function check_homologous_cycles(
    v_chain::Array, 
    w_chain::Array, 
    C::Dict)
    
    # chain of [v] + [w]
    chain_vw = vcat(v_chain, w_chain)

    # only keep simplices that occur odd number of times
    chain_vw = select_odd_count(chain_vw)

    # if chain of v+w is trivial:
    if chain_vw == []
        return true
    else
        # check if [v] + [w] = 0 in H_n(C)
        bounder = bounding_chain(C, chain = chain_vw)
        if bounder == nothing
            return false
        else
            return true
        end
    end
end

function get_vertex_perm(C::Dict)
    # Eirene permutes the vertex indices.
    # Get vertex permutation information from Eirene.
    # note C["nvl2ovl"][i] = k, where i is the Eirene index and k is the original vertex index. 

    """
    --- input ---
    C: (dict) output from Eirene
    --- output ---
    v_perm: (arr) v_perm[k] = i, where k: original vertex index, 
            and i: corresponding Eirene index.
    
            ex) If C = Eirene(distance_matrix), 
            then k corresponds to the kth row/column of distance_matrix
    """

    n_vertices = size(C["nvl2ovl"],1)
    v_perm = zeros(Int64,n_vertices)
    for i=1:n_vertices
        idx = findall(x->x==i, C["nvl2ovl"])[1]
        v_perm[i] = idx
    end
    return v_perm
end




function simplex_to_index(simplex, C)
    # Given a n-simplex [u_0, u_1, ... , u_n], where each u_i is Eirene-indexed,
    # find the index of simplex in C_n according to eirene output C.
    
    """
    --- input ---
    simplex (arr) : [u_0, u_1, ... , u_k], where each u_i is indexed 
                according to Eirene.
                u_0 < u_1 < ... < u_k.
    C (dict): output from Eirene
    --- output ---
    index (int): index of given simplex in C_n according to C.
    --- example use ---
    # 1. get permutation info
    v_perm = get_vertex_perm(C)
    # 2. rewrite simplex in terms of 0-simplices (indexed according to eirene)
    simplex = [96, 183, 188]
    simplex_eirene = sort([v_perm[i] for i in simplex])
    # 3. find simplex index
    simplex_to_index(simplex_eirene, C)
    """ 
    dim = size(simplex,1) - 1

    # given simplex_eirene =[u_0, u_1, ... , u_n], find the portion of row-indices for column u0
    rowvals = Eirene_var.crows(C["firstv"][dim+1], C["farfaces"][dim+1], simplex[1])

    # find location of [u_1, ..., u_n] in rowvals.
    # that is, find the index of [u_1, ..., u_n] among C_{n-1}.
    
    # if simplex is 1-dimensional
    if dim == 1
        rowloc = findall(x->x==simplex[2], rowvals)[1]
        
    # if simplex is higher dimensional
    else
        rowloc0 = simplex_to_index(simplex[2:end], C)
        rowloc = findall(x->x==rowloc0, rowvals)[1]
    end

    # index of simplex according to Eirene
    index = C["firstv"][dim+1][simplex[1]] + rowloc - 1
    
    return index
end



function chain_to_index(
    chain, 
    C::Dict)
    # Given an n-dimensional chain expressed as a list of simplices 
    # chain = [simplex_1, simplex_2, ... , simplex_k],
    # where each simplex_j = [v_0, ... ,v_n], a list of its (n+1) vertices,
    # rewrite the chain as a list of integers [i_1, ..., i_k], 
    # where each $i_j$ is the index of simplex_j in C_n according to Eirene. 

    """
    --- input ---
    chain (arr) : [s1, s2, ... , sn], 
                where each s_j is a n-simplex of form [v_0, v_1, ..., v_n],
                a list of its (n+1) vertices
    C (dict): output from Eirene
    --- output ---
    chain_idx (arr): [i_1, ... ,i_n], where each i_j is the index of s_j according to C
    --- example ---
    test_chain =[[96, 183, 188],[99, 111, 188]]
    chain_idx = chain_to_index(test_chain, C)
    print(chain_idx)
    # check
    dim = size(test_chain[1],1)-1
    for item in chain_idx
        simplex = EireneVar.incidentverts(C["farfaces"], C["firstv"], dim+1, [item])
        simplex = sort(C["nvl2ovl"][simplex])
        print(simplex)
    end
    """
    # get permutation of vertices
    v_perm = get_vertex_perm(C)

    chain_idx = []
    for simplex in chain
        simplex_eirene = sort([v_perm[i] for i in simplex])
        try 
            simplex_idx = simplex_to_index(simplex_eirene, C)
            append!(chain_idx, simplex_idx)
        catch
            # chain doesn't exist
            return nothing
        end
    end

    return chain_idx
end

function prsb_length_Edge(C; requireIntegralSol = false)
    optimized = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, length(C.generators[1]))
    for l in 1: length(C.generators[1])
        uniform_weighted_minimal_gen, uniform_gen_len, uniform_zeroNorm = findLenEdgeOptimalCycles_prs(C, 1, l, optimized, requireIntegralSol, true, true)
        optimized[l] = vcat(uniform_weighted_minimal_gen, zeros(length(C.generators[1][1]) - length(uniform_weighted_minimal_gen.nzval)))
    end
    return optimized
end


"""
    cycle_list_vertices()
Given a cycle of format [[v0, v1], [v1,v2], [v2,v3], ...],
return the list of vertices [v0, v1, v2 ... ]
"""
function cycle_list_vertices(cycle)
    
    return sort(collect(Set(vcat(cycle...))))
end
    

"""
    minimize_cycle_jumps
Given a cycle, decrease the number of jumps between consecutive vertices.
"""
function minimize_cycle_jumps(minimal_cycle, pc, epsilon)
    cycle_v = cycle_list_vertices(minimal_cycle)

    C = Eirene_var.eirene(Array(transpose(pc)), maxrad = epsilon, record = "all", model = "pc");

    end_vertex = cycle_v[end]
    current_idx = 1
    current_vertex = cycle_v[1]
    while current_vertex < end_vertex
    
        # check if there is a jump between current vertex and the next vertex 
        next_vertex = cycle_v[current_idx + 1]
        if next_vertex - current_vertex == 1
            # no jump. update the current index and vertex
            current_idx += 1
            current_vertex = next_vertex
        else
            # a jump has been found

            # find a sequence with smaller jumps (if possible)
            intermediate, start_v, last_v = minimize_jump(cycle_v, current_idx, C)

            if intermediate != nothing
                insert!(cycle_v, current_idx + 1, intermediate)
                
                # also change the minimal cycle
                filter!(x -> x != [start_v,last_v], minimal_cycle)
                push!(minimal_cycle, [start_v, intermediate])
                push!(minimal_cycle, [intermediate, last_v])
                #print("\ndecreasing jump to: ", cycle_v)
            end

            current_idx += 1
            current_vertex = cycle_v[current_idx]

        end
            
    end
    return minimal_cycle, cycle_v
end

"""
    minimize_jump(cycle_v, idx, C)
Given a cycle (as a sorted list of vertices), minimize the jump between cycle[idx] and cycle[idx+1]
"""
function minimize_jump(cycle_v, idx, C)
    start_v = cycle_v[idx]
    last_v = cycle_v[idx + 1]
    
    diff = last_v - start_v
    for i = 1:(diff-1)
        
        # check if there exists a 2-simplex of the form [[start_v, i], [i, last_v], [start_v, last_v]] 
        # check if the cycle tracing the boundary of 2-simplex is null homologous
        mid_v = start_v + i
        test_chain = [[start_v, mid_v], [mid_v, last_v], [start_v, last_v]]
        #print("\ntesting: ", test_chain)
        test_idx = chain_to_index(test_chain, C)
        
        
        if test_idx != nothing
            # check if test_chain is null homologous (just have to check if the 2-simplex exists)
            homologous = check_homologous_cycles(test_idx, [], C)
            if homologous == true
                return start_v + i, start_v, last_v
            end
        end
    end 
    return nothing, nothing, nothing
    
    
end

"""
    get_minimal_cycles_rational(C)
Computes length-minimal cycles over the rationals. 
"""
function get_minimal_cycles_rational(C) 
    print("\nComputing minimal cycles over the rationals")
    # Compute length-optimized generators
    length_optimized = prsb_length_Edge(C)
    C_barcode = permutedims(hcat(C.barCode[1]...));  
    
    n_bars = size(C_barcode,1)
    
    # collect generators
    generators = Dict()
    for i = 1:n_bars
        generators[i] = [j[1] for j in findall(x -> x != 0, length_optimized[i])]
    end
    
    # check for cycles whose coefficients are not +1 or -1
    coefficients_problem = []
    for i = 1:n_bars
        for j in generators[i] 
            if (Float32(length_optimized[i][j]) != 1) & (Float32(length_optimized[i][j]) != -1)
                print("\nCycle representative of bar " * string(i) * " contains coefficients other than 0, 1, -1")
                push!(coefficients_problem, i)
            end
        end
    end
    
    # if there were coefficients other than +1 or -1, then re-run optimization over rationals while requiring integer solutions
    if coefficients_problem != []
       print("\nComputing minimal cycles over the rationals, requiring integer solutions")
       # re-run length minimization, this time requiring Integer solutions  
       length_optimized = prsb_length_Edge(C; requireIntegralSol = true)
        # collect new generators
        generators = Dict()
        for i = 1:n_bars
            generators[i] = [j[1] for j in findall(x -> x != 0, length_optimized[i])]
        end
    
        # check for cycles whose coefficients are not +1 or -1
        coefficients_problem = []
        for i = 1:n_bars
            for j in generators[i] 
                if (Float32(length_optimized[i][j]) != 1) & (Float32(length_optimized[i][j]) != -1)
                    print("\nCycle representative of bar " * string(i) * " contains coefficients other than 0, 1, -1 (required integer solutions")
                    push!(coefficients_problem, i)
                end
            end
        end

        # get vertex representations of the generator
        generators_v = Dict()
        for i = 1:n_bars
            generators_v[i] = [C.permutedlverts[2][:,j] for j in generators[i]]
        end
        
        return generators, generators_v, coefficients_problem, length_optimized
    else
         
        # get vertex representations of the generator
        generators_v = Dict()
        for i = 1:n_bars
            generators_v[i] = [C.permutedlverts[2][:,j] for j in generators[i]]
        end

        return generators, generators_v, coefficients_problem, length_optimized 
    end
end


"""
    get_jump_minimized_cycles(generators_v, C_barcode, pc)
Find homologous cycle with smaller amounts of jumps between consecutive points.
"""
function get_jump_minimized_cycles(generators_v, C_barcode, pc)
    C_eirene = Eirene_var.eirene(Array(transpose(pc)), maxdim = 1, record = "all", model = "pc");
    eirene_barcode = Eirene_var.barcode(C_eirene, dim = 1); # barcode according to Eirene;
    n_bars = size(C_barcode,1)
    
    # minimize cycle jumps
    minimized_generators = Dict()
    for i = 1:n_bars
        # get index of selected bar in Eirene
        eirene_idx = findall(x -> x == C_barcode[i,2], eirene_barcode[:,2])[1]
        # minimize jumps
        epsilon = eirene_barcode[eirene_idx,1]
        jump_minimized_generator, _ = minimize_cycle_jumps(copy(generators_v[i]), pc, epsilon)

        minimized_generators[i] = jump_minimized_generator
    end
    return minimized_generators
end


"""compute_PH_minimal_generators(input_file, output_file; allow_nonbinary_coefficients)
Computes persistent homology in dimension 1 and returns the minimal generators. 

### Inputs
- `input_file`: path to file containing xyz-coordinates of points on a curve. Must end with `.npy` or `.tsv`.
- `output_file`: path to output file. Must end with '.json'
- `allow_nonbinary_coefficients = true`: if set to true, then whenever we encounter a generator with coefficients other than {0, 1, -1}, we consider its support (1-simplices with non-zero coefficients) as the generator. 

### Outputs
- None. Barcode and minimal generators are saved in `output_file``.
"""
function compute_PH_minimal_generators(input_file, output_file; allow_nonbinary_coefficients = true)
    
    # check that input file ends with .npy
    if (split(input_file, ".")[end] != "npy") & (split(input_file, ".")[end] != "tsv")
        print("input file must end with .npy or .tsv")
        throw(error())
    end
    # check that output file ends with .json
    if split(output_file, ".")[end] != "json"
        print("output file must end with .json")
        throw(error())
    end
    
    # read input coordinates
    if split(input_file, ".")[end] == "npy"
        pc = npzread(input_file)
    end

    if split(input_file, ".")[end] == "tsv"
        pc = CSV.read(input_file, DataFrame)
        pc = Array(pc[!,["x", "y", "z"]])
    end
    
    # make sure the input has format (3, n_points)
    if size(pc,1) != 3
        pc = Array(Transpose(pc))
    end
    
    # input should have format ((n_features, n_points))
    C = computeHomology(pc, false, 1)

    # compute barcode
    C_barcode = permutedims(hcat(C.barCode[1]...))

    # check that barcode is nonempty
    if isempty(C_barcode) == true
        print("\nBarcode is trivial")
        dic = Dict()
        dic["barcode"] = C_barcode
        

        open(output_file, "w") do io
            JSON.print(io, dic)
        end

    else
        # find minimal cycles (over rationals)
        generators, generators_v, coefficients_problem, _ = get_minimal_cycles_rational(C)

        if (coefficients_problem == []) | (allow_nonbinary_coefficients == true)
            # find jump-minimized cycles
            minimized_generators = get_jump_minimized_cycles(generators_v, C_barcode, Array(transpose(pc)))

            # save output
            representatives = [minimized_generators[i] for i = 1:length(minimized_generators)]

            dic = Dict()
            dic["barcode"] = C_barcode
            dic["representatives"] = representatives
            dic["non_Z2_coefficients"] = coefficients_problem
            dic["generators_rational"] = generators_v #these are the generators from optimizing over the rationals

            open(output_file, "w") do io
                JSON.print(io, dic)
            end
        # encountered problems with coefficients, and allow_nonbinary_coefficients is false
        else
            print("\nEncountered coefficients other than 0, 1, -1. User specified `allow_nonbinary_coefficients`` parameter to be false.")
            print("\nNot computing jump-minimized generators.")
            dic = Dict()
            dic["barcode"] = C_barcode
            dic["representatives"] = "not computed"
            dic["non_Z2_coefficients"] = coefficients_problem
            dic["generators_rational"] = generators_v #these are the generators from optimizing over the rationals

            open(output_file, "w") do io
                JSON.print(io, dic)
            end

        end
    end
end

end
