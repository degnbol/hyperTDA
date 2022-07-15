module sample_cr
using Eirene
using Plots
using StatsBase
using Random


function get_cyclerep_in_C_epsilon(cyclerep, C_epsilon)
    cyclerep_v = [sort(cyclerep[:,j]) for j = 1:size(cyclerep, 2)]
    cyclerep_Cepsilon = chain_to_index(cyclerep_v, C_epsilon)
    return cyclerep_Cepsilon
end

function plot_PL(P; title = "")
    # P must have size (n, 3), where n is the number of points
    # plot P
    p = scatter3d(P[:,1], P[:,2], P[:,3], 
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

function plot_cycle_3D(P, cycle; title = "")
    # P must have size (n, 3), where n is the number of points
    # plot P
    p = scatter3d(P[:,1], P[:,2], P[:,3], 
                markersize = 1, markercolor = "grey", label = "",
                xaxis = nothing, yaxis = nothing, zaxis = nothing,
                title = title)
    
    # plot cycle (just color the points )
    cycle_v = vec(hcat(cycle...))
    scatter3d!(p, P[cycle_v,1], P[cycle_v,2], P[cycle_v,3], markercolor = "red", markersize = 2, label ="" )  
    
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

function plot_barcode(barcode; 
    color = :grey56, # default bar color
    selected_bars = [], # index of bars to highlight
    epsilon = missing, # if provided, only highlight the portion of bars on the right side of epsilon
    selection_color = :deeppink2,  # highlight color
    v_line = [], # if provided, draw vertical lines at values
    return_perm = false, # whether to return the permutation index or not
    kwargs...)

    # adjust linewidth, depending on number of intervals
    n = size(barcode,1)

    # find ordering according to birth time
    perm = sortperm(barcode[:,1])

    # non-inf maximum death time
    if filter(!isinf,barcode[:,2]) != []
        death_max = maximum(filter(!isinf,barcode[:,2])) 
    else
        death_max = maximum(barcode[:,1]) * 2
    end

    p = plot(framestyle = :box,
            top_margin = 5 * Plots.PlotMeasures.mm, 
            bottom_margin = 5 * Plots.PlotMeasures.mm, 
            yaxis = nothing;
            kwargs...)
    
    # plot all bars
    idx = 1
    for i in perm
        birth = barcode[i,1]
        death = barcode[i,2]
        
        # assign a death time if bar has infinite death time 
        if isinf(death)
            death = death_max * 1.2
        end
        if i in selected_bars
            
            # if epsilon is missing, highlight the entire bar
            if ismissing(epsilon)
                plot!(p,[birth, death], [idx, idx], legend = false, linecolor = selection_color, hover = "class " *string(i); kwargs...)
            
            # if epsilon is provided, only highlight the portion of the bar on the right side of epsilon    
            else 
                if birth <= epsilon
                    plot!(p,[birth, epsilon], [idx, idx], legend = false, linecolor = color, hover = "class " *string(i); kwargs...)
                    plot!(p,[epsilon, death], [idx, idx], legend = false, linecolor = selection_color, hover = "class " *string(i); kwargs...)
                else
                    plot!(p,[birth, death], [idx, idx], legend = false, linecolor = selection_color, hover = "class " *string(i); kwargs...)
                end
            end
        else
            plot!(p,[birth,death],[idx,idx], legend = false, linecolor = color, hover = "class " *string(i); kwargs...)
        end
        idx += 1
    end

    # plot vertical lines 
    if v_line != []
        plot!(v_line, seriestype="vline", linestyle = :dot, linecolor = :red)
    end

    ylims!((-1, n+1))
    
    if return_perm == true
        return p, perm
    else
        return p
    end
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
    rowvals = Eirene.crows(C["firstv"][dim+1], C["farfaces"][dim+1], simplex[1])

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
            print("chain doesn't exist in eirene output")
            print(simplex)
            return nothing
        end
    end

    return chain_idx
end

function select_odd_count(
    orig_list::Array)
    # given an array, return the elements that occur odd number of times.
    """
    --- input ---
    orig_list: (N-element array)
    --- output ---
    new_list: (M-element array) containing only the elements that occur odd number of times
            in orig_list.
    """

    count = countmap(orig_list)
    new_list = [item for item in orig_list if count[item] % 2 != 0]
    return unique(new_list)
end

# given a cycle, find all two simplices that share a boundary 

function sample_cyclereps(initial_cycle, n_difference, C_epsilon)
    cycle = copy(initial_cycle)
    differences = Set()
    rv, cp = Eirene.boundarymatrix(C_epsilon, dim = 2)
    for i=1:n_difference
        # find all two simplices
        two_simplices = find_adjacent_two_simplices(cycle, rv)
     
        # shuffle list of two simplices
        two_simplices = shuffle(two_simplices)
        j = 1
        s = two_simplices[j]

        # re-sample until we get a two-simplex that hasn't been sampled before
        while (s in differences) & (j < length(two_simplices))
            j = j+1
            s = two_simplices[j]
        end

        # check if we iterated through all two_simplices
        if (s in differences) & (j == length(two_simplices))
            print("There are no more two-simplices to add.")
            return cycle, differences
        end


        # find new cycle with added simplex
        rv_idx = cp[s]
        # get the three 1-simplices
        one_simplex1 = rv[rv_idx]
        one_simplex2 = rv[rv_idx+1]
        one_simplex3 = rv[rv_idx+2]


        # update cycle
        cycle = vcat(cycle, [one_simplex1, one_simplex2, one_simplex3])
        cycle = select_odd_count(cycle)

        # update list of differences
        push!(differences, s)

    end
    return cycle, differences
end

function sample_multiple_cyclereps(initial_cycle, n_difference, C_epsilon, n_samples)
    sampled_cycles = Dict()
    sampled_differences = Dict()
    
    for i=1:n_samples
        alt_cycle, differences = sample_cyclereps(initial_cycle, n_difference, C_epsilon);
        sampled_cycles[i] = alt_cycle
        sampled_differences[i] = differences
    end
    return sampled_cycles, sampled_differences
    
end


function chain_to_vertex(chain, C_epsilon; dim = 1)
    chain_v = [Eirene.incidentverts(C_epsilon["farfaces"],C_epsilon["firstv"],dim+1,[item]) for item in chain]
    chain_v = [sort(C_epsilon["nvl2ovl"][item]) for item in chain_v];
    return chain_v
end

# given a cycle, find all two simplices that share a boundary 

function find_adjacent_two_simplices(cycle, rv)
    # cycle: list of indices 
    # rv of boundary matrix from Eirene
    
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


function get_valid_birth_parameter(PH, points, bar_idx; increment = 0.01)
    """Computes the birth parameter of a bar_idx. 
    If this parameter doesn't return a valid boundary matrix, increase the parameter until a valid boundary matrix is returned
    
    --- input ---
    PH: Eirene output
    points: xyz coordiantes of points. Must be of shape (3, n)
    bar_idx: (int) for bar of interest
    
    """
    
    barcode_PH = barcode(PH, dim = 1)
    parameter = barcode_PH[bar_idx, 1]

    # try computing the boundary matrix
    PH_param = eirene(points, maxdim = 1, maxrad = parameter, model = "pc")
    try
        rv, cp = Eirene.boundarymatrix(PH_param, dim = 2)
        return parameter
    catch 
    end

    while true
        # increase parameter by "increment"
        parameter += increment
        PH_param = eirene(points, maxdim = 1, maxrad = parameter, model = "pc")
        
        # try computing the boundary matrix
        try
            rv, cp = Eirene.boundarymatrix(PH_param, dim = 2)
            return parameter
        catch 
        end
    end
end

function get_valid_parameter(PH, points; increment = 0.01)
    """finds the "minimum" parameter at which we can get the boundary matrix from Eirene
    """
    barcode_PH = barcode(PH, dim = 1)
    birth_times = unique(sort(barcode_PH[:,1]))
    
    # try computing the boundary matrix at each birth time 
    for param in birth_times
        PH_param = eirene(points, maxdim = 1, maxrad = param, model = "pc") 
        try
            rv, cp = Eirene.boundarymatrix(PH_param, dim = 2)
            return param
        catch
        end
    end

    
end


end