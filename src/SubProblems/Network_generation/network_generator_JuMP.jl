using JuMP

"""
$(TYPEDEF) 
Holds approximations and approaches to deal with nonconvexities in the objective of the network optimization problem: In particular, approaches related to nonconvexities in LMTD approximation and economies of scale. 
"""
abstract type NetworkObjective end

"""
$(TYPEDEF) 
Minimizes total network area. Uses Arithmetic mean to approximate LMTD. Results in substantial error.
"""
struct AreaArithmeticMean <: NetworkObjective end

"""
$(TYPEDEF) 
Minimizes total network area. Uses Paterson formula to approximate LMTD. 
"""
struct AreaPaterson <: NetworkObjective end

"""
$(TYPEDEF) 
Minimizes total network cost using specified economies of scale. Uses Paterson formula to approximate LMTD. 
"""
struct CostScaledPaterson <: NetworkObjective end


"""
$(TYPEDSIGNATURES)

Generates the Heat Exchanger Network. Define a type of superstructure for each stream. 
"""
function generate_network!(prob::ClassicHENSProblem, EMAT, overall_network::Dict{String, AbstractSuperstructure}; obj_func::NetworkObjective = CostScaledPaterson(), time_limit = 20.0, optimizer, verbose = false, cost_coeff, scaling_coeff, base_cost, save_model = false)
    verbose && @info "Solving the Network Generation subproblem"
    
    haskey(prob.results_dict, :y) || error("Stream match data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :Q) || error("HLD data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :match_list) || error("Match list data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :T_bounds) || generate_temperature_bounds!(prob)

    #y, Q = prob.results_dict[:y],  prob.results_dict[:Q]
    
    model = Model()
    
    # Tuple set with stream and utility s and Edge e
    all_e_tuple_vec = []
    for stream in prob.all_names
        for edge in overall_network[stream].edges
            push!(all_e_tuple_vec, (stream, edge))
        end
    end

    #=
    # Tuple set with only streams (not utilities) and edges e
    stream_e_tuple_vec = []
    for stream in prob.stream_names
        for edge in overall_network[stream].edges
            push!(stream_e_tuple_vec, (stream, edge))
        end
    end
    =#

    # 1. Declaring the stream-wise variables
    @variable(model, 0.0 <= t[all_e_tuple_vec])
    @variable(model, 0.0 <= f[all_e_tuple_vec])

    # 2. Sets stream-wise constraints
    for stream in prob.all_names
        add_stream_constraints!(model, prob.all_dict[stream], overall_network[stream], prob)
    end

    
    # 3. Setting the variables and constraints for the matches
    match_list = Tuple[]
    for hot in prob.hot_names
        for cold in prob.results_dict[:match_list][hot]
            push!(match_list, (hot,cold))
        end
    end
    @variable(model, ΔT_upper[match_list] >= EMAT)
    @variable(model, ΔT_lower[match_list] >= EMAT)
    
    U_dict = Dict()
    for match in match_list
        hot, cold = prob.all_dict[match[1]], prob.all_dict[match[2]]
        add_match_feasibility_constraints!(model, hot, cold, overall_network[hot.name], overall_network[cold.name])
        # Used to avoid `_parse_NL_expr_runtime` error with U as a function.
        push!(U_dict, match => U(hot, cold))
    end
    verbose && println(U_dict)

    
    set_objective_func!(model, match_list, obj_func, prob, EMAT; U_dict = U_dict, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost)

    set_optimizer(model, optimizer)
    set_optimizer_attribute(model, "MaxTime", time_limit)
    set_optimizer_attribute(model, "AbsConFeasTol", 1E-2)

    save_model && push!(prob.results_dict, :network_gen_model => model)

    optimize!(model)
    postprocess_network!(prob, model, match_list)
    return
end


"""
$(TYPEDEF)
Adds constraints for each `stream`.
"""
function add_stream_constraints!(model::AbstractModel, stream::AbstractStream, superstructure::AbstractSplitSuperstructure, prob::ClassicHENSProblem)
    # 1. Add bounds on the `t` and `f` variables
    for edge in superstructure.edges
        #set_lower_bound(model[:t][(stream.name, edge)], prob.results_dict[:T_bounds][stream.name][1])
        #set_upper_bound(model[:t][(stream.name, edge)], prob.results_dict[:T_bounds][stream.name][2])
        #set_upper_bound(model[:f][(stream.name, edge)], stream.mcp)
    end
 
    # 2. Source and Sink constraints on `t` and `f`
    SO_BS_edge = only(out_edges(only(superstructure.source), superstructure))
    @constraint(model, model[:f][(stream.name, SO_BS_edge)] == stream.mcp)
    @constraint(model, model[:t][(stream.name, SO_BS_edge)] == stream.T_in)

    #BM_SK_edge = only(in_edges(only(superstructure.sink), superstructure))
    #@constraint(model, model[:f][(stream.name, BM_SK_edge)] == stream.mcp)
    #@constraint(model, model[:t][(stream.name, BM_SK_edge)] == stream.T_out)  

    # 3. Adding mass balance constraints for all nodes except `Source` and `Sink`
    internal_nodes = filter(superstructure.nodes) do v
        !isa(v, Source) && !isa(v, Sink)
    end
    for node in internal_nodes
        @constraint(model, sum(model[:f][(stream.name, edge)] for edge in in_edges(node, superstructure)) == sum(model[:f][(stream.name, edge)] for edge in out_edges(node, superstructure)))
    end
    
    # 4. Setting equal temperatures for all streams entering and leaving a splitter
    for splitter in superstructure.splitters
        edge_in = only(in_edges(splitter, superstructure)) # Only one in edge allowed
        for out_edge in out_edges(splitter, superstructure)
            @constraint(model, model[:t][(stream.name, out_edge)] == model[:t][(stream.name, edge_in)])
        end
    end

    # 5. Setting mixer energy balance
    for mixer in superstructure.mixers
        edge_out = only(out_edges(mixer, superstructure)) # Only one out edge allowed
        @NLconstraint(model, sum(model[:f][(stream.name, in_edge)]*model[:t][(stream.name, in_edge)] for in_edge in in_edges(mixer, superstructure)) == model[:f][(stream.name, edge_out)]*model[:t][(stream.name, edge_out)])
    end

    # 6. Heat Exchanger energy balance
    Q = prob.results_dict[:Q]
    for hx in superstructure.hxs
        edge_in = only(in_edges(hx, superstructure)) # Only one in edge allowed
        edge_out = only(out_edges(hx, superstructure)) # Only one out edge allowed
        match = hx.match
        if stream isa HotStream # Hot stream gets cooled, match is a ColdStream
            @NLconstraint(model, model[:f][(stream.name, edge_in)]*model[:t][(stream.name, edge_in)] - model[:f][(stream.name, edge_out)]*model[:t][(stream.name, edge_out)]  == Q[match, stream.name])
        elseif stream isa ColdStream # Cold stream gets heat, match is a HotStream
            @NLconstraint(model, model[:f][(stream.name, edge_out)]*model[:t][(stream.name, edge_out)] - model[:f][(stream.name, edge_in)]*model[:t][(stream.name, edge_in)] == Q[stream.name, match])
        else
            error("Stream type $(typeof(stream)) not supported.")
        end
    end
end

function add_stream_constraints!(model::AbstractModel, stream::AbstractUtility, superstructure::AbstractSplitSuperstructure, prob::ClassicHENSProblem)
    #=
    for edge in superstructure.edges
        #set_lower_bound(model[:t][(stream.name, edge)], prob.results_dict[:T_bounds][stream.name][1])
        #set_upper_bound(model[:t][(stream.name, edge)], prob.results_dict[:T_bounds][stream.name][2])
        #set_upper_bound(model[:f][(stream.name, edge)], 0.0)
    end

    # 3. Adding mass balance constraints for all nodes except `Source` and `Sink`
    internal_nodes = filter(superstructure.nodes) do v
        !isa(v, Source) && !isa(v, Sink)
    end
    for node in internal_nodes
        @constraint(model, sum(model[:f][(stream.name, edge)] for edge in in_edges(node, superstructure)) == sum(model[:f][(stream.name, edge)] for edge in out_edges(node, superstructure)))
    end

    # 4. Setting equal temperatures for all streams entering and leaving a splitter
    for splitter in superstructure.splitters
        edge_in = only(in_edges(splitter, superstructure)) # Only one in edge allowed
        for out_edge in out_edges(splitter, superstructure)
            @constraint(model, model[:t][(stream.name, out_edge)] == model[:t][(stream.name, edge_in)])
        end
    end
    =#
    #=
    # 5. Setting mixer energy balance
    for mixer in superstructure.mixers
        edge_out = only(out_edges(mixer, superstructure)) # Only one out edge allowed
        @NLconstraint(model, sum(model[:f][(stream.name, in_edge)]*model[:t][(stream.name, in_edge)] for in_edge in in_edges(mixer, superstructure)) == model[:f][(stream.name, edge_out)]*model[:t][(stream.name, edge_out)])
    end
    =#
    #=
    # Only set edge temperatures for LMTD. All edge temperature before the HX are equal to Source temperature. All after HX are set at Sink temperature.   
    # Hack is to set all edges connected to MajorSplitter at same temperature as stream.T_in, and all to MajorMixer at same temperature as stream.T_out.
    major_splitter = filter(v -> (v isa MajorSplitter), superstructure.nodes)
    major_mixer = filter(v -> (v isa MajorMixer), superstructure.nodes)

    splitter_edges = union(in_edges(major_splitter[1], superstructure), out_edges(major_splitter[1], superstructure))
    mixer_edges = union(in_edges(major_mixer[1], superstructure), out_edges(major_mixer[1], superstructure))
    
    for edge in splitter_edges
        @constraint(model, model[:t][(stream.name, edge)] == stream.T_in)
    end

    for edge in mixer_edges
        @constraint(model, model[:t][(stream.name, edge)] == stream.T_out)
    end
    =#
    return
end

"""
$(TYPEDEF) 
Function used to set the ΔT_upper and ΔT_lower to appropriate temperatures.
"""
function add_match_feasibility_constraints!(model::AbstractModel, hot::Union{HotStream, SimpleHotUtility}, cold::Union{ColdStream, SimpleColdUtility}, hot_superstructure::AbstractSuperstructure, cold_superstructure::AbstractSuperstructure)
    hot_HX = only(filter(hot_superstructure.hxs) do v
         v.match == cold.name
    end) # Should only be 1.
    hot_hx_in_edge, hot_hx_out_edge  = only(in_edges(hot_HX, hot_superstructure)), only(out_edges(hot_HX, hot_superstructure))

    # Match is always a cold stream or cold utility
    cold_HX = only(filter(cold_superstructure.hxs) do v
        v.match == hot.name
   end)
   cold_hx_in_edge, cold_hx_out_edge  = only(in_edges(cold_HX, cold_superstructure)), only(out_edges(cold_HX, cold_superstructure))

   # Setting ΔT_upper
   @constraint(model, model[:ΔT_upper][(hot.name, cold.name)] ==  model[:t][(hot.name, hot_hx_in_edge)] - model[:t][(cold.name, cold_hx_out_edge)])

   # Setting ΔT_lower
   @constraint(model, model[:ΔT_lower][(hot.name, cold.name)] ==  model[:t][(hot.name, hot_hx_out_edge)] - model[:t][(cold.name, cold_hx_in_edge)])
end   

"""
$(TYPEDEF) 
Function used to set the objective of each hot stream problem.
"""
function set_objective_func!(model::AbstractModel, match_list, obj_func::CostScaledPaterson, prob::ClassicHENSProblem, EMAT; U_dict, cost_coeff = 1.0, scaling_coeff = 1, base_cost = 0)
    @variable(model, T_LMTD[match_list] >= EMAT)
    for match in match_list
        @NLconstraint(model, model[:T_LMTD][match] ==  ((2/3)*(model[:ΔT_upper][match]*model[:ΔT_lower][match])^0.5 - (1/6)*(model[:ΔT_upper][match] + model[:ΔT_lower][match])))  
    end
    #@NLobjective(model, Min, sum(base_cost + cost_coeff*((1/(smallest_value+model[:T_LMTD][match]))*(1/(U_dict[match[1], match[2]]))*(prob.results_dict[:Q][match[2], match[1]]))^scaling_coeff for match in match_list))
end

"""
$(TYPEDEF) 
Generates bounds on the temperatures for a stream. 
Creates a `prob.results_dict[:T_bounds]::Dict{String, Tuple}`, where the tuple is `(T_LBD, T_UBD)`
"""
function generate_temperature_bounds!(prob::ClassicHENSProblem)
    T_bounds = Dict{String, Tuple}()
    for (k,v) in prob.all_dict
        push!(T_bounds, k => generate_temperature_bounds!(prob, v))
    end
    prob.results_dict[:T_bounds] = T_bounds
end

function generate_temperature_bounds!(prob::ClassicHENSProblem, stream::HotStream)
    T_UBD = stream.T_in
    T_LBD = minimum([prob.all_dict[match].T_in for match in prob.results_dict[:match_list][stream.name]])
    return (T_LBD, T_UBD)
end

function generate_temperature_bounds!(prob::ClassicHENSProblem, stream::ColdStream)
    T_LBD = stream.T_in
    T_UBD = maximum([prob.all_dict[match].T_in for match in prob.results_dict[:match_list][stream.name]])
    return (T_LBD, T_UBD)
end

function generate_temperature_bounds!(prob::ClassicHENSProblem, stream::SimpleColdUtility)
    T_LBD = stream.T_in
    T_UBD = stream.T_out
    return (T_LBD, T_UBD)
end

function generate_temperature_bounds!(prob::ClassicHENSProblem, stream::SimpleHotUtility)
    T_UBD = stream.T_in
    T_LBD = stream.T_out
    return (T_LBD, T_UBD)
end


#=
From Plasmo file, potentially usefull for other set_objective_func!() methods

"""
$(TYPEDEF) 
Sets the objective of a match.
"""
function set_match_objective!(hot::Union{HotStream, SimpleHotUtility}, cold::Union{ColdStream, SimpleColdUtility}, hot_prob::OptiNode, obj_func::AreaArithmeticMean, prob::ClassicHENSProblem)
    U_ij = U(hot, cold)
    Q_ij = prob.results_dict[:Q][cold.name, hot.name]
    coeff = 2*Q_ij/U_ij
    @NLconstraint(hot_prob, hot_prob[:match_objective][cold.name] == coeff/(hot_prob[:ΔT_upper][cold.name] + hot_prob[:ΔT_lower][cold.name])) # Necessary to avoid the parsing error.
end

"""
$(TYPEDEF) 
Sets the objective of a match.
"""
function set_match_objective!(hot::Union{HotStream, SimpleHotUtility}, cold::Union{ColdStream, SimpleColdUtility}, hot_prob::OptiNode, obj_func::AreaPaterson, prob::ClassicHENSProblem)
    U_ij = U(hot, cold)
    Q_ij = prob.results_dict[:Q][cold.name, hot.name]
    coeff = Q_ij/U_ij
    @NLconstraint(hot_prob, hot_prob[:match_objective][cold.name] == coeff/(((2/3)*(hot_prob[:ΔT_upper][cold.name]*hot_prob[:ΔT_lower][cold.name])^0.5) - (1/6)*(hot_prob[:ΔT_upper][cold.name] + hot_prob[:ΔT_lower][cold.name]))) # Necessary to avoid the parsing error.
end
=#
