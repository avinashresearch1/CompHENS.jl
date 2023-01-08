using Plasmo

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
#function generate_network!(prob::ClassicHENSProblem, EMAT, overall_network::Dict{String, AbstractSuperstructure}; obj_func::NetworkObjective = AreaArithmeticMean(), time_limit = 200.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    #verbose && @info "Solving the Network Generation subproblem"
    obj_func = AreaPaterson()
    haskey(prob.results_dict, :y) || error("Stream match data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :Q) || error("HLD data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :match_list) || error("Match list data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :T_bounds) || generate_temperature_bounds!(prob)

    #y, Q = prob.results_dict[:y],  prob.results_dict[:Q]
    
    HEN = OptiGraph()
    @optinode(HEN, streams[prob.all_names])

    # 1. This generates the stream-wise optimization model i.e., each stream gives rise to an OptiNode
    for stream in prob.all_names
        add_stream_variables!(prob.all_dict[stream], streams[stream], overall_network[stream], prob)
        add_stream_constraints!(prob.all_dict[stream], streams[stream], overall_network[stream], prob)
    end

    #2. Plasmo currently does not allow adding Variables or Objectives to the OptiGraph object (see issue). 
    # The general idea is to treat all the matches as belonging to the hot stream OptiNode. Thus, temperature differences, LMTDs, HX area and HX cost are all treated only on the hot side OptiNode.  
    # The following procedure is followed:
    # i. For every `hot` stream, and for every match: `ΔT_upper` and  `ΔT_lower` variables are declared. hot_node[:ΔT_upper][cold] denotes the ΔT_upper variable for the hot stream on its `cold` match. 
    # ΔT_upper is the temperature difference btn hot stream entrance to HX, and cold stream exit from HX. Both ΔT_upper and ΔT_lower constrained >= EMAT.
    # ii. `add_match_linkconstraints()` function is used to set the ΔT_upper values to be equal to the appropriate temperature difference.
    # iii. `add_stream_objective()` Adds an area objective to ONLY the hot streams and utilities NOT cold streams and objectives. This prevents double-counting area when the overall objective is taken as linear combination of `OptiNode` objectives.  
    # Depending on the type of the `LMTD_approx` kwarg, different approximations may be used for the LMTD and the objective. 
    for hot in prob.hot_names
        cold_match_list = prob.results_dict[:match_list][hot]
        @variable(streams[hot], ΔT_upper[cold_match_list] >= EMAT)
        @variable(streams[hot], ΔT_lower[cold_match_list] >= EMAT)
        @variable(streams[hot], match_objective[cold_match_list])
        for cold in cold_match_list
            add_match_linkconstraints!(prob.hot_dict[hot], prob.cold_dict[cold], HEN, streams[hot], streams[cold], overall_network[hot], overall_network[cold])
            set_match_objective!(prob.hot_dict[hot], prob.cold_dict[cold], streams[hot], obj_func, prob)  # Necessary because of bug in Plasmo in parsing generators.  
        end
        add_stream_objective!(prob.hot_dict[hot], streams[hot], prob)
    end

    @variable(streams["H1"], strm_obj_fn)
    @NLconstraint(streams["H1"], streams["H1"][:strm_obj_fn] == sum(((2/(streams["H1"][:ΔT_upper][cold] + streams["H1"][:ΔT_lower][cold]))*prob.results_dict[:Q][cold, "H1"]) for cold in prob.results_dict[:match_list]["H1"]))
    @objective(streams["H1"], Min, streams["H1"][:strm_obj_fn])    
    
    

    using MadNLP, MadNLPGraph, Plasmo
    MadNLP.optimize!(HEN; print_level=MadNLP.DEBUG, max_iter=100)
    using Ipopt
    #optimize with Ipopt
    set_optimizer(HEN,Ipopt.Optimizer)
    optimize!(HEN)
    
    optimize!(HEN)



    # Also: Adds `ΔT_upper` and `ΔT_lower` for each stream.
    #


# TRY JUNIPER.
# Need to return `streams` for bookkeeping.
#end

"""
$(TYPEDEF)
For a given `stream`, adds the `f` (stream flow rate through edge /mcp) and `t` temperature variables for each edge. `stream_prob` is usually attained from `Plasmo.@optinode`

"""
function add_stream_variables!(stream::AbstractStream, stream_prob::OptiNode, superstructure::AbstractSuperstructure, prob::ClassicHENSProblem, EMAT)
    @variable(stream_prob, prob.results_dict[:T_bounds][stream.name][1] <= t[superstructure.edges] <= prob.results_dict[:T_bounds][stream.name][2])
    @variable(stream_prob, 0 <= f[superstructure.edges] <= stream.mcp)
end

function add_stream_variables!(stream::AbstractUtility, stream_prob::OptiNode, superstructure::AbstractSuperstructure, prob::ClassicHENSProblem)
    # No `f` for utilities.
    @variable(stream_prob, 0 <= t[superstructure.edges])
end

"""
$(TYPEDEF)
Adds constraints for each `stream`.
"""
function add_stream_constraints!(stream::AbstractStream, stream_prob::OptiNode, superstructure::AbstractSplitSuperstructure, prob::ClassicHENSProblem)
    # 1. Source and Sink constraints on `t` and `f`
    SO_BS_edge = out_edges(superstructure.source[1], superstructure)[1]
    @constraint(stream_prob, stream_prob[:f][SO_BS_edge] == stream.mcp)
    @constraint(stream_prob, stream_prob[:t][SO_BS_edge] == stream.T_in)

    BM_SK_edge = in_edges(superstructure.sink[1], superstructure)[1]
    @constraint(stream_prob, stream_prob[:f][BM_SK_edge] == stream.mcp)
    @constraint(stream_prob, stream_prob[:t][BM_SK_edge] == stream.T_out)  

    # 2. Adding mass balance constraints for all nodes except `Source` and `Sink`
    internal_nodes = filter(superstructure.nodes) do v
        !isa(v, Source) && !isa(v, Sink)
    end
    for node in internal_nodes
        @constraint(stream_prob, sum(stream_prob[:f][edge] for edge in in_edges(node, superstructure)) == sum(stream_prob[:f][edge] for edge in out_edges(node, superstructure)))
    end
    
    # 3. Setting equal temperatures for all streams entering and leaving a splitter
    for splitter in superstructure.splitters
        edge_in = in_edges(splitter, superstructure)[1] # Only one in edge allowed
        for out_edge in out_edges(splitter, superstructure)
            @constraint(stream_prob, stream_prob[:t][out_edge] == stream_prob[:t][edge_in])
        end
    end

    # 4. Setting mixer energy balance
    for mixer in superstructure.mixers
        edge_out = out_edges(mixer, superstructure)[1] # Only one out edge allowed
        @NLconstraint(stream_prob, sum(stream_prob[:f][in_edge]*stream_prob[:t][in_edge] for in_edge in in_edges(mixer, superstructure)) == stream_prob[:f][edge_out]*stream_prob[:t][edge_out])
    end

    # 5. Heat Exchanger energy balance
    Q = prob.results_dict[:Q]
    for hx in superstructure.hxs
        edge_in = in_edges(hx, superstructure)[1] # Only one in edge allowed
        edge_out = out_edges(hx, superstructure)[1] # Only one out edge allowed
        match = hx.match
        if stream isa HotStream # Hot stream gets cooled, match is a ColdStream
            @NLconstraint(stream_prob, stream_prob[:f][edge_in]*stream_prob[:t][edge_in] - stream_prob[:f][edge_out]*stream_prob[:t][edge_out]  == Q[match, stream.name])
        elseif stream isa ColdStream # Cold stream gets heat, match is a HotStream
            @NLconstraint(stream_prob, stream_prob[:f][edge_out]*stream_prob[:t][edge_out] - stream_prob[:f][edge_in]*stream_prob[:t][edge_in] == Q[stream.name, match])
        else
            error("Stream type $(typeof(stream)) not supported.")
        end
    end
end

function add_stream_constraints!(stream::AbstractUtility, stream_prob::OptiNode, superstructure::ParallelSplit, prob::ClassicHENSProblem)
    # Only set edge temperatures for LMTD. All edge temperature before the HX are equal to Source temperature. All after HX are set at Sink temperature.   
    # Hack is to set all edges connected to MajorSplitter at same temperature as stream.T_in, and all to MajorMixer at same temperature as stream.T_out.
    major_splitter = filter(v -> (v isa MajorSplitter), superstructure.nodes)
    major_mixer = filter(v -> (v isa MajorMixer), superstructure.nodes)

    splitter_edges = union(in_edges(major_splitter[1], superstructure), out_edges(major_splitter[1], superstructure))
    mixer_edges = union(in_edges(major_mixer[1], superstructure), out_edges(major_mixer[1], superstructure))
    
    for edge in splitter_edges
        @constraint(stream_prob, stream_prob[:t][edge] == stream.T_in)
    end

    for edge in mixer_edges
        @constraint(stream_prob, stream_prob[:t][edge] == stream.T_out)
    end
    return
end

s_n, m_n = "H1", "C1"
hot, cold, overall_prob, hot_prob, cold_prob, hot_superstructure, cold_superstructure = prob.all_dict[s_n], prob.all_dict[m_n], HEN, streams[s_n], streams[m_n], overall_network[s_n], overall_network[m_n]  

"""
$(TYPEDEF) 
Function used to set the ΔT_upper and ΔT_lower for the `hot_prob`
"""
function add_match_linkconstraints!(hot::Union{HotStream, SimpleHotUtility}, cold::Union{ColdStream, SimpleColdUtility}, overall_prob::OptiGraph, hot_prob::OptiNode, cold_prob::OptiNode, hot_superstructure::AbstractSuperstructure, cold_superstructure::AbstractSuperstructure)
    hot_HX = filter(hot_superstructure.hxs) do v
         v.match == cold.name
    end[1] # Should only be 1. Add check?
    hot_hx_in_edge, hot_hx_out_edge  = in_edges(hot_HX, hot_superstructure)[1], out_edges(hot_HX, hot_superstructure)[1]

    # Match is always a cold stream or cold utility
    cold_HX = filter(cold_superstructure.hxs) do v
        v.match == hot.name
   end[1]
   cold_hx_in_edge, cold_hx_out_edge  = in_edges(cold_HX, cold_superstructure)[1], out_edges(cold_HX, cold_superstructure)[1]

   # Setting ΔT_upper
   @linkconstraint(overall_prob, hot_prob[:ΔT_upper][cold.name] ==  hot_prob[:t][hot_hx_in_edge] - cold_prob[:t][cold_hx_out_edge])

    # Setting ΔT_lower
    @linkconstraint(overall_prob, hot_prob[:ΔT_lower][cold.name] ==  hot_prob[:t][hot_hx_out_edge] - cold_prob[:t][cold_hx_in_edge])
end   

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

"""
$(TYPEDEF) 
Function used to set the objective of each hot stream problem. The objective of the network generation problem is a linear combination of objectives of hot stream problems.
"""
function add_stream_objective!(hot::Union{HotStream, SimpleHotUtility}, hot_prob::OptiNode, prob::ClassicHENSProblem)
    @objective(hot_prob, Min, sum(hot_prob[:match_objective]))
    #@NLobjective(hot_prob, Min, sum(((2/(hot_prob[:ΔT_upper][cold] + hot_prob[:ΔT_lower][cold]))*(1/(U(hot, cold)))*prob.results_dict[:Q][cold, hot.name]) for cold in prob.results_dict[:match_list][hot.name]))     
    # @NLobjective(hot_prob, Min, sum(prob.results_dict[:Q][cold, hot.name] for cold in prob.results_dict[:match_list][hot.name]))
    # Currently, with Plasmo need to set as constraint, since @NLobjective is not supported.
    #@variable(hot_prob, stream_obj_fn)
    #@constraint(hot_prob, hot_prob[:stream_obj_fn] == prob.results_dict[:Q][cold.name, hot.name])
    #@NLconstraint(hot_prob, hot_prob[:objective] == 
    # @objective(hot_prob, Min, hot_prob[:objective])
end

#superstructure = overall_network["ST"]
print_node(streams["C1"])
generate_temperature_bounds!(prob)

"""
$(TYPEDEF) 
Generates bounds on the temperatures for a stream. 
Creates a `prob.results_dict[:T_bounds]::Dict{String, Tuple}`, where the tuple is `(T_LBD, T_UBD)`
"""
function generate_temperature_bounds!(prob::ClassicHENSProblem)
    T_bounds = Dict{String, Tuple}()
    for (k,v) in prob.streams_dict
        push!(T_bounds, k => generate_temperature_bounds(prob, v))
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

"""
$(TYPEDEF)
Prints the optimization model for the stream
"""
function print_stream_prob(stream_prob::OptiNode)
    print(jump_model(stream_prob))
end
