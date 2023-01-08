using Plasmo

"""
$(TYPEDSIGNATURES)

Generates the Heat Exchanger Network. Define a type of superstructure for each stream. 
"""
#function generate_network!(prob::ClassicHENSProblem, EMAT, overall_network::Dict{String, AbstractSuperstructure}; time_limit = 200.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    #verbose && @info "Solving the Network Generation subproblem"
    
    haskey(prob.results_dict, :y) || error("Stream match data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :Q) || error("HLD data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :match_list) || error("Match list data not available. Solve corresponding subproblem first.")

    #y, Q = prob.results_dict[:y],  prob.results_dict[:Q]
    
    HEN = OptiGraph()
    @optinode(HEN, streams[prob.all_names])

    for stream in prob.all_names
        add_stream_variables!(prob.all_dict[stream], streams[stream], overall_network[stream])
        add_stream_constraints!(prob.all_dict[stream], streams[stream], overall_network[stream], prob)
    end


# TRY JUNIPER.
# Need to return `streams` for bookkeeping.
#end

"""
$(TYPEDEF)
For a given `stream`, adds the `f` (stream flow rate through edge /mcp) and `t` temperature variables for each edge. `stream_prob` is usually attained from `Plasmo.@optinode`
[TODO:] Add good bounds.
"""
function add_stream_variables!(stream::AbstractStream, stream_prob::OptiNode, superstructure::AbstractSuperstructure)
    @variable(stream_prob, 0 <= t[superstructure.edges])
    @variable(stream_prob, 0 <= f[superstructure.edges] <= stream.mcp)
end

function add_stream_variables!(stream::AbstractUtility, stream_prob::OptiNode, superstructure::AbstractSuperstructure)
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
    @constraint(stream_prob, stream_prob[:t][SO_BS_edge] == stream.T_out)  

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







"""
$(TYPEDEF)
Prints the optimization model for the stream
"""
function print_node(stream_node::OptiNode)
    print(jump_model(stream_node))
end
 
    

#superstructure = overall_network["ST"]
print_node(streams["H1"])
