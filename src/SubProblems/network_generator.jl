using Plasmo

overall_network = construct_superstructure()



"""
$(TYPEDSIGNATURES)

Generates the Heat Exchanger Network. Define a type of superstructure for each stream. 
"""
#function generate_network!(prob::ClassicHENSProblem, EMAT, overall_network::Dict{String, AbstractSuperstructure}; time_limit = 200.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    verbose && @info "Solving the Network Generation subproblem"
    
    haskey(prob.results_dict, :y) || error("Stream match data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :Q) || error("HLD data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :match_list) || error("Match list data not available. Solve corresponding subproblem first.")

    y, Q = prob.results_dict[:y],  prob.results_dict[:Q]
    
    HEN = OptiGraph()
    @optinode(HEN, streams[stream_set])

    for stream in stream_set
        stream_struc =  FloudasCiricGrossmann(stream, prob; verbose = false)
        edges = stream_struc.edges
        push!(overall_network, stream => stream_struc)
        

        # Declare variables stream-wise for all edges
        @variable(streams[stream], 0 <= t[edges]) # `t` denotes the temperatures. [TODO:] Get good bounds!
        if all_stream_set[stream] isa AbstractStream # f only for streams not utilities.
            @variable(streams[stream], 0 <= f[edges]) # `f` denotes the stream flow rate through edge (mcp)
        end

        #add_constraints(all_stream_set[stream], streams, stream_struc)

    end


# TRY JUNIPER.

#end

"""
$(TYPEDEF)
For a given `stream`, adds the `f` (stream flow rate through edge /mcp) and `t` temperature variables for each edge. `stream_prob` is usually attained from `Plasmo.@optinode`
[TODO:] Add good temperature bounds.
"""
function add_stream_variables!(stream::AbstractStream, stream_prob::OptiNode, superstructure::AbstractSuperstructure)
    @variable(stream_prob, 0 <= t[superstructure.edges])
    @variable(stream_prob, 0 <= f[superstructure.edges])
end

function add_stream_variables!(stream::AbstractUtility, stream_prob::OptiNode, superstructure::AbstractSuperstructure)
    # No `f` for utilities.
    @variable(stream_prob, 0 <= t[superstructure.edges])
end


"""
$(TYPEDEF)
Adds constraints for each `stream`. Better not to use MD here. 
"""
function add_stream_constraints!(stream::AbstractStream, streams, stream_struc::AbstractSplitSuperstructure)
    internal_nodes = filter(stream_struc.nodes) do v
        !isa(v, Source) && !isa(v, Sink)
    end

    #=
    # 1. Adding mass balance constraints for all nodes except `Source` and `Sink`
    for node in internal_nodes
        @constraint(streams[stream.name], sum(streams[stream.name][:f][edge] for edge in in_edges(node, stream_struc)) == sum(streams[stream.name][:f][edge] for edge in out_edges(node, stream_struc)))
    end
    
    # 2. Setting equal temperatures for all streams entering and leaving a splitter
    for splitter in values(stream_struc.splitters)
        length(in_edges(splitter, stream_struc)) == 1 || error("This function only works with superstructures where each splitter has a single incoming edge")
        edge_in = in_edges(splitter, stream_struc)[1] # Only one in edge allowed
        for out_edge in out_edges(splitter, stream_struc)
            @constraint(streams[stream.name], streams[stream.name][:t][out_edge] == streams[stream.name][:t][edge_in])
        end
    end
    =#

    # 3. Setting mixer energy balance
    for mixer in values(stream_struc.mixers)
        length(out_edges(mixer, stream_struc)) == 1 || error("This function only works with superstructures where each mixer has a single outgoing edge")
        edge_out = out_edges(mixer, stream_struc)[1] # Only one out edge allowed
        @NLconstraint(streams[stream.name], sum(streams[stream.name][:f][in_edge]*streams[stream.name][:t][in_edge] for in_edge in in_edges(mixer, stream_struc)) == streams[stream.name][:f][edge_out]*streams[stream.name][:t][edge_out])
    end


end


stream_struc = superstructure.streams["H1"]
add_stream_constraints!(all_stream_set["H1"], streams, superstructure.streams["H1"])
print(jump_model(streams["H1"]))
    





#end

#=
"""
 Creates the edge `t` and `f` variables for each `stream`. `streams` is a `DenseAxisArray` usually attained from `Plasmo.@optinode`. 
"""
function create_edge_variables(stream::AbstractStream, streams, superstructure::AbstractSuperstructure)

end
=#


