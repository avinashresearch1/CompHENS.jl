using Plasmo

"""
$(TYPEDSIGNATURES)

Generates the Heat Exchanger Network.

"""
function generate_network!(prob::ClassicHENSProblem, EMAT, superstructure::FloudasCiricGrossmann; level = :quaternary_temperatures, add_units = 1, time_limit = 200.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    verbose && @info "Solving the Network Generation subproblem"
    
    haskey(prob.results_dict, :y) || error("Stream match data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :Q) || error("HLD data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :match_list) || error("Match list data not available. Solve corresponding subproblem first.")

    y, Q = prob.results_dict[:y],  prob.results_dict[:Q]
    
    H_stream_set = merge(prob.hot_utilities_dict, prob.hot_streams_dict)
    H_set = keys(H_stream_set)
    C_stream_set = merge(prob.cold_utilities_dict, prob.cold_streams_dict)
    C_set = keys(C_stream_set)
    stream_set = union(H_set, C_set)

    HEN = OptiGraph()
    @optinode(HEN, streams[stream_set])

    full_superstructure = Dict{String, FloudasCiricGrossmann}()
    # Note: Each element of streams does not have the stream name cleanly. So using axis labels.
    for stream in streams.axes[1]
        push!(full_superstructure, stream => FloudasCiricGrossmann(stream, prob; verbose = false))

        edges = full_superstructure[stream].edges

        # Declare stream-wise variables for all edges
        @variable(streams[stream], 0 <= f[edges]) # `f` denotes the stream flow rate through edge (mcp)
        @variable(streams[stream], 0 <= t[edges]) # `t` denotes the temperatures. [TODO:] Get good bounds!
    end




end


"""
 Creates the edge `t` and `f` variables for each `stream`
"""
function create_edge_variables(stream::AbstractStream, streams, superstructure::AbstractSuperstructure)
end



