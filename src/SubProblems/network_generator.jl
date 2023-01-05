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
    construct_superstructure!(streams, prob, superstructure)



end


"""
$(TYPEDSIGNATURES)

Constucts the streamwise `superstructure::AbstractSuperstructure` for all `streams`. `streams` is usually a `DenseAxisArray`.
"""
function construct_superstructure!(streams, prob::ClassicHENSProblem, superstructure::FloudasCiricGrossmann)
    # Note: Each element of streams does not have the stream name cleanly. So using axis labels; not extensible.
    for stream in streams.axes[1] # Any better option?
        edge_set = [("SO", "BS"), ("BM", "SK")] # Add the first Source -> Big Splitter, and last Big Mixer -> Sink edge 
        for match in prob.results_dict[:match_list][stream]
            push!(edge_set, ("BS", "SM_$(match)")) # Big Splitter to Small Mixer
            push!(edge_set, ("SM_$(match)", "HX_$(match)")) # Small Mixer to HX
            push!(edge_set, ("HX_$(match)", ""))


        add_stream_nodes!()
    end



end

"""
$(TYPEDSIGNATURES)

Add superstructure nodes for each stream
"""
function add_stream_nodes!(stream, prob::ClassicHENSProblem, superstructure::FloudasCiricGrossmann)
    edge_set = [("SO", "BS"), ("BM", "SK")] # Add the first Source -> Big Splitter, and last Big Mixer -> Sink edge 
    for matches in prob.results_dict[:match_list][]
    @variable(stream, f)




- Major stream splitter `BS` (Big Splitter)
- Minor stream mixers `SM` (Small Mixer) - one in each branch corresponding to a match
- Heat Exchanger `HX` - one in each branch corresponding to a match
- Minor stream splitter `SS` (Small Splitter) - one in each branch corresponding to a match
- Major stream mixer `BM` (Big Mixer)
- Sink `SK` 
end
