#[WIP]

#exportall(CompHENS)

#verbose && @info "Solving the Network Generation subproblem"

haskey(prob.results_dict, :y) || error("Stream match data not available. Solve corresponding subproblem first.")
haskey(prob.results_dict, :Q) || error("HLD data not available. Solve corresponding subproblem first.")
haskey(prob.results_dict, :HLD_list) || error("Match list data not available. Solve corresponding subproblem first.")
haskey(prob.results_dict, :T_bounds) || CompHENS.generate_temperature_bounds!(prob)

#y, Q = prob.results_dict[:y],  prob.results_dict[:Q]

model = Model()

# Tuple set with stream and utility s and Edge e
all_e_tuple_vec = []
for stream in prob.all_names
    for edge in overall_network[stream].edges
        push!(all_e_tuple_vec, (stream, edge))
    end
end

# Tuple set with only streams (not utilities) and edges e
stream_e_tuple_vec = []
for stream in prob.stream_names
    for edge in overall_network[stream].edges
        push!(stream_e_tuple_vec, (stream, edge))
    end
end

# 1. Declaring the stream-wise variables
@variable(model, 0.0 <= t[all_e_tuple_vec])
@variable(model, 0.0 <= f[stream_e_tuple_vec])

# Initialize all the Working Terminal Temperatures to T_out: 
for (k,v) in prob.all_dict
    v = prob.all_dict["H2"]
    push!(v.add_user_data, "Working Terminal T" => v.T_out)
end

#TODO:

"""
Test if the `matching_stream` can be feasibly placed as the terminal serial heat exchanger of the `stream`
"""
function test_terminal_match(stream_name::String, matching_stream_name::String)
end

function _test_terminal_match(stream::Union{ColdStream}, matching_stream::Union{HotStream, SimpleHotUtility})
end

stream_name = "C2"
matching_stream_name = "ST"

stream = prob.all_dict[stream_name]
matching_stream = prob.all_dict[matching_stream_name]
 

function get_utility_match_start_vals()
end
    utils_dict = merge(prob.cold_utilities_dict, prob.hot_utilities_dict)
    for (k,v) in utils_dict
        @show k
    end
    (k,v) = first(utils_dict)
    matched_streams = prob.results_dict[:HX_list][k]
    for stream in matched_streams
        @show stream
    end

#function set_start_values!(prob, EMAT, overall_network; verbose = verbose)
    # 1. Initialize the required dictionaries for book-keeping:
    prob.results_dict[:unflagged_matches] = deepcopy(filter(prob.results_dict[:HX_list]) do (k,v)
        k in prob.stream_names
    end
    )

    prob.results_dict[:working_terminal_temp] = Dict(stream => prob.streams_dict[stream].T_out for stream in prob.stream_names)
    prob.results_dict[:working_terminal_node] = Dict(stream => only(overall_network[stream].major_mixer) for stream in prob.stream_names)



stream = prob.all_dict["C1"]
superstructure = overall_network["C1"]
"""
$(TYPEDSIGNATURES)
Used to set the starting values of the network generation NLP.
This specifies that a given `match` HX should be placed serially at the working terminal of a given `stream`.
Initially, the `working_terminal_node` of a stream is its `MajorMixer`. If a `HX` match is set to be at the terminal, then the `working_terminal_node` is changed to the mixer feeding this `HX`.  
The, `working_terminal_temp` is changed from `T_out` and calculated for each added terminal match. 
"""
#function set_terminal_match!(model::AbstractModel, prob, stream::ColdStream, superstructure::AbstractSplitSuperstructure, match ; verbose = true)
    # A. Set starting values for the flows:
    set_start_value(model[:t][(stream.name, only(out_edges(only(superstructure.source), superstructure)))], 69.0)



