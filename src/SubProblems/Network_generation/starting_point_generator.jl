#[WIP]
using CompHENS
pk_dir = pkgdir(CompHENS)
include(joinpath(pk_dir, "Examples", "XLSX_interface", "ClassicHENSProblem", "Gundersen_4_stream", "Gundersen_4_stream.jl"))
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
    push!(v.add_user_data, "Working Terminal T" => v.T_out)
end

#TODO:

function test_terminal_match(hot_T_out, hot_T_in, cold_T_out, cold_T_in, EMAT)
    terminal_feasible = false
    if (hot_T_out >= cold_T_in + EMAT) && (hot_T_in >= cold_T_out + EMAT)
        terminal_feasible = true
    end 
    return (; terminal_feasible, hot_T_out, hot_T_in, cold_T_out, cold_T_in)
end

"""
Test if the `matching_stream` can be feasibly placed as the terminal serial heat exchanger of the `stream`
"""
function test_terminal_match(stream_name::String, matching_stream_name::String, EMAT)
    test_terminal_match(prob.all_dict[stream_name], prob.all_dict[matching_stream_name], EMAT)
end

stream_name = "C2"
matching_stream_name = "ST"

cold_stream = prob.all_dict[stream_name]
hot_matching_stream = prob.all_dict[matching_stream_name]

function test_terminal_match(cold_stream::Union{ColdStream}, hot_matching_stream::Union{HotStream, SimpleHotUtility}, EMAT)

    Q = prob.results_dict[:Q][cold_stream.name, hot_matching_stream.name]

    hot_T_out = hot_matching_stream.add_user_data["Working Terminal T"]
    if hot_matching_stream isa SimpleHotUtility
        hot_T_in = hot_matching_stream.T_in
    else
        hot_T_in = hot_T_out + Q/hot_matching_stream.mcp
    end

    cold_T_out = cold_stream.add_user_data["Working Terminal T"]
    cold_T_in = cold_T_out - Q/cold_stream.mcp

    test_terminal_match(hot_T_out, hot_T_in, cold_T_out, cold_T_in, EMAT)
end

stream_name = "H2"
matching_stream_name = "CW"

hot_stream = prob.all_dict[stream_name]
cold_matching_stream = prob.all_dict[matching_stream_name]

function test_terminal_match(hot_stream::Union{HotStream}, cold_matching_stream::Union{ColdStream, SimpleColdUtility}, EMAT)
    Q = prob.results_dict[:Q][cold_matching_stream.name, hot_stream.name]

    hot_T_out = hot_stream.add_user_data["Working Terminal T"]
    hot_T_in = hot_T_out + Q/hot_stream.mcp
    
    cold_T_out = cold_matching_stream.add_user_data["Working Terminal T"]
    if cold_matching_stream isa SimpleColdUtility
        cold_T_in = cold_matching_stream.T_in
    else
        cold_T_in = cold_T_out - Q/cold_matching_stream.mcp
    end

    test_terminal_match(hot_T_out, hot_T_in, cold_T_out, cold_T_in, EMAT)
end

# This will be mutated in the initialization loop
working_HX_list = deepcopy(prob.results_dict[:HX_list])

# 1. Place utilies at edges
# 1.A Hot utilities: 
hot_utilities = filter(working_HX_list) do (k,v) 
    prob.all_dict[k] isa SimpleHotUtility
end

hot_util = only(keys(hot_utilities))
cold_stream = only(working_HX_list[hot_util]) # Not gen.
#for hot_util in hot_utilities
    for cold_stream in working_HX_list[hot_util]
        @show cold_stream
    end

    (; terminal_feasible, hot_T_out, hot_T_in, cold_T_out, cold_T_in) = test_terminal_match(cold_stream,hot_util, EMAT)
    if terminal_feasible
        set_terminal_match_inits!(model, hot_T_out, hot_T_in, cold_T_out, cold_T_in) 
    end

function set_terminal_match_inits!(model, stream_name::String, matching_stream_name::String; hot_T_out, hot_T_in, cold_T_out, cold_T_in)
end




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
    # A. Set starting values for the temperatures:
    set_start_value(model[:t][(stream.name, only(out_edges(only(superstructure.source), superstructure)))], 69.0)

    # B. Set starting values for the temperatures:
    set_start_value(model[:t][(stream.name, only(out_edges(only(superstructure.source), superstructure)))], 69.0)



