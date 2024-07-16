#[WIP]
using CompHENS
pk_dir = pkgdir(CompHENS)
include(joinpath(pk_dir, "Examples", "XLSX_interface", "ClassicHENSProblem", "Gundersen_4_stream", "Gundersen_4_stream.jl"))
exportall(CompHENS)

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
function test_terminal_match(stream_name::String, matching_stream_name::String, EMAT, overall_network)
    test_terminal_match(prob.all_dict[stream_name], prob.all_dict[matching_stream_name], EMAT, overall_network)
end

stream_name = "C2"
matching_stream_name = "ST"

cold_stream = prob.all_dict[stream_name]
hot_matching_stream = prob.all_dict[matching_stream_name]

function test_terminal_match(cold_stream::Union{ColdStream}, hot_matching_stream::Union{HotStream, SimpleHotUtility}, EMAT, overall_network)

    Q = prob.results_dict[:Q][cold_stream.name, hot_matching_stream.name]

    hot_T_out = overall_network[hot_matching_stream.name].metadata["Working Terminal T"]
    if hot_matching_stream isa SimpleHotUtility
        hot_T_in = hot_matching_stream.T_in
    else
        hot_T_in = hot_T_out + Q/hot_matching_stream.mcp
    end

    cold_T_out = overall_network[cold_stream.name].metadata["Working Terminal T"]
    cold_T_in = cold_T_out - Q/cold_stream.mcp

    test_terminal_match(hot_T_out, hot_T_in, cold_T_out, cold_T_in, EMAT)
end

stream_name = "H2"
matching_stream_name = "CW"

hot_stream = prob.all_dict[stream_name]
cold_matching_stream = prob.all_dict[matching_stream_name]

function test_terminal_match(hot_stream::Union{HotStream}, cold_matching_stream::Union{ColdStream, SimpleColdUtility}, EMAT, overall_network)
    Q = prob.results_dict[:Q][cold_matching_stream.name, hot_stream.name]

    hot_T_out = overall_network[hot_stream.name].metadata["Working Terminal T"]
    hot_T_in = hot_T_out + Q/hot_stream.mcp
    
    cold_T_out = overall_network[cold_matching_stream.name].metadata["Working Terminal T"]
    if cold_matching_stream isa SimpleColdUtility
        cold_T_in = cold_matching_stream.T_in
    else
        cold_T_in = cold_T_out - Q/cold_matching_stream.mcp
    end

    test_terminal_match(hot_T_out, hot_T_in, cold_T_out, cold_T_in, EMAT)
end

## Main while loop for attaining feasible initial values.
# This will be mutated in the initialization loop
working_HX_list = deepcopy(prob.results_dict[:HX_list])

# Initialize all the Working Terminal Temperatures to T_out: 
# Also 
for (k,v) in overall_network
    push!(v.metadata, "Working Terminal Node" => only(v.major_mixer))
    push!(v.metadata, "Working Terminal T" => prob.all_dict[k].T_out)

    # Set the temperatures of the source -> BS edge and the BM -> Sink edge
    superstructure = overall_network[k]
    SO_BS_edge = only(out_edges(superstructure.source[1], superstructure))
    set_start_value(model[:t][(k, SO_BS_edge)], prob.all_dict[k].T_in)
    
    BM_SK_edge = only(in_edges(superstructure.sink[1], superstructure))
    set_start_value(model[:t][(k, BM_SK_edge)], prob.all_dict[k].T_out)

    if prob.all_dict[k] isa Union{HotStream, ColdStream} # Utilities don't have f streams
        set_start_value(model[:f][(k, SO_BS_edge)], prob.all_dict[k].mcp)
        set_start_value(model[:f][(k, BM_SK_edge)], prob.all_dict[k].mcp)
    end
end

only(overall_network["C2"].metadata["Working Terminal Node"] )

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

    (; terminal_feasible, hot_T_out, hot_T_in, cold_T_out, cold_T_in) = test_terminal_match(cold_stream,hot_util, EMAT, overall_network)
    if terminal_feasible
        set_serial_terminal_starting_points!(model, hot_T_out, hot_T_in, cold_T_out, cold_T_in) 
    end



superstructure = overall_network["C2"]

# TODO: 0: Set these inits. Also push all inits to a dict for plotting. 
k = "C2"
prob.all_dict[k].T_in

stream = prob.all_dict[k]
superstructure = overall_network[k]
matching_stream = "ST"
terminal_match_HX = only(filter(k -> k.match == matching_stream, superstructure.hxs))
new_terminal_node = only(get_source_nodes(terminal_match_HX, superstructure))

"""
This function is to be called if a HX is chosen to be placed (by setting starting points) as the last serial HX of the `stream` OR prior to another HX which has also been previously placed as a serial terminal HX.
This function sets starting points such a strictly serial connection is made between the `new_terminal_node` and the `stream`'s current terminal node which can be queried by `superstructure.metadata["Working Terminal Node"]`.
Lastly, the `superstructure.metadata["Working Terminal Node"]` is set to the `new_terminal_node` and the `superstructure.metadata["Working Terminal T"]` is set to the new temperature.
"""
function set_serial_terminal_starting_points!(model, new_terminal_node::Node, stream::ColdStream, superstructure::AbstractSplitSuperstructure; hot_T_out, hot_T_in, cold_T_out, cold_T_in)
end

@assert new_terminal_node isa MinorMixer || error("The new terminal node must be a MinorMixer") 
# 1. Set the starting points for the edge from MinorMixer to HX
MM_HX_edge = only(out_edges(new_terminal_node, superstructure))

function set_serial_terminal_starting_points!(model, new_terminal_node::Node, stream::HotStream, superstructure::AbstractSplitSuperstructure; hot_T_out, hot_T_in, cold_T_out, cold_T_in)
end





#BM_SK_edge = only(in_edges(superstructure.sink[1], superstructure))
#@constraint(model, model[:f][(stream.name, BM_SK_edge)] == stream.mcp)
#@constraint(model, model[:t][(stream.name, BM_SK_edge)] == stream.T_out)  


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



