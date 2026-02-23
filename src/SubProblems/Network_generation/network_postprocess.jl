using DataFrames

function print_stream_results(stream::String, prob::ClassicHENSProblem, model::AbstractModel, superstructure::AbstractSplitSuperstructure; digits = 1)
    for edge in superstructure.edges
        t_val = round(value(model[:t][(stream, edge)]); digits = digits)
        println("$(edge)")
        print("T: $(t_val)      ")
        if prob.all_dict[stream] isa AbstractStream
            f_val = round(value(model[:f][(stream, edge)]); digits = digits)
            print(" f: $(f_val)")
        end
        print("\n \n")
    end
    return
end

"""
$(TYPEDSIGNATURES)

Postprocessing after solving stream generation subproblem. 
Displays the matches and heat load distribution in a 2-D matrix form, maintains stream ordering.
"""
function postprocess_network!(prob::ClassicHENSProblem, model::AbstractModel, HLD_list, overall_network; results_path = nothing, visualize = true, digits = 4, display = true)
    results_df = DataFrame()
    # Persist per-match area values so `get_design_area(prob)` remains valid
    # for all example workflows after network generation.
    areas = Dict{Tuple{String, String}, Float64}()
    match = HLD_list[1]
    for match in HLD_list
        hot, cold = match[1], match[2]
        hot_superstructure, cold_superstructure = overall_network[hot], overall_network[cold]
        hot_HX = only(filter(hot_superstructure.hxs) do v
            v.match == cold
        end)
        hot_hx_in_edge, hot_hx_out_edge  = only(in_edges(hot_HX, hot_superstructure)), only(out_edges(hot_HX, hot_superstructure))
        
        # Match is always a cold stream or cold utility
        cold_HX = only(filter(cold_superstructure.hxs) do v
           v.match == hot
        end)
        cold_hx_in_edge, cold_hx_out_edge  = only(in_edges(cold_HX, cold_superstructure)), only(out_edges(cold_HX, cold_superstructure))
   
        
        ΔT_upper = value(model[:ΔT_upper][match]) 
        ΔT_lower = smallest_value + value(model[:ΔT_lower][match])
        LMTD = (ΔT_upper - ΔT_lower)/(smallest_value + log(ΔT_upper/ΔT_lower)) # Probable numerical issues.
        area = prob.results_dict[:Q][match[2], match[1]]/(LMTD*U(prob.all_dict[hot], prob.all_dict[cold]))
        areas[(hot, cold)] = area
        push!(results_df, (Hot_stream = hot, Q = prob.results_dict[:Q][match[2], match[1]], Tin_Hot = value(model[:t][(hot, hot_hx_in_edge)]), Tout_Hot = value(model[:t][(hot, hot_hx_out_edge)]), Cold_stream = cold, Tin_Cold = value(model[:t][(cold, cold_hx_in_edge)]), Tout_Cold = value(model[:t][(cold, cold_hx_out_edge)]), ΔT_upper = ΔT_upper, ΔT_lower = ΔT_lower, LMTD = LMTD, U = U(prob.all_dict[hot], prob.all_dict[cold]), Area = area))
    end
    prob.results_dict[:areas] = areas
    return results_df
end

"""
$(TYPEDSIGNATURES)

Gets the total designed area. Can only be called after solving `network_generator` subproblem.
"""
function get_design_area(prob::ClassicHENSProblem)
    haskey(prob.results_dict, :areas) || error("Network generation needs to be solved first")
    return sum(values(prob.results_dict[:areas]))
end

function get_design_area(prob::MultiPeriodFlexibleHENSProblem)
    area = Dict()
    
    for t in prob.period_names
        println("PROBLEM $(t)")
        push!(area, t => get_design_area(prob.period_streams_dict[t]))   
        println("Operational Area: $(prob.period_streams_dict[t])")
        print("\n")
    end
    return area
end
