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
function postprocess_network!(prob::ClassicHENSProblem, model::AbstractModel, HLD_list; visualize = true, digits = 4, display = true)
    area_dict = Dict()
    for match in HLD_list
        hot, cold = match[1], match[2]
        ΔT_upper = value(model[:ΔT_upper][match]) 
        ΔT_lower = smallest_value + value(model[:ΔT_lower][match])
        LMTD = (ΔT_upper - ΔT_lower)/(smallest_value + log(ΔT_upper/ΔT_lower)) # Probable numerical issues.
        area = prob.results_dict[:Q][match[2], match[1]]/(LMTD*U(prob.all_dict[hot], prob.all_dict[cold]))
        push!(area_dict, match => area)
    end
    prob.results_dict[:areas] = area_dict
    return
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

