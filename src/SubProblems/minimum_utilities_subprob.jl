using JuMP
using HiGHS

"""
$(TYPEDSIGNATURES)

Constructs and solves the LP transshipment formulation of Papoulias_Grossmann_1983 to determine the minimum utility cost.
Returns:
"""
function solve_minimum_utilities_subproblem!(prob::ClassicHENSProblem; time_limit = 60.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    verbose && @info "Solving the minimum utilities subproblem"
    intervals = CompHENS.generate_heat_cascade_intervals(prob)
    model = Model()
    HU_set = keys(prob.hot_utilities_dict) 
    CU_set = keys(prob.cold_utilities_dict)

    @variable(model, 0 <= Q_in[HU_set])
    @variable(model, 0 <= Q_out[CU_set])
    @variable(model, 0 <= R[intervals]) # Notation: R[interval] is the residual heat exiting a given interval
    JuMP.fix(R[last(intervals)], 0.0; force = true)

    # First interval: Entering == Leaving
    @constraint(model, 
    sum(Q_in[hu] for hu in keys(first(intervals).hot_utilities_contribs)) + first(intervals).total_stream_heat_in == R[first(intervals)] + sum(Q_out[cu] for cu in keys(first(intervals).cold_utilities_contribs)) + first(intervals).total_stream_heat_out)

    # Remaining intervals
    @constraint(model, [k in 2:length(intervals)],
    R[intervals[k-1]] + sum(Q_in[hu] for hu in keys(intervals[k].hot_utilities_contribs)) + intervals[k].total_stream_heat_in == R[intervals[k]] + sum(Q_out[cu] for cu in keys(intervals[k].cold_utilities_contribs)) + intervals[k].total_stream_heat_out)

    # Objective: TODO: Add utility costs.
    @objective(model, Min, sum(Q_in) + sum(Q_out))
    set_optimizer(model, optimizer)
    !verbose && set_silent(model)
    presolve && set_optimizer_attribute(model, "presolve", "on")
    set_optimizer_attribute(model, "time_limit", time_limit)
    optimize!(model)
    if verbose
        @show termination_status(model)
        @show primal_status(model)
        @show dual_status(model)
    end


    # Post-processing
    pinch_points = Tuple[]

    for (k,v) in prob.hot_utilities_dict
        prob.hot_utilities_dict[k].Q = value.(Q_in[k])
    end

    for (k,v) in prob.cold_utilities_dict
        prob.cold_utilities_dict[k].Q = value.(Q_out[k])
    end

    for interval in setdiff(intervals, [last(intervals)])
        if value.(R[interval]) <= smallest_value
            push!(pinch_points, (interval.T_hot_lower, interval.T_cold_lower))
        end
    end
    prob.results_dict[:pinch_points] = pinch_points
    return
end


#=
CODE deprecated in favor of holding results in Problem struct itself.
"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds the solution to the minimum utilities subproblem.
[QN:] Is it better to mutate the `prob` to store the results?
"""
mutable struct MinUtilitiesSolution{T}  <: AbstractSubProblemSolution 
    hot_utilities_consumption::Dict{String, T}
    cold_utilities_consumption::Dict{String, T}
    total_utility_cost::T
    pinch_points::Vector{Tuple}
end

RELEVANT CODE IN SOL MIN UTILITIES PROBLEM
 # Post-processing
    hot_utilities_consumption = Dict{String, Float64}()
    cold_utilities_consumption = Dict{String, Float64}()
    total_utility_cost = objective_value(model)
    pinch_points = Tuple[]
    for hu in HU_set
        push!(hot_utilities_consumption, hu => value.(Q_in[hu]))
    end
    for cu in CU_set
        push!(cold_utilities_consumption, cu => value.(Q_out[cu]))
    end

    for interval in setdiff(intervals, [last(intervals)])
        if value.(R[interval]) <= smallest_value
            push!(pinch_points, (interval.T_hot_lower, interval.T_cold_lower))
        end
    end
    return MinUtilitiesSolution(hot_utilities_consumption, cold_utilities_consumption, total_utility_cost, pinch_points)
=#
