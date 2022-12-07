using JuMP
using HiGHS

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



"""
$(TYPEDSIGNATURES)

Constructs and solves the LP transshipment formulation of Papoulias_Grossmann_1983 to determine the minimum utility cost.
Returns:
- `sol::MinUtilitiesSolution` that contains the utility consumptions, total utility costs and the pinch points.  
[TODO: Utility costs]
"""
function solve_minimum_utilities_subproblem(prob::ClassicHENSProblem; time_limit = 60.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = true)
    @info "Solving the minimum utilities subproblem"
    intervals = CompHENS.generate_heat_cascade_intervals(prob)
    subprob = Model()
    HU_set = keys(prob.hot_utilities_dict) 
    CU_set = keys(prob.cold_utilities_dict)

    @variable(subprob, 0 <= Q_in[HU_set])
    @variable(subprob, 0 <= Q_out[CU_set])
    @variable(subprob, 0 <= R[intervals]) # Notation: R[interval] is the residual heat exiting a given interval
    JuMP.fix(R[last(intervals)], 0.0; force = true)

    # First interval: Entering == Leaving
    @constraint(subprob, 
    sum(Q_in[hu] for hu in keys(first(intervals).hot_utilities_contribs)) + first(intervals).total_stream_heat_in == R[first(intervals)] + sum(Q_out[cu] for cu in keys(first(intervals).cold_utilities_contribs)) + first(intervals).total_stream_heat_out)

    # Remaining intervals
    @constraint(subprob, [i in 2:length(intervals)],
    R[intervals[i-1]] + sum(Q_in[hu] for hu in keys(intervals[i].hot_utilities_contribs)) + intervals[i].total_stream_heat_in == R[intervals[i]] + sum(Q_out[cu] for cu in keys(intervals[i].cold_utilities_contribs)) + intervals[i].total_stream_heat_out)

    # Objective: TODO: Add utility costs.
    @objective(subprob, Min, sum(Q_in) + sum(Q_out))
    set_optimizer(subprob, optimizer)
    presolve && set_optimizer_attribute(subprob, "presolve", "on")
    set_optimizer_attribute(subprob, "time_limit", time_limit)
    optimize!(subprob)
    if verbose
        @show termination_status(subprob)
        @show primal_status(subprob)
        @show dual_status(subprob)
    end


    # Post-processing
    hot_utilities_consumption = Dict{String, Float64}()
    cold_utilities_consumption = Dict{String, Float64}()
    total_utility_cost = objective_value(subprob)
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
end