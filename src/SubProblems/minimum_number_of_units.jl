using JuMP
using HiGHS

"""
$(TYPEDSIGNATURES)

Constructs and solves the MILP transshipment problem of Papoulias_Grossmann_1983 for the minimum number of HX units. 
Returns:
- sol:
"""
function solve_minimum_units_subproblem(prob::ClassicHENSProblem, sol_min_utils::MinUtilitiesSolution; time_limit = 60.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = true)
    @info "Solving the minimum number of units subproblem"
    ΔT_min_absolute = 0.0 # Use the minimum feasible here to get the absolute minimum number of units. 
    intervals = CompHENS.generate_heat_cascade_intervals(prob, ΔT_min_absolute)
    model = Model()
    H_set = keys(prob.hot_utilities_dict) 
    CU_set = keys(prob.cold_utilities_dict)

    @variable(model, 0 <= Q_in[HU_set])
    @variable(model, 0 <= Q_out[CU_set])
    @variable(model, 0 <= R[intervals]) # Notation: R[interval] is the residual heat exiting a given interval
    JuMP.fix(R[last(intervals)], 0.0; force = true)

    # First interval: Entering == Leaving
    @constraint(model, 
    sum(Q_in[hu] for hu in keys(first(intervals).hot_utilities_contribs)) + first(intervals).total_stream_heat_in == R[first(intervals)] + sum(Q_out[cu] for cu in keys(first(intervals).cold_utilities_contribs)) + first(intervals).total_stream_heat_out)

    # Remaining intervals
    @constraint(model, [i in 2:length(intervals)],
    R[intervals[i-1]] + sum(Q_in[hu] for hu in keys(intervals[i].hot_utilities_contribs)) + intervals[i].total_stream_heat_in == R[intervals[i]] + sum(Q_out[cu] for cu in keys(intervals[i].cold_utilities_contribs)) + intervals[i].total_stream_heat_out)

    # Objective: TODO: Add utility costs.
    @objective(model, Min, sum(Q_in) + sum(Q_out))
    set_optimizer(model, optimizer)
    presolve && set_optimizer_attribute(model, "presolve", "on")
    set_optimizer_attribute(model, "time_limit", time_limit)
    optimize!(model)
    if verbose
        @show termination_status(model)
        @show primal_status(model)
        @show dual_status(model)
    end


    # Post-processing
    sol = Dict()
    for hu in HU_set
        push!(sol, hu => value.(Q_in[hu]))
    end
    for cu in CU_set
        push!(sol, cu => value.(Q_out[cu]))
    end
    return sol
end
=#
