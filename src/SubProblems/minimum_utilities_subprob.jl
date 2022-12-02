using JuMP
using HiGHS

"""
$(TYPEDSIGNATURES)

Constructs and solves a minimum utilities optimization problem. 
Returns:
- sol: A dictionary mapping the `keys(prob.hot_utilities_dict)`, `keys(prob.cold_utilities_dict)` to their minimum utility requirements.
Supports multiple utilities if appropriately matched to interval.
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
    sol = Dict()
    for hu in HU_set
        push!(sol, hu => value.(Q_in[hu]))
    end
    for cu in CU_set
        push!(sol, cu => value.(Q_out[cu]))
    end
    return sol
end