using JuMP
using HiGHS

"""
$(TYPEDSIGNATURES)

Constructs and solves the MILP transshipment problem of Papoulias_Grossmann_1983 for the minimum number of HX units. 
Returns:
"""
function solve_minimum_units_subproblem!(prob::ClassicHENSProblem; time_limit = 60.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    verbose && @info "Solving the minimum number of units subproblem"
    # Code design: One alternative here is to define the residual, R variables as `SparseAxisArrays` in JuMP. 
    # However, https://discourse.julialang.org/t/behaviour-of-sparse-array-variables-in-jump/36185/4 suggests instead fixing them to 0.0 and letting the pre-solver take care of them, so doing this.
    # Memory allocation is generally not the issue. Note, that since the stream contributions are known, the presolver should automatically fix them to be 0.0 and get rid of irrelevant residuals. 

    ΔT_min_absolute = 0.0 # Use the minimum feasible here to get the absolute minimum number of units. 
    intervals = generate_transshipment_intervals(prob, ΔT_min_absolute)
    model = Model()
    H_set = union(keys(prob.hot_utilities_dict), keys(prob.hot_streams_dict))
    C_set = union(keys(prob.cold_utilities_dict), keys(prob.cold_streams_dict))
    
    @variable(model, y[H_set, C_set], Bin)
    @variable(model, 0 <= R[H_set, intervals]) # Notation: R[interval] is the residual heat exiting a given interval
    @variable(model, 0 <= Q[H_set, C_set, intervals]) # Heat exchanged by hot stream i to cold stream j in interval k
    for i in H_set
        JuMP.fix(R[i, last(intervals)], 0.0; force = true)
    end
    
    # Getting bounds. TODO: Cross-check this.
    M = Dict()
    for i in H_set
        for j in C_set
            M[(i,j)] = min(sum(get_contribution(i,k.hot_side) for k in intervals), sum(get_contribution(j,k.cold_side) for k in intervals))
        end
    end
    
    # Hot side heat balance
    # First interval: Heat balance for each hot stream. Presolver sets residuals to zero. Entering == Leaving.
    @constraint(model, [i in H_set],
    get_contribution(i, first(intervals).hot_side) == R[i, first(intervals)] + sum(Q[i, j, first(intervals)] for j in C_set))
    
    # Remaining intervals. k used as an int counter here for convenience to avoid issues with reverse iteration.
    @constraint(model, [i in H_set, k in 2:length(intervals)], 
    R[i, intervals[k-1]] + get_contribution(i, intervals[k].hot_side) == R[i, intervals[k]] + sum(Q[i, j, intervals[k]] for j in C_set))
    
    # Cold side heat balance. 
    @constraint(model, [j in C_set, k in intervals],
    get_contribution(j, k.cold_side) == sum(Q[i, j, k] for i in H_set))
    
    # Big-M constraints
    @constraint(model, [i in H_set, j in C_set],
    sum(Q[i, j, k] for k in intervals) <= M[(i, j)]*y[i, j])
    
    @objective(model, Min, sum(y))
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
    min_units = Int(objective_value(model))
    verbose && println("Minimum number of units: $(min_units)")
    prob.results_dict[:min_units] = min_units
return
end
