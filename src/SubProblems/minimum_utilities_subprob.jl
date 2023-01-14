using JuMP
using HiGHS

"""
$(TYPEDSIGNATURES)

Constructs and solves the LP transshipment formulation of Papoulias_Grossmann_1983 to determine the minimum utility cost.
Returns:
"""
function solve_minimum_utilities_subproblem!(prob::ClassicHENSProblem; time_limit = 60.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false, save_model = false)
    verbose && @info "Solving the minimum utilities subproblem"
    intervals = generate_transshipment_intervals(prob)
    model = Model()
    HU_set = keys(prob.hot_utilities_dict) 
    CU_set = keys(prob.cold_utilities_dict)
    
    @variable(model, 0 <= Q_in[HU_set])
    @variable(model, 0 <= Q_out[CU_set])
    @variable(model, 0 <= R[intervals]) # Notation: R[interval] is the residual heat exiting a given interval
    JuMP.fix(R[last(intervals)], 0.0; force = true)
    
    # First interval: Entering == Leaving
    @constraint(model, 
    sum(Q_in[hu] for hu in keys(first(intervals).hot_side.hot_utils)) + first(intervals).hot_side.total_stream_heat_in == R[first(intervals)] + sum(Q_out[cu] for cu in keys(first(intervals).cold_side.cold_utils)) + first(intervals).cold_side.total_stream_heat_out)
    
    # Remaining intervals
    @constraint(model, [k in 2:length(intervals)],
    R[intervals[k-1]] + sum(Q_in[hu] for hu in keys(intervals[k].hot_side.hot_utils)) + intervals[k].hot_side.total_stream_heat_in == R[intervals[k]] + sum(Q_out[cu] for cu in keys(intervals[k].cold_side.cold_utils)) + intervals[k].cold_side.total_stream_heat_out)
    
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

    save_model && push!(prob.results_dict, :min_utils_model => model)

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
            push!(pinch_points, (interval.hot_side.lower.T, interval.cold_side.lower.T))
        end
    end
    prob.results_dict[:pinch_points] = pinch_points
    return
end

"""
$(TYPEDSIGNATURES)

Based on Floudas_Grossmann_1987.
Constructs and solves the LP transshipment formulation of Papoulias_Grossmann_1983 to determine the minimum utility cost in each period of operation. 
The pinch points for each period `<period_name>` are placed in the `prob.period_streams_dict[<period_name>].results_dict[:pinch_points]`, also see 
[TODO:] Implement parallel computing version.
"""
function solve_minimum_utilities_subproblem!(prob::MultiPeriodFlexibleHENSProblem; time_limit = 60.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    for (k,v) in prob.period_streams_dict
        verbose && @info "Problem $(k)"
        solve_minimum_utilities_subproblem!(v; time_limit = time_limit, presolve = presolve, optimizer = optimizer, verbose = verbose)
    end
    return
end

"""
Prints the minimum utility consumption and pinch points for `prob`
Can only be called after the solving `solve_minimum_utilities_subproblem!(prob)`
"""
function print_min_utils_pinch_points(prob::ClassicHENSProblem; drop_infinite = true, digits = 1)
    for (k,v) in prob.hot_utilities_dict
        drop_infinite && v.Q == Inf && continue
        Q_util = round(v.Q; digits = digits)
        println("Hot utility $(k): $(Q_util) ")
    end
    for (k,v) in prob.cold_utilities_dict
        drop_infinite && v.Q == Inf && continue
        Q_util = round(v.Q; digits = digits)
        println("Cold utility $(k): $(Q_util) ")
    end
    for pinch in prob.results_dict[:pinch_points]
        print("Pinch points: $(pinch)    ")
    end
    print("\n")
    return
end

function print_min_utils_pinch_points(prob::MultiPeriodFlexibleHENSProblem; drop_infinite = true, digits = 1)
    for k in prob.period_names
        println("PROBLEM $(k)")
        print_min_utils_pinch_points(prob.period_streams_dict[k]; drop_infinite = drop_infinite, digits = digits)
        print("\n")
    end
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
