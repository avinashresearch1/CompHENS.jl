# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS
using Plots
using JuMP
using HiGHS
using Test

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_ColbergMorari.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0, verbose = true)

# 4. Subdivide into intervals and attain the hot and cold composite curves.

#=
hot_ref_enthalpy, cold_ref_enthalpy = 0.0, 172.596
sorted_intervals = intervals
ylabel = "T [°C or K]"
xlabel = "Heat duty Q"

plt = CompHENS.plot_composite_curve(sorted_intervals; hot_ref_enthalpy, cold_ref_enthalpy, ylabel = "T [°C or K]", xlabel = "Heat duty Q")
ylims!((300,700))
=#

time_limit = 60.0
presolve = true
optimizer = HiGHS.Optimizer
verbose = true

# 5. Solve subproblem 1: minimum utilities. 
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


# Using formulation of Prob. 16.5 Biegler, Grossmann, Westerberg book. Pg. 533.

solve_minimum_utilities_subproblem!(prob)
@test prob.pinch_points == [(517.0, 497.0)]
@test isapprox(prob.hot_utilities_dict["ST"].Q, 244.13; atol = 1)
@test isapprox(prob.cold_utilities_dict["CW"].Q, 172.6; atol = 1)

solve_minimum_units_subproblem!(prob)
@test prob.min_units == 8


