# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS
using Plots
using JuMP
using HiGHS

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_ColbergMorari.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0, verbose = true)

# 4. Subdivide into intervals and attain the hot and cold composite curves.
intervals = CompHENS.generate_heat_cascade_intervals(prob, 0.0)
hot_ref_enthalpy, cold_ref_enthalpy = 0.0, 172.596
sorted_intervals = intervals
ylabel = "T [°C or K]"
xlabel = "Heat duty Q"

plt = CompHENS.plot_composite_curve(sorted_intervals; hot_ref_enthalpy, cold_ref_enthalpy, ylabel = "T [°C or K]", xlabel = "Heat duty Q")
ylims!((300,700))

# 5. Solve subproblem 1: minimum utilities. 

# Using formulation of Prob. 16.5 Biegler, Grossmann, Westerberg book. Pg. 533.

min_utils_soln = solve_minimum_utilities_subproblem(prob)

# Subproblem 2:

ΔT_min_absolute = 0.0 # Use the minimum feasible here to get the absolute minimum number of units. 
intervals = CompHENS.generate_heat_cascade_intervals(prob, ΔT_min_absolute)
model = Model()
H_set = union(keys(prob.hot_utilities_dict), keys(prob.hot_streams_dict))
C_set = union(keys(prob.cold_utilities_dict), keys(prob.cold_streams_dict))

# Code design: One alternative here is to define the residual, R variables as `SparseAxisArrays` in JuMP. 
# However, https://discourse.julialang.org/t/behaviour-of-sparse-array-variables-in-jump/36185/4 suggests instead fixing them to 0.0 and letting the pre-solver take care of them, so doing this.
# Memory allocation is generally not the issue. Note, that since the stream contributions are known, the presolver should automatically fix them to be 0.0 and get rid of irrelevant residuals. 

@variable(model, y[H_set, C_set], Bin)
@variable(model, 0 <= R[H_set, intervals]) # Notation: R[interval] is the residual heat exiting a given interval
@variable(model, 0 <= Q[H_set, C_set, intervals]) # Heat exchanged by hot stream i to cold stream j in interval k
for i in H_set
    JuMP.fix(R[i, last(intervals)], 0.0; force = true)
end

# First interval: Heat balance for each hot stream. Presolver sets residuals to zero.
@constraint(model, 
sum(Q_in[hu] for hu in keys(first(intervals).hot_utilities_contribs)) + first(intervals).total_stream_heat_in == R[first(intervals)] + sum(Q_out[cu] for cu in keys(first(intervals).cold_utilities_contribs)) + first(intervals).total_stream_heat_out)

# Remaining intervals
@constraint(model, [i in 2:length(intervals)],
R[intervals[i-1]] + sum(Q_in[hu] for hu in keys(intervals[i].hot_utilities_contribs)) + intervals[i].total_stream_heat_in == R[intervals[i]] + sum(Q_out[cu] for cu in keys(intervals[i].cold_utilities_contribs)) + intervals[i].total_stream_heat_out)


@constraint(model, 

)

@objective(model, Min, sum(y))






H_set = keys(prob.hot_utilities_dict) 
CU_set = keys(prob.cold_utilities_dict)



