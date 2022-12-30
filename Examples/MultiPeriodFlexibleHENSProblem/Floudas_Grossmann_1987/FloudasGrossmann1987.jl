# Workflow using XLSX input:
# 1. Import necessary packages:
@time using CompHENS
using Plots
using JuMP
using HiGHS
using XLSX
using DataFrames

# 2. Specify path to xlsx file, Construct the appropriate kind of problem: Here it is a `MultiPeriodFlexibleHENSProblem`.
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_FloudasGrossmann.xlsx")
prob = MultiPeriodFlexibleHENSProblem(file_path_xlsx, 3; verbose = true)

# 3. Solve the minimum utilities subproblem:
@time solve_minimum_utilities_subproblem!(prob)
print_min_utils_pinch_points(prob)

# 4. Solve multiperiod minimum units problem
@time solve_minimum_units_subproblem!(prob; verbose = true)


# 5. Solve subproblem 1: minimum utilities. 

# Using formulation of Prob. 16.5 Biegler, Grossmann, Westerberg book. Pg. 533.

min_utils = solve_minimum_utilities_subproblem(prob)



