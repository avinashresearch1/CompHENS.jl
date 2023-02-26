# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS
using Plots
using JuMP
using HiGHS
using XLSX
using DataFrames

using BARON

exportall(CompHENS)

# 2. Specify path to xlsx file, Construct the appropriate kind of problem: Here it is a `MultiPeriodFlexibleHENSProblem`.
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_FloudasGrossmann.xlsx")
prob = MultiPeriodFlexibleHENSProblem(file_path_xlsx, 3; verbose = true)

# 3. Solve the minimum utilities subproblem:
solve_minimum_utilities_subproblem!(prob)
print_min_utils_pinch_points(prob)

# 4. Solve multiperiod minimum units problem
solve_minimum_units_subproblem!(prob; verbose = false)

# 5. Solving stream match generator problem
EMAT = 2.5
prob.results_dict[:add_units] = 1
generate_stream_matches!(prob, EMAT; verbose = false)
print_HLD(prob)

optimizer = optimizer_with_attributes(BARON.Optimizer, "MaxTime" => 20.0, "AbsConFeasTol" => 1)

output_folder = "/home/avinash/Desktop/COMPHENS/CompHENS.jl/Examples/MultiPeriodFlexibleHENSProblem/Floudas_Grossmann_1987/"
generate_network!(prob, EMAT; optimizer = optimizer, verbose = true, output_folder = output_folder)
get_design_area(prob)


print_stream_results