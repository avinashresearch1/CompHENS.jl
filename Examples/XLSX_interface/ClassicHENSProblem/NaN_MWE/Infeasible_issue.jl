# Workflow using XLSX input:
# 1. Import necessary packages:
@time using CompHENS

using Plots
using JuMP
using HiGHS
using Test

using BARON
#using XLSX
#using DataFrames


# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_KangLiu_Period1.xlsx")
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 10, verbose = true)
@time solve_minimum_utilities_subproblem!(prob; verbose = true)
print_min_utils_pinch_points(prob)
(; plt) = plot_composite_curve(prob)
plt |> display

@time solve_minimum_units_subproblem!(prob; verbose = true)

EMAT = 2.5
prob.results_dict[:add_units] = 1

@time generate_stream_matches!(prob, EMAT; digits = 8, verbose = true)
sum(all.(round.(prob.results_dict[:Q] .> 0.0)))

# Network generation:
# Specify which superstructure to use for each stream
# Default is below: 
overall_network = construct_superstructure(prob.all_names, FloudasCiricGrossmann(), prob)
obj_func = CostScaledPaterson()
base_cost, cost_coeff, scaling_coeff =  8333.3,  641.7, 1

optimizer = optimizer_with_attributes(BARON.Optimizer, "MaxTime" => -1, "AbsConFeasTol" => 10)
results_df = generate_network!(prob, EMAT; optimizer = optimizer, obj_func = obj_func, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true)

prob.results_dict[:HLD_list]