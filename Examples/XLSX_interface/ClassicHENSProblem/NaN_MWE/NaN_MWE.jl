@time using CompHENS
using Plots
using JuMP
using HiGHS
using Test

file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_MWE.xlsx")
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 45.0, verbose = true)
@time solve_minimum_utilities_subproblem!(prob; verbose = true)

print_min_utils_pinch_points(prob)
(; plt) = plot_composite_curve(prob)
plt |> display

@time solve_minimum_units_subproblem!(prob; verbose = true)

EMAT = 2.5
prob.results_dict[:add_units] = 1

@time generate_stream_matches!(prob, EMAT; digits = 8, verbose = true)

level = :primary_temperatures
hot_cc, cold_cc = prob.results_dict[level].hot_cc, prob.results_dict[level].cold_cc
print_full(hot_cc)
print_full(cold_cc)

level = :quaternary_temperatures
hot_cc, cold_cc = prob.results_dict[level].hot_cc, prob.results_dict[level].cold_cc
print_full(hot_cc)
print_full(cold_cc)

prob.results_dict[:Q]


# Network generation:
# Specify which superstructure to use for each stream
# Default is below: 
overall_network = construct_superstructure(prob.all_names, FloudasCiricGrossmann(), prob)
obj_func = CostScaledPaterson()
base_cost, cost_coeff, scaling_coeff = 4000, 500, 0.83

optimizer = optimizer_with_attributes(BARON.Optimizer, "MaxTime" => 20.0, "AbsConFeasTol" => 1)
results_df = generate_network!(prob, EMAT; optimizer = optimizer, obj_func = obj_func, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true)

model = prob.results_dict[:network_gen_model]
file_name = "/home/avinash/Desktop/COMPHENS/CompHENS.jl/Examples/XLSX_interface/ClassicHENSProblem/NaN_MWE/NaN_MWE.pdf"

plot_HEN_streamwise(prob, model, overall_network, file_name; digits = 1)
