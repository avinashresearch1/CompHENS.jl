# Workflow using XLSX input:
# 1. Import necessary packages:
@time using CompHENS
using Plots
using JuMP
using HiGHS
using Ipopt
using MathOptInterface
using Test
const MOI = MathOptInterface

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_ColbergMorari.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0, verbose = true)

# 4. Solve minimum utilities problem

@time solve_minimum_utilities_subproblem!(prob)
print_min_utils_pinch_points(prob)
@test prob.pinch_points == [(517.0, 497.0)]
@test isapprox(prob.hot_utilities_dict["ST"].Q, 244.13; atol = 1)
@test isapprox(prob.cold_utilities_dict["CW"].Q, 172.6; atol = 1)

# 5. Solve the minimum number of units subproblem:
@time solve_minimum_units_subproblem!(prob)
@test prob.min_units == 8

# 6. Generate stream matches
EMAT = 2.5
prob.results_dict[:add_units] = 1
@time generate_stream_matches!(prob, EMAT; digits = 8)
prob.results_dict[:Q]

# 7. Network generation:
# Specify which superstructure to use for each stream
obj_func = CostScaledPaterson()
overall_network = merge(construct_superstructure(prob.stream_names, FloudasCiricGrossmann(), prob), construct_superstructure(prob.utility_names, FloudasCiricGrossmann(), prob))
base_cost, cost_coeff, scaling_coeff = 8600, 670, 0.83
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_iter" => 5000, "tol" => 1e-6)

generate_network!(prob, EMAT; overall_network = overall_network, obj_func = obj_func, optimizer = optimizer, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true)
model = prob.results_dict[:network_gen_model]
@test termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL]
@show termination_status(model)
@test primal_status(model) == MOI.FEASIBLE_POINT
# print(model)
value.(model[:ΔT_upper])
value.(model[:ΔT_lower])
value.(model[:T_LMTD])
get_design_area(prob)
