# Workflow using XLSX input:
# 1. Import necessary packages:
@time using CompHENS
using Plots
using JuMP
using HiGHS
using Test

using BARON

exportall(CompHENS)

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "Gundersen_4_stream.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 10.0, verbose = true)

#intervals = CompHENS.generate_transshipment_intervals(prob)
#print_full(intervals)

# 4. Solve minimum utilities problem
@time solve_minimum_utilities_subproblem!(prob; verbose = true)

print_min_utils_pinch_points(prob)

# 5. Solve the minimum number of units subproblem:
@time solve_minimum_units_subproblem!(prob)
@test prob.min_units == 5

# 6. Generate stream matches
EMAT = prob.ΔT_min/8 
add_units = 2
@time generate_stream_matches!(prob, EMAT; add_units = add_units)

# 7. Network generation:
# Specify which superstructure to use for each stream
obj_func = CostScaledPaterson()
overall_network = merge(construct_superstructure(prob.stream_names, FloudasCiricGrossmann(), prob), construct_superstructure(prob.utility_names, FloudasCiricGrossmann(), prob))
base_cost, cost_coeff, scaling_coeff = 4000, 500, 0.83
optimizer = BARON.Optimizer

generate_network!(prob, EMAT, overall_network; obj_func = CostScaledPaterson(), optimizer = optimizer, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true)



