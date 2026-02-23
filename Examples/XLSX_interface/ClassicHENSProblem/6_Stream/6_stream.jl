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
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_SimpleExample.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0, verbose = true)

# 4. Solve minimum utilities problem

@time solve_minimum_utilities_subproblem!(prob)
print_min_utils_pinch_points(prob)

# 5. Solve the minimum number of units subproblem:
@time solve_minimum_units_subproblem!(prob)
@test prob.min_units == 5

# 6. Generate stream matches
EMAT = 2.5
@time generate_stream_matches!(prob, EMAT; add_units = 1)

# 7. Network generation:
# Specify which superstructure to use for each stream
obj_func = CostScaledPaterson()
overall_network = merge(construct_superstructure(prob.stream_names, FloudasCiricGrossmann(), prob), construct_superstructure(prob.utility_names, FloudasCiricGrossmann(), prob))
cost_coeff, scaling_coeff = 670, 0.83
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_iter" => 5000, "tol" => 1e-6)

generate_network!(prob, EMAT; overall_network = overall_network, obj_func = obj_func, optimizer = optimizer, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, save_model = true)
model = prob.results_dict[:network_gen_model]
@test termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL]
@test primal_status(model) == MOI.FEASIBLE_POINT

#=
using Alpine
const alpine = JuMP.optimizer_with_attributes(
    Alpine.Optimizer,
    # "minlp_solver" => minlp_solver,
    "nlp_solver" => JuMP.optimizer_with_attributes(
        Ipopt.Optimizer,
        MOI.Silent() => true,
        "sb" => "yes",
        "max_iter" => Int(1E4),
    ),
    "mip_solver" => JuMP.optimizer_with_attributes(
        HiGHS.Optimizer,
        "presolve" => "on",
        "log_to_console" => false,
    ),
    "presolve_bt" => true,
    "apply_partitioning" => true,
    "partition_scaling_factor" => 10,
)
=#
