# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS
using Plots
using JuMP
using HiGHS
using Ipopt
using MathOptInterface
using XLSX
using DataFrames
const MOI = MathOptInterface

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

optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_iter" => 8000, "tol" => 1e-6)

output_folder = nothing
generate_network!(prob, EMAT; optimizer = optimizer, verbose = true, output_folder = output_folder, save_model = true)
for period in prob.period_names
    model = prob.period_streams_dict[period].results_dict[:network_gen_model]
    status = termination_status(model)
    pstatus = primal_status(model)
    if !(status in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL] && pstatus == MOI.FEASIBLE_POINT)
        error("Period $(period) failed IPOPT solve: termination=$(status), primal=$(pstatus)")
    end
end
get_design_area(prob)
