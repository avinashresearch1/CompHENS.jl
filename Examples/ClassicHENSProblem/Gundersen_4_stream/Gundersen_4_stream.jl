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

# 5. Solve the minimum number of units subproblem:plot_HEN_streamwise(prob::ClassicHENSProblem, model::AbstractModel, overall_network::Dict{String, AbstractSuperstructure}, file_name; digits = 1)
@time solve_minimum_units_subproblem!(prob)
@test prob.min_units == 5

# 6. Generate stream matches
EMAT = prob.ΔT_min/2 
add_units = 2
@time generate_stream_matches!(prob, EMAT; add_units = add_units, digits = 8, verbose  = true)
prob.results_dict[:Q]

# 7. Network generation:
# Specify which superstructure to use for each stream
obj_func = CostScaledPaterson()
overall_network = merge(construct_superstructure(prob.stream_names, FloudasCiricGrossmann(), prob), construct_superstructure(prob.utility_names, FloudasCiricGrossmann(), prob))
base_cost, cost_coeff, scaling_coeff = 4000, 500, 1#0.83

#using Ipopt
#optimizer = Ipopt.Optimizer 
optimizer = BARON.Optimizer

#EMAT = -1
#obj_func = AreaArithmeticMean()
#obj_func = Tupper()

generate_network!(prob, EMAT, overall_network; obj_func = obj_func, optimizer = optimizer, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true)
model = prob.results_dict[:network_gen_model]
print(model)
file_name = "/home/avinash/Desktop/COMPHENS/CompHENS.jl/Result_Plots/Gundersen_4_stream.pdf"

plot_HEN_streamwise(prob, model, overall_network, file_name; digits = 1)
#stream = "C2"
#CompHENS.print_stream_results(stream, prob, model, overall_network[stream])
value.(model[:ΔT_upper])
value.(model[:ΔT_lower])
value.(model[:T_LMTD])
get_design_area(prob)
# 

#print(model)

