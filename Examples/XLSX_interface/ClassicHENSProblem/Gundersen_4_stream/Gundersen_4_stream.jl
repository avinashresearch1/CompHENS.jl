# Workflow using XLSX input:
# 1. Import necessary packages:
@time using CompHENS

using Plots
using JuMP
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
EMAT = prob.ΔT_min/8 
prob.results_dict[:add_units] = 1
@time generate_stream_matches!(prob, EMAT; digits = 8, verbose  = false)
prob.results_dict[:Q]

# 7. Network generation:
# Specify which superstructure to use for each stream
#obj_func = CostScaledPaterson()
#overall_network = merge(construct_superstructure(prob.stream_names, FloudasCiricGrossmann(), prob), construct_superstructure(prob.utility_names, FloudasCiricGrossmann(), prob))
base_cost, cost_coeff, scaling_coeff = 4000, 500, 0.83
obj_func = AreaArithmeticMean()
# using ALPINE for first-pass:
optimizer = SCIP_solver

obj_func = CostScaledPaterson()
#obj_func = Tupper()
#optimizer = optimizer_with_attributes(BARON.Optimizer, "MaxTime" => 20.0, "AbsConFeasTol" => 1)

generate_network!(prob, EMAT; optimizer = optimizer, obj_func = obj_func, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true)
#generate_network!(prob, EMAT, "couenne"; verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true)
model = prob.results_dict[:network_gen_model]
#print(model)
file_name = "/home/avinash/Desktop/COMPHENS/CompHENS.jl/Result_Plots/Gundersen_4_stream.pdf"

plot_HEN_streamwise(prob, model, overall_network, file_name; digits = 1)
#stream = "C2"
#CompHENS.print_stream_results(stream, prob, model, overall_network[stream])
value.(model[:ΔT_upper])
value.(model[:ΔT_lower])
value.(model[:T_LMTD])
get_design_area(prob)
# 

# Setting optimizer for 2nd-pass:
#Get initial values 
var_names = JuMP.all_variables(model)

initial_values = values.(model[JuMP.all_variables(model)])
value(model[Symbol(initial_values[1])])
for var in JuMP.all_variables(model)
    push!(initial_values, var.name => value(model[var]))
end

function get_results_dict(m::JuMP.Model)
    av = JuMP.all_variables(m);
    d = Dict()
    for v in av
        # special handling of time series, sub-indices, etc.
        # e.g. a variable could have three sub-indices in a JuMP.Containers
        v = string(v) # e.g.  "x[a,2,3]"
        k = Symbol(v[1:prevind(v, findfirst('[', v))[1]]) # :x
        vm = value.(m[k])
        d[k] = vm.data
        #  e.g. typeof(vm) is JuMP.Containers.SparseAxisArray{Float64, 3, Tuple{SubString{String}, Int64, Int64}}
    end
    return d
end

res = get_results_dict(model)

set_optimizer(model, IPOPT_solver)


optimize!(model)

variable_refs = []
start_vals  = []

keys(object_dictionary(model))


# Using an ugle hack:
variable_ref = JuMP.all_variables(model)
start_vals = vcat(start_vals, res[:t], res[:f], res[:ΔT_upper], res[:ΔT_lower], res[:T_LMTD]) .+ 0.01

initial_values = (; v_names = variable_ref, v_starts = start_vals)

optimizer = IPOPT_solver
optimizer = JUNIPER_solver
generate_network!(prob, EMAT; optimizer = optimizer, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true, initial_values = initial_values)
