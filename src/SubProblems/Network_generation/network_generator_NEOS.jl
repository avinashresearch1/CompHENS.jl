#=
using JuMP, NEOSServer

function generate_network!(prob::ClassicHENSProblem, EMAT, NEOS_email::String, NEOS_solver::String; overall_network::Dict{String, AbstractSuperstructure} = construct_superstructure(prob.all_names, FloudasCiricGrossmann(), prob), obj_func::NetworkObjective = CostScaledPaterson(), verbose = false, cost_coeff = 100, scaling_coeff = 1, base_cost = 1000, save_model = false, output_file = nothing, set_starting_values = true, initial_values = nothing)
    verbose && @info "Solving the Network Generation subproblem"
    
    haskey(prob.results_dict, :y) || error("Stream match data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :Q) || error("HLD data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :HLD_list) || error("Match list data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :T_bounds) || generate_temperature_bounds!(prob)

    #y, Q = prob.results_dict[:y],  prob.results_dict[:Q]
    
    model = Model() do
        NEOSServer.Optimizer(email=NEOS_email, solver=NEOS_solver) 
    end
    
    # Tuple set with stream and utility s and Edge e
    all_e_tuple_vec = []
    for stream in prob.all_names
        for edge in overall_network[stream].edges
            push!(all_e_tuple_vec, (stream, edge))
        end
    end

    # Tuple set with only streams (not utilities) and edges e
    stream_e_tuple_vec = []
    for stream in prob.stream_names
        for edge in overall_network[stream].edges
            push!(stream_e_tuple_vec, (stream, edge))
        end
    end

    # 1. Declaring the stream-wise variables
    @variable(model, 0.0 <= t[all_e_tuple_vec])
    @variable(model, 0.0 <= f[stream_e_tuple_vec])

    # [WIP]
    # set_starting_values && set_start_values(prob, EMAT, overall_network; verbose = verbose)

    # 2. Sets stream-wise constraints
    for stream in prob.all_names
        add_stream_constraints!(model, prob.all_dict[stream], overall_network[stream], prob)
    end
    
    # 3. Setting the variables and constraints for the matches
    HLD_list = Tuple[]
    for hot in prob.hot_names
        for cold in prob.results_dict[:HLD_list][hot]
            push!(HLD_list, (hot,cold))
        end
    end
    
    @variable(model, ΔT_upper[HLD_list])
    @variable(model, ΔT_lower[HLD_list])
    
    U_dict = Dict()
    for match in HLD_list
        set_lower_bound(model[:ΔT_upper][match], EMAT)
        set_lower_bound(model[:ΔT_lower][match], EMAT)
        hot, cold = prob.all_dict[match[1]], prob.all_dict[match[2]]
        add_match_feasibility_constraints!(model, hot, cold, overall_network[hot.name], overall_network[cold.name])
        # Used to avoid `_parse_NL_expr_runtime` error with U as a function.
        push!(U_dict, match => U(hot, cold))
    end
    
    set_objective_func!(model, HLD_list, obj_func, prob, EMAT; U_dict = U_dict, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost)

    # initial_values are a named tuple of (; v_names::Vector{VariableRef}, v_starts::Vector{Any})
    if !isnothing(initial_values)
        @assert length(initial_values.v_names) == length(initial_values.v_starts)
        for i in 1:length(initial_values.v_names)
            JuMP.set_start_value(initial_values.v_names[i], initial_values.v_starts[i])
        end
    end

    optimize!(model)
    results_df = postprocess_network!(prob, model, HLD_list, overall_network)
    save_model && push!(prob.results_dict, :network_gen_model => model)
    isnothing(output_file) || plot_HEN_streamwise(prob, model, overall_network, output_file; digits = 1)
    return results_df
end

using JuMP, NEOSServer

NEOS_email = "a0088599@gmail.com"
NEOS_solver = "OCTERACT"


model = Model() do
    NEOSServer.Optimizer(email=NEOS_email, solver=NEOS_solver, solver_args = ["octeract_engine.time_limit=10.0"]) 
end

@variable(model, x >= 0)

@variable(model, 0 <= y <= 3)

@objective(model, Min, 12x + 20y)

@constraint(model, c1, 6x + 8y >= 100)

@constraint(model, c2, 7x + 12y >= 120)

print(model)
optimize!(model)

value.(x)


exportall(CompHENS)
prob
EMAT
overall_network = construct_superstructure(prob.all_names, FloudasCiricGrossmann(), prob)
obj_func = CostScaledPaterson()
verbose = false
cost_coeff = 100
scaling_coeff = 1
base_cost = 1000
save_model = false
output_file = nothing
set_starting_values = true
initial_values = nothing

   
    haskey(prob.results_dict, :y) || error("Stream match data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :Q) || error("HLD data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :HLD_list) || error("Match list data not available. Solve corresponding subproblem first.")
    haskey(prob.results_dict, :T_bounds) || generate_temperature_bounds!(prob)

    #y, Q = prob.results_dict[:y],  prob.results_dict[:Q]
    
    model = Model() do
        NEOSServer.Optimizer(email=NEOS_email, solver=NEOS_solver, solver_args = ["octeract_engine.time_limit=10.0"]) 
    end
    
    # Tuple set with stream and utility s and Edge e
    all_e_tuple_vec = []
    for stream in prob.all_names
        for edge in overall_network[stream].edges
            push!(all_e_tuple_vec, (stream, edge))
        end
    end

    # Tuple set with only streams (not utilities) and edges e
    stream_e_tuple_vec = []
    for stream in prob.stream_names
        for edge in overall_network[stream].edges
            push!(stream_e_tuple_vec, (stream, edge))
        end
    end

    # 1. Declaring the stream-wise variables
    @variable(model, 0.0 <= t[all_e_tuple_vec])
    @variable(model, 0.0 <= f[stream_e_tuple_vec])

    # [WIP]
    # set_starting_values && set_start_values(prob, EMAT, overall_network; verbose = verbose)

    # 2. Sets stream-wise constraints
    for stream in prob.all_names
        add_stream_constraints!(model, prob.all_dict[stream], overall_network[stream], prob)
    end
    
    # 3. Setting the variables and constraints for the matches
    HLD_list = Tuple[]
    for hot in prob.hot_names
        for cold in prob.results_dict[:HLD_list][hot]
            push!(HLD_list, (hot,cold))
        end
    end
    
    @variable(model, ΔT_upper[HLD_list])
    @variable(model, ΔT_lower[HLD_list])
    
    U_dict = Dict()
    for match in HLD_list
        set_lower_bound(model[:ΔT_upper][match], EMAT)
        set_lower_bound(model[:ΔT_lower][match], EMAT)
        hot, cold = prob.all_dict[match[1]], prob.all_dict[match[2]]
        add_match_feasibility_constraints!(model, hot, cold, overall_network[hot.name], overall_network[cold.name])
        # Used to avoid `_parse_NL_expr_runtime` error with U as a function.
        push!(U_dict, match => U(hot, cold))
    end
    
    set_objective_func!(model, HLD_list, obj_func, prob, EMAT; U_dict = U_dict, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost)

    # initial_values are a named tuple of (; v_names::Vector{VariableRef}, v_starts::Vector{Any})
    if !isnothing(initial_values)
        @assert length(initial_values.v_names) == length(initial_values.v_starts)
        for i in 1:length(initial_values.v_names)
            JuMP.set_start_value(initial_values.v_names[i], initial_values.v_starts[i])
        end
    end
    #set_time_limit_sec(model, 20.0)
    optimize!(model)
    results_df = postprocess_network!(prob, model, HLD_list, overall_network)
    save_model && push!(prob.results_dict, :network_gen_model => model)
    isnothing(output_file) || plot_HEN_streamwise(prob, model, overall_network, output_file; digits = 1)
    return results_df
    =#