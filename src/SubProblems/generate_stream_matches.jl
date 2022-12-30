using JuMP
using HiGHS
using NamedArrays
using MathOptInterface

"""
$(TYPEDSIGNATURES)

Constructs and solves the MILP transportation problem presented in Anantharaman 2010 P. 132 to generate stream matches. 

"""
function generate_stream_matches!(prob::ClassicHENSProblem, EMAT; level = :quaternary_temperatures, add_units = 1, time_limit = 200.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    verbose && @info "Solving the Stream Match Generator subproblem"

    haskey(prob.results_dict, :min_units) || error("Minimum units data not available. Solve corresponding subproblem first.")
    level == :primary_temperatures && get_primary_temperatures!(prob; verbose = verbose)
    level == :secondary_temperatures && get_secondary_temperatures!(prob; verbose = verbose)
    level == :tertiary_temperatures && get_tertiary_temperatures!(prob, EMAT; verbose = verbose)
    level == :quaternary_temperatures && get_quaternary_temperatures!(prob, EMAT; verbose = verbose)

    hot_cc, cold_cc = prob.results_dict[level].hot_cc, prob.results_dict[level].cold_cc
    
    model = Model()
    H_stream_set = merge(prob.hot_utilities_dict, prob.hot_streams_dict)
    H_set = keys(H_stream_set)
    C_stream_set = merge(prob.cold_utilities_dict, prob.cold_streams_dict)
    C_set = keys(C_stream_set)
    
    @variable(model, y[H_set, C_set], Bin)
    @variable(model, 0 <= Q[H_set, hot_cc, C_set, cold_cc]) # Heat exchanged by hot stream i in hot interval m to cold stream j in cold interval n
        
    for m in hot_cc
        for n in cold_cc
            if is_feasible(m, n; EMAT = EMAT) == 0
                for i in H_set
                    for j in C_set
                        JuMP.fix(Q[i, m, j, n], 0.0; force = true)
                    end
                end
            end
        end
    end
    
    # Getting bounds. TODO: Cross-check this.
    M = Dict()
    for i in H_set
        for j in C_set
            M[(i,j)] = min(sum(get_contribution(i,m) for m in hot_cc), sum(get_contribution(j,n) for n in cold_cc))
        end
    end
    
    # Heat balance: For each hot stream in each hot side temperature interval
    @constraint(model, [i in H_set, m in hot_cc], 
    sum(sum(Q[i, m, j, n] for j in C_set) for n in cold_cc) == get_contribution(i,m))
    
    # Heat balance: For each cold stream in each cold side temperature interval
    @constraint(model, [j in C_set, n in cold_cc], 
    sum(sum(Q[i, m, j, n] for i in H_set) for m in hot_cc) == get_contribution(j,n))
    
    # Big-M constraints
    @constraint(model, [i in H_set, j in C_set],
    sum(sum(Q[i, m, j, n] for m in hot_cc) for n in cold_cc) <= M[(i, j)]*y[i, j])

    # Num HX constraint
    @constraint(model,
    sum(sum(y[i, j] for i in H_set) for j in C_set) == prob.results_dict[:min_units] + add_units)
    
    @objective(model, Min,
    sum(sum(sum(sum((Q[i, m, j, n]* is_feasible(m, n; EMAT =  EMAT)/(U(H_stream_set[i], C_stream_set[j])*LMTD(m, n))) for n in cold_cc) for j in C_set) for m in hot_cc) for i in H_set))
    
    set_optimizer(model, optimizer)
    !verbose && set_silent(model)
    presolve && set_optimizer_attribute(model, "presolve", "on")
    set_optimizer_attribute(model, "time_limit", time_limit)
    optimize!(model)
    if verbose
        @show termination_status(model)
        @show primal_status(model)
        @show dual_status(model)
    end

    # Simple Infeasibility check
    #=
    for i in H_set
        for m in hot_cc
            for j in C_set
                for n in cold_cc
                    if value.(Q[i, m, j, n]) > 0.0 && is_feasible(m,n; EMAT = EMAT) == 0
                        error("Infeasible heat transfer detected")
                    end
                end
            end
        end
    end
    =#

    # Post-processing
    termination_status(model) == MathOptInterface.OPTIMAL || return println("`\n Stream Match Generator problem INFEASIBLE. Try adding more units. \n")
    Q_match, y_match = Dict(), Dict()
    for i in H_set
        for j in C_set
            y_match[(i,j)] = value.(y[i,j])
            Q_match[(i,j)] = round(sum(sum(value.(Q[i,m,j,n]) for m in hot_cc) for n in cold_cc); digits = 4)
        end
    end
    prob.results_dict[:y] = y_match
    prob.results_dict[:Q] = Q_match
return
end

"""
Displays the heat load distribution in a 2-D matrix form. 
Can only be called after the stream match generation subproblem has been solved when `prob.results_dict[:Q]` is available.
"""
function print_HLD(prob::ClassicHENSProblem)
    hot_names = vcat(prob.stream_names[:hot_streams], prob.stream_names[:hot_utilities])
    cold_names = vcat(prob.stream_names[:cold_streams], prob.stream_names[:cold_utilities])
    Q = zeros(Float64, length(cold_names), length(hot_names))
    hld = NamedArray(Q, (cold_names, hot_names))
    for (name,val) in enamerate(hld)
        hld[name[1], name[2]] = round(prob.results_dict[:Q][(name[2], name[1])]; digits = 1) 
    end
    @show hld
    return
end

"""
$(TYPEDSIGNATURES)

Constructs and solves the MILP transportation problem presented in Anantharaman 2010 P. 132 to generate stream matches. 

"""
function generate_stream_matches!(prob::MultiPeriodFlexibleHENSProblem, EMAT; level = :quaternary_temperatures, add_units = 1, time_limit = 200.0, presolve = true, optimizer = HiGHS.Optimizer, verbose = false)
    verbose && @info "Solving the Stream Match Generator subproblem"

    haskey(prob.results_dict, :min_units) || error("Minimum units data not available. Solve corresponding subproblem first.")

    hot_cc, cold_cc = Dict(), Dict()
    for (k,v) in prob.period_streams_dict
        level == :primary_temperatures && get_primary_temperatures!(v; verbose = verbose)
        level == :secondary_temperatures && get_secondary_temperatures!(v; verbose = verbose)
        level == :tertiary_temperatures && get_tertiary_temperatures!(v, EMAT; verbose = verbose)
        level == :quaternary_temperatures && get_quaternary_temperatures!(v, EMAT; verbose = verbose)
        push!(hot_cc, k => v.results_dict[level].hot_cc)
        push!(cold_cc, k => v.results_dict[level].cold_cc)
    end

    periods = keys(prob.period_streams_dict)
    H_stream_set = merge(prob.period_streams_dict[prob.period_names[1]].hot_utilities_dict, prob.period_streams_dict[prob.period_names[1]].hot_streams_dict)
    H_set = keys(H_stream_set)
    C_stream_set = merge(prob.period_streams_dict[prob.period_names[1]].cold_utilities_dict, prob.period_streams_dict[prob.period_names[1]].cold_streams_dict)
    C_set = keys(C_stream_set)

    model = Model()

    @variable(model, y[H_set, C_set], Bin)

    # Constructing tuple set for Q[(i, m, j, n, t)] for i in H_set, m in hot_cc[t], j in C_set, n in cold_cc[t], t in periods
    imjnt_tuple_set = Set()
    for i in H_set
        for j in C_set
            for t in periods
                for m in hot_cc[t]
                    for n in cold_cc[t]
                        push!(imjnt_tuple_set, (i,m,j,n,t))
                    end
                end
            end
        end
    end

    @variable(model, 0 <= Q[imjnt_tuple_set]) # Heat exchanged by hot stream i in hot interval m to cold stream j in cold interval n in period t
    
    for t in periods
        for m in hot_cc[t]
            for n in cold_cc[t]
                if is_feasible(m, n; EMAT = EMAT) == 0
                    for i in H_set
                        for j in C_set
                            JuMP.fix(Q[(i,m,j,n,t)], 0.0; force = true)
                        end
                    end
                end
            end
        end
    end

    @variable(model, 0 <= pseudo_area[periods])
       
    # Heat balance: For each hot stream in each hot side temperature interval
    for t in periods
        @constraint(model, [i in H_set, m in hot_cc[t]], 
        sum(sum(Q[(i, m, j, n, t)] for j in C_set) for n in cold_cc[t]) == get_contribution(i,m))
    end
    
    # Heat balance: For each cold stream in each cold side temperature interval
    for t in periods
        @constraint(model, [j in C_set, n in cold_cc[t]], 
        sum(sum(Q[(i, m, j, n, t)] for i in H_set) for m in hot_cc[t]) == get_contribution(j,n))
    end
    
    # Big-M constraints
    for t in periods
        @constraint(model, [i in H_set, j in C_set],
        sum(sum(Q[(i, m, j, n, t)] for m in hot_cc[t]) for n in cold_cc[t]) <= M(i, j, prob)*y[i, j])
    end

    # Num HX constraint
    @constraint(model,
    sum(sum(y[i, j] for i in H_set) for j in C_set) == prob.results_dict[:min_units] + add_units)

    # Sets the pseudo area for each time period t
    for t in periods
        @constraint(model, 
        pseudo_area[t] == sum(sum(sum(sum((Q[(i, m, j, n, t)]* is_feasible(m, n; EMAT =  EMAT)/(U(H_stream_set[i], C_stream_set[j])*LMTD(m, n))) for n in cold_cc[t]) for j in C_set) for m in hot_cc[t]) for i in H_set)
        )
    end
    
    @objective(model, Min,
    sum(pseudo_area)
    )
    
    set_optimizer(model, optimizer)
    !verbose && set_silent(model)
    presolve && set_optimizer_attribute(model, "presolve", "on")
    set_optimizer_attribute(model, "time_limit", time_limit)
    optimize!(model)
    if verbose
        @show termination_status(model)
        @show primal_status(model)
        @show dual_status(model)
    end

    # Post-processing
    termination_status(model) == MathOptInterface.OPTIMAL || return println("`\n Stream Match Generator problem INFEASIBLE. Try adding more units. \n")
    y_match = Dict()
    for k in periods
        Q_match = Dict()
        for i in H_set
            for j in C_set
                y_match[(i,j)] = value.(y[i,j])
                Q_match[(i,j)] = round(sum(sum(value.(Q[(i,m,j,n,t)]) for m in hot_cc[t]) for n in cold_cc[t]); digits = 4)
            end
        end
        prob.period_streams_dict[k].results_dict[:Q] = Q_match
    end
    prob.results_dict[:y] = y_match
return
end

"""
Displays HLD for Multiperiod problem
"""
function print_HLD(prob::MultiPeriodFlexibleHENSProblem)
    for k in prob.period_names
        println("PROBLEM $(k)")
        print_HLD(prob.period_streams_dict[k])
        print("\n")
    end
    return
end