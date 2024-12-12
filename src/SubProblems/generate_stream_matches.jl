"""
$(TYPEDSIGNATURES)

Constructs and solves the MILP transportation problem presented in Anantharaman 2010 P. 132 to generate stream matches. 

"""
function generate_stream_matches!(prob::ClassicHENSProblem, EMAT; level=:quaternary_temperatures, add_units=prob.results_dict[:add_units], optimizer=HIGHS_solver, verbose=false, digits=4)
    verbose && @info "Solving the Stream Match Generator subproblem"

    haskey(prob.results_dict, :min_units) || error("Minimum units data not available. Solve corresponding subproblem first.")
    level == :primary_temperatures && get_primary_temperatures!(prob; verbose=verbose)
    level == :secondary_temperatures && get_secondary_temperatures!(prob; verbose=verbose)
    level == :tertiary_temperatures && get_tertiary_temperatures!(prob, EMAT; verbose=verbose)
    level == :quaternary_temperatures && get_quaternary_temperatures!(prob, EMAT; verbose=verbose)

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
            if is_feasible(m, n; EMAT=EMAT) == 0
                for i in H_set
                    for j in C_set
                        JuMP.fix(Q[i, m, j, n], 0.0; force=true)
                    end
                end
            end
        end
    end

    # Getting bounds. TODO: Cross-check this.
    M = Dict()
    for i in H_set
        for j in C_set
            M[(i, j)] = min(sum(get_contribution(i, m) for m in hot_cc), sum(get_contribution(j, n) for n in cold_cc))
        end
    end

    # Heat balance: For each hot stream in each hot side temperature interval
    @constraint(model, [i in H_set, m in hot_cc],
        sum(sum(Q[i, m, j, n] for j in C_set) for n in cold_cc) == get_contribution(i, m))

    # Heat balance: For each cold stream in each cold side temperature interval
    @constraint(model, [j in C_set, n in cold_cc],
        sum(sum(Q[i, m, j, n] for i in H_set) for m in hot_cc) == get_contribution(j, n))

    # Big-M constraints
    @constraint(model, [i in H_set, j in C_set],
        sum(sum(Q[i, m, j, n] for m in hot_cc) for n in cold_cc) <= M[(i, j)] * y[i, j])

    # Num HX constraint
    @constraint(model,
        sum(sum(y[i, j] for i in H_set) for j in C_set) == prob.results_dict[:min_units] + add_units)

    @objective(model, Min,
        sum(sum(sum(sum((Q[i, m, j, n] * is_feasible(m, n; EMAT=EMAT) / (U(H_stream_set[i], C_stream_set[j]) * LMTD(m, n; verbose=verbose))) for n in cold_cc) for j in C_set) for m in hot_cc) for i in H_set))

    set_optimizer(model, optimizer)
    !verbose && set_silent(model)
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
    termination_status(model) == MathOptInterface.OPTIMAL || return println("`\n Stream Match Generator problem INFEASIBLE. Try changing number of units. \n")
    post_HLD_matches!(prob, model, level; digits=digits, display=verbose)
    return
end

"""
$(TYPEDSIGNATURES)

Postprocessing after solving stream generation subproblem. 
Displays the matches and heat load distribution in a 2-D matrix form, maintains stream ordering.
"""
function post_HLD_matches!(prob::ClassicHENSProblem, model::AbstractModel, level=:quaternary_temperatures; digits=4, display=true)
    H_set = prob.hot_names # ordered sets
    C_set = prob.cold_names
    hot_cc, cold_cc = prob.results_dict[level].hot_cc, prob.results_dict[level].cold_cc

    Q_match = zeros(Float64, length(C_set), length(H_set))
    y_match = zeros(Int, length(C_set), length(H_set))
    Q = NamedArray(Q_match, (C_set, H_set))
    y = NamedArray(y_match, (C_set, H_set))

    for i in H_set
        for j in C_set
            y[j, i] = round(value(model[:y][i, j]); digits=0)
            Q[j, i] = sum(sum(value(model[:Q][i, m, j, n]) for m in hot_cc) for n in cold_cc)
        end
    end

    println("\n") # Double line after optimizer output.

    display && @show y
    display && @show Q


    # Important to distinguish between HX and Heat loads, can have a HX with Q = 0, especially for multiperiod cases.
    HX_list = Dict{String,Vector{String}}()
    HLD_list = Dict{String,Vector{String}}()
    for i in H_set
        matches = filter(j -> y[j, i] == 1, C_set)
        push!(HX_list, i => matches)
        HLDs = filter(C_set) do j
            round(Q[j, i]; digits=2) > 0.0 # Rounding should only be for the filtering
        end
        push!(HLD_list, i => HLDs)
    end
    for j in C_set
        matches = filter(i -> y[j, i] == 1, H_set)
        push!(HX_list, j => matches)
        HLDs = filter(H_set) do i
            round(Q[j, i]; digits=2) > 0.0 # Rounding should only be for the filtering
        end
        push!(HLD_list, j => HLDs)
    end

    prob.results_dict[:y] = y
    prob.results_dict[:Q] = Q
    prob.results_dict[:HX_list] = HX_list
    prob.results_dict[:HLD_list] = HLD_list

    # Print
    display && println("Num HXs: $(sum(all.(y .== 1))), Num HLDs: $(sum(all.(round.(Q; digits = 2) .> 0.0)))")
    return
end

#= Deprecated
"""
$(TYPEDSIGNATURES)

Displays the matches and heat load distribution in a 2-D matrix form. 
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
    return hld
end
=#

"""
$(TYPEDSIGNATURES)

Constructs and solves the MILP transportation problem presented in Anantharaman 2010 P. 132 to generate stream matches. 

"""
function generate_stream_matches!(prob::MultiPeriodFlexibleHENSProblem, EMAT; level=:quaternary_temperatures, add_units=prob.results_dict[:add_units], optimizer=HIGHS_solver, verbose=false, digits=4)
    verbose && @info "Solving the Stream Match Generator subproblem"

    haskey(prob.results_dict, :min_units) || error("Minimum units data not available. Solve corresponding subproblem first.")

    hot_cc, cold_cc = Dict(), Dict()
    for (k, v) in prob.period_streams_dict
        level == :primary_temperatures && get_primary_temperatures!(v; verbose=verbose)
        level == :secondary_temperatures && get_secondary_temperatures!(v; verbose=verbose)
        level == :tertiary_temperatures && get_tertiary_temperatures!(v, EMAT; verbose=verbose)
        level == :quaternary_temperatures && get_quaternary_temperatures!(v, EMAT; verbose=verbose)
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
                        push!(imjnt_tuple_set, (i, m, j, n, t))
                    end
                end
            end
        end
    end

    @variable(model, 0 <= Q[imjnt_tuple_set]) # Heat exchanged by hot stream i in hot interval m to cold stream j in cold interval n in period t

    for t in periods
        for m in hot_cc[t]
            for n in cold_cc[t]
                if is_feasible(m, n; EMAT=EMAT) == 0
                    for i in H_set
                        for j in C_set
                            JuMP.fix(Q[(i, m, j, n, t)], 0.0; force=true)
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
            sum(sum(Q[(i, m, j, n, t)] for j in C_set) for n in cold_cc[t]) == get_contribution(i, m))
    end

    # Heat balance: For each cold stream in each cold side temperature interval
    for t in periods
        @constraint(model, [j in C_set, n in cold_cc[t]],
            sum(sum(Q[(i, m, j, n, t)] for i in H_set) for m in hot_cc[t]) == get_contribution(j, n))
    end

    # Big-M constraints
    for t in periods
        @constraint(model, [i in H_set, j in C_set],
            sum(sum(Q[(i, m, j, n, t)] for m in hot_cc[t]) for n in cold_cc[t]) <= M(i, j, prob) * y[i, j])
    end

    # Num HX constraint
    @constraint(model,
        sum(sum(y[i, j] for i in H_set) for j in C_set) == prob.results_dict[:min_units] + add_units)

    # Sets the pseudo area for each time period t
    for t in periods
        @constraint(model,
            pseudo_area[t] == sum(sum(sum(sum((Q[(i, m, j, n, t)] * is_feasible(m, n; EMAT=EMAT) / (U(H_stream_set[i], C_stream_set[j]) * LMTD(m, n; verbose=verbose))) for n in cold_cc[t]) for j in C_set) for m in hot_cc[t]) for i in H_set)
        )
    end

    @objective(model, Min,
        sum(pseudo_area)
    )

    set_optimizer(model, optimizer)
    !verbose && set_silent(model)
    optimize!(model)
    if verbose
        @show termination_status(model)
        @show primal_status(model)
        @show dual_status(model)
    end

    # Post-processing
    termination_status(model) == MathOptInterface.OPTIMAL || return println("`\n Stream Match Generator problem INFEASIBLE. Try adding more units. \n")
    post_HLD_matches!(prob, model, level; digits=digits, display=verbose)
    #=
    y_match = Dict()
    for t in periods
        Q_match = Dict()
        for i in H_set
            for j in C_set
                y_match[(i,j)] = value.(y[i,j])
                Q_match[(i,j)] = round(sum(sum(value.(Q[(i,m,j,n,t)]) for m in hot_cc[t]) for n in cold_cc[t]); digits = 4)
            end
        end
        prob.period_streams_dict[t].results_dict[:Q] = Q_match
    end
    prob.results_dict[:y] = y_match
    =#
    return
end

"""
$(TYPEDSIGNATURES)

Postprocessing after solving stream generation subproblem. 
Displays the matches and heat load distribution in a 2-D matrix form, maintains stream ordering.
"""
function post_HLD_matches!(prob::MultiPeriodFlexibleHENSProblem, model::AbstractModel, level=:quaternary_temperatures; digits=4, display=true)
    H_set = prob.period_streams_dict[prob.period_names[1]].hot_names # ordered sets
    C_set = prob.period_streams_dict[prob.period_names[1]].cold_names

    hot_cc, cold_cc = Dict(), Dict()
    for (k, v) in prob.period_streams_dict
        push!(hot_cc, k => v.results_dict[level].hot_cc)
        push!(cold_cc, k => v.results_dict[level].cold_cc)
    end

    y_match = zeros(Int, length(C_set), length(H_set))
    y = NamedArray(y_match, (C_set, H_set))
    for t in prob.period_names
        Q_match = zeros(Float64, length(C_set), length(H_set))
        Q = NamedArray(Q_match, (C_set, H_set))
        for i in H_set
            for j in C_set
                y[j, i] = round(value(model[:y][i, j]); digits=0)
                Q[j, i] = sum(sum(value.(model[:Q][(i, m, j, n, t)]) for m in hot_cc[t]) for n in cold_cc[t])
            end
        end
        # All periods must have the same HX list but may have different HLD_lists. The same matrix y is copied to each period's `results_dict` in order to allow reuse of `ClassicHENSProblem` code. 
        # y is also placed in the `results_dict` of the overall prob.
        # The `Q` values are usually different from period to period. 
        prob.period_streams_dict[t].results_dict[:y] = y
        prob.period_streams_dict[t].results_dict[:Q] = Q
    end
    prob.results_dict[:y] = y # Remeber assignment not deep copy. Changing one changes all.
    println("\n") # Double line after optimizer output.

    for t in prob.period_names
        HX_list = Dict{String,Vector{String}}()
        HLD_list = Dict{String,Vector{String}}()
        for i in H_set
            matches = filter(j -> prob.period_streams_dict[t].results_dict[:y][j, i] == 1, C_set)
            HLDs = filter(j -> (round(prob.period_streams_dict[t].results_dict[:Q][j, i]; digits=0) > 0.0), C_set)
            push!(HX_list, i => matches)
            push!(HLD_list, i => HLDs)
        end
        for j in C_set
            matches = filter(i -> prob.period_streams_dict[t].results_dict[:y][j, i] == 1, H_set)
            HLDs = filter(i -> (round(prob.period_streams_dict[t].results_dict[:Q][j, i]; digits=0) > 0.0), H_set)
            push!(HX_list, j => matches)
            push!(HLD_list, j => HLDs)
        end
        prob.period_streams_dict[t].results_dict[:HX_list] = HX_list
        prob.period_streams_dict[t].results_dict[:HLD_list] = HLD_list

        # Print
        display && println("$t : Num HXs: $(sum(all.(prob.period_streams_dict[t].results_dict[:y] .== 1))), Num HLDs: $(sum(all.(round.(prob.period_streams_dict[t].results_dict[:Q]; digits = 0) .> 0.0)))")
    end
    return
end

"""
Displays HLD for Multiperiod problem
"""
function print_HLD(prob::MultiPeriodFlexibleHENSProblem)
    println("Overall Units")
    @show prob.results_dict[:y]
    println()
    for t in prob.period_names
        println("PROBLEM $(t)")
        println("HLD \n")
        @show prob.period_streams_dict[t].results_dict[:Q]
        #println("HLD list")
        #@show prob.period_streams_dict[t].results_dict[:HLD_list]
        print("\n")
    end
    return
end

