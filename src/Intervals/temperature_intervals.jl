using Plots

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds a single temperature interval. Note are attained from either the hot or cold composite curves, not both. 
For double-sided intervals covering both hot and cold side composite curves, use `TransshipmentIntervals`. 
Usually this single-sided `TemperatureInterval` is sufficient for transportation problems.
"""
mutable struct TemperatureInterval{R <: Real}
    index::Int64
    upper::Point{R}
    lower::Point{R}
    """Participating streams in interval"""
    streams::Dict{String, Union{HotStream, ColdStream, SimpleHotUtility, SimpleColdUtility}}
    """Contributions of each of the participating streams. Hot adding, Cold removing heat."""
    contributions::Dict{String, R}
    @add_kwonly function TemperatureInterval{R}(index, upper, lower, streams = Dict{String, Union{HotStream, ColdStream, SimpleHotUtility, SimpleColdUtility}}(), contributions = Dict{String, R}()) where {R}
        new(index, upper, lower, streams, contributions)
    end
end

Base.show(io::IO, interval::TemperatureInterval) = print(io, "itv_$(interval.index)")

# Code design choice: Instead of having fields in `TemperatureInterval` for each stream type,
# idea is that it is better/more extensible to compute (lazily) as below.
function Base.getproperty(interval::TemperatureInterval, sym::Symbol)   
    if sym == :hot_utils
        return filter(interval.contributions) do (k,v)
            interval.streams[k] isa SimpleHotUtility
        end
    elseif sym == :cold_utils
        return filter(interval.contributions) do (k,v)
            interval.streams[k] isa SimpleColdUtility
        end
    elseif sym == :hot_streams
        return filter(interval.contributions) do (k,v)
            interval.streams[k] isa HotStream
        end
    elseif sym == :cold_streams
        return filter(interval.contributions) do (k,v)
            interval.streams[k] isa ColdStream
        end
    elseif sym == :total_stream_heat_in
        return sum(values(interval.hot_streams))
    elseif sym == :total_stream_heat_out
        return sum(values(interval.cold_streams))
    elseif sym == :total_heat_in # Both stream and utility heat in
        return sum(values(interval.hot_streams)) + sum(values(interval.hot_utils))
    elseif sym == :total_heat_out # Both stream and utility heat out
        return sum(values(interval.cold_streams)) + sum(values(interval.cold_utils))
    else # fallback to getfield
        return getfield(interval, sym)
    end
end

"""
$(TYPEDSIGNATURES)

Gets the contribution of first argument to the temperature interval.
Returns 0.0 if no contribution.
"""
function get_contribution(stream::String, interval::TemperatureInterval)
    # [QN: Any advantages of multiple dispatch here? Or just use basic string comparison?
    # TODO: Potential bug returning 0.0 for nonexistent streams.].
    stream in keys(interval.contributions) && return interval.contributions[stream]
    return 0.0
end

"""
$(TYPEDSIGNATURES)

Get contribution of `stream_type` to interval. `stream_type` is a vector of stream types. 
"""
function get_contribution(stream_type::Vector{DataType}, interval::TemperatureInterval)
    streams_dict = filter(interval.contributions) do (k,v)
        typeof(interval.streams[k]) ∈ stream_type
    end
    !isempty(streams_dict) && return sum(values(streams_dict))
    return 0.0
end

function print_full(intervals::Vector{TemperatureInterval}; digits = 1)
    for interval in intervals
        print("itv_", interval.index, ": [", interval.upper.T, ", ", interval.lower.T, "]")
        for (k,v) in interval.contributions
            Q = round(v; digits)
            print(" $k: $Q")
        end
        print("\n")
    end
end


"""
$(TYPEDSIGNATURES)

Returns a vector of sorted `TemperatureInterval`s with the `upper` and `lower` temperatures populated but with the `streams` and `contributions` empty. 
"""
function initialize_temperature_intervals(temp_vec::Vector)
    # QN: Does this mess up type inference? Put Float64 for now, should find how to generalize 
    sort!(unique!(temp_vec), rev = true)
    intervals = TemperatureInterval[]
    for i in 1:length(temp_vec)-1
        upper = Point{Float64}(T = temp_vec[i])
        lower = Point{Float64}(T = temp_vec[i+1])
        push!(intervals, TemperatureInterval{Float64}(index = i, upper = upper, lower = lower))
    end
    return intervals
end

"""
$(TYPEDSIGNATURES)

Mutates the `streams` and `contributions`. Method dispatched depends on type of stream. Only for streams not utilities.
"""
function assign_stream!(interval::TemperatureInterval, stream::HotStream)
    if (stream.T_in >= interval.upper.T) && (stream.T_out <= interval.lower.T)
        push!(interval.streams, stream.name => stream) 
        Q_contrib = (interval.upper.T - interval.lower.T)*stream.mcp
        push!(interval.contributions, stream.name => Q_contrib)
    end
end

function assign_stream!(interval::TemperatureInterval, stream::ColdStream)
    if (stream.T_in <= interval.lower.T) && (stream.T_out >= interval.upper.T)
        push!(interval.streams, stream.name => stream) 
        Q_contrib = (interval.upper.T - interval.lower.T)*stream.mcp
        push!(interval.contributions, stream.name => Q_contrib)
    end
end

"""
$(TYPEDSIGNATURES)

Assigns utility based on cheapest utility principle. 
Thus each hot utility is assigned only to the hottest interval that can accept heat from it.
Each cold utility is assigned to the coldest interval that can reject heat to it.
Note: `sorted_intervals` should be a vector of the entire vector of `TemperatureInterval`s sorted from hot to cold.
- `drop_infinite` does not assign the `hot/cold_utility` if `hot/cold_utility.Q == Inf`
"""
function assign_utility!(sorted_intervals::Vector{TemperatureInterval}, hot_utility::SimpleHotUtility; drop_infinite = false)
    drop_infinite && hot_utility.Q == Inf && return
    for interval in sorted_intervals # Has to be higher than upper? [POTENTIAL BUG]
        if hot_utility.T_out >= interval.lower.T
            push!(interval.streams, hot_utility.name => hot_utility)
            push!(interval.contributions, hot_utility.name => hot_utility.Q)
            break
        end
    end
end

function assign_utility!(sorted_intervals::Vector{TemperatureInterval}, cold_utility::SimpleColdUtility; drop_infinite = false)
    drop_infinite && cold_utility.Q == Inf && return
    for interval in reverse(sorted_intervals) # Get first interval that is hotter than cu.
        if interval.upper.T >= cold_utility.T_out
            push!(interval.streams, cold_utility.name => cold_utility)
            push!(interval.contributions, cold_utility.name => cold_utility.Q)
            break
        end
    end
end

"""
$(TYPEDSIGNATURES)

Assigns all streams and utilities specified in `prob` to the intervals `hot_sorted_intervals` and cold_sorted_intervals as appropriate. Calls `assign_utility!` and `assign_stream!` 
- `drop_infinite` does not assign the `hot/cold_utility` if `hot/cold_utility.Q == Inf`
"""
function assign_all_streams_and_utilities!(prob::ClassicHENSProblem, hot_sorted_intervals::Vector{TemperatureInterval}, cold_sorted_intervals::Vector{TemperatureInterval}; drop_infinite = false)
    for interval in hot_sorted_intervals
        for stream in values(prob.hot_streams_dict)
            assign_stream!(interval, stream)
        end
    end
    
    for interval in cold_sorted_intervals
        for stream in values(prob.cold_streams_dict)
            assign_stream!(interval, stream)
        end
    end
    
    for utility in values(prob.hot_utilities_dict)
        assign_utility!(hot_sorted_intervals, utility; drop_infinite = drop_infinite)
    end
    
    for utility in values(prob.cold_utilities_dict)
        assign_utility!(cold_sorted_intervals, utility; drop_infinite = drop_infinite)
    end
end


"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds two vectors of `TemperatureInterval`s one for the hot side, one for the cold side.
.
"""
mutable struct TransshipmentInterval{R <: Real}
    index::Int64
    hot_side::TemperatureInterval{R}
    cold_side::TemperatureInterval{R}
    ΔT_min::Float64
    @add_kwonly function TransshipmentInterval{R}(hot_side, cold_side, ΔT_min = 10.0) where {R}
        index = hot_side.index
        index == cold_side.index || error("Mismatched hot/cold indexes")
        (hot_side.upper.T - cold_side.upper.T == ΔT_min) && (hot_side.lower.T - cold_side.lower.T == ΔT_min) || error("Inconsistent ΔT_min")
        for (k,v) in hot_side.streams
            v isa Union{HotStream, SimpleHotUtility} || error("Non-hot stream assigned to hot_side of transshipment interval.")
        end
        for (k,v) in cold_side.streams
            v isa Union{ColdStream, SimpleColdUtility} || error("Non-cold stream assigned to cold_side of transshipment interval.")
        end
        new(index, hot_side, cold_side, ΔT_min)
    end
end

Base.show(io::IO, interval_tship::TransshipmentInterval) = print(io, "itv_$(interval_tship.index)")

function print_full(intervals::Vector{TransshipmentInterval}; digits = 1)
    for interval in intervals
        println("itv_", interval.index, " H: [", interval.hot_side.upper.T, ", ", interval.hot_side.lower.T, "]", " C: [", interval.cold_side.upper.T, ", ", interval.cold_side.lower.T, "]")
    end
end



"""
$(TYPEDSIGNATURES)

Generates the double-sided intervals necessary for a transshipment problem (as visualized in a heat cascade diagram).
"""
function generate_transshipment_intervals(prob::ClassicHENSProblem, ΔT_min = prob.ΔT_min)
    # Getting the hot and cold temps are tailored for each subproblem type. 
    hot_side_temps, cold_side_temps = Float64[], Float64[]
    for (k,v) in prob.hot_streams_dict
        push!(hot_side_temps, v.T_in, v.T_out)
    end
    
    for (k,v) in prob.cold_streams_dict
        push!(hot_side_temps, v.T_in + ΔT_min, v.T_out + ΔT_min)
    end
    
    for (k,v) in prob.hot_streams_dict
        push!(cold_side_temps, v.T_in - ΔT_min, v.T_out - ΔT_min)
    end
    
    for (k,v) in prob.cold_streams_dict
        push!(cold_side_temps, v.T_in, v.T_out)
    end
    
    hot_side = initialize_temperature_intervals(hot_side_temps)
    cold_side = initialize_temperature_intervals(cold_side_temps)
    
    length(hot_side) == length(cold_side) || error("Inconsistent lengths")

    assign_all_streams_and_utilities!(prob, hot_side, cold_side; drop_infinite = false)
    intervals_tship = TransshipmentInterval[]
    for i in 1:length(hot_side)
        push!(intervals_tship, TransshipmentInterval{Float64}(hot_side[i], cold_side[i], ΔT_min))
    end
    return intervals_tship
end

"""
$(TYPEDSIGNATURES)

Get the primary temperatures for the hot and cold side composite curves.
Adds `prob.results_dict[primary_temperatures] = (; hot_cc, cold_cc, hot_temps, cold_temps)`. `hot_cc` and `cold_cc` are ::Vector{TemperatureInterval} with both the `T` and `H` value for `interval.upper` and `interval.lower` fixed.
`hot_temps` and `cold_temps` gives the vector of primary temperatures only.
Note: If the amount of utility consumption in `prob` is `Inf` (e.g., prior to solution of the Minimum Utilities subproblem), then the corresponding utility stream is ignored.
"""
function get_primary_temperatures!(prob::ClassicHENSProblem; verbose = false)
    # Getting the hot and cold temps are tailored for each subproblem type. 
    hot_temps, cold_temps = Float64[], Float64[]
    for (k,v) in prob.hot_streams_dict
        push!(hot_temps, v.T_in, v.T_out)
    end
    
    for (k,v) in prob.cold_streams_dict
        push!(cold_temps, v.T_in, v.T_out)
    end
    
    for (k,v) in prob.hot_utilities_dict
        if v.Q != Inf
            push!(hot_temps, v.T_in, v.T_out)
        end
    end
    
    for (k,v) in prob.cold_utilities_dict
        if v.Q != Inf
            push!(cold_temps, v.T_in, v.T_out)
        end
    end

    hot_cc = initialize_temperature_intervals(hot_temps)
    cold_cc = initialize_temperature_intervals(cold_temps)

    verbose && "Hot intervals: $(length(hot_cc)), Cold intervals $(length(cold_cc)) \n"

    assign_all_streams_and_utilities!(prob, hot_cc, cold_cc; drop_infinite = true)
    calculate_enthalpies!([HotStream, SimpleHotUtility], hot_cc)
    calculate_enthalpies!([ColdStream, SimpleColdUtility], cold_cc)

    isapprox(first(hot_cc).upper.H, first(cold_cc).upper.H; atol = 1.0) || error("Primary temperature composite curve not balanced")
    prob.results_dict[:primary_temperatures] =  (; hot_cc, cold_cc, hot_temps, cold_temps)
    return 
end

"""
$(TYPEDSIGNATURES)

For a given sorted set of intervals `sorted_intervals`, populates the `upper.H` and `lower.H` values starting from `ref_enthalpy = 0.0`.
Only the contributions of streams of type listed in stream_type are counted.
"""
function calculate_enthalpies!(stream_type::Vector{DataType}, sorted_intervals::Vector{TemperatureInterval}; ref_enthalpy = 0.0)
    total_sum = ref_enthalpy
    for interval in reverse(sorted_intervals)
        interval.lower.H = total_sum
        total_sum += get_contribution(stream_type, interval)
        interval.upper.H = total_sum
    end
end

"""
$(TYPEDSIGNATURES)

Plots a single composite curve given `sorted_intervals`
"""
function plot_composite_curve(sorted_intervals::Vector{TemperatureInterval}; verbose = false, color = :blue, shape = :circle, ylabel = "T [°C or K]", xlabel = "Heat duty Q", kwargs...)
    T_vals, H_vals = [], []
    for interval in reverse(sorted_intervals)
        push!(T_vals, interval.lower.T, interval.upper.T)
        push!(H_vals, interval.lower.H, interval.upper.H)
    end
    verbose && println("T_vals: $(T_vals), H_vals = $(H_vals)")
    plot(H_vals, T_vals, ylabel = ylabel, xlabel = xlabel, color = color, shape = shape, legend = false, kwargs...)
end

"""
$(TYPEDSIGNATURES)

Get the secondary temperatures for the hot and cold side composite curves.
Loads `prob.results_dict[:primary_temperatures]` if available.
Adds `prob.results_dict[:secondary_temperatures] = (; hot_cc, cold_cc, hot_temps, cold_temps)`
"""
function get_secondary_temperatures!(prob::ClassicHENSProblem; verbose = false)
    hot_temps, merged_hot_temps, cold_temps, merged_cold_temps = Float64[], Float64[], Float64[], Float64[]
    if !haskey(prob.results_dict, :primary_temperatures)
        get_primary_temperatures!(prob)
    end

    hot_cc = prob.results_dict[:primary_temperatures].hot_cc
    cold_cc = prob.results_dict[:primary_temperatures].cold_cc

    # Secondary temps
    for cold_interval in cold_cc
        lower_match = get_enthalpy_match(cold_interval.lower, hot_cc)
        upper_match = get_enthalpy_match(cold_interval.upper, hot_cc)
        push!(hot_temps, lower_match.T, upper_match.T)
    end

    for hot_interval in hot_cc
        lower_match = get_enthalpy_match(hot_interval.lower, cold_cc)
        upper_match = get_enthalpy_match(hot_interval.upper, cold_cc)
        push!(cold_temps, lower_match.T, upper_match.T)
    end

    sort!(unique!(hot_temps), rev = true)
    sort!(unique!(cold_temps), rev = true)

    merged_hot_temps = vcat(prob.results_dict[:primary_temperatures].hot_temps, hot_temps) # Primary and secondary
    merged_cold_temps = vcat(prob.results_dict[:primary_temperatures].cold_temps, cold_temps)

    hot_cc = initialize_temperature_intervals(merged_hot_temps)
    cold_cc = initialize_temperature_intervals(merged_cold_temps)

    verbose && "Hot intervals: $(length(hot_cc)), Cold intervals $(length(cold_cc)) \n"

    assign_all_streams_and_utilities!(prob, hot_cc, cold_cc; drop_infinite = true)
    calculate_enthalpies!([HotStream, SimpleHotUtility], hot_cc)
    calculate_enthalpies!([ColdStream, SimpleColdUtility], cold_cc)

    isapprox(first(hot_cc).upper.H, first(cold_cc).upper.H; atol = 1.0) || error("Secondary temperature composite curve not balanced")

    prob.results_dict[:secondary_temperatures] = (; hot_cc, cold_cc, hot_temps, cold_temps)
    return
end

"""
$(TYPEDSIGNATURES)

Given a `point::Point`, returns the point `match_point::Point` that has the same enthalpy value interpolated from the `sorted_intervals` vector.
"""
function get_enthalpy_match(point::Point, sorted_intervals::Vector{TemperatureInterval}; verbose = false)
    verbose && println("Getting match for point $(point)")
    for interval in sorted_intervals # Can implement binary search here.
        if point.H >= floor(interval.lower.H) && point.H <= ceil(interval.upper.H) # Use floor and ceil to avoid precision issues. 
            match_T = interval.lower.T + (interval.upper.T - interval.lower.T)*(point.H - interval.lower.H)/(interval.upper.H - interval.lower.H)
            verbose && println("Match T: $(match_T), Match interval: T: [$(interval.lower.T), $(interval.upper.T)], H: [$(interval.lower.H), $(interval.upper.H)]")
            return Point{Float64}(T = match_T, H = point.H)
        end
    end
    error("Match point for $(point) not found")
end

"""
$(TYPEDSIGNATURES)

Get the tertiary temperatures for the hot and cold side composite curves for a given `prob`.
Args:
- EMAT: Specified Exchanger Minimum Approach Temperature.
Loads `prob.results_dict[:primary_temperatures]` and `prob.results_dict[:secondary_temperatures]` if available.   
Adds `prob.results_dict[:tertiary_temperatures] = (; hot_cc, cold_cc)`
"""
function get_tertiary_temperatures!(prob::ClassicHENSProblem, EMAT; verbose = false)
    # Load secondary temperatures.
    hot_temps, merged_hot_temps, cold_temps, merged_cold_temps = Float64[], Float64[], Float64[], Float64[]
    if !haskey(prob.results_dict, :secondary_temperatures)
        get_secondary_temperatures!(prob)
    end

    # Get tertiary temperatures:
    hot_targets, cold_targets = Float64[], Float64[]
    for hot_stream in values(merge(prob.hot_streams_dict, prob.hot_utilities_dict))
        push!(hot_targets, hot_stream.T_out)
    end
    for cold_stream in values(merge(prob.cold_streams_dict, prob.cold_utilities_dict))
        push!(cold_targets, cold_stream.T_out)
    end

    coldest_hot_target = minimum(hot_targets)
    hottest_cold_target = maximum(cold_targets)

    for cold_stream in values(merge(prob.cold_streams_dict, prob.cold_utilities_dict))
        T_tertiary = cold_stream.T_in + EMAT
        if T_tertiary >= coldest_hot_target 
            push!(hot_temps, T_tertiary)
        end
    end

    for hot_stream in values(merge(prob.hot_streams_dict, prob.hot_utilities_dict))
        T_tertiary = hot_stream.T_in - EMAT
        if T_tertiary <= hottest_cold_target 
            push!(cold_temps, T_tertiary)
        end
    end

    merged_hot_temps = vcat(prob.results_dict[:primary_temperatures].hot_temps, prob.results_dict[:secondary_temperatures].hot_temps, hot_temps) # Primary, secondary and tertiary
    merged_cold_temps = vcat(prob.results_dict[:primary_temperatures].cold_temps, prob.results_dict[:secondary_temperatures].cold_temps , cold_temps)

    hot_cc = initialize_temperature_intervals(merged_hot_temps)
    cold_cc = initialize_temperature_intervals(merged_cold_temps)

    verbose && "Hot intervals: $(length(hot_cc)), Cold intervals $(length(cold_cc)) \n"

    assign_all_streams_and_utilities!(prob, hot_cc, cold_cc; drop_infinite = true)
    calculate_enthalpies!([HotStream, SimpleHotUtility], hot_cc)
    calculate_enthalpies!([ColdStream, SimpleColdUtility], cold_cc)

    isapprox(first(hot_cc).upper.H, first(cold_cc).upper.H; atol = 1.0) || error("Tertiary temperature composite curve not balanced")
    prob.results_dict[:tertiary_temperatures] = (; hot_cc, cold_cc, hot_temps, cold_temps)
    return
end

"""
$(TYPEDSIGNATURES)

Get the quaternary temperatures for the hot and cold side composite curves for a given `prob`.
Args:
- EMAT: Specified Exchanger Minimum Approach Temperature.
Loads `prob.results_dict[:primary_temperatures]`, `prob.results_dict[:secondary_temperatures]` and `prob.results_dict[:tertiary_temperatures]` if available.   
Adds `prob.results_dict[:quaternary_temperatures] = (; hot_cc, cold_cc)`
"""
function get_quaternary_temperatures!(prob::ClassicHENSProblem, EMAT; verbose = false)
    # Load tertiary temperatures.
    merged_hot_temps, merged_cold_temps = Float64[], Float64[]
    if !haskey(prob.results_dict, :tertiary_temperatures)
        get_tertiary_temperatures!(prob, EMAT)
    end

    # Get quaternary temperatures:
    hot_temps = prob.results_dict[:secondary_temperatures].cold_temps .+ EMAT
    cold_temps = prob.results_dict[:secondary_temperatures].hot_temps .- EMAT

    # Merging 
    merged_hot_temps = vcat(prob.results_dict[:primary_temperatures].hot_temps, prob.results_dict[:secondary_temperatures].hot_temps, prob.results_dict[:tertiary_temperatures].hot_temps, hot_temps) # Primary, secondary, tertiary and quaternary
    merged_cold_temps = vcat(prob.results_dict[:primary_temperatures].cold_temps, prob.results_dict[:secondary_temperatures].cold_temps, prob.results_dict[:tertiary_temperatures].cold_temps, cold_temps)

    hot_cc = initialize_temperature_intervals(merged_hot_temps)
    cold_cc = initialize_temperature_intervals(merged_cold_temps)

    verbose && "Hot intervals: $(length(hot_cc)), Cold intervals $(length(cold_cc)) \n"

    assign_all_streams_and_utilities!(prob, hot_cc, cold_cc; drop_infinite = true)
    calculate_enthalpies!([HotStream, SimpleHotUtility], hot_cc)
    calculate_enthalpies!([ColdStream, SimpleColdUtility], cold_cc)

    isapprox(first(hot_cc).upper.H, first(cold_cc).upper.H; atol = 1.0) || error("Tertiary temperature composite curve not balanced")
    prob.results_dict[:quaternary_temperatures] = (; hot_cc, cold_cc, hot_temps, cold_temps)
    return
end

"""
$(TYPEDSIGNATURES)

Get LMTD for heat transfer between `hot_interval` and `cold_interval`. Return `smallest_value` with verbose warning if infeasible.
"""
function LMTD(hot_interval::TemperatureInterval, cold_interval::TemperatureInterval; verbose = false)
    ΔT_upper = hot_interval.upper.T - cold_interval.upper.T
    ΔT_lower = hot_interval.lower.T - cold_interval.lower.T
    if (ΔT_upper < 0.0) || (ΔT_lower < 0.0)
        verbose && @warn "Infeasible heat transfer between hot $(hot_interval) and $(cold_interval)"
        return smallest_value
    end

    if isapprox(ΔT_upper, ΔT_lower; atol = smallest_value) # For numerical robustness when temp differences are similar.
        return ΔT_upper
    end

    return (ΔT_upper - ΔT_lower)/log(ΔT_upper/ΔT_lower)
end

"""
$(TYPEDSIGNATURES)

Returns 1 if heat transfer between `hot_interval` and  `cold_interval` is feasible, 0 otherwise.

"""
function is_feasible(hot_interval::TemperatureInterval, cold_interval::TemperatureInterval; EMAT = 0.0)
    ΔT_upper = hot_interval.upper.T - cold_interval.upper.T
    ΔT_lower = hot_interval.lower.T - cold_interval.lower.T
    if (ΔT_upper >= EMAT) && (ΔT_lower >= EMAT)
        return true
    end
    return false
end





            







    #=
"""
$(TYPEDSIGNATURES)

Plots the hot-side composite curve. Assume intervals are already sorted e.g., attained from the `generate_heat_cascade_intervals` function. 
"""
function plot_hot_composite_curve(prob::ClassicHENSProblem; ref_enthalpy = 0.0, ylabel = "T [°C or K]", xlabel = "Heat duty Q", verbose = false)
    temperatures = Float64[]
    for 
    sorted_intervals = initialize_temperature_intervals([prob.hot_streams_dict])
    T_vals = Float64[last(sorted_intervals).T_hot_lower]
    Q_vals = Float64[ref_enthalpy]
    for interval in reverse(sorted_intervals)
        !verbose || println(interval.T_hot_upper, " ", interval.T_hot_lower, " ", interval.total_stream_heat_in)
        push!(T_vals, interval.T_hot_upper)
        ref_enthalpy += interval.total_stream_heat_in
        push!(Q_vals, ref_enthalpy)
    end
    plt = plot(Q_vals, T_vals, ylabel = ylabel, xlabel = xlabel, color = :red, shape = :circle, legend = false)
    return plt, Q_vals, T_vals
end

"""
$(TYPEDSIGNATURES)

Plots the cold-side composite curve. Assume intervals are already sorted e.g., attained from the `generate_heat_cascade_intervals` function. 
"""
function plot_cold_composite_curve(sorted_intervals::Vector{TemperatureInterval}; ref_enthalpy = 0.0, ylabel = "T [°C or K]", xlabel = "Heat duty Q", verbose = false)
    T_vals = Float64[last(sorted_intervals).T_cold_lower]
    Q_vals = Float64[ref_enthalpy]
    for interval in reverse(sorted_intervals)
        !verbose || println(interval.T_cold_upper, " ", interval.T_cold_lower, " ", interval.total_stream_heat_out)
        push!(T_vals, interval.T_cold_upper)
        ref_enthalpy += interval.total_stream_heat_out
        push!(Q_vals, ref_enthalpy)
    end
    plt = plot(Q_vals, T_vals, ylabel = ylabel, xlabel = xlabel, color = :blue, shape = :circle, legend = false)
    return plt, Q_vals, T_vals
end

"""
$(TYPEDSIGNATURES)

Plots the cold-side composite curve. Assume intervals are already sorted e.g., attained from the `generate_heat_cascade_intervals` function. 
"""
function plot_composite_curve(sorted_intervals::Vector{TemperatureInterval}; hot_ref_enthalpy = 0.0, cold_ref_enthalpy = 0.0, ylabel = "T [°C or K]", xlabel = "Heat duty Q")
    plt_hot, Q_hot_ints, T_hot_ints = CompHENS.plot_hot_composite_curve(sorted_intervals; ref_enthalpy = hot_ref_enthalpy); 
    plt_cold, Q_cold_int, T_cold_int = CompHENS.plot_cold_composite_curve(sorted_intervals; ref_enthalpy = cold_ref_enthalpy); 
    
    # No point manipulating the plots. plot!(plt1, plt2) unpacks the data anyway. 
    x_limits = (floor(Int64, min(hot_ref_enthalpy, cold_ref_enthalpy, minimum(Q_cold_int), minimum(Q_hot_ints))), ceil(Int64, round(max(maximum(Q_hot_ints), maximum(Q_cold_int)), sigdigits = 2)))
    y_limits = (floor(Int64, min(minimum(T_hot_ints), minimum(T_cold_int))), ceil(Int64, max(maximum(T_hot_ints), maximum(T_cold_int))))
    plot(Q_hot_ints, T_hot_ints, ylabel = ylabel, xlabel = xlabel, color = :red, shape = :circle, legend = false, xlims = x_limits, ylims = y_limits)
    plot!(Q_cold_int, T_cold_int, color = :blue, shape = :circle, legend = false)
end


#=
"""
$(TYPEDSIGNATURES)

Gets the total heat into interval from participating heat streams.
"""
function total_stream_heat_in(interval::TemperatureInterval)
    return sum(values(hot_streams_contribs)) 


end
=#



"""
$(TYPEDSIGNATURES)

Gets the contribution of a stream to the interval. Returns 0.0 if no contribution.
"""
function get_contribution(stream::String, interval::TemperatureInterval,)
    # [QN: Any advantages of multiple dispatch here? Or just use basic string comparison?
    # TODO: Potential bug returning 0.0 for nonexistent streams.].
    stream in keys(interval.hot_streams_contribs) && return interval.hot_streams_contribs[stream]
    stream in keys(interval.cold_streams_contribs) && return interval.cold_streams_contribs[stream]
    stream in keys(interval.hot_utilities_contribs) && return interval.hot_utilities_contribs[stream]
    stream in keys(interval.cold_utilities_contribs) && return interval.cold_utilities_contribs[stream]
    return 0.0
end
=#

