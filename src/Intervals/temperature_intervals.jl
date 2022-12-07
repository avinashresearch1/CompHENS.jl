using Plots

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds a single temperature interval. Each temperature interval is defined by a box containing the upper and lower temperatures on the hot and cold sides. The `*_contribs` contain only the streams that participate in the interval.
"""
mutable struct TemperatureInterval{R}
    index::Int64
    T_hot_upper::point{R}
    T_hot_lower::point{R}
    T_cold_upper::point{R}
    T_cold_lower::point{R}
    """The heat contribution of the hot stream to the interval"""
    hot_streams_contribs::Dict{String, R}
    """The heat removed by cold stream from the interval"""
    cold_streams_contribs::Dict{String, R}
    """Hot utility entering system at interval. Note: Each HU can only enter at one interval"""
    hot_utilities_contribs::Dict{String, R}
    """Cold utility entering system at interval. Note: Each CU can only enter at one interval"""
    cold_utilities_contribs::Dict{String, R}
    @add_kwonly function TemperatureInterval{R}(index, T_hot_upper, T_hot_lower, T_cold_upper, T_cold_lower, hot_streams_contribs = Dict{String, Float64}(), cold_streams_contribs = Dict{String, Float64}(), hot_utilities_contribs = Dict{String, Float64}(), cold_utilities_contribs = Dict{String, Float64}()) where {R}
        new(index, T_hot_upper, T_hot_lower, T_cold_upper, T_cold_lower, hot_streams_contribs, cold_streams_contribs, hot_utilities_contribs, cold_utilities_contribs)
    end
end

Base.show(io::IO, interval::TemperatureInterval) = print(io, "itv_$(interval.index)")

#= DEPRECATED gives strange bugs with JuMP set iteration to have members printed differently from collection.
function Base.show(io::IO, intervals::Vector{TemperatureInterval})
    for interval in intervals
        println(io, "itv_", interval.index, ":  Hot side: [", interval.T_hot_upper, ", ", interval.T_hot_lower, "]  Cold side: [", interval.T_cold_upper, ", ", interval.T_cold_lower, "]")
    end
end
=#
#=
function print_full(intervals::Vector{TemperatureInterval})
    for interval in intervals
        println("itv_", interval.index, ":  Hot side: [", interval.T_hot_upper.T, ", ", interval.T_hot_lower.T, "]  Cold side: [", interval.T_cold_upper.T, ", ", interval.T_cold_lower.T, "]")
    end
end
=#


"""
$(TYPEDSIGNATURES)

Generates the intervals necessary for the heat cascade diagram. Currently only works with one hot and one cold utility.
TODO: Multiple utilities.

Potential refactor: 
1) Do the sorting and first four lines of four loop that generates `Vector{TemperatureInterval}`. 
2) Write four different _get_stream_contribution(interval::TemperatureInterval, stream::<VariousStreamTypes>) here
"""
function generate_heat_cascade_intervals(prob::ClassicHENSProblem, ΔT_min = prob.ΔT_min)
    intervals = TemperatureInterval[]
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
    sort!(unique!(hot_side_temps), rev = true)
    sort!(unique!(cold_side_temps), rev = true)
    
    length(hot_side_temps) == length(cold_side_temps) && all(hot_side_temps .- cold_side_temps .== ΔT_min) || error("Inconsistency in attaining sorted temperature intervals.")
    
    for i in 1:length(hot_side_temps)-1
        hot_streams_contribs, cold_streams_contribs, hot_utilities_contribs, cold_utilities_contribs = Dict{String, Float64}(), Dict{String, Float64}(), Dict{String, Float64}(), Dict{String, Float64}()
        # hot_utilities_contribs, cold_utilities_contribs = Dict{String, Float64}(keys(prob.hot_utilities_dict) .=> 0.0), Dict{String, Float64}(keys(prob.cold_utilities_dict) .=> 0.0)
        T_hot_upper = hot_side_temps[i]
        T_hot_lower = hot_side_temps[i+1]
        T_cold_upper = cold_side_temps[i]
        T_cold_lower = cold_side_temps[i+1]
        
        # Filtering the participating streams. Would change for forbidden matches. No easy way to avoid inner for loop:
        for (k,v) in prob.hot_streams_dict
            if (v.T_in >= T_hot_upper) && (v.T_out <= T_hot_lower)
                Q_contrib = (T_hot_upper - T_hot_lower)*v.mcp
                push!(hot_streams_contribs, k => Q_contrib)
            end
        end

        for (k,v) in prob.cold_streams_dict
            if (v.T_in <= T_cold_lower) && (v.T_out >= T_cold_upper)
                Q_contrib = (T_hot_upper - T_hot_lower)*v.mcp
                push!(cold_streams_contribs, k => Q_contrib)
            end
        end

        if i == 1 # TODO: Extend for multiple utilities. 
            for (k,v) in prob.hot_utilities_dict
                push!(hot_utilities_contribs, k => v.Q)
            end
        end

        if i == length(hot_side_temps)-1
            for (k,v) in prob.cold_utilities_dict
                push!(cold_utilities_contribs, k => v.Q)
            end
        end

        push!(intervals, TemperatureInterval(i, T_hot_upper, T_hot_lower, T_cold_upper, T_cold_lower, hot_streams_contribs, cold_streams_contribs, hot_utilities_contribs, cold_utilities_contribs))
    end
    return intervals
end

"""
$(TYPEDSIGNATURES)

Plots the hot-side composite curve. Assume intervals are already sorted e.g., attained from the `generate_heat_cascade_intervals` function. 
"""
function plot_hot_composite_curve(sorted_intervals::Vector{TemperatureInterval}; ref_enthalpy = 0.0, ylabel = "T [°C or K]", xlabel = "Heat duty Q", verbose = false)
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

function Base.getproperty(interval::TemperatureInterval, sym::Symbol)
    if sym == :total_stream_heat_in
        return sum(values(interval.hot_streams_contribs))
    elseif sym == :total_stream_heat_out
        return sum(values(interval.cold_streams_contribs))
    else # fallback to getfield
        return getfield(interval, sym)
    end
end

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

