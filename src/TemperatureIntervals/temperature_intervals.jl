using Plots

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds a single temperature interval. Each temperature interval is defined by a box containing the upper and lower temperatures on the hot and cold sides. The `*_dict` contain only the streams that participate in the interval.

QN: Make this immutable?
"""
mutable struct TemperatureInterval
    index::Int64
    T_hot_upper::Float64
    T_hot_lower::Float64
    T_cold_upper::Float64
    T_cold_lower::Float64
    """The residual heat supplied to the interval by the adjacent hotter `TemperatureInterval`. Defined after solving optimization problem."""
    R_in::Float64
    """The residual heat supplied from the interval to the adjacent colder `TemperatureInterval`. Defined after solving optimization problem."""
    R_out::Float64
    """The heat contribution of the hot stream to the interval"""
    hot_streams_contribs::Dict{String, Float64}
    """The heat removed by cold stream from the interval"""
    cold_streams_contribs::Dict{String, Float64}
    """Hot utility entering system at interval. Note: Each HU can only enter at one interval"""
    hot_utilities_contribs::Dict{String, Float64}
    """Cold utility entering system at interval. Note: Each CU can only enter at one interval"""
    cold_utilities_contribs::Dict{String, Float64}
    @add_kwonly function TemperatureInterval(index, T_hot_upper, T_hot_lower, T_cold_upper, T_cold_lower, R_in = 0.0, R_out = 0.0, hot_streams_contribs = Dict{String, Float64}(), cold_streams_contribs = Dict{String, Float64}(), hot_utilities_contribs = Dict{String, Float64}(), cold_utilities_contribs = Dict{String, Float64}())
        R_in >= 0.0 && R_out >= 0.0 || error("Residuals to temperature interval negative")
        new(index, T_hot_upper, T_hot_lower, T_cold_upper, T_cold_lower, R_in, R_out, hot_streams_contribs, cold_streams_contribs, hot_utilities_contribs, cold_utilities_contribs)
    end
end

Base.show(io::IO, interval::TemperatureInterval) = print(io, "itv_$(interval.index)")
function Base.show(io::IO, intervals::Vector{TemperatureInterval})
    for interval in intervals
        println(io, "itv_", interval.index, ":  Hot side: [", interval.T_hot_upper, ", ", interval.T_hot_lower, "]  Cold side: [", interval.T_cold_upper, ", ", interval.T_cold_lower, "]")
    end
end

"""
$(TYPEDSIGNATURES)

Generates the intervals necessary for the heat cascade diagram. Currently only works with one hot and one cold utility.
TODO: Multiple utilities.
"""
function generate_heat_cascade_intervals(prob::ClassicHENSProblem)
    intervals = TemperatureInterval[]
    hot_side_temps, cold_side_temps = Float64[], Float64[]
    for (k,v) in prob.hot_streams_dict
        push!(hot_side_temps, v.T_in, v.T_out)
    end

    for (k,v) in prob.cold_streams_dict
        push!(hot_side_temps, v.T_in + prob.ΔT_min, v.T_out + prob.ΔT_min)
    end

    for (k,v) in prob.hot_streams_dict
        push!(cold_side_temps, v.T_in - prob.ΔT_min, v.T_out - prob.ΔT_min)
    end

    for (k,v) in prob.cold_streams_dict
        push!(cold_side_temps, v.T_in, v.T_out)
    end
    sort!(unique!(hot_side_temps), rev = true)
    sort!(unique!(cold_side_temps), rev = true)
    
    length(hot_side_temps) == length(cold_side_temps) && all(hot_side_temps .- cold_side_temps .== prob.ΔT_min) || error("Inconsistency in attaining sorted temperature intervals.")
    
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
                push!(hot_utilities_contribs, k => Inf)
            end
        end

        if i == length(hot_side_temps)-1
            for (k,v) in prob.cold_utilities_dict
                push!(cold_utilities_contribs, k => Inf)
            end
        end

        R_in, R_out = 0.0, 0.0 # Place holders. To be attained from optimizer.
        push!(intervals, TemperatureInterval(i, T_hot_upper, T_hot_lower, T_cold_upper, T_cold_lower, R_in, R_out, hot_streams_contribs, cold_streams_contribs, hot_utilities_contribs, cold_utilities_contribs))
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


