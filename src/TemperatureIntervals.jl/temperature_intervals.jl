"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds a single temperature interval. Each temperature interval is defined by a box containing the upper and lower temperatures on the hot and cold sides. The `*_dict` contain only the streams that participate in the interval.

TODO: Make this immutable, replace dictionaries with vectors. 
"""
mutable struct TemperatureInterval
    T_hot_upper::Float64
    T_hot_lower::Float64
    T_cold_upper::Float64
    T_cold_lower::Float64
    """The residual heat supplied to the interval by the adjacent hotter `TemperatureInterval`"""
    R_in::Float64
    """The residual heat supplied from the interval to the adjacent colder `TemperatureInterval`"""
    R_out::Float64
    hot_streams_dict = Dict{String, HotStream}()
    cold_streams_dict = Dict{String, ColdStream}()
    hot_utilities_dict = Dict{String, SimpleHotUtility}()
    cold_utilities_dict = Dict{String, SimpleColdUtility}()
    @add_kwonly function TemperatureInterval(T_hot_upper, T_hot_lower, T_cold_upper, T_cold_lower, R_in = 0.0, R_out = 0.0, hot_streams_dict = Dict{String, HotStream}(), cold_streams_dict = Dict{String, ColdStream}(), hot_utilities_dict = Dict{String, SimpleHotUtility}(), cold_utilities_dict = Dict{String, SimpleColdUtility}())
        R_in >= 0.0 && R_out >= 0.0 || error("Residuals to temperature interval negative")
        new(T_hot_upper, T_hot_lower, T_cold_upper, T_cold_lower, R_in, R_out, hot_streams_dict, cold_streams_dict, hot_utilities_dict, cold_utilities_dict)
    end
end

"""
$(TYPEDSIGNATURES)

Generates the intervals necessary for the heat cascade diagram. 
"""
function generate_heat_cascade_intervals(prob::ClassicHENSProblem)
    hot_side_temps, cold_side_temps = Set{Float64}(), Set{Float64}()
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
    return hot_side_temps, cold_side_temps
end

# generate_intervals_from_T_set