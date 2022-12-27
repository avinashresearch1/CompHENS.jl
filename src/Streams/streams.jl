abstract type AbstractStream end


"""
$(TYPEDEF)
$(TYPEDFIELDS)

Single hot stream
"""
mutable struct HotStream <: AbstractStream
    name::String # May eventually need parametric types when these can be variables.
    T_in::Float64
    T_out::Float64 
    mcp::Float64   
    h::Float64 # Film heat transfer coefficient
    add_user_data::Dict{String, Any} # Used to hold useful information moving forwards. `Any` here may ruin type inference?
    calc::Dict{String, Float64} # New data calculated from input data
    @add_kwonly function HotStream(name, T_in, T_out, mcp, h, add_user_data = Dict{String, Any}(), calc = Dict{String, Float64}())
        T_in isa Real && T_out isa Real && mcp isa Real && h isa Real || error("Input data contains a non-real number")
        T_in >= T_out || error("Supply and Target temperature don't match stream type")
        mcp > smallest_value && h > smallest_value || error("mcp or h values infeasible")
        new(name, Float64(T_in), Float64(T_out), Float64(mcp), Float64(h), add_user_data, calc)
    end 
end


"""
$(TYPEDEF)
$(TYPEDFIELDS)

Single cold stream
"""
mutable struct ColdStream <: AbstractStream
    name::String
    T_in::Float64
    T_out::Float64 
    mcp::Float64   
    h::Float64
    add_user_data::Dict{String, Any} 
    calc::Dict{String, Float64}
    @add_kwonly function ColdStream(name, T_in, T_out, mcp, h, add_user_data = Dict{String, Any}(), calc = Dict{String, Float64}())
        T_in isa Real && T_out isa Real && mcp isa Real && h isa Real || error("Input data contains a non-real number")
        T_in <= T_out || error("Supply and Target temperature don't match stream type")
        mcp > smallest_value && h > smallest_value || error("mcp or h values infeasible")
        new(name, Float64(T_in), Float64(T_out), Float64(mcp), Float64(h), add_user_data, calc)
    end 
end

abstract type AbstractUtility end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Simple hot utility stream. Single stream, no configuration information.
"""
mutable struct SimpleHotUtility <: AbstractUtility
    name::String
    T_in::Float64  
    T_out::Float64   
    h::Float64 
    """Duty initialized to `Inf`. Fixed during solution."""
    Q::Float64 
    add_user_data::Dict{String, Any}
    calc::Dict{String, Float64}
    @add_kwonly function SimpleHotUtility(name, T_in, T_out, h, Q = Inf, add_user_data = Dict{String, Any}(), calc = Dict{String, Float64}())
        T_in isa Real && T_out isa Real && h isa Real || error("Input data contains a non-real number")
        h > smallest_value || error("h value infeasible")
        new(name, Float64(T_in), Float64(T_out), Float64(h), Q, add_user_data, calc)
    end
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Simple cold utility stream. Single stream, no configuration information.
"""
mutable struct SimpleColdUtility <: AbstractUtility
    name::String
    T_in::Float64 
    T_out::Float64   
    h::Float64
    """Duty initialized to `Inf`. Fixed during solution."""
    Q::Float64 
    add_user_data::Dict{String, Any}
    calc::Dict{String, Float64}
    @add_kwonly function SimpleColdUtility(name, T_in, T_out, h, Q = Inf, add_user_data = Dict{String, Any}(), calc = Dict{String, Float64}())
        T_in isa Real && T_out isa Real && h isa Real || error("Input data contains a non-real number")
        h > smallest_value || error("h value infeasible")
        new(name, Float64(T_in), Float64(T_out), Float64(h), Q, add_user_data, calc)
    end 
end

"""
$(TYPEDSIGNATURES)

Returns the overall heat transfer coefficient between the `hot_stream` and the  `cold_stream`.

"""
function U(hot_stream::Union{HotStream, SimpleHotUtility}, cold_stream::Union{ColdStream, SimpleColdUtility})
    u_ij = hot_stream.h*cold_stream.h/(hot_stream.h + cold_stream.h)
    return max(u_ij, smallest_value) # Numerically don't want to return a 0.
end
