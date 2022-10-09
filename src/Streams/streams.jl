abstract type AbstractStream end

"""
$(TYPEDEF)

Single hot stream
"""
mutable struct HotStream <: AbstractStream
    name::String
    T_in::Float64 # May eventually need parametric types when these can be variables. 
    T_out::Float64 
    mcp::Float64   
    h::Float64 # Film heat transfer coefficient
    info::Dict{String, Any} # Used to hold useful info moving forwards. `Any` here may ruin type inference?
end


"""
$(TYPEDEF)

Single cold stream
"""
mutable struct ColdStream <: AbstractStream
    name::String
    T_in::Float64 # May eventually need parametric types when these can be variables. 
    T_out::Float64 
    mcp::Float64   
    h::Float64 # Film heat transfer coefficient
    info::Dict{String, Any} # Used to hold useful info moving forwards. `Any` here may ruin type inference?
end

abstract type AbstractUtility end

"""
$(TYPEDEF)

Simple hot utility stream. Single stream, no configuration information.
"""
mutable struct SimpleHotUtility <: AbstractUtility
    name::String
    T_in::Float64 # May eventually need parametric types when these can be variables. 
    T_out::Float64   
    h::Float64 # Film heat transfer coefficient
    info::Dict{String, Any} # Used to hold useful info moving forwards. `Any` here may ruin type inference?
end

"""
$(TYPEDEF)

Simple cold utility stream. Single stream, no configuration information.
"""
mutable struct SimpleColdUtility <: AbstractUtility
    name::String
    T_in::Float64 # May eventually need parametric types when these can be variables. 
    T_out::Float64   
    h::Float64 # Film heat transfer coefficient
    info::Dict{String, Any} # Used to hold useful info moving forwards. `Any` here may ruin type inference?
end
