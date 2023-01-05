"""
$(TYPEDEF)

Holds an abstract superstructure 
"""
abstract type AbstractSuperstructure end

abstract type Node end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Source node for a single stream
"""
struct Source <: Node
    name::String
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Mixer node for a single stream
"""
struct Mixer <: Node
    name::String
    """match is the stream the HX immediately after the Mixer matches with. Set to `nothing` for Big Mixer"""
    match::{Nothing, String}
    function Mixer(name, match = nothing)
        new(name, match)
    end
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Heat Exchanger node for a single stream
"""
struct HX <: Node
    name::String
    match::String
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Mixer node for a single stream
"""
struct Splitter <: Node
    name::String
    """match is the stream the HX immediately before the Splitter matches with. Set to `nothing` for Big Splitter"""
    match::{Nothing, String}
    function Splitter(name, match = nothing)
        new(name, match)
    end
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Sink node for a single stream
"""
struct Sink <: Node
    name::String
end

struct Edge
    in::Node
    out::Node
end 



"""
$(TYPEDEF)
For each stream, the following nodes are defined in the superstructure with abbreviations provided:
#  `nodes`:
- Source `SO`
- Major stream splitter `BS` (Big Splitter)
- Minor stream mixers `SM` (Small Mixer) - one in each branch corresponding to a match
- Heat Exchanger `HX` - one in each branch corresponding to a match
- Minor stream splitter `SS` (Small Splitter) - one in each branch corresponding to a match
- Major stream mixer `BM` (Big Mixer)
- Sink `SK` 
"""
struct FloudasCiricGrossmann <: AbstractSuperstructure
    nodes::Vector{Node}
    edges::Vector{Edge} 
end

function FloudasCiricGrossmann(stream::String, prob::ClassicHENSProblem; verbose = true)
    verbose && @info "Using Superstructure: Floudas, C.A., Ciric, A.R. and Grossmann, I.E., Automatic synthesis of optimum heat exchanger network configurations. AIChE Journal. 1986." 
    
    new([],[])
end



