"""
$(TYPEDEF)
This denotes the type of the superstructure used. By definition, the superstructure is defined for a single stream.
The general procedure to define a superstructure is as follows:
1. Define the concrete structure and constructor which typically has two fields: `nodes` and  `edges` and `<: AbstractSuperstructure`
2. For each type of `Node`, it is generally required to specify methods for `define_out_nodes`
3. Define a `construct_superstructure` method that in turn calls the constructor of 1. This allows a user to specify the superstructure for several streams at once. 
"""
abstract type AbstractSuperstructure end

"""
$(TYPEDEF)
Holds stream superstructure types for which the stream has splits.
"""
abstract type AbstractSplitSuperstructure <: AbstractSuperstructure end 


abstract type Node end
Base.show(io::IO, node::Node) = print(io, "$(node.name)")

struct Source <: Node 
    name::String
end

abstract type Mixer <: Node end
struct MajorMixer <: Mixer 
    name::String
end

struct HX <: Node
    name::String
    match::String
end

abstract type Splitter <: Node end
struct MajorSplitter <: Splitter
    name::String
end

struct Sink <: Node
    name::String
end

struct Edge
    in::Node
    out::Node
end 

Base.show(io::IO, edge::Edge) = print(io, "($(edge.in.name), $(edge.out.name))")

# Helpful functions to make superstructure construction easier:
"""
$(TYPEDSIGNATURES)

Used to construct the same type of superstructure for all elements of the `streams` vector.
Returns: A `Dict{String, AbstractSuperstructure}`.
"""
function construct_superstructure(streams::Vector{String}, superstructure::AbstractSuperstructure, prob::ClassicHENSProblem; verbose = false)
    overall_network = Dict{String, AbstractSuperstructure}() # Can potentially use broadcast here. 
    for stream in streams 
        push!(overall_network, stream => construct_superstructure(stream, superstructure, prob; verbose = verbose))
    end
    return overall_network
end
        

"""
$(TYPEDEF)
    Returns a `Vector{Node}` for all `out` nodes connected to by in `node`. Needs to be defined for all node types for any new superstructure.
"""
function define_out_nodes(node::Source, nodes::Vector{Node}, superstructure::AbstractSplitSuperstructure)
    out_nodes = filter(nodes) do v
        v isa MajorSplitter
    end

    length(out_nodes) == 1 || error("Not compatible with FloudasCiricGrossmann")
    return out_nodes
end

function define_out_nodes(node::MajorMixer, nodes::Vector{Node}, superstructure::AbstractSplitSuperstructure) 
    out_nodes = filter(nodes) do v # A MajorMixer only connects to a Sink
        v isa Sink
    end
    length(out_nodes) == 1 || error("Not compatible with FloudasCiricGrossmann")
    return out_nodes
end

function define_out_nodes(node::Sink, nodes::Vector{Node}, superstructure::AbstractSplitSuperstructure) 
    return Node[]
end

function Base.getproperty(superstructure::AbstractSplitSuperstructure, sym::Symbol)   
    if sym == :splitters
        return filter(superstructure.nodes) do v
            v isa Splitter
        end  
    elseif sym == :mixers
        return filter(superstructure.nodes) do v
            v isa Mixer
        end  
    elseif sym == :hxs
        return filter(superstructure.nodes) do v
            v isa HX
        end
    elseif sym == :source
        return filter(superstructure.nodes) do v
            v isa Source
        end
    elseif sym == :sink
        return filter(superstructure.nodes) do v
            v isa Sink
        end
    else # fallback to getfield
        return getfield(superstructure, sym)
    end
end

"""
$(TYPEDEF)
Given a `node::Node`, get all the outgoing edges.
""" 
function out_edges(node::Node, superstructure::AbstractSplitSuperstructure)
    return filter(edge -> edge.in == node, superstructure.edges)
end 

"""
$(TYPEDEF)
Given a `node::Node`, get all the incoming edges.
""" 
function in_edges(node::Node, superstructure::AbstractSplitSuperstructure)
    return filter(edge -> edge.out == node, superstructure.edges)
end

"""
$(TYPEDEF)
Given a `node::Node`, get all the source nodes that have edges connecting to it.
""" 
function get_source_nodes(node::Node, superstructure::AbstractSplitSuperstructure)
    source_nodes = filter(superstructure.nodes) do v
        Edge(v,node) ∈ superstructure.edges
    end
    return source_nodes
end

"""
$(TYPEDEF)
Given a `node::Node`, get all the destination nodes that it connects to through an edge.
""" 
function get_destination_nodes(node::Node, superstructure::AbstractSplitSuperstructure)
    destination_nodes = filter(superstructure.nodes) do v
        Edge(node, v) ∈ superstructure.edges
    end
    return destination_nodes
end
