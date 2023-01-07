"""
$(TYPEDEF)
This denotes the type of the superstructure used. By definition, the superstructure is defined for a single stream.
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
struct MinorMixer <: Mixer
    name::String
    """match is the stream the HX immediately after the Mixer matches with. Set to `nothing` for Big Mixer"""
    match::String
end

struct HX <: Node
    name::String
    match::String
end

abstract type Splitter <: Node end
struct MajorSplitter <: Splitter
    name::String
end
struct MinorSplitter <: Splitter
    name::String
    """match is the stream the HX immediately before the Splitter matches with. Set to `nothing` for Big Splitter"""
    match::String
end

struct Sink <: Node
    name::String
end

struct Edge
    in::Node
    out::Node
end 

Base.show(io::IO, edge::Edge) = print(io, "($(edge.in.name), $(edge.out.name)")

"""
$(TYPEDEF)
For each stream, a parallel split superstructure is defined. A single (major) split is made from the `Source` and a single (major) mix is made prior to the  `Sink`.

Naming convention for `nodes`. Type of the node (see below) `<type>_<stream_name>_<branch>` :
<type> names are used for nodes.
- Source `SO`
- Major stream splitter `BS` (Big Splitter)
- Heat Exchanger `HX` - one in each branch corresponding to a match
- Major stream mixer `BM` (Big Mixer)
- Sink `SK` 
"""
struct ParallelSplit <: AbstractSplitSuperstructure
    nodes::Union{Nothing, Vector{Node}}
    edges::Union{Nothing, Vector{Edge}} 
end

function ParallelSplit(; verbose = true)
    verbose && @info "Using Superstructure: Parallel split of streams once at start and merge once at end" 
    ParallelSplit(nothing, nothing)
end

function ParallelSplit(stream::String, prob::ClassicHENSProblem; verbose = true)
    nodes = Node[]
    push!(nodes, Source("SO_$(stream)"))
    push!(nodes, MajorSplitter("BS_$(stream)"))
    for match in prob.results_dict[:match_list][stream]
        push!(nodes,  HX("HX_$(stream)_$(match)", match))
    end
    push!(nodes,  MajorMixer("BM_$(stream)"))
    push!(nodes,  Sink("SK_$(stream)"))

    edges = Edge[]
    for v in nodes
        for node in define_out_nodes(v, nodes, ParallelSplit(; verbose = false))
            push!(edges, Edge(v, node))
        end
    end
    length(nodes) == length(prob.results_dict[:match_list][stream]) + 4 || error("Incorrect number of nodes generated!")
    return ParallelSplit(nodes, edges)
end


"""
$(TYPEDEF)
For each stream, the `FloudasCiricGrossmann` superstructure is defined.

Naming convention for `nodes`. Type of the node (see below) `<type>_<stream_name>_<branch>` :
<type> names are used for nodes.
- Source `SO`
- Major stream splitter `BS` (Big Splitter)
- Minor stream mixers `SM` (Small Mixer) - one in each branch corresponding to a match
- Heat Exchanger `HX` - one in each branch corresponding to a match
- Minor stream splitter `SS` (Small Splitter) - one in each branch corresponding to a match
- Major stream mixer `BM` (Big Mixer)
- Sink `SK` 
"""
struct FloudasCiricGrossmann <: AbstractSplitSuperstructure
    nodes::Union{Nothing, Vector{Node}}
    edges::Union{Nothing, Vector{Edge}} 
end

function FloudasCiricGrossmann(; verbose = true)
    verbose && @info "Using Superstructure: Floudas, C.A., Ciric, A.R. and Grossmann, I.E., Automatic synthesis of optimum heat exchanger network configurations. AIChE Journal. 1986." 
    FloudasCiricGrossmann(nothing, nothing)
end

function FloudasCiricGrossmann(stream::String, prob::ClassicHENSProblem; verbose = true)
    nodes = Node[]
    push!(nodes, Source("SO_$(stream)"))
    push!(nodes, MajorSplitter("BS_$(stream)"))
    for match in prob.results_dict[:match_list][stream]
        push!(nodes,  MinorMixer("SM_$(stream)_$(match)", match))
        push!(nodes,  HX("HX_$(stream)_$(match)", match))
        push!(nodes, MinorSplitter("SS_$(stream)_$(match)", match))
    end
    push!(nodes,  MajorMixer("BM_$(stream)"))
    push!(nodes,  Sink("SK_$(stream)"))

    edges = Edge[]
    for v in nodes
        for node in define_out_nodes(v, nodes, FloudasCiricGrossmann(; verbose = false))
            push!(edges, Edge(v, node))
        end
    end
    length(nodes) == 3*length(prob.results_dict[:match_list][stream]) + 4 || error("Incorrect number of nodes generated!")
    return FloudasCiricGrossmann(nodes, edges)
end

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
        

# Macro to do this?
function construct_superstructure(stream::String, superstructure::ParallelSplit, prob::ClassicHENSProblem; verbose = true)
    ParallelSplit(stream, prob; verbose =  verbose)
end

function construct_superstructure(stream::String, superstructure::FloudasCiricGrossmann, prob::ClassicHENSProblem; verbose = true)
    FloudasCiricGrossmann(stream, prob; verbose =  verbose)
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

function define_out_nodes(node::MajorSplitter, nodes::Vector{Node}, superstructure::AbstractSplitSuperstructure) 
    # For FloudasCiricGrossmann superstructure. Returns a `Vector{Node}` for all `out` nodes connected to by  in `node`.
    out_nodes = filter(nodes) do v # MajorSplitter connects to all MinorMixer nodes
        v isa MinorMixer
    end
    return out_nodes
end

function define_out_nodes(node::MinorMixer, nodes::Vector{Node}, superstructure::AbstractSplitSuperstructure) 
    out_nodes = filter(nodes) do v # MinorMixer connects to HX with the same match. 
        v isa HX && node.match == v.match
    end
    return out_nodes
end

function define_out_nodes(node::HX, nodes::Vector{Node}, superstructure::AbstractSplitSuperstructure) 
    out_nodes = filter(nodes) do v # HXs connects to MinorSplitter with the same match. 
        v isa MinorSplitter && node.match == v.match
    end
    return out_nodes
end

function define_out_nodes(node::MinorSplitter, nodes::Vector{Node}, superstructure::AbstractSplitSuperstructure) 
    out_nodes = filter(nodes) do v # A MinorSplitter only connects to MajorMixer
        v isa MajorMixer
    end
    return out_nodes
end

function define_out_nodes(node::MinorSplitter, nodes::Vector{Node}, superstructure::FloudasCiricGrossmann) 
    out_nodes = filter(nodes) do v # A MinorSplitter connects to every MinorMixer that does not have the same match. Also connects to MajorMixer
        (v isa MajorMixer) || (v isa MinorMixer && node.match != v.match)
    end
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
