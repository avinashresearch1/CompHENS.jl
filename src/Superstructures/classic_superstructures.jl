"""
$(TYPEDEF)

Holds an abstract superstructure 
"""
abstract type AbstractSuperstructure end

abstract type Node end

struct Source <: Node end

abstract type Mixer <: Node end
struct MajorMixer <: Mixer end
struct MinorMixer <: Mixer
    """match is the stream the HX immediately after the Mixer matches with. Set to `nothing` for Big Mixer"""
    match::String
end

struct HX <: Node
    match::String
end

abstract type Splitter <: Node end
struct MajorSplitter <: Splitter end
struct MinorSplitter <: Splitter
    """match is the stream the HX immediately before the Splitter matches with. Set to `nothing` for Big Splitter"""
    match::String
end

struct Sink <: Node end

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
    nodes::Dict{String, Node}
    edges::Vector{Edge} 
end

function FloudasCiricGrossmann(stream::String, prob::ClassicHENSProblem; verbose = true)
    verbose && @info "Using Superstructure: Floudas, C.A., Ciric, A.R. and Grossmann, I.E., Automatic synthesis of optimum heat exchanger network configurations. AIChE Journal. 1986." 
    nodes = Dict{String, Node}()
    push!(nodes, "SO" => Source())
    push!(nodes, "BS" => MajorSplitter())
    for match in prob.results_dict[:match_list][stream]
        push!(nodes, "SM_$(match)" => MinorMixer(match))
        push!(nodes, "HX_$(match)" => HX(match))
        push!(nodes, "SS_$(match)" => MinorSplitter(match))
    end
    push!(nodes, "BM" => MajorMixer())
    push!(nodes, "SK" => Sink())

    edges = Edge[]
    for (k,v) in nodes
        for node in _fcg_out_nodes(v, nodes)
            push!(edges, Edge(v, node))
        end
    end
    length(nodes) == 3*length(prob.results_dict[:match_list][stream]) + 4 || error("Incorrect number of nodes generated!")
    return nodes, edges
end

# [TODO:] More elegant way to do this? Every other approach doesn't seem manual and not extensible to other superstructures.
function _fcg_out_nodes(node::Source, nodes::Dict{String, Node})
    # For FloudasCiricGrossmann superstructure. Returns a `Vector{Node}` for all `out` nodes connected to by  in `node`.
    out_nodes = filter(nodes) do (k,v)
        v isa MajorSplitter
    end

    length(out_nodes) == 1 || error("Not compatible with FloudasCiricGrossmann")
    return values(out_nodes)
end

function _fcg_out_nodes(node::MajorSplitter, nodes::Dict{String, Node}) 
    # For FloudasCiricGrossmann superstructure. Returns a `Vector{Node}` for all `out` nodes connected to by  in `node`.
    out_nodes = filter(nodes) do (k,v) # MajorSplitter connects to all MinorMixer nodes
        v isa MinorMixer
    end
    return values(out_nodes)
end

function _fcg_out_nodes(node::MinorMixer, nodes::Dict{String, Node}) 
    out_nodes = filter(nodes) do (k,v) # MinorMixer connects to HX with the same match. 
        v isa HX && node.match == v.match
    end
    return values(out_nodes)
end

function _fcg_out_nodes(node::HX, nodes::Dict{String, Node}) 
    out_nodes = filter(nodes) do (k,v) # HXs connects to MinorSplitter with the same match. 
        v isa MinorSplitter && node.match == v.match
    end
    return values(out_nodes)
end

function _fcg_out_nodes(node::MinorSplitter, nodes::Dict{String, Node}) 
    out_nodes = filter(nodes) do (k,v) # A MinorSplitter connects to every MinorMixer that does not have the same match. Also connects to MajorMixer
        (v isa MajorMixer) || (v isa MinorMixer && node.match != v.match)
    end
    return values(out_nodes)
end

function _fcg_out_nodes(node::MajorMixer, nodes::Dict{String, Node}) 
    out_nodes = filter(nodes) do (k,v) # A MajorMixer only connects to a Sink
        v isa Sink
    end
    length(out_nodes) == 1 || error("Not compatible with FloudasCiricGrossmann")
    return values(out_nodes)
end

function _fcg_out_nodes(node::Sink, nodes::Dict{String, Node}) 
    return Node[]
end

