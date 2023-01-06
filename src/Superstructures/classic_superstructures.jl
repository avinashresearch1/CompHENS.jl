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
    #new([],[])
    return nodes
end


function _fcg_out_nodes(node::Source, nodes::Dict{String, Node})
    # For FloudasCiricGrossmann superstructure. Returns a `Vector{Node}` for all `out` nodes connected to by  in `node`.
    out_nodes = filter(nodes) do (k,v)
        v isa MajorMixer
    end

    length(out_nodes) == 1 || error("Not compatible with FloudasCiricGrossmann")
    return values(out_nodes)
end

function _fcg_out_nodes(node::MajorMixer, nodes::Dict{String, Node})
    # For FloudasCiricGrossmann superstructure. Returns a `Vector{Node}` for all `out` nodes connected to by  in `node`.
    out_nodes = filter(nodes) do (k,v)
        v isa MajorMixer
    end

    length(out_nodes) == 1 || "Not compatible with FloudasCiricGrossmann"
    return values(out_nodes)
end

