struct MinorMixer <: Mixer
    name::String
    """match is the stream the HX immediately after the Mixer matches with. Set to `nothing` for Big Mixer"""
    match::String
end

struct MinorSplitter <: Splitter
    name::String
    """match is the stream the HX immediately before the Splitter matches with. Set to `nothing` for Big Splitter"""
    match::String
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

    # Sanity tests
    length(nodes) == 3*length(prob.results_dict[:match_list][stream]) + 4 || error("Incorrect number of nodes generated!")
    length(edges) == (2 + 4*length(prob.results_dict[:match_list][stream]) + length(prob.results_dict[:match_list][stream])*(length(prob.results_dict[:match_list][stream])-1)) || error("Incorrect number of edges generated!")
    
    for splitter in filter(v -> (v isa Splitter), nodes)
        length(filter(edge -> edge.out == splitter, edges)) == 1 || error("Each splitter only has single incoming edge")
    end

    for mixer in filter(v -> (v isa Mixer), nodes)
        length(filter(edge -> edge.in == mixer, edges)) == 1 || error("Each mixer only has a single outgoing edge")
    end
    
    return FloudasCiricGrossmann(nodes, edges)
end

function construct_superstructure(stream::String, superstructure::FloudasCiricGrossmann, prob::ClassicHENSProblem; verbose = true)
    FloudasCiricGrossmann(stream, prob; verbose =  verbose)
end

function define_out_nodes(node::MajorSplitter, nodes::Vector{Node}, superstructure::FloudasCiricGrossmann) 
    # For FloudasCiricGrossmann superstructure. Returns a `Vector{Node}` for all `out` nodes connected to by  in `node`.
    out_nodes = filter(nodes) do v # MajorSplitter connects to all MinorMixer nodes
        v isa MinorMixer
    end
    return out_nodes
end

function define_out_nodes(node::MinorMixer, nodes::Vector{Node}, superstructure::FloudasCiricGrossmann) 
    out_nodes = filter(nodes) do v # MinorMixer connects to HX with the same match. 
        v isa HX && node.match == v.match
    end
    return out_nodes
end

function define_out_nodes(node::HX, nodes::Vector{Node}, superstructure::FloudasCiricGrossmann) 
    out_nodes = filter(nodes) do v # HXs connects to MinorSplitter with the same match. 
        v isa MinorSplitter && node.match == v.match
    end
    return out_nodes
end

function define_out_nodes(node::MinorSplitter, nodes::Vector{Node}, superstructure::FloudasCiricGrossmann) 
    out_nodes = filter(nodes) do v # A MinorSplitter only connects to MajorMixer
        v isa MajorMixer || (v isa MinorMixer && node.match != v.match)
    end
    return out_nodes
end
