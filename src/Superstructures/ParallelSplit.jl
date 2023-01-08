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

    # Sanity tests
    length(nodes) == length(prob.results_dict[:match_list][stream]) + 4 || error("Incorrect number of nodes generated!")
    for splitter in superstructure.splitters
        length(in_edges(splitter, superstructure)) == 1 || error("This function only works with superstructures where each splitter has a single incoming edge")
    end
    for mixer in superstructure.mixers
        length(out_edges(mixer, superstructure)) == 1 || error("This function only works with superstructures where each mixer has a single outgoing edge")
    end
    return ParallelSplit(nodes, edges)
end

# Macro to do this?
function construct_superstructure(stream::String, superstructure::ParallelSplit, prob::ClassicHENSProblem; verbose = true)
    ParallelSplit(stream, prob; verbose =  verbose)
end

function define_out_nodes(node::MajorSplitter, nodes::Vector{Node}, superstructure::ParallelSplit) 
    # For FloudasCiricGrossmann superstructure. Returns a `Vector{Node}` for all `out` nodes connected to by  in `node`.
    out_nodes = filter(nodes) do v # MajorSplitter connects to all MinorMixer nodes
        v isa HX
    end
    return out_nodes
end

function define_out_nodes(node::HX, nodes::Vector{Node}, superstructure::ParallelSplit) 
    out_nodes = filter(nodes) do v # HXs connects to MinorSplitter with the same match. 
        v isa MajorMixer
    end
    return out_nodes
end