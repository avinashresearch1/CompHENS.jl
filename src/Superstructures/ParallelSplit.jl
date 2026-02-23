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

"""
Throw a clear warning + error for unsupported bounded-parallel superstructures.

The current NLP workflow is validated only with `FloudasCiricGrossmann()`.
"""
function _parallel_split_unsupported()
    @warn "BoundedParallel/ParallelSplit is unsupported for network generation. Use FloudasCiricGrossmann() for all streams and utilities."
    error("ParallelSplit is disabled. Use FloudasCiricGrossmann().")
end

function ParallelSplit(; verbose = true)
    _parallel_split_unsupported()
end

function ParallelSplit(stream::String, prob::ClassicHENSProblem; verbose = true)
    _parallel_split_unsupported()
end

# Macro to do this?
function construct_superstructure(stream::String, superstructure::ParallelSplit, prob::ClassicHENSProblem; verbose = true)
    _parallel_split_unsupported()
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
