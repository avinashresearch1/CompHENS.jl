"""
$(TYPEDSIGNATURES)
The visualization tools for the HEN uses the Python `NetworkX` and `MatPlotlib` packages. 
By default, OS-X and Windows users should PyCall configured to use the Miniconda environment installed at `.julia/conda/` (this can be attained by calling `Conda.ROOTENV`).
For Linux users: By default the default system installation is used, and it is necessary to over-ride this as follows:
- Set 
```
ENV["PYTHON"]=""
using Pkg
Pkg.build("PyCall")
```
Then **restart** the Julia process. 
"""
function plot_HEN_streamwise(args...; kwargs...)
    # Concrete plotting implementation was moved to optional extension
    # `CompHENSNetworkXPlotsExt` to avoid hard dependency on PyCall.
    error("Plotting requires the optional extension `CompHENSNetworkXPlotsExt` (install/load PyCall).")
end

"""


Gets a NetworkX Digraph for each stream
    Returns: (g, edge_labels, node_labels, position, node_size)
"""
function get_stream_graph(args...; kwargs...)
    # Concrete graph builder lives in optional extension
    # `CompHENSNetworkXPlotsExt`.
    error("Graph plotting requires the optional extension `CompHENSNetworkXPlotsExt` (install/load PyCall).")
end
