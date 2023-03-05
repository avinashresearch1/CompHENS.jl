using Conda, PyCall

const nx = PyNULL()
const plt = PyNULL()
const back_pdf = PyNULL()


function __init__()
    copy!(nx, pyimport_conda("networkx", "networkx"))
    copy!(plt, pyimport_conda("matplotlib.pyplot", "matplotlib"))
    copy!(back_pdf, pyimport_conda("matplotlib.backends.backend_pdf", "matplotlib"))
end

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
function plot_HEN_streamwise(prob::ClassicHENSProblem, model::AbstractModel, overall_network::Dict{String, AbstractSuperstructure}, file_name; digits = 1)
    __init__()
    pdf = CompHENS.back_pdf.PdfPages(file_name)
    for stream in prob.all_names
        plt.close()
        (g, edge_labels, node_labels, position, node_size) = get_stream_graph(prob.all_dict[stream], prob, model, overall_network[stream]; digits = digits)
        nx.draw_networkx(g, position, node_size = node_size, with_labels=false, arrowstyle="->", node_shape = "s")
        nx.draw_networkx_labels(g, position, labels = node_labels, font_size = 3, horizontalalignment = "left", verticalalignment = "top")
        nx.draw_networkx_edge_labels(g, position, edge_labels = edge_labels, horizontalalignment = "center", verticalalignment = "center", font_size = 2, rotate=false)
        pdf.savefig()
    end
    pdf.close()
end

"""


Gets a NetworkX Digraph for each stream
    Returns: (g, edge_labels, node_labels, position, node_size)
"""
function get_stream_graph(stream::AbstractStream, prob::ClassicHENSProblem, model::AbstractModel, superstructure::AbstractSplitSuperstructure; digits = 1)
    g = nx.DiGraph()
    edge_labels = Dict()
    node_labels = Dict()

    # Add nodes to NetworkX graph
    for node in superstructure.nodes
        g.add_node(node.name)
        if node isa HX
            match = node.match
            if stream isa HotStream
                node_labels[node.name] = "match: $(match), Q = $(round((prob.results_dict[:Q][match, stream.name]), digits = digits))"  
            elseif stream isa ColdStream
                node_labels[node.name] = "match: $(match), Q = $(round((prob.results_dict[:Q][stream.name, match]), digits = digits))"
            end
        end
    end

    # Only add edges if f_edge is greater than 0
    for edge in superstructure.edges
        if value(model[:f][(stream.name, edge)]) > 0.0
            g.add_edge(edge.in.name, edge.out.name)
            edge_labels[(edge.in.name, edge.out.name)] = "m: $(round(value(model[:f][(stream.name, edge)]); digits = digits)) T: $(round(value(model[:t][(stream.name, edge)]); digits = digits))"
        end
    end

    num_matches = length(prob.results_dict[:HLD_list][stream.name])    
    # Copied from SeqHENS.jl
    # Setting the coordinates of the nodes
    position = Dict()
    coordinates = [] # Array to keep pushing paired coordinates to, will then go through and assign to dictionary
    bfix = 0.4 # Determines the vertical spacing

    push!(coordinates, (0,0))
    push!(coordinates, (1,0))
    for e in 1:num_matches
        push!(coordinates,(2,0-(bfix*(e-1))))
        push!(coordinates,(5,0-(bfix*(e-1))))
        push!(coordinates,(9,0-(bfix*(e-1))))
    end
    push!(coordinates,(10,0))
    push!(coordinates,(11,0))

    i = 0 # hack for now
    for node in superstructure.nodes
        i += 1
        position[node.name] = coordinates[i]
    end

    node_size = fill(1.0, length(superstructure.nodes))
    return(g, edge_labels, node_labels, position, node_size)
end

function get_stream_graph(stream::AbstractUtility, prob::ClassicHENSProblem, model::AbstractModel, superstructure::AbstractSplitSuperstructure; digits = 1)
    g = nx.DiGraph()
    edge_labels = Dict()
    node_labels = Dict()

    # Add nodes to NetworkX graph
    for node in superstructure.nodes
        g.add_node(node.name)
        if node isa HX
            match = node.match
            if stream isa HotStream
                node_labels[node.name] = "o: $(match), Q = $(round((prob.results_dict[:Q][match, stream.name]), digits = digits))"  
            elseif stream isa ColdStream
                node_labels[node.name] = "o: $(match), Q = $(round((prob.results_dict[:Q][stream.name, match]), digits = digits))"
            end
        end
    end

    # Only add edges if f_edge is greater than 0
    for edge in superstructure.edges
        if edge.in isa HX || edge.out isa HX # Only care about edges attached to HX
            g.add_edge(edge.in.name, edge.out.name)
            edge_labels[(edge.in.name, edge.out.name)] = "T: $(round(value(model[:t][(stream.name, edge)]); digits = digits))"
        end
    end

    num_matches = length(prob.results_dict[:HLD_list][stream.name])    
    # Copied from SeqHENS.jl
    # Setting the coordinates of the nodes
    position = Dict()
    coordinates = [] # Array to keep pushing paired coordinates to, will then go through and assign to dictionary
    bfix = 0.4 # Determines the vertical spacing

    push!(coordinates, (0,0))
    push!(coordinates, (1,0))
    for e in 1:num_matches
        push!(coordinates,(2,0-(bfix*(e-1))))
        push!(coordinates,(5,0-(bfix*(e-1))))
        push!(coordinates,(9,0-(bfix*(e-1))))
    end
    push!(coordinates,(10,0))
    push!(coordinates,(11,0))

    i = 0 # hack for now
    for node in superstructure.nodes
        i += 1
        position[node.name] = coordinates[i]
    end

    node_size = fill(1.0, length(superstructure.nodes))
    return(g, edge_labels, node_labels, position, node_size)
end
