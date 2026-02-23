using JuMP

"""
Heuristic warm-start builder for network generation NLPs.

Implementation notes:
1. Mirrors the combinatorial mCp flow at a practical level: prioritize utility-quality matches,
   attempt serial placement when feasible, and use paraT/alpha weighting for fallback parallel splits.
2. Produces starts for stream-edge temperatures/flows and match-level ΔT/LMTD terms.
3. Keeps all logic data-driven from `prob.results_dict[:HLD_list]` and `prob.results_dict[:Q]`.
"""

"""
Return heat load `Q` for a hot stream to one matched cold stream/utility.
"""
function _stream_heat_load(prob::ClassicHENSProblem, stream::HotStream, match_name::String)
    return prob.results_dict[:Q][match_name, stream.name]
end

"""
Return heat load `Q` for a cold stream to one matched hot stream/utility.
"""
function _stream_heat_load(prob::ClassicHENSProblem, stream::ColdStream, match_name::String)
    return prob.results_dict[:Q][stream.name, match_name]
end

"""
Return mCp for streams and a safe proxy for utilities.
Utilities do not own mCp in this model, so a positive proxy keeps split formulas stable.
"""
function _mcp_or_proxy(stream)
    if stream isa Union{HotStream, ColdStream}
        return max(stream.mcp, smallest_value)
    end
    return 1.0
end

"""
Order matches for a hot stream by utility preference and heat quality.

Rationale:
1. Cold utilities are prioritized at stream terminal positions.
2. Then sort by descending cold `T_out` to place tighter-quality cold matches earlier.
"""
function _ordered_matches(prob::ClassicHENSProblem, stream::HotStream)
    matches = copy(prob.results_dict[:HLD_list][stream.name])
    sort!(matches; by = m -> begin
        s = prob.all_dict[m]
        util_priority = s isa SimpleColdUtility ? 1 : 0
        quality = -s.T_out
        (util_priority, quality)
    end)
    return matches
end

"""
Order matches for a cold stream by utility preference and heat quality.

Rationale:
1. Hot utilities are prioritized at stream terminal positions.
2. Then sort by ascending hot `T_out` so lower-quality hot matches are allocated first.
"""
function _ordered_matches(prob::ClassicHENSProblem, stream::ColdStream)
    matches = copy(prob.results_dict[:HLD_list][stream.name])
    sort!(matches; by = m -> begin
        s = prob.all_dict[m]
        util_priority = s isa SimpleHotUtility ? 1 : 0
        quality = s.T_out
        (util_priority, quality)
    end)
    return matches
end

function _ordered_matches(prob::ClassicHENSProblem, stream::AbstractUtility)
    return copy(prob.results_dict[:HLD_list][stream.name])
end

"""
Return the unique HX node in one stream superstructure corresponding to one match.
"""
function _find_hx(superstructure::AbstractSplitSuperstructure, match_name::String)
    return only(filter(v -> v.match == match_name, superstructure.hxs))
end

"""
Return the unique edge from `in_node` to `out_node` inside a superstructure.
"""
function _find_edge(superstructure::AbstractSplitSuperstructure, in_node::Node, out_node::Node)
    return only(filter(e -> e.in == in_node && e.out == out_node, superstructure.edges))
end

"""
Propagate source-side stream temperature across outgoing splitter edges.

Why:
To keep physically consistent starts before nonlinear mixer constraints are enforced.
"""
function _set_splitter_temperatures!(t_start::Dict{Edge, Float64}, stream_name::String, superstructure::AbstractSplitSuperstructure)
    source_edge = only(out_edges(only(superstructure.source), superstructure))
    source_temp = t_start[source_edge]
    for splitter in superstructure.splitters
        edge_in = only(in_edges(splitter, superstructure))
        source_temp = get(t_start, edge_in, source_temp)
        for edge_out in out_edges(splitter, superstructure)
            t_start[edge_out] = source_temp
        end
    end
    return
end

"""
Check serial-feasibility for one hot stream across an ordered match list.

The check uses both upper and lower approach temperature constraints with EMAT.
"""
function _serial_feasible(prob::ClassicHENSProblem, stream::HotStream, matches::Vector{String}, EMAT)
    mcp = max(stream.mcp, smallest_value)
    current_T = stream.T_in
    for match in matches
        q = _stream_heat_load(prob, stream, match)
        other = prob.all_dict[match]
        hot_in = current_T
        hot_out = current_T - q / mcp

        if other isa ColdStream
            other_mcp = max(other.mcp, smallest_value)
            cold_in = other.T_in
            cold_out = cold_in + q / other_mcp
        elseif other isa SimpleColdUtility
            cold_in = other.T_in
            cold_out = other.T_out
        else
            return false
        end

        if !(hot_in >= cold_out + EMAT && hot_out >= cold_in + EMAT)
            return false
        end
        current_T = hot_out
    end
    return true
end

"""
Check serial-feasibility for one cold stream across an ordered match list.

The check uses both upper and lower approach temperature constraints with EMAT.
"""
function _serial_feasible(prob::ClassicHENSProblem, stream::ColdStream, matches::Vector{String}, EMAT)
    mcp = max(stream.mcp, smallest_value)
    current_T = stream.T_in
    for match in matches
        q = _stream_heat_load(prob, stream, match)
        other = prob.all_dict[match]
        cold_in = current_T
        cold_out = current_T + q / mcp

        if other isa HotStream
            other_mcp = max(other.mcp, smallest_value)
            hot_in = other.T_in
            hot_out = hot_in - q / other_mcp
        elseif other isa SimpleHotUtility
            hot_in = other.T_in
            hot_out = other.T_out
        else
            return false
        end

        if !(hot_in >= cold_out + EMAT && hot_out >= cold_in + EMAT)
            return false
        end
        current_T = cold_out
    end
    return true
end

function _serial_feasible(prob::ClassicHENSProblem, stream::AbstractUtility, matches::Vector{String}, EMAT)
    return false
end

"""
Compute combinatorial-style parallel split weights for one stream.

Formula implemented from extracted method:
1. `paraT = quality_dT + Q/mCp_other`
2. `alpha = Q/paraT`
3. `split = alpha / sum(alpha)`
"""
function _parallel_split_weights(prob::ClassicHENSProblem, stream::Union{HotStream, ColdStream}, matches::Vector{String}, EMAT)
    alpha = Float64[]
    for match in matches
        q = _stream_heat_load(prob, stream, match)
        other = prob.all_dict[match]
        other_mcp = _mcp_or_proxy(other)

        quality_dT = if stream isa HotStream
            stream.T_in - other.T_out
        else
            other.T_in - stream.T_in
        end
        paraT = max(quality_dT, EMAT / 10) + q / other_mcp
        push!(alpha, q / max(paraT, smallest_value))
    end

    total = sum(alpha)
    if total <= smallest_value
        qvals = [_stream_heat_load(prob, stream, m) for m in matches]
        total_q = max(sum(qvals), smallest_value)
        return [q / total_q for q in qvals]
    end
    return [a / total for a in alpha]
end

"""
Uniform split fallback for utility-only branches.
"""
function _parallel_split_weights(prob::ClassicHENSProblem, stream::AbstractUtility, matches::Vector{String}, EMAT)
    n = max(length(matches), 1)
    return fill(1.0 / n, n)
end

"""
Initialize per-edge temperature/flow starts with source/sink anchors.
"""
function _initialize_stream_defaults!(prob::ClassicHENSProblem, stream_name::String, superstructure::AbstractSplitSuperstructure)
    stream = prob.all_dict[stream_name]
    t_start = Dict{Edge, Float64}(edge => stream.T_in for edge in superstructure.edges)
    f_start = Dict{Edge, Float64}(edge => 0.0 for edge in superstructure.edges)

    source_edge = only(out_edges(only(superstructure.source), superstructure))
    sink_edge = only(in_edges(only(superstructure.sink), superstructure))
    t_start[source_edge] = stream.T_in
    t_start[sink_edge] = stream.T_out

    if stream isa Union{HotStream, ColdStream}
        f_start[source_edge] = stream.mcp
        f_start[sink_edge] = stream.mcp
    end

    return t_start, f_start
end

"""
Apply serial warm-starts over ordered matches for split superstructures.

Effect:
1. Sends full stream mCp through one active branch at a time.
2. Connects branch exits serially toward subsequent matches.
3. Updates edge temperatures through each exchanged heat load.
"""
function _apply_serial_starts!(prob::ClassicHENSProblem, stream::Union{HotStream, ColdStream}, superstructure::AbstractSplitSuperstructure, matches::Vector{String}, t_start::Dict{Edge, Float64}, f_start::Dict{Edge, Float64})
    mcp = max(stream.mcp, smallest_value)
    source_node = only(superstructure.source)
    bs_node = only(get_destination_nodes(source_node, superstructure))

    current_T = stream.T_in
    source_edge = only(out_edges(source_node, superstructure))
    t_start[source_edge] = stream.T_in

    for edge in out_edges(bs_node, superstructure)
        t_start[edge] = stream.T_in
    end

    for (i, match) in enumerate(matches)
        hx = _find_hx(superstructure, match)
        hx_in_edge = only(in_edges(hx, superstructure))
        hx_out_edge = only(out_edges(hx, superstructure))

        sm_node = only(get_source_nodes(hx, superstructure))
        bs_to_sm = _find_edge(superstructure, bs_node, sm_node)

        # In serial mode, only the first match receives flow directly from the
        # major splitter. Downstream matches are fed by upstream minor splitters.
        if i == 1
            f_start[bs_to_sm] = stream.mcp
        else
            f_start[bs_to_sm] = 0.0
        end
        f_start[hx_in_edge] = stream.mcp
        f_start[hx_out_edge] = stream.mcp

        t_start[hx_in_edge] = current_T

        q = _stream_heat_load(prob, stream, match)
        next_T = if stream isa HotStream
            current_T - q / mcp
        else
            current_T + q / mcp
        end
        t_start[hx_out_edge] = next_T

        ss_node = only(get_destination_nodes(hx, superstructure))
        if i < length(matches)
            next_hx = _find_hx(superstructure, matches[i + 1])
            next_sm = only(get_source_nodes(next_hx, superstructure))
            link_edge = _find_edge(superstructure, ss_node, next_sm)
            f_start[link_edge] = stream.mcp
            t_start[link_edge] = next_T
        else
            bm_node = only(superstructure.major_mixer)
            ss_to_bm = _find_edge(superstructure, ss_node, bm_node)
            f_start[ss_to_bm] = stream.mcp
            t_start[ss_to_bm] = next_T
        end

        current_T = next_T
    end

    sink_edge = only(in_edges(only(superstructure.sink), superstructure))
    t_start[sink_edge] = current_T
    f_start[sink_edge] = stream.mcp

    _set_splitter_temperatures!(t_start, stream.name, superstructure)
    return
end

"""
Apply parallel warm-starts over ordered matches using computed split fractions.

Effect:
1. Splits mCp across branches according to paraT/alpha weighting.
2. Computes each branch outlet temperature from assigned branch flow and heat load.
3. Mixes branch outlets to seed sink temperature.
"""
function _apply_parallel_starts!(prob::ClassicHENSProblem, stream::Union{HotStream, ColdStream}, superstructure::AbstractSplitSuperstructure, matches::Vector{String}, t_start::Dict{Edge, Float64}, f_start::Dict{Edge, Float64}, EMAT)
    split = _parallel_split_weights(prob, stream, matches, EMAT)
    source_node = only(superstructure.source)
    bs_node = only(get_destination_nodes(source_node, superstructure))
    bm_node = only(superstructure.major_mixer)

    source_edge = only(out_edges(source_node, superstructure))
    t_start[source_edge] = stream.T_in

    for edge in out_edges(bs_node, superstructure)
        t_start[edge] = stream.T_in
    end

    branch_T_out = Float64[]
    branch_f = Float64[]

    for (i, match) in enumerate(matches)
        hx = _find_hx(superstructure, match)
        hx_in_edge = only(in_edges(hx, superstructure))
        hx_out_edge = only(out_edges(hx, superstructure))

        # Keep branch-flow conservation exact at the major mixer.
        # A nonzero denominator safeguard is applied only in temperature update.
        f_branch = split[i] * stream.mcp
        q = _stream_heat_load(prob, stream, match)
        t_out = if stream isa HotStream
            stream.T_in - q / max(f_branch, smallest_value)
        else
            stream.T_in + q / max(f_branch, smallest_value)
        end

        if hx_in_edge.in == bs_node
            f_start[hx_in_edge] = f_branch
            t_start[hx_in_edge] = stream.T_in
        else
            sm_node = hx_in_edge.in
            bs_to_sm = _find_edge(superstructure, bs_node, sm_node)
            f_start[bs_to_sm] = f_branch
            t_start[bs_to_sm] = stream.T_in
            f_start[hx_in_edge] = f_branch
            t_start[hx_in_edge] = stream.T_in
        end

        f_start[hx_out_edge] = f_branch
        t_start[hx_out_edge] = t_out

        if hx_out_edge.out == bm_node
            f_start[hx_out_edge] = f_branch
            t_start[hx_out_edge] = t_out
        else
            ss_node = hx_out_edge.out
            ss_to_bm = _find_edge(superstructure, ss_node, bm_node)
            f_start[ss_to_bm] = f_branch
            t_start[ss_to_bm] = t_out
        end

        push!(branch_T_out, t_out)
        push!(branch_f, f_branch)
    end

    sink_edge = only(in_edges(only(superstructure.sink), superstructure))
    mixed_T = sum(branch_T_out[i] * branch_f[i] for i in eachindex(branch_f)) / max(sum(branch_f), smallest_value)
    t_start[sink_edge] = mixed_T
    f_start[sink_edge] = sum(branch_f)

    _set_splitter_temperatures!(t_start, stream.name, superstructure)
    return
end

"""
Seed utility stream temperatures around HX-adjacent edges.

Utilities do not carry `f` variables in this formulation.
"""
function _apply_utility_starts!(stream::AbstractUtility, superstructure::AbstractSplitSuperstructure, t_start::Dict{Edge, Float64})
    source_edge = only(out_edges(only(superstructure.source), superstructure))
    sink_edge = only(in_edges(only(superstructure.sink), superstructure))
    t_start[source_edge] = stream.T_in
    t_start[sink_edge] = stream.T_out

    for edge in superstructure.edges
        if edge.in isa HX
            t_start[edge] = stream.T_out
        elseif edge.out isa HX
            t_start[edge] = stream.T_in
        end
    end
    return
end

"""
Validate that candidate stream starts stay within declared model bounds.

This is used as a guardrail for parallel warm starts: if a split allocation
drives branch temperatures/flows out of bounds, we fall back to serial starts.
"""
function _starts_within_bounds(prob::ClassicHENSProblem, stream::Union{HotStream, ColdStream}, superstructure::AbstractSplitSuperstructure, t_start::Dict{Edge, Float64}, f_start::Dict{Edge, Float64}; atol = 1e-6)
    t_low, t_high = prob.results_dict[:T_bounds][stream.name]
    f_high = max(stream.mcp, 0.0)

    for edge in superstructure.edges
        t_val = get(t_start, edge, stream.T_in)
        if !isfinite(t_val) || t_val < t_low - atol || t_val > t_high + atol
            return false
        end

        f_val = get(f_start, edge, 0.0)
        if !isfinite(f_val) || f_val < -atol || f_val > f_high + atol
            return false
        end
    end

    return true
end

"""
Build edge-level warm starts for one stream.

Selection rule:
1. Utilities use utility-specific temperature seeding.
2. Process streams use serial starts when feasible.
3. Otherwise, process streams use parallel paraT/alpha starts.
"""
function _starts_for_stream!(prob::ClassicHENSProblem, stream_name::String, EMAT, superstructure::AbstractSplitSuperstructure)
    stream = prob.all_dict[stream_name]
    t_start, f_start = _initialize_stream_defaults!(prob, stream_name, superstructure)

    if stream isa AbstractUtility
        _apply_utility_starts!(stream, superstructure, t_start)
        return t_start, f_start
    end

    matches = _ordered_matches(prob, stream)
    if isempty(matches)
        return t_start, f_start
    end

    can_serial = superstructure isa FloudasCiricGrossmann && _serial_feasible(prob, stream, matches, EMAT)
    if can_serial
        _apply_serial_starts!(prob, stream, superstructure, matches, t_start, f_start)
    else
        _apply_parallel_starts!(prob, stream, superstructure, matches, t_start, f_start, EMAT)
        if !_starts_within_bounds(prob, stream, superstructure, t_start, f_start)
            _apply_serial_starts!(prob, stream, superstructure, matches, t_start, f_start)
        end
    end

    return t_start, f_start
end

"""
Create starts for match-level ΔT variables and optional LMTD variables.

Why:
The NLP objective and feasibility constraints depend on these match variables; supplying
consistent starts reduces first-iteration infeasibility and improves IPOPT stability.
"""
function _match_temperature_starts!(prob::ClassicHENSProblem, model::AbstractModel, HLD_list, overall_network, EMAT)
    start_vals = Dict{VariableRef, Float64}()

    has_T_lmtd = haskey(JuMP.object_dictionary(model), :T_LMTD)
    for match in HLD_list
        hot_name, cold_name = match
        hot_super = overall_network[hot_name]
        cold_super = overall_network[cold_name]
        hot_hx = _find_hx(hot_super, cold_name)
        cold_hx = _find_hx(cold_super, hot_name)

        hot_hx_in = only(in_edges(hot_hx, hot_super))
        hot_hx_out = only(out_edges(hot_hx, hot_super))
        cold_hx_in = only(in_edges(cold_hx, cold_super))
        cold_hx_out = only(out_edges(cold_hx, cold_super))

        t_hot_in = JuMP.start_value(model[:t][(hot_name, hot_hx_in)])
        t_hot_out = JuMP.start_value(model[:t][(hot_name, hot_hx_out)])
        t_cold_in = JuMP.start_value(model[:t][(cold_name, cold_hx_in)])
        t_cold_out = JuMP.start_value(model[:t][(cold_name, cold_hx_out)])

        dt_u = max(t_hot_in - t_cold_out, EMAT + 1e-3)
        dt_l = max(t_hot_out - t_cold_in, EMAT + 1e-3)

        start_vals[model[:ΔT_upper][match]] = dt_u
        start_vals[model[:ΔT_lower][match]] = dt_l

        if has_T_lmtd
            lmtd = smallest_value + (2 / 3) * sqrt(dt_u * dt_l) + (1 / 6) * (dt_u + dt_l)
            start_vals[model[:T_LMTD][match]] = max(lmtd, smallest_value)
        end
    end

    return start_vals
end

"""
Public entry-point: generate warm starts for network-generation NLP.

Returned structure:
`(; v_names::Vector{VariableRef}, v_starts::Vector{Float64})`

This function sets starts directly on model variables and also returns them to preserve
the existing `initial_values` interface used by `generate_network!`.
"""
function build_network_start_values(prob::ClassicHENSProblem, model::AbstractModel, EMAT, overall_network::Dict{String, AbstractSuperstructure}, HLD_list; verbose = false)
    v_names = VariableRef[]
    v_starts = Float64[]

    for stream_name in prob.all_names
        superstructure = overall_network[stream_name]
        superstructure isa AbstractSplitSuperstructure || continue

        t_start, f_start = _starts_for_stream!(prob, stream_name, EMAT, superstructure)

        for edge in superstructure.edges
            t_var = model[:t][(stream_name, edge)]
            t_val = get(t_start, edge, prob.all_dict[stream_name].T_in)
            push!(v_names, t_var)
            push!(v_starts, t_val)
            JuMP.set_start_value(t_var, t_val)

            if prob.all_dict[stream_name] isa AbstractStream
                f_var = model[:f][(stream_name, edge)]
                f_val = get(f_start, edge, 0.0)
                push!(v_names, f_var)
                push!(v_starts, f_val)
                JuMP.set_start_value(f_var, f_val)
            end
        end
    end

    match_starts = _match_temperature_starts!(prob, model, HLD_list, overall_network, EMAT)
    for (v, s) in match_starts
        push!(v_names, v)
        push!(v_starts, s)
        JuMP.set_start_value(v, s)
    end

    verbose && @info "Generated starting values from combinatorial mCp heuristics" n_vars = length(v_names)
    return (; v_names, v_starts)
end
