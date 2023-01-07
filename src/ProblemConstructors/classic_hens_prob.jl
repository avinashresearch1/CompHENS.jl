using XLSX
using DataFrames
"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds the classical HENS problem with fixed stream data. 
[Refactor ideas:]
- Add the overall algorithm as a field (just like JuMP does)
"""
mutable struct ClassicHENSProblem  <: AbstractSynthesisProblem 
    names::Dict{Symbol, Vector{String}}
    hot_streams_dict::Dict{String, HotStream}
    cold_streams_dict::Dict{String, ColdStream}
    hot_utilities_dict::Dict{String, SimpleHotUtility}
    cold_utilities_dict::Dict{String, SimpleColdUtility}
    ΔT_min::Float64
    results_dict::Dict{Symbol,Any}
    @add_kwonly function ClassicHENSProblem(names, hot_streams_dict, cold_streams_dict, hot_utilities_dict = Dict{String, SimpleHotUtility}(), cold_utilities_dict = Dict{String, SimpleColdUtility}(), results_dict = Dict{Symbol,Any}(); ΔT_min = 10)
        new(names, hot_streams_dict, cold_streams_dict, hot_utilities_dict, cold_utilities_dict, ΔT_min, results_dict)
    end
end


function Base.show(io::IO, prob::ClassicHENSProblem)
    println(io, "Classic HEN synthesis problem:")
    println(io, "$(length(prob.hot_streams_dict)) hot streams, $(length(prob.cold_streams_dict)) cold streams, $(length(prob.hot_utilities_dict)) hot utilities, $(length(prob.cold_utilities_dict)) cold utilities.")
end

function Base.getproperty(prob::ClassicHENSProblem, sym::Symbol)
    if sym == :pinch_points || sym == :min_units # Put all fields I expect to `results_dict` add later here.
        sym in keys(prob.results_dict) && return prob.results_dict[sym]
        error("$(sym) not yet added to `prob.results_dict`")
    elseif sym == :hot_dict # Hot streams and Hot utilities
        return merge(prob.hot_utilities_dict, prob.hot_streams_dict)
    elseif sym == :hot_names # Hot streams and Hot utilities
        return union(prob.names[:hot_streams], prob.names[:hot_utilities])
    elseif sym == :cold_dict # Cold streams and Cold utilities
        return merge(prob.cold_utilities_dict, prob.cold_streams_dict)
    elseif sym == :cold_names # Cold streams and Cold utilities
        return union(prob.names[:cold_streams], prob.names[:cold_utilities])
    elseif sym == :streams_dict # All streams
        return merge(prob.hot_streams_dict, prob.cold_streams_dict)
    elseif sym == :stream_names # All streams
        return union(prob.names[:hot_streams], prob.names[:cold_streams])
    elseif sym == :utilities_dict # All utilities
        return merge(prob.hot_utilities_dict, prob.cold_utilities_dict)
    elseif sym == :utility_names # All utilities
        return union(prob.names[:hot_utilities], prob.names[:cold_utilities])
    elseif sym == :all_dict # All streams and utilities
        return merge(prob.hot_dict, prob.cold_dict)
    elseif sym == :all_names # All streams and utilities
        return union(prob.hot_names, prob.cold_names)   
    else # fallback to getfield
        return getfield(prob, sym)
    end
end

"""
$(TYPEDSIGNATURES)

Reads data from an XSLX file in `file_path_xlsx` and constructs a `ClassicHENSProblem`.

- **`file_path_xlsx`** needs to be a string that ends in .xlsx
- By default `ΔT_min` is set to 10 °C.
"""
function ClassicHENSProblem(file_path_xlsx::String; ΔT_min, verbose = false)
    stream_data_dfs = DataFrame[]
    XLSX.openxlsx(file_path_xlsx) do xf
        verbose && println("Num worksheets imported: $(XLSX.sheetcount(xf))")
        for sheet in XLSX.sheetnames(xf)
            verbose && println("Importing sheet: $(sheet)")
            push!(stream_data_dfs, DataFrame(XLSX.gettable(xf[sheet]; infer_eltypes=true)))
        end
    end

    length(stream_data_dfs) == 1 || error("Only 1 stream data worksheet can be imported for `ClassicHENSProblem`")
    return ClassicHENSProblem(stream_data_dfs[1]; ΔT_min)
end

"""
$(TYPEDSIGNATURES)

Given a `DataFrame` object (from a single excel worksheet), returns a `ClassicHENSProblem`
- By default `ΔT_min` is set to 10 °C.

Note that this function can also called for other problem types e.g., `MultiPeriodFlexibleHENSProblem`. 
"""
function ClassicHENSProblem(stream_data_df::DataFrame; ΔT_min)
    # Column names of XLSX interface:
    stream_label, type_label, t_in_label, t_out_label, heat_cap_label, heat_coeff_label, cost_label, forbidden_label, compulsory_label = "Stream", "Type [H, C, HU or CU]", "Supply Temperature T_in [C or K]", "Target Temperature T_out [C or K]", "Heat Capacity mCp [kW/C or kW/K]", "Heat transfer coefficient h [kW/m2C or kW/m2K]", "Cost [\$/kW]", "Forbidden Matches", "Compulsory Matches"
    additional_user_fields = Set{String}(["Cost [\$/kW]", "Forbidden Matches", "Compulsory Matches", "Maximum Temperature", "Mininum Temperature"])

    s_names = Dict{Symbol, Vector{String}}()
    hot_streams_dict = Dict{String, HotStream}()
    cold_streams_dict = Dict{String, ColdStream}()
    hot_utilities_dict = Dict{String, SimpleHotUtility}()
    cold_utilities_dict = Dict{String, SimpleColdUtility}()

    names_stream_df = names(stream_data_df)

    hot_stream_names, cold_stream_names, hot_utility_names, cold_utility_names = String[], String[], String[], String[]

    for row in eachrow(stream_data_df)
        # add_user_data field in all stream types.
        add_user_data = Dict{String, Any}(k => row[k] for k in additional_user_fields if (k in names_stream_df && !ismissing(row[k])))
        if row[type_label] == "H" # Hot stream
            push!(hot_stream_names, row[stream_label])
            hot_streams_dict[row[stream_label]] = HotStream(name = row[stream_label], T_in = row[t_in_label], T_out = row[t_out_label], mcp = row[heat_cap_label], h = row[heat_coeff_label], add_user_data = add_user_data)
        elseif row[type_label] == "C" # Cold stream
            push!(cold_stream_names, row[stream_label])
            cold_streams_dict[row[stream_label]] = ColdStream(name = row[stream_label], T_in = row[t_in_label], T_out = row[t_out_label], mcp = row[heat_cap_label], h = row[heat_coeff_label], add_user_data = add_user_data)  
        elseif row[type_label] == "HU" # Hot utility stream
            push!(hot_utility_names, row[stream_label])
            hot_utilities_dict[row[stream_label]] = SimpleHotUtility(name = row[stream_label], T_in = row[t_in_label], T_out = row[t_out_label], h = row[heat_coeff_label], add_user_data = add_user_data)
        elseif row[type_label] == "CU" # Cold utility stream
            push!(cold_utility_names, row[stream_label])
            cold_utilities_dict[row[stream_label]] = SimpleColdUtility(name = row[stream_label], T_in = row[t_in_label], T_out = row[t_out_label], h = row[heat_coeff_label], add_user_data = add_user_data)
        end
    end
    s_names[:hot_streams] = hot_stream_names
    s_names[:cold_streams] = cold_stream_names
    s_names[:hot_utilities] = hot_utility_names
    s_names[:cold_utilities] = cold_utility_names
    ## TODO: Logic to get `ΔT_min` from XLSX sheet.  
    return ClassicHENSProblem(s_names, hot_streams_dict, cold_streams_dict, hot_utilities_dict, cold_utilities_dict; ΔT_min)
end


function M(hot_stream::String, cold_stream::String, prob::ClassicHENSProblem) 
    merged_hot = merge(prob.hot_streams_dict, prob.hot_utilities_dict)
    merged_cold = merge(prob.cold_streams_dict, prob.cold_utilities_dict)
    return min(merged_hot[hot_stream].Q, merged_cold[cold_stream].Q)
end




