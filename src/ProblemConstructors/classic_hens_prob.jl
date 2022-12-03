using XLSX
using DataFrames
"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds the classical HENS problem with fixed stream data. 
"""
struct ClassicHENSProblem  <: AbstractSynthesisProblem 
    hot_streams_dict::Dict{String, HotStream}
    cold_streams_dict::Dict{String, ColdStream}
    hot_utilities_dict::Dict{String, SimpleHotUtility}
    cold_utilities_dict::Dict{String, SimpleColdUtility}
    ΔT_min::Float64
    @add_kwonly function ClassicHENSProblem(hot_streams_dict, cold_streams_dict, hot_utilities_dict = Dict{String, SimpleHotUtility}(), cold_utilities_dict = Dict{String, SimpleColdUtility}(); ΔT_min = 10)
        new(hot_streams_dict, cold_streams_dict, hot_utilities_dict, cold_utilities_dict, ΔT_min)
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

    hot_streams_dict = Dict{String, HotStream}()
    cold_streams_dict = Dict{String, ColdStream}()
    hot_utilities_dict = Dict{String, SimpleHotUtility}()
    cold_utilities_dict = Dict{String, SimpleColdUtility}()

    names_stream_df = names(stream_data_df)

    for row in eachrow(stream_data_df)
        # add_user_data field in all stream types.
        add_user_data = Dict{String, Any}(k => row[k] for k in additional_user_fields if (k in names_stream_df && !ismissing(row[k])))
        if row[type_label] == "H" # Hot stream
            hot_streams_dict[row[stream_label]] = HotStream(row[stream_label], row[t_in_label], row[t_out_label], row[heat_cap_label], row[heat_coeff_label], add_user_data)
        elseif row[type_label] == "C" # Cold stream
            cold_streams_dict[row[stream_label]] = ColdStream(row[stream_label], row[t_in_label], row[t_out_label], row[heat_cap_label], row[heat_coeff_label], add_user_data)  
        elseif row[type_label] == "HU" # Hot utility stream
            hot_utilities_dict[row[stream_label]] = SimpleHotUtility(row[stream_label], row[t_in_label], row[t_out_label], row[heat_coeff_label], add_user_data)
        elseif row[type_label] == "CU" # Cold utility stream
            cold_utilities_dict[row[stream_label]] = SimpleColdUtility(row[stream_label], row[t_in_label], row[t_out_label], row[heat_coeff_label], add_user_data)
        end
    end

    ## TODO: Logic to get `ΔT_min` from XLSX sheet.  
    return ClassicHENSProblem(hot_streams_dict, cold_streams_dict, hot_utilities_dict, cold_utilities_dict; ΔT_min)
end



