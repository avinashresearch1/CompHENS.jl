using XLSX
using DataFrames
using DocStringExtensions

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds the overall multiperiod HENS problem. For each period a `ClassicHENSProblem` is defined.
The `period_streams_dict` is a mapping from the `period_name` to `ClassicHENSProblem`.
"""
mutable struct MultiPeriodFlexibleHENSProblem{I<:Integer}  <: AbstractSynthesisProblem 
    num_periods::I
    period_names::Vector{String}
    period_streams_dict::Dict{String, ClassicHENSProblem}
    results_dict::Dict{Symbol,Any}
end


"""
$(TYPEDSIGNATURES)

Reads data from an XSLX file in `file_path_xlsx` and constructs a `MultiPeriodFlexibleHENSProblem`.

- **`file_path_xlsx`** needs to be a string that ends in .xlsx
- By default `ΔT_min` is set to 10 °C.
"""
function MultiPeriodFlexibleHENSProblem(file_path_xlsx::String, num_periods::Integer; verbose = false, ΔT_min = 10)
    stream_data_dfs = Dict{String, DataFrame}()
    period_streams_dict = Dict{String, ClassicHENSProblem}()
    period_names = String[]
    results_dict = Dict{Symbol,Any}()
    
    XLSX.openxlsx(file_path_xlsx) do xf
        verbose && println("Num worksheets imported: $(XLSX.sheetcount(xf))")
        for sheet in XLSX.sheetnames(xf)
            verbose && println("Importing sheet: $(sheet)")
            push!(period_names, sheet)
            push!(stream_data_dfs, sheet => DataFrame(XLSX.gettable(xf[sheet]; infer_eltypes=true)))
        end
    end

    length(stream_data_dfs) == num_periods || error("Inconsistent number of imported worksheets and `num_periods`")
    for (k,v) in stream_data_dfs
        push!(period_streams_dict, k => ClassicHENSProblem(v; ΔT_min))
    end
    return MultiPeriodFlexibleHENSProblem(num_periods, period_names, period_streams_dict, results_dict)
end

function M(hot_stream::String, cold_stream::String, prob::MultiPeriodFlexibleHENSProblem) 
    M_all = []
    for (k,v) in prob.period_streams_dict
        push!(M_all, M(hot_stream, cold_stream, v))
    end
    return maximum(M_all)
end
