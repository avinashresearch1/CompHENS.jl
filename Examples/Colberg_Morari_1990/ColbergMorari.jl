# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_ColbergMorari.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0)

    hot_side_temps, cold_side_temps = Float64[], Float64[]
    for (k,v) in prob.hot_streams_dict
        push!(hot_side_temps, v.T_in, v.T_out)
    end

    for (k,v) in prob.cold_streams_dict
        push!(hot_side_temps, v.T_in + prob.ΔT_min, v.T_out + prob.ΔT_min)
    end

    for (k,v) in prob.hot_streams_dict
        push!(cold_side_temps, v.T_in - prob.ΔT_min, v.T_out - prob.ΔT_min)
    end

    for (k,v) in prob.cold_streams_dict
        push!(cold_side_temps, v.T_in, v.T_out)
    end
    unique!(hot_side_temps)
    unique!(cold_side_temps)
    length(hot_side_temps) == length(cold_side_temps) || error("Inconsistent number of hot side and cold side temperatures.")
    
    return hot_side_temps, cold_side_temps
end
