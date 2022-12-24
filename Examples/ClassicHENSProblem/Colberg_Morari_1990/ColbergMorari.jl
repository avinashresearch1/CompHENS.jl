# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS
using Plots
using JuMP
using HiGHS
using Test

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_ColbergMorari.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0, verbose = true)

# 4. Subdivide into intervals and attain the hot and cold composite curves.

#=
hot_ref_enthalpy, cold_ref_enthalpy = 0.0, 172.596
sorted_intervals = intervals
ylabel = "T [°C or K]"
xlabel = "Heat duty Q"

plt = CompHENS.plot_composite_curve(sorted_intervals; hot_ref_enthalpy, cold_ref_enthalpy, ylabel = "T [°C or K]", xlabel = "Heat duty Q")
ylims!((300,700))
=#

# 5. Solve minimum utilities problem

solve_minimum_utilities_subproblem!(prob)
@test prob.pinch_points == [(517.0, 497.0)]
@test isapprox(prob.hot_utilities_dict["ST"].Q, 244.13; atol = 1)
@test isapprox(prob.cold_utilities_dict["CW"].Q, 172.6; atol = 1)

# 6. Solve the minimum number of units subproblem:
solve_minimum_units_subproblem!(prob)
@test prob.min_units == 8

#get_primary_temperatures(prob)
# TODO: Further refactor the extraction of temperatures and assignment of streams to intervals using multiple dispatch. Worth it?
hot_temps, cold_temps = Float64[], Float64[]
for (k,v) in prob.hot_streams_dict
    push!(hot_temps, v.T_in, v.T_out)
end

for (k,v) in prob.cold_streams_dict
    push!(cold_temps, v.T_in, v.T_out)
end

for (k,v) in prob.hot_utilities_dict
    if v.Q != Inf
        push!(hot_temps, v.T_in, v.T_out)
    end
end

for (k,v) in prob.cold_utilities_dict
    if v.Q != Inf
        push!(cold_temps, v.T_in, v.T_out)
    end
end

hot_cc = initialize_temperature_intervals(hot_temps)
cold_cc = initialize_temperature_intervals(cold_temps)

for interval in hot_cc
    for stream in values(prob.hot_streams_dict)
        assign_stream!(interval, stream)
    end
end

for interval in cold_cc
    for stream in values(prob.cold_streams_dict)
        assign_stream!(interval, stream)
    end
end

for utility in values(prob.hot_utilities_dict)
    assign_utility!(hot_cc, utility)
end

for utility in values(prob.cold_utilities_dict)
    assign_utility!(cold_cc, utility)
end

@run assign_utility!(hot_cc, prob.hot_utilities_dict["ST"])

