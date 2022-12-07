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
ΔT_min = 20.0
hot_side_temps, cold_side_temps = Float64[], Float64[]
for (k,v) in prob.hot_streams_dict
    push!(hot_side_temps, v.T_in, v.T_out)
end

for (k,v) in prob.cold_streams_dict
    push!(hot_side_temps, v.T_in + ΔT_min, v.T_out + ΔT_min)
end

for (k,v) in prob.hot_streams_dict
    push!(cold_side_temps, v.T_in - ΔT_min, v.T_out - ΔT_min)
end

for (k,v) in prob.cold_streams_dict
    push!(cold_side_temps, v.T_in, v.T_out)
end

hot_side = initialize_temperature_intervals(hot_side_temps)
cold_side = initialize_temperature_intervals(cold_side_temps)

for interval in hot_side
    for stream in values(prob.hot_streams_dict)
        assign_stream!(interval, stream)
    end
end

for interval in cold_side
    for stream in values(prob.cold_streams_dict)
        assign_stream!(interval, stream)
    end
end

for utility in values(prob.hot_utilities_dict)
    assign_utility!(hot_side, utility)
end

for utility in values(prob.cold_utilities_dict)
    assign_utility!(cold_side, utility)
end

intervals_tship = TransshipmentIntervals{Float64}(hot_side, cold_side, ΔT_min)


#=
hot_ref_enthalpy, cold_ref_enthalpy = 0.0, 172.596
sorted_intervals = intervals
ylabel = "T [°C or K]"
xlabel = "Heat duty Q"

plt = CompHENS.plot_composite_curve(sorted_intervals; hot_ref_enthalpy, cold_ref_enthalpy, ylabel = "T [°C or K]", xlabel = "Heat duty Q")
ylims!((300,700))
=#

# 5. Solve subproblem 1: minimum utilities. 

# Using formulation of Prob. 16.5 Biegler, Grossmann, Westerberg book. Pg. 533.

solve_minimum_utilities_subproblem!(prob)
@test prob.pinch_points == [(517.0, 497.0)]
@test isapprox(prob.hot_utilities_dict["ST"].Q, 244.13; atol = 1)
@test isapprox(prob.cold_utilities_dict["CW"].Q, 172.6; atol = 1)

solve_minimum_units_subproblem!(prob)
@test prob.min_units == 8


