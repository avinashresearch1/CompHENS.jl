# Workflow using XLSX input:
# 1. Import necessary packages:
@time using CompHENS
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

@time solve_minimum_utilities_subproblem!(prob)
print_min_utils_pinch_points(prob)
@test prob.pinch_points == [(517.0, 497.0)]
@test isapprox(prob.hot_utilities_dict["ST"].Q, 244.13; atol = 1)
@test isapprox(prob.cold_utilities_dict["CW"].Q, 172.6; atol = 1)

# 6. Solve the minimum number of units subproblem:
@time solve_minimum_units_subproblem!(prob)
@test prob.min_units == 8

# 7.
EMAT = 2.5
@time generate_stream_matches!(prob, EMAT; add_units = 1, verbose = true)
print_HLD(prob)

#=
#i) get primary temperatures.

CompHENS.get_primary_temperatures!(prob)
print_full(prob.results_dict[:primary_temperatures].hot_cc)
print_full(prob.results_dict[:primary_temperatures].cold_cc)
plt1 = plot_composite_curve(prob.results_dict[:primary_temperatures].hot_cc; verbose = true, color = :red)
plt2 = plot_composite_curve(prob.results_dict[:primary_temperatures].cold_cc; color = :blue)

# ii. get secondary temperatures.
#CompHENS.get_enthalpy_match(pr_temps.cold_cc[5].upper, pr_temps.hot_cc; verbose = true)

CompHENS.get_secondary_temperatures!(prob; verbose = true)
print_full(prob.results_dict[:secondary_temperatures].hot_cc)
print_full(prob.results_dict[:secondary_temperatures].cold_cc)
#plt1 = plot_composite_curve(prob.results_dict[:secondary_temperatures].hot_cc; verbose = true, color = :red)
#plt2 = plot_composite_curve(prob.results_dict[:secondary_temperatures].cold_cc; color = :blue)

# iii. get tertiary temperatures
EMAT = 2.5
CompHENS.get_tertiary_temperatures!(prob, EMAT; verbose = true)
print_full(prob.results_dict[:tertiary_temperatures].hot_cc)
print_full(prob.results_dict[:tertiary_temperatures].cold_cc)
#plt1 = plot_composite_curve(prob.results_dict[:tertiary_temperatures].hot_cc; verbose = true, color = :red)
#plt2 = plot_composite_curve(prob.results_dict[:tertiary_temperatures].cold_cc; color = :blue)

# iii. get quaternary temperatures
EMAT = 2.5
CompHENS.get_quaternary_temperatures!(prob, EMAT; verbose = true)
print_full(prob.results_dict[:quaternary_temperatures].hot_cc)
print_full(prob.results_dict[:quaternary_temperatures].cold_cc)
println("Hot - Pri: $(prob.results_dict[:primary_temperatures].hot_temps) \n \n Sec: $(prob.results_dict[:secondary_temperatures].hot_temps) \n \n Ter: $(prob.results_dict[:tertiary_temperatures].hot_temps) \n \n Q: $(prob.results_dict[:quaternary_temperatures].hot_temps)") 
println("Cold - Pri: $(prob.results_dict[:primary_temperatures].cold_temps) \n \n Sec: $(prob.results_dict[:secondary_temperatures].cold_temps) \n \n Ter: $(prob.results_dict[:tertiary_temperatures].cold_temps) \n \n Q: $(prob.results_dict[:quaternary_temperatures].cold_temps)") 

#plt1 = plot_composite_curve(prob.results_dict[:tertiary_temperatures].hot_cc; verbose = true, color = :red)
#plt2 = plot_composite_curve(prob.results_dict[:tertiary_temperatures].cold_cc; color = :blue)

# 7.b. Generate stream matches:
@time generate_stream_matches!(prob::ClassicHENSProblem)
CompHENS.print_HLD(prob)
=#