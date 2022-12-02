# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS
using Plots

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_ColbergMorari.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0)
intervals = CompHENS.generate_heat_cascade_intervals(prob)
hot_ref_enthalpy, cold_ref_enthalpy = 0.0, 200.0
sorted_intervals = intervals
ylabel = "T [°C or K]"
xlabel = "Heat duty Q"

plt = CompHENS.plot_composite_curve(sorted_intervals; hot_ref_enthalpy, cold_ref_enthalpy, ylabel = "T [°C or K]", xlabel = "Heat duty Q")


plt_hot, Q_hot_ints, T_hot_ints = CompHENS.plot_hot_composite_curve(sorted_intervals; ref_enthalpy = hot_ref_enthalpy); 
plt_cold, Q_cold_int, T_cold_int = CompHENS.plot_cold_composite_curve(sorted_intervals; ref_enthalpy = cold_ref_enthalpy); 
x_limits = (floor(Int64, (min(hot_ref_enthalpy, cold_ref_enthalpy, minimum(Q_cold_int), minimum(Q_hot_ints)))), ceil(Int64, round(max(maximum(Q_hot_ints), maximum(Q_cold_int)), sigdigits = 2)))
plot(Q_hot_ints, T_hot_ints, ylabel = ylabel, xlabel = xlabel, color = :red, shape = :circle, legend = false, xlims = x_limits)
plot!(Q_cold_int, T_cold_int, color = :blue, shape = :circle, legend = false)

x_limits = (floor(Int64, (min(hot_ref_enthalpy, cold_ref_enthalpy, minimum(Q_cold_int), minimum(Q_hot_ints)))), ceil(Int64, round(max(maximum(Q_hot_ints), maximum(Q_cold_int)), sigdigits = 2)))
xlims!(plt, (0,1800))
