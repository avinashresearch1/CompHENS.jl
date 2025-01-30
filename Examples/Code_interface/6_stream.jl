# 1. Import necessary packages:
using CompHENS
using DataFrames, Plots, Test
using JuMP, HiGHS
using BenchmarkTools

# 2. Specifyx the dataFrame
df = DataFrame("Stream" => ["C1", "H2", "C3", "H4", "S", "CW"],
    "Type [H, C, HU or CU]" => ["C", "H", "C", "H", "HU", "CU"],
    "Supply Temperature T_in [C or K]" => [60, 160, 116, 249, 270, 38],
    "Target Temperature T_out [C or K]" => [160, 93, 260, 138, 269, 82],
    "Heat Capacity mCp [kW/C or kW/K]" => [7.62, 8.791, 6.0833, 10.54954, NaN, NaN],
    "Heat transfer coefficient h [kW/m2C or kW/m2K]" => [0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
    "Cost [\$/kW]" => [0, 0, 0, 0, 200.0, 20.0],
    "Forbidden Matches" => fill(nothing, 6),
    "Compulsory Matches" => fill(nothing, 6));
# df = DataFrame("Stream" => ["H1", "H2", "C1", "ST", "CW"],
#     "Type [H, C, HU or CU]" => ["H", "H", "C", "HU", "CU"],
#     "Supply Temperature T_in [C or K]" => [45, 50, 50, 210, 30],
#     "Target Temperature T_out [C or K]" => [40, 45, 60, 209, 40],
#     "Heat Capacity mCp [kW/C or kW/K]" => [10, 10, 0.1, NaN, NaN],
#     "Heat transfer coefficient h [kW/m2C or kW/m2K]" => [0.2, 0.2, 0.2, 0.2, 0.2],
#     "Cost [\$/kW]" => [0, 0, 0, 200.0, 20.0],
#     "Forbidden Matches" => fill(nothing, 5),
#     "Compulsory Matches" => fill(nothing, 5));

# 3. Construct a classic HENS problem
prob = ClassicHENSProblem(df; ΔT_min=5.0)
res = plot_composite_curve(prob; balanced=false, cold_ref_enthalpy=0, background_color=:transparent, foreground_color_text=RGB(0.95), foreground_color_axis=RGB(0.95), foreground_color_border=RGB(0.95), foreground_color_guide=RGB(0.95), linewidth=1.5, verbose=true, ylimit=(0, 220));
# balanced: get a balanced composite curve including utilities

# 4. Solve minimum utilities problem
# solve_minimum_utilities_subproblem!(prob) # minimize_utility_costs!()
@time solve_minimum_utilities_subproblem!(prob, verbose=true) # default to HIGHS_solver
# GLPK  232.000 μs (3277 allocations: 124.55 KiB)
# Clp   320.300 μs (3539 allocations: 144.77 KiB)
# HiGHS 644.700 μs (3358 allocations: 129.66 KiB)
# performance enhancements to 624.600 μs (3091 allocations: 111.84 KiB)
# Cbc   2.637 ms (3524 allocations: 145.45 KiB) 

res = plot_composite_curve(prob; balanced=true, background_color=:transparent, foreground_color_text=RGB(0.95), foreground_color_axis=RGB(0.95), foreground_color_border=RGB(0.95), foreground_color_guide=RGB(0.95), linewidth=1.5, verbose=true, ylimit=(0, 220));

println("Objective value = ", objective_value(prob.results_dict[:min_utils_model]), " USD")
print_min_utils_pinch_points(prob)
@test prob.pinch_points == [(125.0, 115.0)]
@test isapprox(prob.hot_utilities_dict["ST"].Q, 300.0; atol=1)
@test isapprox(prob.cold_utilities_dict["CW"].Q, 220.0; atol=1)

# 5. Solve the minimum number of units subproblem:
@time solve_minimum_units_subproblem!(prob, verbose=true)
# GLPK  816.100 μs (8842 allocations: 475.00 KiB)
# HiGHS 2.020 ms (8817 allocations: 483.30 KiB)
# Cbc   4.549 ms (8716 allocations: 514.66 KiB)
@test prob.min_units == 5

# 6. Generate stream matches
EMAT = 2.5 # Exchanger Minimum Approach Temperature.
@time generate_stream_matches!(prob, EMAT; add_units=1, verbose=true)
# GLPK  10.810 ms (191113 allocations: 12.25 MiB)
# HiGHS 17.093 ms (190990 allocations: 12.34 MiB)
# Cbc   25.547 ms (197424 allocations: 12.88 MiB)
prob.results_dict[:y]
prob.results_dict[:Q]

## 7. Network generation:
# Specify which superstructure to use for each stream
obj_func = CostScaledPaterson()
overall_network = merge(construct_superstructure(prob.stream_names, FloudasCiricGrossmann(), prob), construct_superstructure(prob.utility_names, ParallelSplit(), prob))
cost_coeff, scaling_coeff = 670, 0.83
optimizer = BARON.Optimizer

generate_network!(prob, EMAT; optimizer, overall_network, obj_func, verbose=true, cost_coeff, scaling_coeff)