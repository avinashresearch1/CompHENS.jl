# 1. Import necessary packages:
using CompHENS
using DataFrames, Plots, Test
using JuMP, HiGHS, BARON
using BenchmarkTools

# 2. Specify the dataFrame
df = DataFrame("Stream" => ["H1", "H2", "C1", "C2", "ST", "CW"],
    "Type [H, C, HU or CU]" => ["H", "H", "C", "C", "HU", "CU"],
    "Supply Temperature T_in [C or K]" => [175, 125, 20, 40, 180, 15],
    "Target Temperature T_out [C or K]" => [45, 65, 155, 112, 179, 25],
    "Heat Capacity mCp [kW/C or kW/K]" => [10, 40, 20, 15, NaN, NaN],
    "Heat transfer coefficient h [kW/m2C or kW/m2K]" => [0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
    "Cost [\$/kW]" => [0, 0, 0, 0, 200.0, 20.0],
    "Forbidden Matches" => [nothing, nothing, nothing, nothing, nothing, nothing],
    "Compulsory Matches" => [nothing, nothing, nothing, nothing, nothing, nothing]);

# 3. Construct a classic HENS problem
prob = ClassicHENSProblem(df; ΔT_min=10.0)
res = plot_composite_curve(prob; balanced=false, cold_ref_enthalpy=0.0);
display(res.plt)

# 4. Solve minimum utilities problem
solve_minimum_utilities_subproblem!(prob)
@btime solve_minimum_utilities_subproblem!(prob) # default to HIGHS_solver
# GLPK  232.000 μs (3277 allocations: 124.55 KiB)
# Clp   320.300 μs (3539 allocations: 144.77 KiB)
# HiGHS 644.700 μs (3358 allocations: 129.66 KiB)
# performance enhancements to 624.600 μs (3091 allocations: 111.84 KiB)
# Cbc   2.637 ms (3524 allocations: 145.45 KiB) 
res = plot_composite_curve(prob; balanced=true);
display(res.plt)
println("Objective value = ", objective_value(prob.results_dict[:min_utils_model]), " USD")
print_min_utils_pinch_points(prob)
@test prob.pinch_points == [(125.0, 115.0)]
@test isapprox(prob.hot_utilities_dict["ST"].Q, 300.0; atol=1)
@test isapprox(prob.cold_utilities_dict["CW"].Q, 220.0; atol=1)

# 5. Solve the minimum number of units subproblem:
@time solve_minimum_units_subproblem!(prob)
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

## 7. Network generation:
# Specify which superstructure to use for each stream
obj_func = CostScaledPaterson()
overall_network = merge(construct_superstructure(prob.stream_names, FloudasCiricGrossmann(), prob), construct_superstructure(prob.utility_names, ParallelSplit(), prob))
cost_coeff, scaling_coeff = 670, 0.83
optimizer = BARON.Optimizer

generate_network!(prob, EMAT; optimizer, overall_network, obj_func, verbose=true, cost_coeff, scaling_coeff)