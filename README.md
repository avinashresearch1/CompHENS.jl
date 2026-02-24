# CompHENS

[![Build Status](https://github.com/avinashresearch1/CompHENS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/avinashresearch1/CompHENS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![DOI](https://zenodo.org/badge/511286580.svg)](https://zenodo.org/badge/latestdoi/511286580)

This software provides a Julia-based toolkit for synthesis of Heat Exchanger Networks (HENs) using a mathematical programming framework. Currently, a sequential algorithm is implemented whereby an LP is formulated to determine the minimum utility consumption, MILP to determine the minimum number of units and the stream matches. Finally, an NLP is formulated to generate the network and calculate the HEN area. 

If you use this toolkit, please cite:
Avinash Subramanian, Flemming Holtorf, Rahul Anantharaman, Truls Gundersen. (2023). CompHENS: Computational Tools for Heat Exchanger Network Synthesis (Version v0.1.0) [Computer software]
DOI: 10.5281/zenodo.7545869

## XLSX interface usage (recommended): 
1. Clone this repo. Open the Julia repl, then `] activate .` and `]instantiate`. This will install all required dependencies.
2. Test that this works on one of the examples such as by running the code in `Examples/XLSX_interface/ClassicHENSProblem/Colberg_Morari_1990/ColbergMorari.jl`
3. For your own problem:
   - Fill in the stream data in an xlsx file in a similar way as `CompHENS_interface_ColbergMorari.xlsx`
   - The structure of the code that runs the problem is as follows:
  
```
# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS
using Plots
using JuMP
using HiGHS
using Ipopt
using MathOptInterface
using Test
const MOI = MathOptInterface

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_ColbergMorari.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`. Specify the ΔT_min.
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0, verbose = true)

# 4. Solve minimum utilities problem

solve_minimum_utilities_subproblem!(prob)
print_min_utils_pinch_points(prob)

# 5. Solve the minimum number of units subproblem:
solve_minimum_units_subproblem!(prob)
@show prob.min_units

# 6. Generate stream matches. This requires specifying the EMAT and the number of additional units. 
EMAT = 2.5
prob.results_dict[:add_units] = 1
@time generate_stream_matches!(prob, EMAT; digits = 8)

# View the results. 
prob.results_dict[:Q]

# 7. Network generation:
# Specify which superstructure to use for each stream
obj_func = CostScaledPaterson()
overall_network = merge(construct_superstructure(prob.stream_names, FloudasCiricGrossmann(), prob), construct_superstructure(prob.utility_names, FloudasCiricGrossmann(), prob))

# Specify base costs, scaling coefficinets. 
base_cost, cost_coeff, scaling_coeff = 8600, 670, 0.83
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "max_iter" => 5000, "tol" => 1e-6)

generate_network!(prob, EMAT; overall_network = overall_network, obj_func = obj_func, optimizer = optimizer, verbose = true, cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true)
model = prob.results_dict[:network_gen_model]

# If this fails, the optimizer has failed to find a local minimum
@test termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL]
@show termination_status(model)
@test primal_status(model) == MOI.FEASIBLE_POINT
# print(model)
value.(model[:ΔT_upper])
value.(model[:ΔT_lower])
value.(model[:T_LMTD])
get_design_area(prob)
```

## No code usage:
1. Download the 2 interface files from: https://github.com/avinashresearch1/CompHENS.jl/tree/main/NoCode_Interface and put them in a folder of your choice (same files from email). 
2. Type in the stream data in the `InputData.xlsx` file. Note that all streams must have a sensible temperature difference (say use a 1 C temperature difference for condensing steam).
3. Download and Install Julia from [Julia website](https://julialang.org/downloads/). The No-Code interface only requires the Julia REPL.
4. From the Julia REPL, access the Package Manager by typing: `]`. Install `CompHENS` and `Pluto`
![image](https://user-images.githubusercontent.com/90404321/217259675-2c48f58c-bd7a-4a86-9d76-1da82989c559.png)
4. Once everything is installed, exit the package manager by typing backspace. Type `using Pluto; Pluto.run()`. This will launch the browser.
5. Navigate to your directory to the `Interface.jl file`. Pluto will launch. 
6. The slider can be used to move the composite curves. The curves update automatically with changing the `DT_min`.

**Note:** It may be necessary to re-run the Pluto notebook in case of errors: Click out of any cell, and `Ctrl+A`followed by `Shift+Enter`.

## Slide Deck with basic overview:



![image](https://user-images.githubusercontent.com/90404321/223507752-45147c98-860f-4f94-b874-5fade560d9b3.png)
![image](https://user-images.githubusercontent.com/90404321/223507848-016b7c1b-30a3-4af3-bb89-267f091c505d.png)
![image](https://user-images.githubusercontent.com/90404321/223507941-fa4c6579-b921-4ae7-88bc-0c857a012db6.png)
![image](https://user-images.githubusercontent.com/90404321/223507975-b3a07e28-918f-4893-9049-1de8a112a07c.png)
![image](https://user-images.githubusercontent.com/90404321/223508019-8e85c1f4-aed0-46c9-938c-62d54be467d0.png)
![image](https://user-images.githubusercontent.com/90404321/223508183-7346b624-3e2f-4630-8be1-af921c9725f6.png)
![image](https://user-images.githubusercontent.com/90404321/223508221-e8f7cd81-4a5c-4a2f-9204-6f68b10f820a.png)
![image](https://user-images.githubusercontent.com/90404321/223508267-e60e9ef9-1f78-4afb-a14d-a06f31e4c3a7.png)
![image](https://user-images.githubusercontent.com/90404321/223508335-6e79647b-4a57-42a9-a4b6-952948ea2678.png)
![image](https://user-images.githubusercontent.com/90404321/223508377-b4225c39-ec43-4f3f-921d-da74b6dd65b2.png)
![image](https://user-images.githubusercontent.com/90404321/223508415-fe8c050f-37db-449e-9587-61a6a9c74cdc.png)
![image](https://user-images.githubusercontent.com/90404321/223508455-4899d253-423e-4ca2-bb05-75232f256d1a.png)













