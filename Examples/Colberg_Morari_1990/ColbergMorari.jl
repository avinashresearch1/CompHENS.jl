# Workflow using XLSX input:
# 1. Import necessary packages:
using CompHENS

# 2. Specify path to xlsx file
file_path_xlsx = joinpath(@__DIR__, "CompHENS_interface_ColbergMorari.xlsx")

# 3. Construct the appropriate kind of problem: Here it is a `ClassicHENSProblem`
prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = 20.0)
intervals = CompHENS.generate_heat_cascade_intervals(prob)