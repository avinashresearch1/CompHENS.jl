EMAT = 2.5
optimizer = BARON.Optimizer
verbose = true
save_model =  true
time_limit = 100.0
output_folder = "/home/avinash/Desktop/COMPHENS/CompHENS.jl/Examples/MultiPeriodFlexibleHENSProblem/"
overall_network::Dict{String, AbstractSuperstructure} = construct_superstructure(prob.all_names, FloudasCiricGrossmann(), prob)

