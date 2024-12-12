# Default solver configurations:
# MILP
const HIGHS_solver = JuMP.optimizer_with_attributes(
    HiGHS.Optimizer,
    "presolve" => "on",
    "log_to_console" => false,
    "time_limit" => 20.0
)

#=
# NLP
const IPOPT_solver = JuMP.optimizer_with_attributes(
        Ipopt.Optimizer,
        #MOI.Silent() => true,
        "sb" => "yes",
        "max_iter" => Int(1E4),
        "max_wall_time" => 20.0
    )

const JUNIPER_solver = JuMP.optimizer_with_attributes(
    Juniper.Optimizer,
    "nl_solver"=>optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
) 

#    optimizer_with_attributes(optimizer, "nl_solver" => optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))

#const COUENNE_solver = JuMP.optimizer_with_attributes(
#    AmplNLWriter.Optimizer(Couenne_jll.amplexe)
#)



const SCIP_solver = JuMP.optimizer_with_attributes(
    SCIP.Optimizer,
    "display/verblevel" => 1,
    "limits/time" => 100.0
)
=#

#=
const ALPINE_solver = JuMP.optimizer_with_attributes(
    Alpine.Optimizer,
    "nlp_solver" => JuMP.optimizer_with_attributes(
        Ipopt.Optimizer,
        MOI.Silent() => true,
        "sb" => "yes",
        #"max_iter" => Int(1E4),
    ),
    "mip_solver" => JuMP.optimizer_with_attributes(
        HiGHS.Optimizer,
        "presolve" => "on",
        "log_to_console" => false),
    "presolve_bt" => true,
    "apply_partitioning" => true,
    "partition_scaling_factor" => 10,
    "max_iter" => 1,
    #"time_limit" => 100.0
)
=#