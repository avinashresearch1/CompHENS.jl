module CompHENS

using DocStringExtensions

"""
$(TYPEDEF)

Umbrella for all problems where `CompHENS.jl` is relevant
"""
abstract type AbstractCompHENSProblem end

"""
$(TYPEDEF)

The classical Heat Exchanger Network Synthesis (HENS) problem.
"""
abstract type AbstractHENSProblem  <: AbstractCompHENSProblem end

"""
$(TYPEDEF)

Subproblem formulated while solving an `AbstractCompHENSProblem` 
"""
abstract type AbstractSubProblem end


"""
$(TYPEDEF)

Type for technique to solve one or more `AbstractCompHENSProblem` types.
"""
abstract type AbstractSynthesisAlgorithm end

"""
$(TYPEDEF)

Algorithm to solve an `AbstractSubProblem`
"""
abstract type AbstractSubProblemAlgorithm end

"""
$(TYPEDEF)

Holds the solution of an `AbstractCompHENSProblem`
"""
abstract type AbstractSolution end

"""
$(TYPEDEF)

Holds the solution of an `AbstractSubProblem`
"""
abstract type AbstractSubProblemSolution end



end
