module CompHENS

using DocStringExtensions
using Kwonly


"""
$(TYPEDEF)

Umbrella for all problems where `CompHENS.jl` is relevant
"""
abstract type AbstractSynthesisProblem end

"""
$(TYPEDEF)

The classical Heat Exchanger Network Synthesis (HENS) problem.
"""
abstract type AbstractHENSProblem  <: AbstractSynthesisProblem end

"""
$(TYPEDEF)

Subproblem formulated while solving an `AbstractSynthesisProblem` 
"""
abstract type AbstractSubProblem end


"""
$(TYPEDEF)

Type for technique to solve one or more `AbstractSynthesisProblem` types.
"""
abstract type AbstractSynthesisAlgorithm end

"""
$(TYPEDEF)

Algorithm to solve an `AbstractSubProblem`
"""
abstract type AbstractSubProblemAlgorithm end

"""
$(TYPEDEF)

Holds the solution of an `AbstractSynthesisProblem`
"""
abstract type AbstractSolution end

"""
$(TYPEDEF)

Holds the solution of an `AbstractSubProblem`
"""
abstract type AbstractSubProblemSolution end

const smallest_value = 1e-4

# Holds structures for streams
export AbstractStream, HotStream, ColdStream, AbstractUtility, SimpleHotUtility, SimpleColdUtility 
include("Streams/streams.jl")

# Hold structures of problem types
export ClassicHENSProblem
include("ProblemConstructors/classic_hens_prob.jl")

end
