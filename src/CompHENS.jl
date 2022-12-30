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

export AbstractSynthesisProblem, AbstractSynthesisAlgorithm, AbstractSolution  

# Define own suitable for HENS.
const smallest_value = 1e-8
export smallest_value

# Holds structures for streams
export AbstractStream, HotStream, ColdStream, AbstractUtility, SimpleHotUtility, SimpleColdUtility, U, M
include("Streams/streams.jl")

# Hold structures of problem types
export ClassicHENSProblem, MultiPeriodFlexibleHENSProblem
include("ProblemConstructors/classic_hens_prob.jl")
include("ProblemConstructors/multiperiod_flexible_hens_prob.jl")

# Holds structures for processing the composite curve e.g., kink points
export Point
include("Intervals/curve_points.jl")

# Holds all kinds of temperature intervals
export TemperatureInterval, TransshipmentInterval, 
generate_transshipment_intervals, plot_hot_composite_curve, plot_cold_composite_curve, plot_composite_curve, 
get_contribution, print_full, initialize_temperature_intervals, assign_stream!, assign_utility!, assign_all_streams_and_utilities!,
get_primary_temperatures!, calculate_enthalpies!, get_secondary_temperatures!, get_tertiary_temperatures!, get_quaternary_temperatures!,
LMTD, is_feasible

include("Intervals/temperature_intervals.jl")

export solve_minimum_utilities_subproblem!, print_min_utils_pinch_points
include("SubProblems/minimum_utilities_subprob.jl")

export solve_minimum_units_subproblem!
include("SubProblems/minimum_number_of_units.jl")

export generate_stream_matches!, print_HLD
include("SubProblems/generate_stream_matches.jl")

end
