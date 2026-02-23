using CompHENS
using Test

@testset "CompHENS.jl" begin
    @testset "Problem Construction" begin
        include("classic_hens_constructor.jl")
    end

    @testset "XLSX Interface Examples" begin
        # Ensure plotting backends run headless in CI.
        ENV["GKSwstype"] = "100"

        examples_root = joinpath(@__DIR__, "..", "Examples", "XLSX_interface")
        example_files = [
            joinpath(examples_root, "ClassicHENSProblem", "6_Stream", "6_stream.jl"),
            joinpath(examples_root, "ClassicHENSProblem", "Colberg_Morari_1990", "ColbergMorari.jl"),
            joinpath(examples_root, "ClassicHENSProblem", "TEP4125_Examples", "Assignment_2.xlsx", "Assignment_2.jl"),
            joinpath(examples_root, "ClassicHENSProblem", "NaN_MWE", "NaN_MWE.jl"),
            joinpath(examples_root, "ClassicHENSProblem", "NaN_MWE", "Infeasible_issue.jl"),
            joinpath(examples_root, "ClassicHENSProblem", "Gundersen_4_stream", "Gundersen_4_stream.jl"),
            joinpath(examples_root, "MultiPeriodFlexibleHENSProblem", "Floudas_Grossmann_1987", "FloudasGrossmann1987.jl"),
        ]

        for file in example_files
            @testset "$(basename(file))" begin
                include(file)
            end
        end
    end
end
