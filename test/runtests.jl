using CompHENS
using Test

@testset "CompHENS.jl" begin
    @testset "Problem Construction" begin
        include("classic_hens_constructor.jl")
    end
end
