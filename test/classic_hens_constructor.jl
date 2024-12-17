@testset "Colberg Morari 1990" begin
    test_file = joinpath(@__DIR__, "..", "Examples", "XLSX_interface", "ClassicHENSProblem", "Colberg_Morari_1990", "CompHENS_interface_ColbergMorari.xlsx")
    isfile(test_file) || error("Test file not found: $test_file")

    prob = ClassicHENSProblem(test_file)
    @test length(prob.hot_streams_dict) == 3 && length(prob.cold_streams_dict) == 4 && length(prob.hot_utilities_dict) == 1 && length(prob.cold_utilities_dict) == 1
    @test prob.cold_utilities_dict["CW"].add_user_data["Cost [\$/kW]"] == 20 && prob.hot_utilities_dict["ST"].add_user_data["Cost [\$/kW]"] == 200

    solve_minimum_utilities_subproblem!(prob)
    @test prob.pinch_points[] == (507.0, 497.0)
    @test prob.hot_utilities_dict["ST"].Q ≈ 182.521
    @test prob.cold_utilities_dict["CW"].Q ≈ 110.986

    solve_minimum_units_subproblem!(prob)
    @test prob.min_units == 9

    EMAT = 2.5
    generate_stream_matches!(prob, EMAT; add_units=1)
    @test prob.results_dict[:y] == [1 1 1 1; 0 1 1 0; 0 0 1 0; 0 1 1 0; 0 0 1 0]
    @test all(prob.results_dict[:Q] .≈ [392.08 164.272 93.891 182.521; 0.0 63.314 56.553 0.0; 0.0 0.0 457.62 0.0; 0.0 68.445 359.125 0.0; 0.0 0.0 110.986 0.0])
    @test count(==(1), prob.results_dict[:y]) == 10 && count(>(0), prob.results_dict[:Q]) == 10
end;

@testset "6 stream" begin
    test_file = joinpath(@__DIR__, "..", "Examples", "XLSX_interface", "ClassicHENSProblem", "6_Stream", "CompHENS_interface_SimpleExample.xlsx")
    isfile(test_file) || error("Test file not found: $test_file")

    prob = ClassicHENSProblem(test_file)
    @test length(prob.hot_streams_dict) == 2 && length(prob.cold_streams_dict) == 2 && length(prob.hot_utilities_dict) == 1 && length(prob.cold_utilities_dict) == 1
    @test prob.cold_utilities_dict["CW"].add_user_data["Cost [\$/kW]"] == 20 && prob.hot_utilities_dict["ST"].add_user_data["Cost [\$/kW]"] == 200

    solve_minimum_utilities_subproblem!(prob)
    @test prob.pinch_points[] == (125.0, 115.0)
    @test prob.hot_utilities_dict["ST"].Q == 300.0
    @test prob.cold_utilities_dict["CW"].Q == 220.0

    solve_minimum_units_subproblem!(prob)
    @test prob.min_units == 5

    EMAT = 2.5
    generate_stream_matches!(prob, EMAT; add_units=1)
    @test prob.results_dict[:y] == [1 1 1; 0 1 0; 1 1 0]
    @test all(prob.results_dict[:Q] .≈ [1176.0 1224.0 300.0; -7.815970093361102e-14 1080.0 0.0; 124.0 96.0 0.0])
    @test count(==(1), prob.results_dict[:y]) == 6 && count(>(0), prob.results_dict[:Q]) == 6
end;

@testset "Gundersen 4-stream" begin
    test_file = joinpath(@__DIR__, "..", "Examples", "XLSX_interface", "ClassicHENSProblem", "Gundersen_4_stream", "Gundersen_4_stream.xlsx")
    isfile(test_file) || error("Test file not found: $test_file")

    prob = ClassicHENSProblem(test_file)
    @test length(prob.hot_streams_dict) == 2 && length(prob.cold_streams_dict) == 2 && length(prob.hot_utilities_dict) == 1 && length(prob.cold_utilities_dict) == 1
    @test prob.cold_utilities_dict["CW"].add_user_data["Cost [\$/kW]"] == 200 && prob.hot_utilities_dict["ST"].add_user_data["Cost [\$/kW]"] == 20

    solve_minimum_utilities_subproblem!(prob)
    @test prob.pinch_points[] == (170.0, 160.0)
    @test prob.hot_utilities_dict["ST"].Q == 600.0
    @test prob.cold_utilities_dict["CW"].Q == 400.0

    solve_minimum_units_subproblem!(prob)
    @test prob.min_units == 5

    EMAT = 2.5
    generate_stream_matches!(prob, EMAT; add_units=1)
    @test prob.results_dict[:y] == [1 1 0; 1 1 1; 0 1 0]
    @test all(prob.results_dict[:Q] .≈ [1015.9428571428573 2184.057142857143 0.0; 964.0571428571427 935.9428571428573 600.0; 0.0 400.0 0.0])
    @test count(==(1), prob.results_dict[:y]) == 6 && count(>(0), prob.results_dict[:Q]) == 6
end;