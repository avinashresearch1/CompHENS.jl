@info "Colberg Morari 1990"

test_file = joinpath(pwd(), "Examples", "XLSX_interface", "ClassicHENSProblem", "Colberg_Morari_1990", "CompHENS_interface_ColbergMorari.xlsx")
isfile(test_file) || error("Test file not found: $test_file")

prob = ClassicHENSProblem(test_file)
@test length(prob.hot_streams_dict) == 3 && length(prob.cold_streams_dict) == 4 && length(prob.hot_utilities_dict) == 1 && length(prob.cold_utilities_dict) == 1
@test prob.cold_utilities_dict["CW"].add_user_data["Cost [\$/kW]"] == 20 && prob.hot_utilities_dict["ST"].add_user_data["Cost [\$/kW]"] == 200
