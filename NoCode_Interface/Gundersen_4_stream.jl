### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a7355458-a4ec-11ed-3d18-c3794e1eca79
# ╠═╡ show_logs = false
begin
	using Pkg
	curr_dir = pwd()
	Pkg.activate(curr_dir);
	Pkg.add("Plots"); Pkg.add("PlutoUI"); Pkg.add("HiGHS"); Pkg.add("Interpolations")
	#Pkg.activate(abspath(joinpath(curr_dir, "../../..")));
	using CompHENS, Plots, HiGHS, PlutoUI, Interpolations, JuMP, BARON
end

# ╔═╡ 2d3f9546-4a1a-48b7-b53f-574ef42a633f
md""" # CompHENS: Computational Tools for Heat Exchanger Network Synthesis"""

# ╔═╡ bbec22e9-2362-4e42-b82a-fdf19a73cbdb
@bind inits PlutoUI.combine() do Child
	md"""
	**Insert stream data path below:**
	
	$(
		Child(TextField((70,1); default = "Gundersen_4_stream.xlsx"))
	)
	

	\
	
	**ΔT_min :**  $(
		Child(NumberField(0.0:200; default = 10.0))
	)
	"""
end

# ╔═╡ 11737402-8a05-46ac-805c-b24ca75a7af5
md""" ### Graphical Approach:"""

# ╔═╡ f0867f5d-d516-4414-b028-bf3daa2f1613
	@bind cold_ref_enthalpy Slider(0.0:10000)#; show_value = true)

# ╔═╡ 0bc8889e-7d4b-41a8-bc13-8addfc0a105b
md""" ### Heat Cascade"""

# ╔═╡ d959aa8a-a4b3-4abb-ae4c-66768e5a67c7
md"""### Mathematical Programming Approach"""

# ╔═╡ 5615c739-0077-48ea-9993-71878cc4a3a4
md"""**1. Minimum utility consumption**"""

# ╔═╡ 3ff5e8f4-914c-47ca-90c8-86a916db4245
md"""**Balanced Composite Curve**"""

# ╔═╡ adeb338f-2fbd-4d40-97d5-0137bb219196
md"""**2. Minimum number of units:**"""

# ╔═╡ 0e6c28ca-5371-4d9b-b081-ca7e1c3a1b39
md"""**EMAT:**"""

# ╔═╡ cdb4978e-ea10-48a5-897d-a2b409d898ea
@bind EMAT NumberField(0.0:0.01:20.0; default = inits[2]/8)

# ╔═╡ 4b9fb493-8426-4381-8f8c-64f06b4fbf3c
md"""**Additional Units:**"""

# ╔═╡ a10e4bb1-9daa-4df3-a6ed-8b4a43651e19
@bind add_units NumberField(0:10; default = 1)

# ╔═╡ 173b873a-9924-48bc-b1e2-28323724246d
md"""**Heat Load Distributions:**"""

# ╔═╡ d24d5dcb-6c32-4436-9cc8-f71d53b2448b
md"""**HX Base cost:**"""

# ╔═╡ c34244fb-bc88-4355-a651-5dc0a5b4ed8e
@bind base_cost NumberField(0.0:10000; default = 4000.0)


# ╔═╡ 3ff77685-a92b-4401-9f1d-ef070cb0c1e4
md"""**Cost Coefficient:**"""

# ╔═╡ f2e50f32-bb24-4391-a6ce-2285a8136c07
@bind cost_coeff NumberField(0.0:10000; default = 500.0)

# ╔═╡ e2e089c0-8d8f-4ee0-b162-db3b4632a99b
md"""**Scaling Factor:**"""

# ╔═╡ 322c8c7a-0f82-4eb2-ae8d-e65cc0bde709
@bind scaling_coeff NumberField(0.0:0.01:1.0; default = 0.83)

# ╔═╡ b9eb8c3c-0051-45ee-898b-dd9ff27effe2
md"""**Solver time:**"""

# ╔═╡ e654296e-927f-40d5-8934-d898dc2359c1
@bind solver_time NumberField(0.0:500; default = 20)

# ╔═╡ 14a8fe11-0df8-4b59-9ebc-14d7aea12cfc
md"""**Results:**"""

# ╔═╡ 3589836e-c720-426f-8580-dcb617fda150
md"""**Output file path:**"""

# ╔═╡ ef27eeed-4b5f-4110-8cc7-890c2ef6a99a
@bind output_path TextField((70,1); default = "Gundersen_4_stream.pdf")

# ╔═╡ 72964326-b9e0-4c88-afc8-cb15d05879b3
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(5px, 5%);
    	padding-right: max(5px, 5%);
	}
</style>
"""

# ╔═╡ d313bd2d-ddba-4244-802f-c833b16a9944
# ╠═╡ show_logs = false
begin 
	file_path_xlsx = abspath(joinpath(curr_dir, inits[1]));
	prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = inits[2], verbose = true);
end;

# ╔═╡ 12ba84b1-36e6-4c54-b1fc-6e7b3a0054e6
with_terminal() do
	print_full(generate_transshipment_intervals(prob))
end

# ╔═╡ fc3a40a2-6c4c-4256-977c-9b61c7bfab1b
begin
	prob2 = deepcopy(prob)
	a1 = solve_minimum_utilities_subproblem!(prob2; verbose = false)
	
end

# ╔═╡ 223c191e-5d0d-45c6-a46c-e9398c3fd2b8
with_terminal() do
	print_min_utils_pinch_points(prob2)
end


# ╔═╡ 19b2497c-c4a0-43cc-84bb-11f50de4fc6d
# ╠═╡ show_logs = false
begin
	a1;
	a2 = solve_minimum_units_subproblem!(prob2)
	push!(prob2.results_dict, :add_units => convert(Int64, add_units))
	a3 = nothing;
	
end

# ╔═╡ 02a29d27-1cc9-4e9d-8a60-b5a96d428495
begin
	a2;
	prob2.results_dict[:min_units]
end

# ╔═╡ 2798fb3f-5c53-4d5a-8028-283001d170b4
# ╠═╡ show_logs = false
begin
	a3;
	a4 = generate_stream_matches!(prob2, convert(Float64, EMAT); digits = 8, verbose  = false)
end

# ╔═╡ 6c50557f-f964-4f77-af51-5086c4a517d4
begin
	a4;
	prob2.results_dict[:Q]
end

# ╔═╡ 4d970d06-318a-4959-a40c-8b9d457f2973
begin
	a4;
	output_file = abspath(joinpath(curr_dir, output_path));
	results_df = generate_network!(prob2, EMAT; optimizer = optimizer_with_attributes(BARON.Optimizer, "MaxTime" => solver_time, "AbsConFeasTol" => 1), cost_coeff = cost_coeff, scaling_coeff = scaling_coeff, base_cost = base_cost, save_model = true, output_file = output_file);
end; 

# ╔═╡ 1d3f6f45-b00c-49f2-b1bc-9e63bf7f8387
results_df

# ╔═╡ 833ae097-cdd4-412f-a695-69ce0e9e164a
res2 = plot_composite_curve(prob2; balanced = true);

# ╔═╡ 12072330-5790-4421-89c0-dbd8e1990c19
res2.plt


# ╔═╡ aeee17c8-cc3d-4540-a87d-82e19d1536a5
res = plot_composite_curve(prob; cold_ref_enthalpy = cold_ref_enthalpy);

# ╔═╡ 8af04aea-5378-4a42-9f3d-e94335458ac7
res.plt

# ╔═╡ 1248fe03-c1a6-4536-a955-2e5e0ebd861b
# ╠═╡ show_logs = false
begin
	hot_H_mod = deepcopy(res.H_vals_hot)
	hot_T_mod = deepcopy(res.T_vals_hot)
	cold_H_mod = deepcopy(res.H_vals_cold)
	cold_T_mod = deepcopy(res.T_vals_cold)
	min_overlap = max(first(res.H_vals_cold), first(res.H_vals_hot))
	max_overlap = min(last(res.H_vals_cold), last(res.H_vals_hot))
	#push!(hot_H_mod, maximum(res.H_vals_cold)) # To allow interpolation of all cold H
	#push!(hot_T_mod, maximum(res.T_vals_hot)+1)
	hot_interp = linear_interpolation(hot_H_mod, hot_T_mod)
	cold_interp = linear_interpolation(res.H_vals_cold, res.T_vals_cold)
	Dt_vec = [hot_interp(h) - cold_interp(h) for h in union(res.H_vals_hot, res.H_vals_cold) if (h > min_overlap && h < max_overlap)]# Only kinks matter
	utils = (Q_cold = cold_ref_enthalpy, Q_hot=max(0.0, maximum(res.H_vals_cold)-maximum(res.H_vals_hot)), ΔT = minimum(Dt_vec))
end;

# ╔═╡ 1004e6f3-2e94-4561-8ce5-6b58f6af7110
utils

# ╔═╡ Cell order:
# ╟─2d3f9546-4a1a-48b7-b53f-574ef42a633f
# ╟─bbec22e9-2362-4e42-b82a-fdf19a73cbdb
# ╟─11737402-8a05-46ac-805c-b24ca75a7af5
# ╟─8af04aea-5378-4a42-9f3d-e94335458ac7
# ╟─f0867f5d-d516-4414-b028-bf3daa2f1613
# ╟─1004e6f3-2e94-4561-8ce5-6b58f6af7110
# ╟─0bc8889e-7d4b-41a8-bc13-8addfc0a105b
# ╟─12ba84b1-36e6-4c54-b1fc-6e7b3a0054e6
# ╟─d959aa8a-a4b3-4abb-ae4c-66768e5a67c7
# ╟─5615c739-0077-48ea-9993-71878cc4a3a4
# ╟─3ff5e8f4-914c-47ca-90c8-86a916db4245
# ╟─12072330-5790-4421-89c0-dbd8e1990c19
# ╟─223c191e-5d0d-45c6-a46c-e9398c3fd2b8
# ╟─adeb338f-2fbd-4d40-97d5-0137bb219196
# ╟─02a29d27-1cc9-4e9d-8a60-b5a96d428495
# ╟─0e6c28ca-5371-4d9b-b081-ca7e1c3a1b39
# ╟─cdb4978e-ea10-48a5-897d-a2b409d898ea
# ╟─4b9fb493-8426-4381-8f8c-64f06b4fbf3c
# ╟─a10e4bb1-9daa-4df3-a6ed-8b4a43651e19
# ╟─173b873a-9924-48bc-b1e2-28323724246d
# ╟─6c50557f-f964-4f77-af51-5086c4a517d4
# ╟─d24d5dcb-6c32-4436-9cc8-f71d53b2448b
# ╟─c34244fb-bc88-4355-a651-5dc0a5b4ed8e
# ╟─3ff77685-a92b-4401-9f1d-ef070cb0c1e4
# ╟─f2e50f32-bb24-4391-a6ce-2285a8136c07
# ╟─e2e089c0-8d8f-4ee0-b162-db3b4632a99b
# ╟─322c8c7a-0f82-4eb2-ae8d-e65cc0bde709
# ╟─b9eb8c3c-0051-45ee-898b-dd9ff27effe2
# ╟─e654296e-927f-40d5-8934-d898dc2359c1
# ╟─4d970d06-318a-4959-a40c-8b9d457f2973
# ╟─14a8fe11-0df8-4b59-9ebc-14d7aea12cfc
# ╟─1d3f6f45-b00c-49f2-b1bc-9e63bf7f8387
# ╟─3589836e-c720-426f-8580-dcb617fda150
# ╟─ef27eeed-4b5f-4110-8cc7-890c2ef6a99a
# ╟─fc3a40a2-6c4c-4256-977c-9b61c7bfab1b
# ╟─19b2497c-c4a0-43cc-84bb-11f50de4fc6d
# ╟─2798fb3f-5c53-4d5a-8028-283001d170b4
# ╟─1248fe03-c1a6-4536-a955-2e5e0ebd861b
# ╟─aeee17c8-cc3d-4540-a87d-82e19d1536a5
# ╟─833ae097-cdd4-412f-a695-69ce0e9e164a
# ╟─a7355458-a4ec-11ed-3d18-c3794e1eca79
# ╟─72964326-b9e0-4c88-afc8-cb15d05879b3
# ╟─d313bd2d-ddba-4244-802f-c833b16a9944
