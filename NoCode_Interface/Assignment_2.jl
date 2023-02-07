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
	using CompHENS, Plots, HiGHS, PlutoUI, Interpolations
end

# ╔═╡ 2d3f9546-4a1a-48b7-b53f-574ef42a633f
md""" # CompHENS: Computational Tools for Heat Exchanger Network Synthesis"""

# ╔═╡ bbec22e9-2362-4e42-b82a-fdf19a73cbdb
@bind inits PlutoUI.combine() do Child
	md"""
	**Insert stream data path below:**
	
	$(
		Child(TextField((70,1); default = "InputData.xlsx"))
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

# ╔═╡ 19b2497c-c4a0-43cc-84bb-11f50de4fc6d
begin 
	prob2 = deepcopy(prob)
	solve_minimum_utilities_subproblem!(prob2; verbose = false)
end

# ╔═╡ 223c191e-5d0d-45c6-a46c-e9398c3fd2b8
with_terminal() do
	print_min_utils_pinch_points(prob2)
end


# ╔═╡ 02a29d27-1cc9-4e9d-8a60-b5a96d428495
prob2.results_dict[:min_units]

# ╔═╡ 833ae097-cdd4-412f-a695-69ce0e9e164a
res2 = plot_composite_curve(prob2; balanced = true);

# ╔═╡ 12072330-5790-4421-89c0-dbd8e1990c19
res2.plt

# ╔═╡ e6924ad2-ca1b-4ea7-a318-76b0ba402aab
solve_minimum_units_subproblem!(prob2)

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
# ╟─833ae097-cdd4-412f-a695-69ce0e9e164a
# ╟─e6924ad2-ca1b-4ea7-a318-76b0ba402aab
# ╟─19b2497c-c4a0-43cc-84bb-11f50de4fc6d
# ╟─1248fe03-c1a6-4536-a955-2e5e0ebd861b
# ╟─aeee17c8-cc3d-4540-a87d-82e19d1536a5
# ╟─a7355458-a4ec-11ed-3d18-c3794e1eca79
# ╟─72964326-b9e0-4c88-afc8-cb15d05879b3
# ╟─d313bd2d-ddba-4244-802f-c833b16a9944
