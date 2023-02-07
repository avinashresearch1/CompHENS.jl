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

# ╔═╡ 210ec76a-a488-11ed-3f57-5552a8e1a4fa
# ╠═╡ show_logs = false
begin
	using Pkg
	curr_dir = pwd()
	Pkg.activate(curr_dir);
	Pkg.add("CompHENS")
	#Pkg.activate(abspath(joinpath(curr_dir, "../../..")));
	using CompHENS
	using JuMP, Plots, HiGHS, PlutoUI
end

# ╔═╡ 2d0fc1e8-def8-42a3-8bd7-22e3b93ac3a8
md""" ## CompHENS: Computational Tools for Heat Exchanger Network Synthesis"""

# ╔═╡ 6dfeb381-4276-4987-b054-b821b57fb4c0
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

# ╔═╡ a8a93e4d-007c-4626-84a5-68afc446d490
file_path_xlsx = abspath(joinpath(curr_dir, inits[1]));


# ╔═╡ 07e341d8-87cd-4f26-86d7-1303ea337323
inits[2]

# ╔═╡ 41b25b11-8dea-4cbb-9ff5-a8d22afe2b56
# ╠═╡ show_logs = false
begin
	prob = ClassicHENSProblem(file_path_xlsx; ΔT_min = inits[2], verbose = true);
	intervals = CompHENS.generate_transshipment_intervals(prob)
end

# ╔═╡ 5d424122-2303-43d9-9bfd-8486d78f5e87
get_primary_temperatures!(prob; verbose = true, balanced = false)
#plot_composite_curve(intervals.hot_side)

# ╔═╡ 2f4bbf85-4ed3-464d-9c5b-77b491158eed
plot_composite_curve(prob.results_dict[:primary_temperatures].hot_cc)

# ╔═╡ d61de861-22ce-4eed-a769-54ff7ffde106
plot_composite_curve(prob; cold_ref_enthalpy = 1000)

# ╔═╡ 8944db57-2600-492c-a5e5-eff66132c370
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

# ╔═╡ Cell order:
# ╟─2d0fc1e8-def8-42a3-8bd7-22e3b93ac3a8
# ╠═210ec76a-a488-11ed-3f57-5552a8e1a4fa
# ╠═6dfeb381-4276-4987-b054-b821b57fb4c0
# ╠═a8a93e4d-007c-4626-84a5-68afc446d490
# ╠═07e341d8-87cd-4f26-86d7-1303ea337323
# ╠═41b25b11-8dea-4cbb-9ff5-a8d22afe2b56
# ╠═5d424122-2303-43d9-9bfd-8486d78f5e87
# ╠═2f4bbf85-4ed3-464d-9c5b-77b491158eed
# ╠═d61de861-22ce-4eed-a769-54ff7ffde106
# ╟─8944db57-2600-492c-a5e5-eff66132c370
