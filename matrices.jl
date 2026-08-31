### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 2fdacc9a-a541-11f1-bc99-9dead9035b56
vec1=Float64[2, 1.0,3]

# ╔═╡ 1deaa66e-764b-4cb8-bf15-5c1fa35e82d9
typeof(vec1)

# ╔═╡ 629e5fa2-17af-4d6a-843b-101f7d37d2e7
size(vec1)

# ╔═╡ 65be2151-50d7-45dc-862a-f37b4067f398
vec2 = Float64[2 1.0 3]

# ╔═╡ 07a471c5-2858-4b62-b74e-458579e796b7
mat1 = [1.0 2.3 3.1; 
		3.0 1.0 -7]

# ╔═╡ 7d0d19ed-c60c-47b2-b7de-ecc41b3e1884
mat1[1,2]

# ╔═╡ db37f73a-51cd-48ba-b7fc-8d8a4c9be118
mat1[end, end]

# ╔═╡ a6d6817d-e158-4eca-97a7-63f41cf131d3
mat1[begin,end-1]

# ╔═╡ 6709316a-4f67-4a35-a196-26bf52e5f7ad
mat1[: , [1,3]]

# ╔═╡ bc555db8-5cf2-4f2e-bd89-8f0cc1b6f066
size(vec3)

# ╔═╡ cfff7317-778e-4cbb-bf45-1ae3629cfbdb
rand(Int, 5,5)

# ╔═╡ f6f925a8-8cbf-4237-86f0-260d49a17379
zeros(2,3)

# ╔═╡ 96d4dedf-a27f-40f9-8dd0-b3dd2b73cbf0
mat1 .+ 10

# ╔═╡ 081663cc-b1c4-4a6b-bffb-6b74be2cc672
mat1

# ╔═╡ 0493962a-a8b5-4f50-a7db-f3348413b1aa
mat1[1,:]

# ╔═╡ 6b838126-fa50-4221-ba3f-e8539e032e5d
[mat1[1,:] -3*mat1[1,:]+mat1[2,:]]'

# ╔═╡ 21c5d075-7578-44cb-b29b-17dbf57dc05e
# ╠═╡ disabled = true
#=╠═╡
vec3 = vec2'
  ╠═╡ =#

# ╔═╡ 15050ff9-f783-4ada-a531-46dd21b8ac18
vec3 = [rand()*5-2.5 for i in 1:20]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╠═2fdacc9a-a541-11f1-bc99-9dead9035b56
# ╠═1deaa66e-764b-4cb8-bf15-5c1fa35e82d9
# ╠═629e5fa2-17af-4d6a-843b-101f7d37d2e7
# ╠═65be2151-50d7-45dc-862a-f37b4067f398
# ╠═21c5d075-7578-44cb-b29b-17dbf57dc05e
# ╠═07a471c5-2858-4b62-b74e-458579e796b7
# ╠═7d0d19ed-c60c-47b2-b7de-ecc41b3e1884
# ╠═db37f73a-51cd-48ba-b7fc-8d8a4c9be118
# ╠═a6d6817d-e158-4eca-97a7-63f41cf131d3
# ╠═6709316a-4f67-4a35-a196-26bf52e5f7ad
# ╠═15050ff9-f783-4ada-a531-46dd21b8ac18
# ╠═bc555db8-5cf2-4f2e-bd89-8f0cc1b6f066
# ╠═cfff7317-778e-4cbb-bf45-1ae3629cfbdb
# ╠═f6f925a8-8cbf-4237-86f0-260d49a17379
# ╠═96d4dedf-a27f-40f9-8dd0-b3dd2b73cbf0
# ╠═081663cc-b1c4-4a6b-bffb-6b74be2cc672
# ╠═0493962a-a8b5-4f50-a7db-f3348413b1aa
# ╠═6b838126-fa50-4221-ba3f-e8539e032e5d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
