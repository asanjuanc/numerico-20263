### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ d35d9fb1-a54d-43d9-9f11-ea18f06d9099
using Pkg

# ╔═╡ 5849a8a8-b2ad-11f1-8961-8fc246ed10d9
begin 
	using LinearAlgebra
	using PlutoUI
	using PrettyTables
end

# ╔═╡ c28f1ae2-ae22-4053-a351-7b4b3dbeebdf
md"""
# Factorización LU, norma de matrices y complejidad
"""

# ╔═╡ b2cab27e-c727-488b-a884-a55983a31f89
md" ## Factoización LU "

# ╔═╡ 59d2dc68-f91b-4a1a-8a32-fc73b1a70c1c
A = rand(-1:3, 3,3)

# ╔═╡ 2a7cb052-288e-40fa-b028-a48dff2f0f16
U = zeros(3,3)

# ╔═╡ 4ea73341-a655-4cb5-8a6f-429d2e09dd2a
L = diagm(ones(3))

# ╔═╡ 95208661-1c1b-4772-ac03-9a7153a0201a
L

# ╔═╡ 725f5c1e-1ecb-4823-850f-00fef5257967
U[1,:]=A[1,:]

# ╔═╡ 67534eb6-e116-4a8e-8942-2f5bf0b722bb
L[:,1] = A[:,1]/U[1,1]

# ╔═╡ 413e5563-940b-43f2-931e-cebb41746f43
A₂=A-L[:,1]*U[1,:]'

# ╔═╡ 8fcb735a-c549-4978-ac41-459a8f239b85
U[2,:] = A₂[2,:]

# ╔═╡ 5035b670-abd3-4a25-acf5-10b57b8591d1
L[:,2] = A₂[:,2]/U[2,2]

# ╔═╡ bb691682-58fe-475a-a7bb-dbbffc47efb0
L

# ╔═╡ ee57155f-79f8-4605-aecf-47eb99d90efa
A₃ = A₂ - L[:,2]*U[2,:]'

# ╔═╡ 5e856e08-72c0-4c25-9554-e990dec8c62f
U[3,:]=A₃[3,:]

# ╔═╡ 87e6ddd3-305d-43fa-90b7-cdebc78b5307
U

# ╔═╡ 8878bd1d-e49d-4510-b400-c3bd3a35e160
L*U-A

# ╔═╡ e45f900e-e6ba-4766-9dbf-4c1384b360be
md" ## Pivoteando"

# ╔═╡ 23ca78f0-fa45-4fd1-9d40-3f200b6c3e7f
@bind ϵ  Slider(0.00000001:0.00000001:0.1, show_value=true, default=10e-2)

# ╔═╡ 23500ee1-0b3e-4ad5-95d2-0fd443868cad
Aϵ = [-ϵ 1; 1 -1]

# ╔═╡ 10f6f546-b39d-4442-9348-23f4044677cc
LUϵ_sin_pivote = lu(Aϵ, NoPivot())

# ╔═╡ f1f41562-ab5c-4454-b804-1c03ea26ec22
cond(Aϵ)

# ╔═╡ cdc5176c-9adb-4572-921d-931a720bee78
cond(LUϵ_sin_pivote.U)

# ╔═╡ 2a4f90cd-f34f-4467-911f-5b9944d55ef2
LUϵ_con_pivote = lu(Aϵ)

# ╔═╡ 30de9ef9-b0cd-4761-af6d-b2e297209fe6
cond(LUϵ_con_pivote.U)

# ╔═╡ 8d13be97-411c-4e8c-ac81-6b4da1d0823c
cond(LUϵ_con_pivote.L)

# ╔═╡ 0cf52eaf-97d2-438f-97c8-771f5718470d
opnorm(LUϵ_con_pivote.L)*opnorm(inv(LUϵ_con_pivote.L)) 
# la norma de A por la norma de A inversa con la norma inducida con p=2 es igual al numero de condicionamiento de la matriz.

# ╔═╡ 4a8c29fc-bd4d-4b38-ac50-9c7e9afaf1f5
md" ## Numero de condicinamiento en la matriz problemática"

# ╔═╡ 57e5694f-f080-47aa-9468-c0bcace76bb4
begin
	α = 0.1
	@bind n Slider(1:16, show_value=true)
end

# ╔═╡ ed464406-038c-4c04-b47d-b323b8aabc50
β = 10^n

# ╔═╡ dd82c3d8-ee24-4d58-a0ff-8c8995739d82
begin B = [1 -1 0 α-β β;
	0 1 -1 0 0;
	 	0 0 1 -1 0;
	 	0 0 0 1 -1;
	 	0 0 0 0 1
	]
	b = [α, 0 ,0, 0, 1]
	B
end

# ╔═╡ 5b25a053-65eb-43f0-adcb-bfe73fc468da
x = B\b

# ╔═╡ 07bf7f95-54bf-45cd-9704-0f128e3db9af
x[1]

# ╔═╡ ca3342c1-0b13-4926-9667-5e71d27911c6
begin
	print([β, abs(x[1]-1), cond(B)  ])
end

# ╔═╡ ec330705-6838-4ce4-93f1-cb8e725148d2
md" esta es una matriz altamanente mal condicionada"

# ╔═╡ 9201a3b5-c75a-4b20-b38b-1bf9545f1b08
md" ## Complejidad y conteo de flops"

# ╔═╡ e0828aba-38e5-492b-abcb-aa34651ce4d3
begin
    m = 1000:1000:5000
    t = []
    for m in m
        C = randn(m, m)  
        y = randn(m)
        time = @elapsed for j in 1:80; C * y; end
        push!(t, time)
    end
end

# ╔═╡ 047566da-b0fd-4cf9-bc72-0ae8c730815d
pretty_table([m t];
     column_labels=["size", "time (sec.)"])

# ╔═╡ 133c3b47-d977-4a43-9802-90b2e8a4c00b
 t[m.==4000] ./ t[m.==2000];

# ╔═╡ 92c7ad08-4164-4641-8bd9-2777cbc9f010
t[m.==2000]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[compat]
PlutoUI = "~0.7.83"
PrettyTables = "~3.4.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "09440ddf0ab3eea5b5f397fa468846f65c3527aa"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.1+2"

[[deps.Crayons]]
git-tree-sha1 = "54b76cbb40d9a0f5368c880725b2f141da77c94f"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.2.0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05f45c2e0de6259db764adbfd2f1dc6d3f8de13c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "5005266de4bfe50e53ff44a5cb5c540b6e47a254"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1b8aa19f229b1cea7fc93874a52e49db6a854450"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.4.8"

    [deps.PrettyTables.extensions]
    PrettyTablesExcelExt = "XLSX"
    PrettyTablesTypstryExt = "Typstry"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"
    XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "773065c6e0e903924a9d838259be74338422aef2"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.5.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─c28f1ae2-ae22-4053-a351-7b4b3dbeebdf
# ╟─b2cab27e-c727-488b-a884-a55983a31f89
# ╠═d35d9fb1-a54d-43d9-9f11-ea18f06d9099
# ╟─5849a8a8-b2ad-11f1-8961-8fc246ed10d9
# ╠═59d2dc68-f91b-4a1a-8a32-fc73b1a70c1c
# ╠═2a7cb052-288e-40fa-b028-a48dff2f0f16
# ╠═4ea73341-a655-4cb5-8a6f-429d2e09dd2a
# ╠═95208661-1c1b-4772-ac03-9a7153a0201a
# ╠═725f5c1e-1ecb-4823-850f-00fef5257967
# ╠═67534eb6-e116-4a8e-8942-2f5bf0b722bb
# ╠═413e5563-940b-43f2-931e-cebb41746f43
# ╠═8fcb735a-c549-4978-ac41-459a8f239b85
# ╠═5035b670-abd3-4a25-acf5-10b57b8591d1
# ╠═bb691682-58fe-475a-a7bb-dbbffc47efb0
# ╠═ee57155f-79f8-4605-aecf-47eb99d90efa
# ╠═5e856e08-72c0-4c25-9554-e990dec8c62f
# ╠═87e6ddd3-305d-43fa-90b7-cdebc78b5307
# ╠═8878bd1d-e49d-4510-b400-c3bd3a35e160
# ╟─e45f900e-e6ba-4766-9dbf-4c1384b360be
# ╠═23ca78f0-fa45-4fd1-9d40-3f200b6c3e7f
# ╠═23500ee1-0b3e-4ad5-95d2-0fd443868cad
# ╠═10f6f546-b39d-4442-9348-23f4044677cc
# ╠═f1f41562-ab5c-4454-b804-1c03ea26ec22
# ╠═cdc5176c-9adb-4572-921d-931a720bee78
# ╠═2a4f90cd-f34f-4467-911f-5b9944d55ef2
# ╠═30de9ef9-b0cd-4761-af6d-b2e297209fe6
# ╠═8d13be97-411c-4e8c-ac81-6b4da1d0823c
# ╠═0cf52eaf-97d2-438f-97c8-771f5718470d
# ╟─4a8c29fc-bd4d-4b38-ac50-9c7e9afaf1f5
# ╟─57e5694f-f080-47aa-9468-c0bcace76bb4
# ╟─ed464406-038c-4c04-b47d-b323b8aabc50
# ╟─dd82c3d8-ee24-4d58-a0ff-8c8995739d82
# ╟─5b25a053-65eb-43f0-adcb-bfe73fc468da
# ╟─07bf7f95-54bf-45cd-9704-0f128e3db9af
# ╟─ca3342c1-0b13-4926-9667-5e71d27911c6
# ╠═ec330705-6838-4ce4-93f1-cb8e725148d2
# ╟─9201a3b5-c75a-4b20-b38b-1bf9545f1b08
# ╠═e0828aba-38e5-492b-abcb-aa34651ce4d3
# ╠═047566da-b0fd-4cf9-bc72-0ae8c730815d
# ╠═133c3b47-d977-4a43-9802-90b2e8a4c00b
# ╠═92c7ad08-4164-4641-8bd9-2777cbc9f010
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
