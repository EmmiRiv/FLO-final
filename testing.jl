### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 8ebeb3b0-6510-11f0-1a05-a9bb7c54b706
using CSV, Printf

# ╔═╡ 6dd7b362-9df1-497e-ab1c-50fe1b0a15c8
fC = CSV.File(open("subsets/adult-236-36.csv"))

# ╔═╡ 81441eca-f2f0-4dcc-8d4f-6be700139a54
fF = CSV.File(open("subsets/adult-236-36-fac.csv"))

# ╔═╡ 489e49bd-c92d-4f4a-b6d2-bde7828d3e33
function l2(v1, v2)
	d = 0
	for i in 1:size(v1)[1]
		d += (v1[i]-v2[i])^2
	end
	return sqrt(d)
end

# ╔═╡ 5f2e41b1-7f61-4f23-9756-988f6c34b96b
function process()
	fC = CSV.File(open("subsets/adult-236-36.csv"))
	fF = CSV.File(open("subsets/adult-236-36-fac.csv"))
	Fpos = []
	dist = []
	distMatr = []
	gps = Dict()
	gC = []
	groupsF = Dict()
	dmin = Inf
	dmax = 0
	om = 0

	iind = 0
	for rwi in fF
		iind += 1
		push!(Fpos, [rwi.d1, rwi.d2, rwi.d3, rwi.d4, rwi.d5, rwi.d6])
	end

	jind = 0
	for rwj in fC
		jind += 1
		jdist = []
		jdistM = []
		jpos = [rwj.d1, rwj.d2, rwj.d3, rwj.d4, rwj.d5, rwj.d6]
		for ind in 1:size(Fpos)[1]
			ijdist = l2(jpos, Fpos[ind])
			push!(jdistM, ijdist)
			push!(jdist, [ind,ijdist])
			if ijdist < dmin
				dmin = ijdist
			end
			if ijdist > dmax
				dmax = ijdist
			end
		end
		sort!(jdist, by = x -> x[2])
		
		push!(dist, jdist)
		push!(distMatr, jdistM)

		if !haskey(gps, rwj.c2)
			om += 1
			gps[rwj.c2] = om
			groupsF[jind] = om
			append!(gC, 1)
		else
			groupsF[jind] = gps[rwj.c2]
			gC[gps[rwj.c2]] += 1
		end
	
	end
	return dist, distMatr, jind, size(Fpos)[1], gC, groupsF, om, dmin, dmax
		
end

# ╔═╡ d5b1c885-4dbd-4702-b33d-26c95dc9435b
function processSynth()
	fC = CSV.File(open("datasets-synth/synth-90-10.csv"))
	fF = CSV.File(open("datasets-synth/synth-90-10-fac.csv"))
	Fpos = []
	dist = []
	distMatr = []
	gps = Dict()
	gC = []
	groupsF = Dict()
	groupsO = Dict()
	dmin = Inf
	dmax = 0
	om = 0
	costs = Dict()

	iind = 0
	for rwi in fF
		iind += 1
		push!(Fpos, [rwi.d1, rwi.d2])
		costs[iind] = rwi.c
	end

	jind = 0
	for rwj in fC
		jind += 1
		jdist = []
		jdistM = []
		jpos = [rwj.d1, rwj.d2]
		for ind in 1:size(Fpos)[1]
			ijdist = l2(jpos, Fpos[ind])
			push!(jdistM, ijdist)
			push!(jdist, [ind,ijdist])
			if ijdist < dmin
				dmin = ijdist
			end
			if ijdist > dmax
				dmax = ijdist
			end
		end
		sort!(jdist, by = x -> x[2])
		
		push!(distMatr, jdistM)
		push!(dist,jdist)

		if !haskey(gps, rwj.c1)
			om += 1
			gps[rwj.c1] = om
			groupsF[jind] = om
			append!(gC, 1)
		else
			groupsF[jind] = gps[rwj.c1]
			gC[gps[rwj.c1]] += 1
		end
		groupsO[jind] = 1
	
	end
	return dist, distMatr, jind, iind, gC, groupsF, om, dmin, dmax
		
end

# ╔═╡ 5c4f8425-f831-4bca-8fcf-624b28d6cd69
dist, distMatr, n, m, gC, groupsF, om, dmin, dmax = process()

# ╔═╡ b7375671-d7ec-4f02-a5e0-136dac5bab1d
function cliDist(ctrs, j)
	for fac in dist[j]
		if in(fac[1],ctrs)
			return fac
		end
	end
end

# ╔═╡ 080e7d37-073d-4f5e-8738-d8e8530eea23
function cliSort(ctrs)
	cli = [j for j in 1:n]
	sort!(cli, by = x -> cliDist(ctrs, x)[2])
	return cli
end

# ╔═╡ 7085de6f-b838-4368-bea1-5a8d148545c2
function takeCensus(outliers)
	cen = zeros(om)
	for j in outliers
		cen[groupsF[j]] += 1
	end
	return cen
end

# ╔═╡ 47c6a83b-ab45-4e64-85b5-84ef462cf2f3
function takeRealCensus(assignment)
	cen = zeros(om)
	cost = 0
	for j in 1:n
		if assignment[j] == -1
			cen[groupsF[j]] += 1
		else
			cost += distMatr[j][assignment[j]]
		end
	end
	return cen, cost
end

# ╔═╡ 27bd7ab6-bf9e-4cdd-a7cf-ffb200d1292a
function computeCost(ctrs, served)
	cost = 0
	for j in served
		cost += cliDist(ctrs, j)[2]
	end
	return cost
end

# ╔═╡ 75494f47-3f51-4142-b6ae-748d67342cc3
function controlOut(out, ctrs)
	jOrder = cliSort(ctrs)
	outliers = jOrder[n-out+1:end]
	served = jOrder[1:n-out]
	cens = takeCensus(outliers)
	cost = computeCost(ctrs, served)
	return cens, cost
end

# ╔═╡ 98da20cf-2ceb-4ede-bb35-8ec6d9c33d6a
k = 5

# ╔═╡ 41ec7430-120a-4403-95a1-c7c7953c296b
function initLS(groups, pens)
	assignment = Dict()
	orig = [i for i in 1:m]
	init = []
	cost = 0
	
	while size(init)[1] < k
		fnd = rand(1:size(orig)[1])
		push!(init, orig[fnd])
		deleteat!(orig, fnd)
	end

	for j in 1:n
		lmin = pens[groups[j]]
		fmin = -1
		for i in init
			if distMatr[j][i] < lmin
				lmin = distMatr[j][i]
				fmin = i
			end
		end
		assignment[j] = fmin
		cost += lmin
	end

	return init, orig, assignment, cost
end

# ╔═╡ 1a6c3d6e-a8dc-48d6-af03-26c2ca76aaec
function localSearch(X0, FX0, assignment0, cost0, groups, pens)
	red = 0.99
	update = true
	for fnd in 1:k 
		fout = X0[fnd]
		xS = copy(X0)
		deleteat!(xS, fnd)
		for fndp in 1:m-k
			costN = 0
			assignmentN = Dict()
			fin = FX0[fndp]
			XN = copy(xS)
			push!(XN, fin)
			for j in 1:n
				if assignment0[j] == -1
					if distMatr[j][fin] < pens[groups[j]]
						assignmentN[j] = fin
						costN += distMatr[j][fin]
					else
						assignmentN[j] = -1
						costN += pens[groups[j]]
					end
				elseif assignment0[j] == fout
					lmin = pens[groups[j]]
					fmin = -1
					for i in XN
						if distMatr[j][i] < lmin
							lmin = distMatr[j][i]
							fmin = i
						end
					end
					assignmentN[j] = fmin
					costN += lmin
				else
					if distMatr[j][fin] < distMatr[j][assignment0[j]]
						assignmentN[j] = fin
						costN += distMatr[j][fin]
					else
						assignmentN[j] = assignment0[j]
						costN += distMatr[j][assignment0[j]]
					end
				end
			end
			if costN < red*cost0
				FXN = copy(FX0)
				deleteat!(FXN, fndp)
				push!(FXN, fout)
				update = false
				return update, XN, FXN, assignmentN, costN
			end
		end
	end
	
	return update, X0, FX0, assignment0, cost0
			
end

# ╔═╡ 9518c643-2d40-47b8-b38c-f03e61989947
function loopSearch(groups, guess, lG, gam)
	pens = [guess/(l*gam) for l in lG]
	X, FX, assignment, cost = initLS(groups, pens)
	fl = false
	while !fl
		fl, X, FX, assignment, cost = localSearch(X, FX, assignment, cost, groups, pens)
	end
	return X, FX, assignment, cost
end

# ╔═╡ 98fe5a47-f953-48da-b680-0da57e1940ab
function loopGuess(groups, lG, gam)
	eps = gam/om
	minCost = (n-sum(lG))*dmax
	minC = []
	minAssignment = Dict()
	for j in 1:n
		minAssignment[j] = -1
	end
	if dmin > 0
		guess = (n-sum(lG))*dmin
	else
		guess = 1
	end

	while guess < (n-sum(lG))*dmax
		C, FC, assignment, cost = loopSearch(groups, guess, lG, gam)
		if cost < minCost
			minCost = cost
			minC = C
			minAssignment = assignment
		end
		guess = guess*(1+eps)
	end
	return minC, minAssignment
end

# ╔═╡ 74717bb6-4af7-44d8-8251-bda4d089d926
function tstFair(lG, gam, out)
	ctrsF, assignmentF = loopGuess(groupsF, lG, gam)
	realCensus, realCost = takeRealCensus(assignmentF)
	censusF, costF = controlOut(out, ctrsF)
	return censusF, costF, realCensus, realCost
end

# ╔═╡ 5f32345e-5706-4027-8f75-7b643f8423b7
function tstOut(gam, out)
	groupsO = Dict()
	for j in 1:n
		groupsO[j] = 1
	end
	ctrsO, assignmentO = loopGuess(groupsO, [out], gam)
	realCensus, realCost = takeRealCensus(assignmentO)
	censusO, costO = controlOut(out, ctrsO)
	return censusO, costO, realCensus, realCost
end

# ╔═╡ 7ef27f26-1357-4c8a-bf39-af6d68f329db
function loopSearchBasic(groups, guess)
	pens = [Inf for l in 1:om]
	X, FX, assignment, cost = initLS(groups, pens)
	fl = false
	while !fl
		fl, X, FX, assignment, cost = localSearch(X, FX, assignment, cost, groups, pens)
	end
	return X, FX, assignment, cost
end

# ╔═╡ 9fa0610b-c971-44d6-b185-9cdd93c330fe
function loopGuessBasic(groups, gam)
	eps = om/gam
	minCost = n*dmax
	minC = []
	if dmin > 0
		guess = n*dmin
	else
		guess = 1
	end

	while guess < n*dmax
		C, FC, assignment, cost = loopSearchBasic(groups, guess)
		if cost < minCost
			minCost = cost
			minC = C
		end
		guess = guess*(1+eps)
	end
	return minC
end

# ╔═╡ b4af48ab-9f94-4c51-8816-c816fb73c4ef
function tstBasic(gam, out)
	ctrsB = loopGuessBasic(groupsF, gam)
	censusB, costB = controlOut(out, ctrsB)
	return ctrsB, censusB, costB
end

# ╔═╡ f182b5e4-caec-4c22-8bef-6d418cb693ed
function randCenters()
	orig = [i for i in 1:m]
	init = []
	
	while size(init)[1] < k
		fnd = rand(1:size(orig)[1])
		push!(init, orig[fnd])
		deleteat!(orig, fnd)
	end
	return init
end

# ╔═╡ 2ad3c215-d24d-48fb-9513-ee8115cc5977
function tstRand(out)
	ctrsR = randCenters()
	censusR, costR = controlOut(out, ctrsR)
	return ctrsR, censusR, costR
end

# ╔═╡ 953c5eb2-2a00-4c24-9061-ce9a6e79c509
function changeRat()
	costF = []
	costO = []
	costRF = []
	costRO = []
	dispR = []
	dispF = []
	dispO = []
	dispB = []
	dispRF = []
	dispRO = []
	costB = []
	for p in 1:3
		r = p/10
		lG = [floor(l*r) for l in gC]
		print(lG)
		out = trunc(Int, sum(lG))
		gam = r/10


		for gg in 1:5
			gam = (gg+p)/10
			cenF, cstF, cenRF, cstRF = tstFair(lG, gam, out)
			cenO, cstO, cenRO, cstRO = tstOut(gam/om, out)
			cB, cenB, cstB = tstBasic(gam, out)
	
			push!(costF, cstF)
			push!(costRF, cstRF)
			push!(costO, cstO)
			push!(costRO, cstRO)
			push!(costB, cstB)
	
			ratF = [cenF[g]/lG[g] for g in 1:om]
			ratO = [cenO[g]/lG[g] for g in 1:om]
			ratB = [cenB[g]/lG[g] for g in 1:om]
			ratRF = [cenRF[g]/lG[g] for g in 1:om]
			ratRO = [cenRO[g]/lG[g] for g in 1:om]
	
			push!(dispF, maximum(ratF))
			push!(dispO, maximum(ratO))
			push!(dispB, maximum(ratB))
			push!(dispRF, maximum(ratRF))
			push!(dispRO, maximum(ratRO))
		end
	end
	print(costRF, costF, costRO, costO, costB, dispRF, dispF, dispRO, dispO, dispB)
end

# ╔═╡ 2d97ac84-8a33-4e59-901b-59b642c895ca
changeRat()


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
CSV = "~0.10.15"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "1a81d4ab6a19fa608c3e1c1dad1cdcac849d9478"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

    [deps.Compat.weakdeps]
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

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

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

    [deps.FilePathsBase.weakdeps]
    Mmap = "a63ad114-7e13-5084-954f-fe012c677804"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.InlineStrings]]
git-tree-sha1 = "6a9fde685a7ac1eb3495f8e812c5a7c3711c2d5e"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.3"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "44f6c1f38f77cafef9450ff93946c53bd9ca16ff"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

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
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"
"""

# ╔═╡ Cell order:
# ╠═8ebeb3b0-6510-11f0-1a05-a9bb7c54b706
# ╠═6dd7b362-9df1-497e-ab1c-50fe1b0a15c8
# ╠═81441eca-f2f0-4dcc-8d4f-6be700139a54
# ╠═489e49bd-c92d-4f4a-b6d2-bde7828d3e33
# ╠═5f2e41b1-7f61-4f23-9756-988f6c34b96b
# ╠═d5b1c885-4dbd-4702-b33d-26c95dc9435b
# ╠═5c4f8425-f831-4bca-8fcf-624b28d6cd69
# ╠═41ec7430-120a-4403-95a1-c7c7953c296b
# ╠═1a6c3d6e-a8dc-48d6-af03-26c2ca76aaec
# ╠═9518c643-2d40-47b8-b38c-f03e61989947
# ╠═98fe5a47-f953-48da-b680-0da57e1940ab
# ╠═7ef27f26-1357-4c8a-bf39-af6d68f329db
# ╠═9fa0610b-c971-44d6-b185-9cdd93c330fe
# ╠═b7375671-d7ec-4f02-a5e0-136dac5bab1d
# ╠═080e7d37-073d-4f5e-8738-d8e8530eea23
# ╠═7085de6f-b838-4368-bea1-5a8d148545c2
# ╠═47c6a83b-ab45-4e64-85b5-84ef462cf2f3
# ╠═27bd7ab6-bf9e-4cdd-a7cf-ffb200d1292a
# ╠═f182b5e4-caec-4c22-8bef-6d418cb693ed
# ╠═75494f47-3f51-4142-b6ae-748d67342cc3
# ╠═2ad3c215-d24d-48fb-9513-ee8115cc5977
# ╠═74717bb6-4af7-44d8-8251-bda4d089d926
# ╠═b4af48ab-9f94-4c51-8816-c816fb73c4ef
# ╠═5f32345e-5706-4027-8f75-7b643f8423b7
# ╠═98da20cf-2ceb-4ede-bb35-8ec6d9c33d6a
# ╠═953c5eb2-2a00-4c24-9061-ce9a6e79c509
# ╠═2d97ac84-8a33-4e59-901b-59b642c895ca
# ╠═359ee077-1209-41ac-8d3e-ab5f4823e5bc
# ╠═991f3042-7e2e-4239-8da2-0513e8980b5a
# ╠═2d95c10f-5ce5-4c48-947c-091c6f0c1e52
# ╠═ea98cd79-f4e9-4051-b4ae-8fbf0fb9a4d0
# ╠═b9960958-2ffc-4537-8ac0-292208a3df0a
# ╠═42757bd3-1eb3-4c65-8402-6032b2e29663
# ╠═8aa3d11d-2e81-4057-8b77-5fb261b2b6f5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
