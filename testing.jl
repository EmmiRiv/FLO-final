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
	minC = cost0
	minX = X0
	minFX = FX0
	minAssignment = assignment0
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
			if costN < red*cost0 && costN < minC
				FXN = copy(FX0)
				deleteat!(FXN, fndp)
				push!(FXN, fout)
				update = false
				minC = costN
				minX = XN
				minFX = FXN
				minAssignment = assignmentN
			end
		end
	end
	
	return update, minX, minFX, minAssignment, minC
			
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

# ╔═╡ 359ee077-1209-41ac-8d3e-ab5f4823e5bc
"""
[17.0, 6.0][34.0, 13.0][51.0, 19.0]


[10.366888337270446, 2.7662911860695725, 0.40820779221242043, 1.10270681903917, 0, 0.0, 0, 0, 0, 0.0, 0, 0.0, 0, 0, 0]Any[298.26850955033854, 316.3216626808096, 339.11626308505424, 306.86459163763, 304.8409546807165, 328.48924090824477, 285.53492366025165, 282.3772869861773, 314.7033129926668, 256.6636823648396, 213.77663991369454, 237.20557138546152, 231.83019690637218, 216.59817581494943, 203.6884144096667]Any[0.17507920645939315, 1.0697135143869057, 0.0, 0.0, 0.10033388370098634, 0.0, 0, 0, 0.0, 0, 0.0, 0, 0.0, 0.0, 0]Any[383.0957360999701, 304.19343152388456, 335.19512684920176, 370.0509358216024, 297.4139056660129, 309.1561484187586, 243.00897409744545, 257.52462958292085, 313.87926472511947, 274.53361345890335, 219.88352506436675, 242.60156024838182, 213.06138303377227, 201.62922361036988, 230.32635828512113]Any[290.9074832559619, 290.2511347602454, 290.9074832559619, 290.9074832559619, 290.9074832559619, 232.9935200326406, 232.9935200326406, 232.9935200326406, 232.9935200326406, 232.9935200326406, 187.3037926213792, 187.3037926213792, 187.3037926213792, 187.3037926213792, 187.40637342395178]Any[9.764705882352942, 9.941176470588236, 10.833333333333334, 10.333333333333334, 11.0, 5.0, 5.076923076923077, 5.076923076923077, 5.076923076923077, 5.076923076923077, 3.473684210526316, 3.4210526315789473, 3.473684210526316, 3.473684210526316, 3.473684210526316]Any[1.1764705882352942, 1.0588235294117647, 1.1176470588235294, 1.1764705882352942, 1.2352941176470589, 1.0294117647058822, 1.0, 1.2058823529411764, 1.1538461538461537, 1.0588235294117647, 1.0526315789473684, 1.2105263157894737, 1.2105263157894737, 1.0196078431372548, 1.1176470588235294]Any[10.833333333333334, 10.666666666666666, 11.0, 10.833333333333334, 11.0, 5.076923076923077, 5.076923076923077, 5.076923076923077, 5.076923076923077, 5.076923076923077, 3.473684210526316, 3.473684210526316, 3.4210526315789473, 3.4210526315789473, 3.473684210526316]Any[1.1764705882352942, 1.2352941176470589, 1.1764705882352942, 1.1176470588235294, 1.1764705882352942, 1.0588235294117647, 1.0588235294117647, 1.0294117647058822, 1.088235294117647, 1.1538461538461537, 1.2105263157894737, 1.1372549019607843, 1.1764705882352942, 1.1372549019607843, 1.2105263157894737]Any[1.2352941176470589, 1.1764705882352942, 1.2352941176470589, 1.2352941176470589, 1.2352941176470589, 1.1176470588235294, 1.1176470588235294, 1.1176470588235294, 1.1176470588235294, 1.1176470588235294, 1.0196078431372548, 1.0196078431372548, 1.0196078431372548, 1.0196078431372548, 1.0]
"""

# ╔═╡ 991f3042-7e2e-4239-8da2-0513e8980b5a
"""
[9.0, 1.0][18.0, 2.0][27.0, 3.0]

[971.3899527981814, 971.3899527981814, 971.3899527981814, 971.3899527981814, 809.1488871595795, 
784.1083730561743, 631.1572327130186, 470.2632801574731, 411.7231425289309, 333.4886427053707, 
340.0870540530564, 305.1756515801352, 224.30722814923107, 200.3467437781408, 168.32451023505934]


[717.09461008427, 717.09461008427, 717.09461008427, 717.09461008427, 717.09461008427, 
552.3612851579877, 549.1730261045398, 575.2383859662033, 566.9413248875428, 566.9413248875428, 
430.4150261763669, 430.4150261763669, 430.4150261763669, 430.85263879831155, 480.311845730163]


[971.3899527981814, 971.3899527981814, 971.3899527981814, 924.2391446858551, 890.4044729418845, 
894.0601205056287, 822.8101843484418, 609.4855249735381, 535.1455813287658, 433.81630270649487, 
431.1856210439554, 285.24852651565084, 218.4950228974721, 154.9539389078236, 128.81802826277485]



[717.09461008427, 717.09461008427, 717.09461008427, 717.09461008427, 717.09461008427, 
552.3612851579877, 557.8319615239002, 552.3612851579877, 549.1730261045398, 549.1730261045398, 
431.18562104395534, 422.33352195232686, 430.2517747715022, 426.7256357093375, 426.7256357093375]


[720.9314524976058, 717.09461008427, 717.09461008427, 726.2055174662953, 717.09461008427, 
549.1730261045398, 549.1730261045398, 549.1730261045398, 549.1730261045398, 549.1730261045398, 
422.4639233138159, 422.4639233138159, 430.85263879831155, 422.4639233138159, 422.4639233138159]


[0.0, 0.0, 0.0, 0.0, 0.8888888888888888, 
0.5555555555555556, 1.1666666666666667, 2.111111111111111, 2.4444444444444446, 3.0555555555555554, 
2.0, 2.2222222222222223, 2.4814814814814814, 2.740740740740741, 2.740740740740741]


[5.0, 5.0, 5.0, 5.0, 5.0, 
3.5, 4.0, 3.5, 3.5, 3.5, 
2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665]


[0.0, 0.0, 0.0, 1.0, 2.0, 
1.0, 2.0, 3.0, 4.0, 4.0, 
2.6666666666666665, 2.6666666666666665, 3.0, 3.3333333333333335, 3.3333333333333335]


[5.0, 5.0, 5.0, 5.0, 5.0, 
3.5, 3.5, 3.5, 4.0, 4.0, 
2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665]


[5.0, 5.0, 5.0, 5.0, 5.0, 
4.0, 4.0, 4.0, 4.0, 4.0, 
2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665]
"""

# ╔═╡ 2d95c10f-5ce5-4c48-947c-091c6f0c1e52
"""
[9.0, 1.0][18.0, 2.0][27.0, 3.0]

[971.3899527981814, 971.3899527981814, 971.3899527981814, 925.0679477763908, 822.9492834498217, 
760.6330040721729, 616.5141018654821, 463.5814965413115, 397.0225816928832, 333.4886427053707, 
340.0870540530564, 313.63103810463963, 221.92964036329636, 215.3937564378681, 158.96921995766908]

[717.09461008427, 717.09461008427, 717.09461008427, 717.09461008427, 720.9314524976058, 
549.1730261045398, 558.283933486565, 566.9413248875428, 582.2628025380975, 566.9413248875428, 
430.4150261763669, 433.72472073319705, 422.4639233138159, 481.7981172808462, 478.70263274603286]

[971.3899527981814, 890.4044729418845, 840.5454750859809, 702.1150352658117, 609.4855249735381, 
535.1455813287658, 355.60547887034755, 213.15477750780968, 137.9048243505141, 127.37265530005028, 
128.81802826277485, 79.37276963764742, 58.184270426087686, 43.70624017701824, 34.16514678751104]

[717.09461008427, 717.09461008427, 717.09461008427, 720.9314524976058, 715.4607761316934, 
549.1730261045398, 553.3793728913448, 568.2477393213968, 553.3793728913448, 582.2628025380975, 
426.7256357093375, 426.7256357093375, 472.30092614751067, 472.30092614751067, 472.30092614751067]

[715.4607761316934, 717.09461008427, 717.09461008427, 717.09461008427, 717.09461008427, 
552.3612851579877, 552.3612851579877, 549.1730261045398, 549.1730261045398, 549.1730261045398, 
422.4639233138159, 422.4639233138159, 422.4639233138159, 422.4639233138159, 422.4639233138159]

[0.0, 0.0, 0.0, 0.2222222222222222, 0.8888888888888888, 
0.6111111111111112, 1.2777777777777777, 2.111111111111111, 2.5555555555555554, 3.0555555555555554, 
2.0, 2.185185185185185, 2.5185185185185186, 2.5925925925925926, 2.8518518518518516]

[5.0, 5.0, 5.0, 5.0, 5.0, 
4.0, 4.0, 3.5, 3.5, 3.5, 
2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665]

[0.0, 2.0, 3.0, 6.0, 6.0, 
4.0, 4.0, 4.5, 5.0, 4.5, 
3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335, 3.3333333333333335]

[5.0, 5.0, 5.0, 5.0, 5.0, 
4.0, 4.0, 3.5, 4.0, 3.5, 
2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665]

[5.0, 5.0, 5.0, 5.0, 5.0, 
3.5, 3.5, 4.0, 4.0, 4.0, 
2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665, 2.6666666666666665]
"""

# ╔═╡ ea98cd79-f4e9-4051-b4ae-8fbf0fb9a4d0
"""
236-36 k=5
realfair = [367.5821006090546, 367.5821006090546, 339.9308564556597, 324.91550864571605, 272.3059567385554, 182.20278761737984, 124.58216223023297, 90.55475239111873, 73.66575579289916, 63.96865202348378]

fair = [362.02928691652636, 348.82720700301877, 330.45394121210575, 321.09102372307365, 319.99407605919225, 312.1936318708818, 312.97246784726815, 299.7131669252102, 290.7727818404495, 279.0044179356801]

realout = [367.5821006090546, 367.5821006090546, 339.9308564556597, 270.5758110641916, 195.77250100589112, 71.00484796415243, 45.77768750558939, 21.117561447561904, 11.79367017992354, 1.412367146349064]

out = [362.02928691652636, 348.82720700301877, 330.45394121210575, 333.9627353166192, 319.99407605919225, 332.88764369251544, 328.9492666648642, 297.9080953211539, 302.414742794426, 380.3514376379319]

basic = [362.02928691652636, 352.9083335819541, 341.4465680640389, 334.7959774145329, 325.0522513255055, 318.8398686249275, 312.918078468806, 304.44127503990313, 298.9740745522122, 290.9074832559619]


goal = [[1.0, 0.0][3.0, 1.0][5.0, 1.0][6.0, 2.0][8.0, 3.0][10.0, 3.0][11.0, 4.0][13.0, 5.0][15.0, 5.0][17.0, 6.0]]

realfair = [NaN, NaN, Inf, Inf, Inf, Inf, 8.909090909090908, 7.884615384615384, 6.761904761904762, 5.223529411764705]

fair = [NaN, Inf, Inf, Inf, 3.75, 3.5999999999999996, 5.090909090909091, 3.076923076923077, 1.888888888888889, 1.6764705882352942]

realout = [NaN, NaN, Inf, 1.5, 1.725, 1.1447811447811447, 1.0043478260869567, 1.0207100591715976, 1.0784313725490196, 1.1094377510040159]

out = [NaN, Inf, Inf, Inf, 3.75, 3.5999999999999996, 5.090909090909091, 1.9230769230769231, 3.0, 1.2395833333333335]

basic = [NaN, Inf, 1.0, 2.3333333333333335, 3.75, 1.6500000000000001, 2.3636363636363638, 3.076923076923077, 3.0, 3.7058823529411766]
"""

# ╔═╡ b9960958-2ffc-4537-8ac0-292208a3df0a
# ╠═╡ disabled = true
#=╠═╡
begin
	#4520-678
	costR = [9790.402511961718, 8459.326506499316, 7995.062633673838, 7758.692455648049, 8857.577546470315, 6953.4211597505655]
		
	costF = [6645.566847089055, 6661.844155741353, 6518.948386942639, 6340.6134522438215, 6195.423414557967, 5812.331255832305]
		
	costO = [6616.538004995221, 6546.531086061736, 6529.493555819822, 6331.768985944946, 6162.136599424345, 5905.885523801596]
		
	costB = [6559.942668286338, 6551.4484637407895, 6518.5426538994, 6282.47152045906, 6229.78941401983, 5560.64067723848]
		
	dispR = [NaN, NaN, Inf, 1.5866666666666667, 1.8148148148148149, 2.155327342747112]
		
	dispF = [NaN, NaN, Inf, 1.5866666666666667, 1.8148148148148149, 2.290485829959514]
		
	dispO = [NaN, NaN, Inf, 1.5866666666666667, 2.1, 2.2212171052631575]
		
	dispB = [NaN, NaN, Inf, 1.5866666666666667, 1.8148148148148149, 1.4844497607655502]
end
  ╠═╡ =#

# ╔═╡ 42757bd3-1eb3-4c65-8402-6032b2e29663
"""
2266 k=5
Any[3959.5528061240675, 4597.539354641413, 5051.377905971037, 3536.537580456031, 4088.251669766518, 3272.9768268686807, 3319.3615064326023, 3569.4879885146765, 3115.064837094923, 3140.0167828371023]Any[3193.851884482705, 3136.281938483185, 3008.651579073222, 2954.885819612603, 2963.630458592958, 2762.5185655416526, 2750.263914113105, 2742.2574730726546, 2700.7204048229178, 2599.0642536686955]Any[3203.1083429131013, 3101.6963632383804, 3038.60098771188, 2964.4603453471186, 3006.4813394176654, 2870.6534154393553, 2897.071116487048, 2849.874032163973, 3072.709112393235, 2971.8166564439307]Any[3178.446453928473, 3103.587291000106, 3035.4052225626374, 2949.539766375231, 2886.3795578638915, 2844.0533226667867, 2741.8360966546725, 2691.573109437083, 2620.5907996927067, 2542.784451701784]Any[2.1, 1.75, 2.488888888888889, 2.714285714285714, 2.576296296296296, 2.1999999999999997, 1.7262585034013607, 1.4986225895316805, 1.5075000000000003, 1.2302598064187467]Any[2.1, 1.2307692307692306, 1.2350877192982457, 1.2307692307692306, 2.576296296296296, 2.3125, 2.0965079365079364, 1.7479338842975207, 1.4305712669683257, 1.313821832941679]Any[2.1, 1.5454545454545454, 1.2350877192982457, 2.9615384615384612, 2.576296296296296, 2.0961538461538463, 1.7262585034013607, 1.5439519158527422, 1.394284128745838, 1.582045621780721]Any[2.1, 1.5454545454545454, 1.070899470899471, 1.2999999999999998, 1.1317647058823528, 1.034090909090909, 1.0290310650887573, 1.0515816471929325, 1.0, 1.1540708998831322]
"""

# ╔═╡ 8aa3d11d-2e81-4057-8b77-5fb261b2b6f5
"""
2266-340-sex
rand_cost_k5 = [4105.520220477397, 3892.3363388245643, 3353.567130223865, 3659.0115625799763, 4023.6243101511272, 3744.6094368562494, 3144.309427571429, 3048.2181539263415, 3119.9101520001796, 3975.0223992696415]

fair_cost_k5 = [3216.775552090925, 3101.322565455392, 3025.661682257725, 2949.514843539642, 2856.773170730138, 2781.6704111094896, 2776.3155048297112, 2662.425534966781, 2648.964928577097, 2599.8703763511603]

outlier_cost_k5 = [3212.790435513854, 3124.148053823756, 3025.9508794517296, 2954.7207609590255, 2953.9545921105714, 2917.5283992734753, 3015.7258791695176, 3000.948015653535, 3167.459853991453, 3907.643590170008]

basic_cost_k5 = [3210.5396356611536, 3106.4791142409017, 3019.1993934625957, 2965.4990905964805, 2880.6332407301634, 2791.862362516795, 2750.127345406527, 2734.253452307608, 2642.0034835970464, 2567.574760936451]

rand_disp_k5 = [2.1, 2.714285714285714, 2.786666666666667, 3.25, 1.2333333333333334, 2.0, 1.6645502645502646, 1.8054672600127146, 1.5075000000000003, 1.661996943453897]

fair_disp_k5 = [2.1, 1.5454545454545454, 1.2350877192982457, 1.1666666666666667, 2.576296296296296, 2.1999999999999997, 2.0965079365079364, 1.8054672600127146, 1.5483193277310925, 1.4058816926703332]

outlier_disp_k5 = [2.1, 1.2307692307692306, 1.3308641975308644, 3.25, 2.576296296296296, 2.0961538461538463, 1.8608946608946608, 1.4986225895316805, 1.7297385620915033, 1.9435646415202994]

basic_disp_k5 = [2.1, 1.75, 1.2350877192982457, 1.2307692307692306, 1.2333333333333334, 1.034090909090909, 1.0915451895043733, 1.0766949152542373, 1.0699678308823528, 1.1301468471062481]

rand_fair_outlier_basic_cens_k5 = [[[18.0, 4.0], [18.0, 4.0], [18.0, 4.0], [18.0, 4.0]], [[38.0, 7.0], [34.0, 11.0], [32.0, 13.0], [35.0, 10.0]], [[57.0, 10.0], [48.0, 19.0], [49.0, 18.0], [48.0, 19.0]], [[78.0, 12.0], [63.0, 27.0], [78.0, 12.0], [64.0, 26.0]], [[80.0, 32.0], [94.0, 18.0], [94.0, 18.0], [80.0, 32.0]], [[108.0, 27.0], [110.0, 25.0], [109.0, 26.0], [91.0, 44.0]], [[121.0, 36.0], [127.0, 30.0], [124.0, 33.0], [108.0, 49.0]], [[142.0, 39.0], [142.0, 39.0], [136.0, 45.0], [118.0, 63.0]], [[153.0, 50.0], [154.0, 49.0], [158.0, 45.0], [139.0, 64.0]], [[174.0, 52.0], [167.0, 59.0], [180.0, 46.0], [157.0, 69.0]]]

p_out = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
"""

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
