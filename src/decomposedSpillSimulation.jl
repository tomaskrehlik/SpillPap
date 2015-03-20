# Searching for some process that has nicely distributed variance
using SpillPap

burnout = 100
sims=1000

A = [0.045, 0.045]

B = [-0.1 0.00;
	 0.00 -0.1]

C = [1.00 0.00;
	 0.00 1.00]

R = [0.98 0.00;
	 0.00 0.98]

P = [0.01, 0.01]

O = [0.01, 0.03]

(qdata, sdata) = simBivar(A, B, C, O, P, R)

(qdata, sdata) = (qdata[:,burnout:end], sdata[:,burnout:end])

levels = 4

dat = modwt((qdata+sdata)',"la8",levels,"periodic")

getEnergies(dat, 4, 2)

# Simulations

addprocs(3)

@everywhere using SpillPap

function pointwiseArrayMap(f::Function, a)
	s1, s2 = size(a[1])
	return [f(convert(Vector{Float64},map((x)->x[i,j], a))) for i=1:s1, j=1:s2]
end

@everywhere begin
	burnout = 100
	sims = 1000
	length = 2000
	function getSimulation(A, B, C, O, P, R, length, burnout, filter, restricted, egls)		

		(qdata, sdata) = simBivar(A, B, C, O, P, R, length)

		(qdata, sdata) = (qdata[:,burnout:end], sdata[:,burnout:end])

		return spilloverFinal((qdata+sdata)', 4, 10, filter, restricted, egls)
	end
end

		# The original values
	# A = [0.01, 0.01]

	# B = [0.75 0.00;
	# 	 0.00 0.75]

	# C = [1.00 0.00;
	# 	 0.00 1.00]

	# R = [0.99 0.00;
	# 	 0.00 0.99]

	# P = [0.01, 0.01]

	# O = [0.01, 0.03]

# Set the default values

sims=100
@everywhere begin
	A = [0.01, 0.01]

	B = [-0.2 0.00;
		 0.00 -0.2]

	C = [1.00 0.00;
		 0.00 1.00]

	R = [0.70 0.00;
		 0.00 0.70]

	P = [0.01, 0.01]

	O = [0.01, 0.03]
end

@time withoutSpills = pmap((x) -> getSimulation(A, B, C, O, P, R, length, burnout, "la16", true, true), 1:sims);
pointwiseArrayMap(mean, [withoutSpills[i][1] for i=1:sims])
pointwiseArrayMap(std, [withoutSpills[i][1] for i=1:sims])
pointwiseArrayMap(mean, [withoutSpills[i][3] for i=1:sims])
pointwiseArrayMap(std, [withoutSpills[i][3] for i=1:sims])
mean([withoutSpills[i][2] for i=1:sims])
std(convert(Vector{Float64},[withoutSpills[i][2] for i=1:sims]))

@everywhere begin
	C = [1.00 0.9;
		 0.9 1.00]
end

@time covarianceSpillsRes = pmap((x) -> getSimulation(A, B, C, O, P, R, length, burnout, "la16", true, true), 1:sims);
pointwiseArrayMap(mean, [covarianceSpillsRes[i][1] for i=1:sims])
pointwiseArrayMap(std, [covarianceSpillsRes[i][1] for i=1:sims])
pointwiseArrayMap(mean, [covarianceSpillsRes[i][3] for i=1:sims])
pointwiseArrayMap(std, [covarianceSpillsRes[i][3] for i=1:sims])
mean([covarianceSpillsRes[i][2] for i=1:sims])
std(convert(Vector{Float64},[covarianceSpillsRes[i][2] for i=1:sims]))



# Set the default values
@everywhere begin
	A = [0.01, 0.01]

	B = [-0.2 0.00;
		 0.00 -0.2]

	C = [1.00 0.00;
		 0.00 1.00]

	R = [0.70 0.00;
		 0.00 0.70]

	P = [0.01, 0.01]

	O = [0.01, 0.03]
end
@everywhere begin
	B = [-0.2 0.8;
		  0.8 -0.2]
end
@time shortTermSpills = pmap((x) -> getSimulation(A, B, C, O, P, R, length, burnout, "la16", true, true), 1:sims);
pointwiseArrayMap(mean, [shortTermSpills[i][1] for i=1:sims])
pointwiseArrayMap(std, [shortTermSpills[i][1] for i=1:sims])
pointwiseArrayMap(mean, [shortTermSpills[i][3] for i=1:sims])
pointwiseArrayMap(std, [shortTermSpills[i][3] for i=1:sims])
mean([shortTermSpills[i][2] for i=1:sims])
std(convert(Vector{Float64},[shortTermSpills[i][2] for i=1:sims]))


@everywhere begin
	A = [0.01, 0.01]

	B = [-0.2 0.00;
		 0.00 -0.2]

	C = [1.00 0.00;
		 0.00 1.00]

	R = [0.70 0.00;
		 0.00 0.70]

	P = [0.01, 0.01]

	O = [0.01, 0.03]
end
@everywhere begin
	R = [0.70 0.29;
		 0.29 0.70]
end

@time longTermSpills = pmap((x) -> getSimulation(A, B, C, O, P, R, length, burnout, "la8", true, true), 1:sims);
pointwiseArrayMap(mean, [longTermSpills[i][1] for i=1:sims])
pointwiseArrayMap(std, [longTermSpills[i][1] for i=1:sims])
pointwiseArrayMap(mean, [longTermSpills[i][3] for i=1:sims])
pointwiseArrayMap(std, [longTermSpills[i][3] for i=1:sims])
mean([longTermSpills[i][2] for i=1:sims])
std(convert(Vector{Float64},[longTermSpills[i][2] for i=1:sims]))

[diag(pointwiseArrayMap(mean, [withoutSpills[i][1] for i=1:sims])) diag(pointwiseArrayMap(std, [withoutSpills[i][1] for i=1:sims])) diag(pointwiseArrayMap(mean, [covarianceSpillsRes[i][1] for i=1:sims])) diag(pointwiseArrayMap(std, [covarianceSpillsRes[i][1] for i=1:sims])) diag(pointwiseArrayMap(mean, [shortTermSpills[i][1] for i=1:sims])) diag(pointwiseArrayMap(std, [shortTermSpills[i][1] for i=1:sims])) diag(pointwiseArrayMap(mean, [longTermSpills[i][1] for i=1:sims])) diag(pointwiseArrayMap(std, [withoutSpills[i][1] for i=1:sims]))]'
[mean([withoutSpills[i][2] for i=1:sims]) std(convert(Vector{Float64},[withoutSpills[i][2] for i=1:sims])) mean([covarianceSpillsRes[i][2] for i=1:sims]) std(convert(Vector{Float64},[covarianceSpillsRes[i][2] for i=1:sims])) mean([shortTermSpills[i][2] for i=1:sims]) std(convert(Vector{Float64},[shortTermSpills[i][2] for i=1:sims])) mean([longTermSpills[i][2] for i=1:sims]) std(convert(Vector{Float64},[longTermSpills[i][2] for i=1:sims]))]'






@everywhere begin
	B = [-0.7 0.29;
		  0.29 -0.7]
end
@time shortTermSpills12 = pmap((x) -> getSimulation(A, B, C, O, P, R, length, burnout, "la16", true, true), 1:sims);
pointwiseArrayMap(mean, [shortTermSpills12[i][1] for i=1:sims])
pointwiseArrayMap(std, [shortTermSpills12[i][1] for i=1:sims])
pointwiseArrayMap(mean, [shortTermSpills12[i][3] for i=1:sims])
pointwiseArrayMap(std, [shortTermSpills12[i][3] for i=1:sims])
mean([shortTermSpills12[i][2] for i=1:sims])
std(convert(Vector{Float64},[shortTermSpills12[i][2] for i=1:sims]))
