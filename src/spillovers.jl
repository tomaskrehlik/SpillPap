using Wavelets
using Distributions

function testNonOverlapping(data, levels::Int, lags::Int)	
	(obs, series) = size(data)
	dat = modwt(data, "la8", levels, "periodic");
	(D, S) = mra(dat);
	ΣDs = [VAREST(apply(hcat, [D[:,j,i] for i=1:series]), lags, "Const").Σ for j=1:levels]
	ΣSs = VAREST(apply(hcat, [S[:,levels,i] for i=1:series]), lags, "Const").Σ

	da = hcat(apply(hcat, [apply(hcat, [D[:,j,i] for i=1:series]) for j=1:levels]), apply(hcat, [S[:,levels,i] for i=1:series]))
	
	LLik = obs*(log(prod(map(det, ΣDs)) * det(ΣSs)) - log(det(VAREST(da, lags, "Const").Σ)))

	rest = (series*levels)^2 - series^2*levels
	
	return ccdf(Chisq(rest), LLik)
end

function geomMean(a::Vector{Float64})
	return prod(a)^(1/length(a))
end

function getEnergies(dat, levels, series)
	v=mapslices(var, dat.W, [1])
	weights = [geomMean(convert(Array{Float64,1},[v[1,i,j] for j=1:series])) for i = 1:levels]
	s=mapslices(var, dat.V, [1])
	weights = [weights, geomMean(convert(Array{Float64,1},[s[1,levels,j] for j=1:series]))]
	weights = weights / sum(weights)
	return weights
end


# function spilloverFinal(data, levels::Int, H::Int, filter::ASCIIString, restricted::Bool, egls::Bool)
# 	(obs, series) = size(data)
# 	dat = modwt(data, filter, levels, "periodic")
	
# 	(D, S) = mra(dat)

# 	# Make the threeway system

# 	L = deepcopy(S[:,end,:])
# 	M = [mapslices(sum, D[:,2:end,i], [2]) for i=1:series]
# 	S = deepcopy(D[:,1,:])

# 	system = [apply(hcat, [S[:,:,i] for i=1:series]) apply(hcat,M) apply(hcat, [L[:,:,i] for i=1:series])]

# 	if restricted
# 		t = fill(true, series, series)
# 		f = fill(false, series, series)
# 		c = fill(true, series)
# 		mat = [[c t f f], [c f t f], [c f f t]]
# 		return generalisedWavSpillover(restrictVAR2(VAREST(system, 1, "Const"), mat, egls), H, 2)
# 	else
# 		return generalisedWavSpillover(VAREST(system, 1, "Const"), H, 2)
# 	end
# end

function generalisedWavDecompSpillover(data, levels::Int, lags::Int, H::Int, filter::ASCIIString, method::ASCIIString, restricted::Bool)
	(obs, series) = size(data)
	dat = modwt(data, filter, levels, "periodic");
	v=mapslices(var, dat.W, [1])
	weights = [geomMean(convert(Array{Float64,1},[v[1,i,j] for j=1:series])) for i = 1:levels]
	s=mapslices(var, dat.V, [1])
	weights = [weights, geomMean(convert(Array{Float64,1},[s[1,levels,j] for j=1:series]))]
	weights = weights / sum(weights)
	if method=="mra"
		(D, S) = mra(dat);
	else
		(D, S) = (dat.W, dat.V);
	end
	if restricted
		spillovers = [generalisedSpillover(restrictVAR(VAREST(apply(hcat, [D[:,j,i] for i=1:series]), lags, "Const"), 2.0), H) for j=1:levels]
		spillovers = [spillovers, generalisedSpillover(restrictVAR(VAREST(apply(hcat, [S[:,levels,i] for i=1:series]), lags, "Const"), 2.0), H)]
	else
		spillovers = [generalisedSpillover(VAREST(apply(hcat, [D[:,j,i] for i=1:series]), lags, "Const"), H) for j=1:levels]
		spillovers = [spillovers, generalisedSpillover(VAREST(apply(hcat, [S[:,levels,i] for i=1:series]), lags, "Const"), H)]
	end
	return (spillovers, weights)
end

