"""
The `spilloverTable` function computes the the spillover tables as in Diebold & Yilmaz (2009, 2012), DY,
and Barunik & Krehlik (2016), BK from now on. The function takes the following compulsory arguments:
	
	- `est` which is supposed to be an estimate of VAR or VECM from the VARmodels package,
	- `H` the horizon for which you want to compute the spillover,

The following are the named arguments that have default values.

	- `fevd_function` is the function to compute the spillover.
		- `fevd` corresponds to DY (2009)
		- `genFEVD` (*default*) corresponds to DY (2012)
		- `fftFEVD` corresponds to the frequency version of DY (2009) from BK (2016)
		- `fftGenFEVD` corresponds to the frequency version of DY (2012) from BK (2016)
	- `bounds` the frequency bounds how to decompose the table. Only used for 
	BK (2016) versions of the spillovers. The default is `nothing`.
	- `nocorr` if you want to nullify the collerations in the system. For motivation
	as why to do it, look into BK (2016). *Default* value is `false`, ie. the standard
	estimate.

The function return `KxKxB-1` array, where `K` is the number of variables within the system
and `B` is the size of the `bounds` vector. In the case for DY (2009, 2012), the `B = 2`

Use the function such as
````
using VARmodels
usign SpillPap

# Assuming that data hold some data you want to compute the spillover for.

lags = 2
type = \"Const\"
bounds = [π + 0.0001; π/5; π/20; π/60; 0]
spillTab = spilloverTable(varEstimate(data, lags, typ), 600; fevd_function = fftGenFEVD, bounds = bounds, nocorr = false)
````
"""
function spilloverTable(est::VARmodels.VARRepresentation, H; fevd_function = genFEVD, bounds = nothing, nocorr = false)
	@assert fevd_function in [fftGenFEVD, fftFEVD, fevd, genFEVD] "The only supported types of FEVD computation are [:fftGenFEVD, :fftFEVD, :fevd, :genFEVD]"
	

	if fevd_function in [fftGenFEVD, fftFEVD]
		@assert typeof(bounds) <: Vector "Bounds should be an instance of Vector containing bounds for decomposition"
		levels = size(bounds)[1]-2
		r = collect(1:H)[getIndexes(H, max(bounds...), min(bounds...))]
		r = (min(r...), max(r...))

		decomp = fevd_function(est, H; nocorr = nocorr, range = r)
		k = size(decomp)[1]
		
		output = zeros(Float64, k, k, levels + 1)
		for i = 1:(levels+1)
			temp = decomp[:,:,getIndexes(H, bounds[i], bounds[i+1])]
			output[:,:,i] = mapslices(sum, temp, collect(3))
		end
	else
		decomp = fevd_function(est, H; nocorr = nocorr)
		output = reshape(decomp, size(decomp)[1], size(decomp)[2], 1)
	end

	return 100*output
end

"""
The `spilloverMeasure` returns composite measures of spillovers over the spillover tables.
In the case of the non-frequency domain spillovers, see documentation for `spilloverTable` for
details, the possible measures are `overall, from, to, pairwise, net`. In case of the frequency
domain spillover tables, you have to additionally specify whether you want `absolute` or `within`
spillover.

The function takes in two arguments:

	- `table` output of the `spilloverTable` function.
	- `which` is an array specifying which spillovers you want to compute. In case of the non-frequency
	spillovers it is an array of strings, eg. `[\"overall\", \"from\"]`. In case of the frequency domain
	spillovers it is an array of tuples, eg. `[(\"overall\", \"absolute\"), (\"to\", \"within\")]`. Additionally
	for convenience, one can supply to the function string \"all\" saying that you want to compute all
	possible measures.

The function returns array of arrays where the ordering corresponds to the following:

	- `overall, from, to, net, pairwise` for non-frequency spillovers
	- `(overall, absolute), (from, absolute), (to, absolute), (net, absolute), (pairwise, absolute), (overall, within), (from, within), (to, within), (net, within), (pairwise, within)`

If you want multiple measures, it is better to specify them within one run of the function than to call the function
multiple times.

Use the function as
````
# For spillTab see the help for spilloverTable.
spilloverMeasure(spillTab, [\"overall\", \"to\", \"from\"])
````
"""
function spilloverMeasure(table, which)
	@assert typeof(table) <: Array "The first argument should be a spillover table produced by function spilloverTable."

	function get_pairwise(tab)
		mat = (triu(tab) - tril(tab)')
		k = size(tab)[1]
		return map((x) -> getindex(mat, x...), collect(combinations(collect(1:k), 2)))
	end

	if which=="all"&&size(table)[3]==1
		which = ["overall", "from", "to", "pairwise", "net"]
	elseif which=="all"&&size(table)[3]>1
		which = [("overall", "absolute"), ("from", "absolute"), ("to", "absolute"), ("pairwise", "absolute"), ("net", "absolute"), ("overall", "within"), ("from", "within"), ("to", "within"), ("pairwise", "within"), ("net", "within")]
	end

	output = []

	if size(table)[3]==1
		if "overall" in which
			append!(output, Any[mapslices((x) -> 1 - sum(diag(x))/sum(x), table, (1,2))])
		end
		if "from" in which
			append!(output, Any[mapslices((x) -> (collect(sum(x, 1)) - diag(x))/sum(x), table, (1,2))])
		end
		if "to" in which
			append!(output, Any[mapslices((x) -> (collect(sum(x, 2)) - diag(x))/sum(x), table, (1,2))])
		end
		if "net" in which
			t = spilloverMeasure(table, ["from"])[1]
			f = spilloverMeasure(table, ["to"])[1]
			append!(output, Any[reshape(vcat(map(x->f[:,:,x] - t[:,:,x], 1:size(table)[3])...), size(table)[1], 1, size(table)[3])]/100)
		end
		if "pairwise" in which
			append!(output, Any[mapslices((x) -> get_pairwise(x)/sum(x), table, (1,2))])
		end
	elseif size(table)[3]>1
		if ("overall", "absolute") in which
			append!(output, Any[mapslices((x) -> mapslices((j) -> (sum(j) - sum(diag(j)))/sum(x), x, (1,2)), table, (1,2,3))])
		end
		if ("from", "absolute") in which
			append!(output, Any[mapslices((x) -> mapslices((j) -> (collect(sum(j, 1)') - diag(j))/sum(x), x, (1,2)), table, (1,2,3))])
		end
		if ("to", "absolute") in which
			append!(output, Any[mapslices((x) -> mapslices((j) -> (collect(sum(j, 2)) - diag(j))/sum(x), x, (1,2)), table, (1,2,3))])
		end
		if ("net", "absolute") in which
			f = spilloverMeasure(table, [("from", "absolute")] )[1]
			t = spilloverMeasure(table, [("to", "absolute")] )[1]
			append!(output, Any[reshape(vcat(map(x->f[:,:,x] - t[:,:,x], 1:size(table)[3])...), size(table)[1], 1, size(table)[3])]/100)
		end
		if ("pairwise", "absolute") in which
			append!(output, Any[mapslices((x) -> mapslices((j) -> (get_pairwise(j))/sum(x), x, (1,2)), table, (1,2,3))])
		end
		if ("overall", "within") in which
			append!(output, Any[mapslices((x) -> 1 - sum(diag(x))/sum(x), table, (1,2))])
		end
		if ("from", "within") in which
			append!(output, Any[mapslices((x) -> (collect(sum(x, 1)') - diag(x))/sum(x), table, (1,2))])
		end
		if ("to", "within") in which
			append!(output, Any[mapslices((x) -> (collect(sum(x, 2)) - diag(x))/sum(x), table, (1,2))])
		end
		if ("net", "within") in which
			f = spilloverMeasure(table, [("from", "within")])[1]
			t = spilloverMeasure(table, [("to", "within")])[1]
			append!(output, Any[reshape(vcat(map(x->f[:,:,x] - t[:,:,x], 1:size(table)[3])...), size(table)[1], 1, size(table)[3])]/100)
		end
		if ("pairwise", "within") in which
			append!(output, Any[mapslices((x) -> get_pairwise(x)/sum(x), table, (1,2))])
		end
	end

	return 100*output
end

"""
The function `bootSpillovers` uses bootstrap to approximate the standard errors of spillovers. The usage should be self-evident
given the help to other spillover function. Additionally, other parameters that you need to specify are:

	- `reps` number of repetitions within the bootstrap

Named arguments
	- `burnout` the burnout for the simulation, *default* is 100.
	- `quantiles` the quantiles of empirical distribution that you
	 want to get, *default* is `[0.05, 0.95]`. It should be array of floats.

The function returns array of size 2 of arrays. The first one holds the estimated values of spillvovers. The second element
holds the quantiles for the respective estimates from element one.

Use as

````
window = 250
H = 600
p = 2
typ = \"Const\"
bounds = [π + 0.0001; π/5; π/20; π/60; 0]

bootSpillovers(varEstimate(data[(1:window)+1,:], p, typ), H, 100, \"all\"; fevd_function = fftGenFEVD, bounds = bounds, nocorr = false)
````
"""
function bootSpillovers(est::VARmodels.VARRepresentation, H, reps, which; fevd_function = genFEVD, bounds = nothing, nocorr = false, burnout = 100, quantiles = [0.05, 0.95], boot_sample_length = 0)
	# Compute spillovers
	if boot_sample_length==0
		boot_sample_length=est.obs
	end

	spillTable = spilloverTable(est::VARmodels.VARRepresentation, H; fevd_function = fevd_function, bounds = bounds, nocorr = nocorr)
	measures = spilloverMeasure(spillTable, which)
	if typeof(est)==VARmodels.varEstimate
		spillTableBooted = map(i->spilloverTable(varEstimate(simulateVAR(est, boot_sample_length, burnout), est.lags, est.typ), H; fevd_function = fevd_function, bounds = bounds, nocorr = nocorr), 1:reps)
	elseif typeof(est)==VARmodels.varEstimateShrink
		spillTableBooted = map(i->spilloverTable(varEstimateShrink(simulateVAR(est, boot_sample_length, burnout), est.lags, est.typ, est.λ), H; fevd_function = fevd_function, bounds = bounds, nocorr = nocorr), 1:reps)
	elseif typeof(est)==VARmodels.varRepresentationVECM
		spillTableBooted = map(i->spilloverTable(varRepresentationVECM(CointegrationEstimate(simulateVAR(est, boot_sample_length, burnout), "none", est.lags, "longrun"; season = "no", dumvar = "no"), est.r), H; fevd_function = fevd_function, bounds = bounds, nocorr = nocorr), 1:reps)
	else
		error("Unsupported type of estimate for bootstrap.")
	end	
	spillTableBooted = cat(length(size(spillTableBooted[1])) + 1, spillTableBooted...)
	
	sds_spillTable = map((x) ->  quantile(x, quantiles), mapslices(EmpiricalUnivariateDistribution, spillTableBooted, [4]))

	measuresBooted = mapslices((x) -> spilloverMeasure(x, which), spillTableBooted, 1:(length(size(spillTableBooted))-1))
	measuresBooted = [[hcat(map(x->x[:,:,i], collect(measuresBooted[j,1,1,:]))...) for i=1:size(measuresBooted[1,1,1,1])[3]] for j=1:size(measuresBooted)[1]]
	sds_measures = [[mapslices((x) -> quantile(EmpiricalUnivariateDistribution(x), quantiles), measuresBooted[i][j], [2]) for j=1:size(measuresBooted[1])[1]] for i=1:size(measures)[1]]
	return ((spillTable, measures), (sds_spillTable, sds_measures))
end