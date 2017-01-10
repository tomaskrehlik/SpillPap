
function fftSpillover09(data, lags, H, typ)
	(n, k) = size(data)
	est = varEstimate(data, lags, typ)
	decomp = fftFEVD(est, H)
	return 1-sum(mapslices((x)->sum(diag(x))/k, decomp, [1, 2]))
end

function fftSpillover12(data, lags, H, typ)
	(n, k) = size(data)
	est = varEstimate(data, lags, typ)
	decomp = fftGenFEVD(est, H)
	return 1-sum(mapslices((x)->sum(diag(x))/k, decomp, [1, 2]))
end

function getIndexes(len, up, down)
	space = collect(0:fld(len,2))./(fld(len,2)).*π
	lb = space.>=down
	ub = space.<up
	output = (lb & ub)
	if (len%2 == 0)
		output = [output; reverse(output[2:(end-1)])]
	else
		output = [output; reverse(output[2:end])]
	end
	return output
end

function testIndexes()
	bounds = [pi+0.00001, [pi/j for j = [5, 25, 75]], 0]
	
	covers = all(mapslices(any, apply(hcat, map((i) -> getIndexes(401, bounds[i], bounds[i+1]), [1:(size(bounds)[1]-1)])), [2]))
	covers2 = all(mapslices(any, apply(hcat, map((i) -> getIndexes(400, bounds[i], bounds[i+1]), [1:(size(bounds)[1]-1)])), [2]))

	return (covers & covers2)
end


function fftSpilloverDec09(data, lags, H, typ, bounds, proportions)
	est = varEstimate(data, lags, typ)
	decomp = fftFEVD(est, H)
	(k, k, H) = size(decomp)
	levels = size(bounds)[1]-2
	output = zeros(Float64, levels + 1)
	# bounds = [pi, [pi/j for j = [5, 25, 75]],0]
	if proportions
		for i = 1:(levels+1)
			ind = getIndexes(H, bounds[i], bounds[i+1])
			temp = decomp[:,:,ind]
			output[i] = (sum(temp) - sum(mapslices((x)->sum(diag(x)), temp, [1, 2])))./(k*sum(temp))
		end
	else
		for i = 1:(levels+1)
			ind = getIndexes(H, bounds[i], bounds[i+1])
			temp = decomp[:,:,ind]
			output[i] = (sum(temp) - sum(mapslices((x)->sum(diag(x)), temp, [1, 2])))./k
		end
	end
	return output
end

function fftSpilloverDec12(data, lags, H, typ, bounds, proportions; nocorr = false)
	est = varEstimate(data, lags, typ)
	decomp = fftGenFEVD(est, H; nocorr = nocorr)
	(k, k, H) = size(decomp)
	levels = size(bounds)[1]-2
	output = zeros(Float64, levels + 1)
	
	ind = [getIndexes(H, bounds[i], bounds[i+1]) for i = 1:(levels+1)]
	all_ind = mapslices(any, apply(hcat, ind), [1])

	if proportions
		for i = 1:(levels+1)
			# ind = getIndexes(H, bounds[i], bounds[i+1])
			temp = decomp[:,:,ind[i]]
			output[i] = (sum(temp) - sum(mapslices((x) -> sum(diag(x)), temp, [1, 2])))./(sum(temp))
		end
	else
		for i = 1:(levels+1)
			# ind = getIndexes(H, bounds[i], bounds[i+1])
			temp = decomp[:,:,ind[i]]
			output[i] = (sum(temp) - sum(mapslices((x) -> sum(diag(x)), temp, [1, 2])))./k
		end
	end
	return output
end

function fftSpilloverTableDec09(data, lags, H, typ, bounds)
	est = varEstimate(data, lags, typ)
	decomp = fftFEVD(est, H)
	(k, k, H) = size(decomp)
	levels = size(bounds)[1]-2
	output = zeros(Float64, k, k, levels + 1)
	# bounds = [pi, [pi/j for j = [5, 25, 75]],0]
	for i = 1:(levels+1)
		ind = getIndexes(H, bounds[i], bounds[i+1])
		temp = decomp[:,:,ind]
		output[:,:,i] = mapslices(sum, temp, [3])
		# output[i] = (sum(temp) - sum(mapslices((x)->sum(diag(x)), temp, [1, 2])))./k
	end
	return output
end

function from(table, typ)
	if typ=="within"
		return mapslices((x) -> (collect(sum(x, 1)) - diag(x))/sum(x), table, (1,2))
	elseif typ=="absolute"
		return mapslices((x) -> mapslices((j) -> (collect(sum(j, 1)) - diag(j))/sum(x), x, (1,2)), table, (1,2,3))
	else
		error("Wrong definition of which spillover you want, say either absolute or within.")
	end
end

function to(table, typ)
	if typ=="within"
		return mapslices((x) -> (collect(sum(x, 2)) - diag(x))/sum(x), table, (1,2))
	elseif typ=="absolute"
		return mapslices((x) -> mapslices((j) -> (collect(sum(j, 2)) - diag(j))/sum(x), x, (1,2)), table, (1,2,3))
	else
		error("Wrong definition of which spillover you want, say either absolute or within.")
	end
end

function pairwise(table, typ)
	function get_pairwise(tab)
		mat = (triu(tab) - tril(tab)')
		k = size(tab)[1]
		return map((x) -> getindex(mat, x...), collect(combinations(collect(1:k), 2)))
	end
	if typ=="within"
		return mapslices((x) -> get_pairwise(x)/sum(x), table, (1,2))
	elseif typ=="absolute"
		return mapslices((x) -> mapslices((j) -> (get_pairwise(j))/sum(x), x, (1,2)), table, (1,2,3))
	else
		error("Wrong definition of which spillover you want, say either absolute or within.")
	end
end

function overall(table, typ)
	if typ=="within"
		return mapslices((x) -> 1 - sum(diag(x))/sum(x), table, (1,2))
	elseif typ=="absolute"
		return mapslices((x) -> mapslices((j) -> (sum(j) - sum(diag(j)))/sum(x), x, (1,2)), table, (1,2,3))
	else
		error("Wrong definition of which spillover you want, say either absolute or within.")
	end
end

function fftSpilloverTableDec12(data, lags, H, typ, bounds, nocorr = false; boot = false, quantiles = (0.95, 0.05), reps = 1000)
	est = varEstimate(data, lags, typ)
	decomp = fftGenFEVD(est, H, nocorr)
	(k, k, H) = size(decomp)
	levels = size(bounds)[1]-2
	output = zeros(Float64, k, k, levels + 1)
	
	for i = 1:(levels+1)
		ind = getIndexes(H, bounds[i], bounds[i+1])
		temp = decomp[:,:,ind]
		output[:,:,i] = mapslices(sum, temp, collect(3))
	end

	# output = mapslices(augment, output, (1,2))

	if boot
		function get_sd(object, quantiles)
			ecdfs = mapslices(EmpiricalUnivariateDistribution, object, [4])
			return map((x) ->  quantile(x, quantiles), ecdfs)
		end			
		a = [fftSpilloverTableDec12(simulateVAR(est, est.obs), lags, 100, typ, bounds)[1] for i=1:reps]
		a = cat(4, a...)
		table_sd = get_sd(a, [down, up])
		fw_sd = get_sd(from(a, "within"), [down, up])
		fa_sd = get_sd(from(a, "absolute"), [down, up])
		tw_sd = get_sd(to(a, "within"), [down, up])
		ta_sd = get_sd(to(a, "absolute"), [down, up])
		ow_sd = get_sd(overall(a, "within"), [down, up])
		oa_sd = get_sd(overall(a, "absolute"), [down, up])
		pairwisew_sd = get_sd(pairwise(a, "within"), [down, up])
		pairwisea_sd = get_sd(pairwise(a, "absolute"), [down, up])

		return ((output, (overall(output, "within"), overall(output, "absolute")), (to(output, "within"), to(output, "absolute")), (from(output, "within"), from(output, "absolute")), (pairwise(output, "within"), pairwise(output, "absolute"))), (table_sd, (ow_sd, oa_sd), (tw_sd, ta_sd), (fw_sd, fa_sd), (pairwisew_sd, pairwisea_sd)))
	end

	return (output, (overall(output, "within"), overall(output, "absolute")), (to(output, "within"), to(output, "absolute")), (from(output, "within"), from(output, "absolute")), (pairwise(output, "within"), pairwise(output, "absolute")))
end