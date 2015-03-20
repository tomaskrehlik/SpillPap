
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
	space = [(0:fld(len,2))/fld(len,2)]*π
	lb = space.>down
	ub = space.<=up
	output = (lb & ub)
	if (len%2 == 0)
		output = [output, reverse(output)][2:(end-1)]	
	else
		output = [output, reverse(output)][1:(end-1)]
	end
	return output
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

function fftSpilloverDec12(data, lags, H, typ, bounds, proportions)
	est = varEstimate(data, lags, typ)
	decomp = fftGenFEVD(est, H)
	(k, k, H) = size(decomp)
	levels = size(bounds)[1]-2
	output = zeros(Float64, levels + 1)
	# bounds = [[pi/2^j for j = 0:levels],0]
	# bounds = [pi, [pi/j for j = [5, 25, 75]],0]
	if proportions
		for i = 1:(levels+1)
			ind = getIndexes(H, bounds[i], bounds[i+1])
			temp = decomp[:,:,ind]
			output[i] = (sum(temp) - sum(mapslices((x) -> sum(diag(x)), temp, [1, 2])))./(sum(temp))
		end
	else
		for i = 1:(levels+1)
			ind = getIndexes(H, bounds[i], bounds[i+1])
			temp = decomp[:,:,ind]
			output[i] = (sum(temp) - sum(mapslices((x) -> sum(diag(x)), temp, [1, 2])))./k
		end
	end
	return output
end

function fftSpilloverDec12(data, lags, H, typ, bounds, proportions, nocorr = false)
	est = varEstimate(data, lags, typ)
	decomp = fftGenFEVD(est, H, nocorr)
	(k, k, H) = size(decomp)
	levels = size(bounds)[1]-2
	output = zeros(Float64, levels + 1)
	# bounds = [[pi/2^j for j = 0:levels],0]
	# bounds = [pi, [pi/j for j = [5, 25, 75]],0]
	if proportions
		for i = 1:(levels+1)
			ind = getIndexes(H, bounds[i], bounds[i+1])
			temp = decomp[:,:,ind]
			output[i] = (sum(temp) - sum(mapslices((x) -> sum(diag(x)), temp, [1, 2])))./(sum(temp))
		end
	else
		for i = 1:(levels+1)
			ind = getIndexes(H, bounds[i], bounds[i+1])
			temp = decomp[:,:,ind]
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

function fftSpilloverTableDec12(data, lags, H, typ, bounds, nocorr = false)
	est = varEstimate(data, lags, typ)
	decomp = fftGenFEVD(est, H, nocorr)
	(k, k, H) = size(decomp)
	levels = size(bounds)[1]-2
	output = zeros(Float64, k, k, levels + 1)
	
	for i = 1:(levels+1)
		ind = getIndexes(H, bounds[i], bounds[i+1])
		temp = decomp[:,:,ind]
		output[:,:,i] = mapslices(sum, temp, [3])
	end
	return output
end