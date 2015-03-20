function fftFEVD(data, lags, H, typ)
	(n, k) = size(data)
	a = VAREST(data, lags, typ)

	Phi(a, H)

	P = chol(a.Σ)

	ft = mapslices(fft, a.Phi, [3])

	decomp = [(abs(ft[:,i,z]'*P'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
	denom = mapslices(sum, decomp, [3, 2])

	return  mapslices((x) -> x./denom[:,:,1], decomp, [1,2])
end

function fftGenFEVD(data, lags, H, typ, nocorr = false)
	(n, k) = size(data)
	a = VAREST(data, lags, typ)

	Phi(a, H)

	if nocorr
		Σ = diagm(diag(a.Σ))
	else
		Σ = a.Σ
	end

	P = chol(Σ)

	ft = mapslices(fft, a.Phi, [3])

	decompNew = [(abs(ft[:,i,z]'*Σ'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
	decomp = [(abs(ft[:,i,z]'*P'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
	denom = mapslices(sum, decomp, [3, 2])

	θ = zeros(k,k,H)

	for i = 1:k
		for j = 1:k
			for h = 1:H
				θ[i,j,h] = decompNew[i,j,h]/(denom[i,1,1]*Σ[i,i])
			end
		end
	end

	div = mapslices(sum, θ, [1,3])

	for i = 1:k
		for j = 1:k
			for h = 1:H
				θ[i,j,h] = θ[i,j,h] / div[1,j,1]
			end
		end
	end

	return θ
end

# function fftGenFEVD(data, lags, H, typ)
# 	(n, k) = size(data)
# 	a = VAREST(data, lags, typ)

# 	Phi(a, H)

# 	P = chol(a.Σ)

# 	ft = mapslices(fft, a.Phi, [3])

# 	decompNew = [(abs(ft[:,i,z]'*a.Σ'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
# 	decomp = [(abs(ft[:,i,z]'*P'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
# 	denom = mapslices(sum, decomp, [3, 2])

# 	θ = zeros(k,k,H)

# 	for i = 1:k
# 		for j = 1:k
# 			for h = 1:H
# 				θ[i,j,h] = decompNew[i,j,h]/(denom[i,1,1]*a.Σ[i,i])
# 			end
# 		end
# 	end

# 	div = mapslices(sum, θ, [1,3])

# 	for i = 1:k
# 		for j = 1:k
# 			for h = 1:H
# 				θ[i,j,h] = θ[i,j,h] / div[1,j,1]
# 			end
# 		end
# 	end

# 	# θ = [decompNew[i,j,h]/(denom[i,1,1]*a.Σ[i,i]) for i=1:k, j=1:k, h = 1:H]
# 	# θ = [θ[i,j,h]/sum(θ[:,j,:]) for i=1:k, j=1:k, h = 1:H]
# 	# θ = convert(Array{Float64}, θ)

# 	return θ
# end

function fftSpillover09(data, lags, H, typ)
	(n, k) = size(data)
	decomp = fftFEVD(data, lags, H, typ)
	return 1-sum(mapslices((x)->sum(diag(x))/k, decomp, [1, 2]))
end

function fftSpillover12(data, lags, H, typ)
	(n, k) = size(data)
	decomp = fftGenFEVD(data, lags, H, typ)
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
	decomp = fftFEVD(data, lags, H, typ)
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
	decomp = fftGenFEVD(data, lags, H, typ)
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
	decomp = fftGenFEVD(data, lags, H, typ, nocorr)
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
	decomp = fftFEVD(data, lags, H, typ)
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
	decomp = fftGenFEVD(data, lags, H, typ, nocorr)
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