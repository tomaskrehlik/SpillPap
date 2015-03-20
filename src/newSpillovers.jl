function generalisedWavSpillover(estimate::varEstimate, H::Int, levels::Int)
	series = convert(Int, estimate.vars/(levels + 1))
	θ = genFEVD(estimate, H)
	diags = [sum(diag(θ[(1:series)+series*i, (1:series)+series*j])) for i=0:(levels), j=0:(levels)]./series
	tot = [sum(θ[(1:series)+series*i, (1:series)+series*j]) for i=0:(levels), j=0:(levels)]./series
	return ((tot-diags)/(levels+1), 1-sum(diags)/sum(tot), tot)
end

function spilloverOurs(data, levels, H, filter, restricted, egls)
	(obs, series) = size(data)
	dat = modwt(data, filter, levels, "periodic")
	
	(D, S) = mra(dat)

	# Make the threeway system

	L = deepcopy(S[:,end,:])
	M = [mapslices(sum, D[:,2:end,i], [2]) for i=1:series]
	S = deepcopy(D[:,1,:])

	system = [apply(hcat, [S[:,:,i] for i=1:series]) apply(hcat,M) apply(hcat, [L[:,:,i] for i=1:series])]

	if restricted
		t = fill(true, series, series)
		f = fill(false, series, series)
		c = fill(true, series)
		mat = [[c t f f], [c f t f], [c f f t]]
		return generalisedWavSpillover(restrictVAR2(varEstimate(system, 5, "Const"), mat, egls), H, 2)
	else
		return generalisedWavSpillover(varEstimate(system, 5, "Const"), H, 2)
	end
end

function spilloverOursRolling(data, levels, H, filter, restricted, egls, window)
	(obs, series) = size(data)
	dat = modwt(data, filter, levels, "periodic")
	
	(D, S) = mra(dat)

	# Make the threeway system

	L = deepcopy(S[:,end,:])
	M = [mapslices(sum, D[:,2:end,i], [2]) for i=1:series]
	S = deepcopy(D[:,1,:])

	system = [apply(hcat, [S[:,:,i] for i=1:series]) apply(hcat,M) apply(hcat, [L[:,:,i] for i=1:series])]

	if restricted
		t = fill(true, series, series)
		f = fill(false, series, series)
		c = fill(true, series)
		mat = [[c t f f], [c f t f], [c f f t]]
		out = @parallel (vcat) for i=0:(obs-window)
			generalisedWavSpillover(restrictVAR2(varEstimate(system[(1:window)+i,:], 1, "Const"), mat, egls), H, 2)
		end
		return out
	else
		return [generalisedWavSpillover(varEstimate(system[(1:window)+i,:], 2, "Const"), H, 2) for i=1:(obs-window)]
	end
end