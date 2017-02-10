type SpilloverTableBooted
	estimate::SpilloverTable
	mean::SpilloverTable
	up::SpilloverTable
	down::SpilloverTable
	quantiles::Vector{Float64}
end

function show(io::IO, tab::SpilloverTableBooted)
	print("The object contains bootstrapped confidence intervals not shown in here.\n")
	printTable(tab.estimate)
end

function SpilloverTableBooted(est::VARmodels.VARRepresentation, H, reps; fevd_function = genFEVD, nocorr = false, burnout = 100, quantiles = [0.05, 0.95], boot_sample_length = 0, names = nothing)
	# Compute spillovers
	if boot_sample_length==0
		boot_sample_length=est.obs
	end

	if names==nothing
		names = ["v" * j for j=[sprintf1("%d", i) for i=1:est.vars]]
	end

	spillsEstimate = SpilloverTable(spilloverTable(est::VARmodels.VARRepresentation, H; fevd_function = fevd_function, nocorr = nocorr), names)

	if typeof(est)==VARmodels.varEstimate
		spillTableBooted = map(i->SpilloverTable(spilloverTable(varEstimate(simulateVAR(est, boot_sample_length, burnout), est.lags, est.typ), H; fevd_function = fevd_function, nocorr = nocorr)), 1:reps)
	elseif typeof(est)==VARmodels.varEstimateShrink
		spillTableBooted = map(i->SpilloverTable(spilloverTable(varEstimateShrink(simulateVAR(est, boot_sample_length, burnout), est.lags, est.typ, est.λ), H; fevd_function = fevd_function, nocorr = nocorr)), 1:reps)
	elseif typeof(est)==VARmodels.varRepresentationVECM
		spillTableBooted = map(i->SpilloverTable(spilloverTable(varRepresentationVECM(CointegrationEstimate(simulateVAR(est, boot_sample_length, burnout), "none", est.lags, "longrun"; season = "no", dumvar = "no"), est.r), H; fevd_function = fevd_function, nocorr = nocorr)), 1:reps)
	else
		error("Unsupported type of estimate for bootstrap.")
	end

	down_ind = convert(Int, floor(reps*quantiles[1]))
	up_ind = convert(Int, floor(reps*quantiles[2]))

	(mean, up, down) = get_quantiles_from_SpilloverTable(spillTableBooted, down_ind, up_ind)

	SpilloverTableBooted(spillsEstimate, mean, up, down, quantiles)
end


type SpilloverTableFrequencyBooted
	estimate::SpilloverTableFrequency
	mean::SpilloverTableFrequency
	up::SpilloverTableFrequency
	down::SpilloverTableFrequency
	quantiles::Vector{Float64}
end


function show(io::IO, tab::SpilloverTableFrequencyBooted)
	print("The object contains bootstrapped confidence intervals not shown in here.\n")
	printTable(tab.estimate)
end

function SpilloverTableFrequencyBooted(est::VARmodels.VARRepresentation, H, reps; bounds = nothing, fevd_function = fftGenFEVD, nocorr = false, burnout = 100, quantiles = [0.05, 0.95], boot_sample_length = 0)
	# Compute spillovers
	if boot_sample_length==0
		boot_sample_length=est.obs
	end

	spillsEstimate = SpilloverTableFrequency(spilloverTable(est::VARmodels.VARRepresentation, H; fevd_function = fevd_function, bounds = bounds, nocorr = nocorr), bounds)

	if typeof(est)==VARmodels.varEstimate
		spillTableBooted = map(i->SpilloverTableFrequency(spilloverTable(varEstimate(simulateVAR(est, boot_sample_length, burnout), est.lags, est.typ), H; fevd_function = fevd_function, bounds = bounds, nocorr = nocorr), bounds), 1:reps)
	elseif typeof(est)==VARmodels.varEstimateShrink
		spillTableBooted = map(i->SpilloverTableFrequency(spilloverTable(varEstimateShrink(simulateVAR(est, boot_sample_length, burnout), est.lags, est.typ, est.λ), H; fevd_function = fevd_function, bounds = bounds, nocorr = nocorr), bounds), 1:reps)
	elseif typeof(est)==VARmodels.varRepresentationVECM
		spillTableBooted = map(i->SpilloverTableFrequency(spilloverTable(varRepresentationVECM(CointegrationEstimate(simulateVAR(est, boot_sample_length, burnout), "none", est.lags, "longrun"; season = "no", dumvar = "no"), est.r), H; fevd_function = fevd_function, bounds = bounds, nocorr = nocorr), bounds), 1:reps)
	else
		error("Unsupported type of estimate for bootstrap.")
	end

	down_ind = convert(Int, floor(reps*quantiles[1]))
	up_ind = convert(Int, floor(reps*quantiles[2]))

	tabs_within = mapslices(q->get_quantiles_from_SpilloverTable(q, down_ind, up_ind), [i.tables_within[j] for i=spillTableBooted, j=1:(size(bounds)[1]-1)], (1))
	tabs_abs = mapslices(q->get_quantiles_from_SpilloverTable(q, down_ind, up_ind), [i.tables_absolute[j] for i=spillTableBooted, j=1:(size(bounds)[1]-1)], (1))

	SpilloverTableFrequencyBooted(spillsEstimate, SpilloverTableFrequency(spillsEstimate.bands, vec([i[1] for i=tabs_within]), vec([i[1] for i=tabs_abs])), SpilloverTableFrequency(spillsEstimate.bands, vec([i[2] for i=tabs_within]), vec([i[2] for i=tabs_abs])), SpilloverTableFrequency(spillsEstimate.bands, vec([i[3] for i=tabs_within]), vec([i[3] for i=tabs_abs])), quantiles)
end

# This function takes a list of objects SpilloverTable and selects [ind_down, ind_up] from ordered entries.
# So, if I want a 95% quantile, I just need to select properly the ind_up.
#
# This is only a helper function all the machinery is taken care of by the constructor functions of the ...Booted objects.

function get_quantiles_from_SpilloverTable(stl, ind_down, ind_up)
	names = stl[1].names
	reps = size(stl)[1]
	# table
	bb = hcat(map((i)->vec(stl[i].table), 1:reps)...)
	res = zeros(size(bb)[1],2)
	res_mean = zeros(size(bb)[1],1)
	for i=1:size(bb)[1]
	    res[i, :] = sort(bb[i,:])[[ind_down, ind_up]]
			res_mean[i, 1] = mean(bb[i,:])
	end
	table_down=reshape(res[:,1], size(stl[1].table))
	table_up=reshape(res[:,2], size(stl[1].table))
	table_mean=reshape(res_mean[:,1], size(stl[1].table))

	# pairwise
	bb = hcat(map((i)->vec(stl[i].pairwise), 1:reps)...)
	res = zeros(size(bb)[1],2)
	res_mean = zeros(size(bb)[1],1)
	for i=1:size(bb)[1]
	    res[i, :] = sort(bb[i,:])[[ind_down, ind_up]]
			res_mean[i, 1] = mean(bb[i,:])
	end
	pairwise_down=reshape(res[:,1], size(stl[1].pairwise))
	pairwise_up=reshape(res[:,2], size(stl[1].pairwise))
	pairwise_mean=reshape(res_mean[:,1], size(stl[1].pairwise))

	# from
	bb = hcat(map((i)->stl[i].from, 1:reps)...)
	res = zeros(size(bb)[1],2)
	res_mean = zeros(size(bb)[1],1)
	for i=1:size(bb)[1]
	    res[i, :] = sort(bb[i,:])[[ind_down, ind_up]]
			res_mean[i, 1] = mean(bb[i,:])
	end
	from_down=reshape(res[:,1], size(stl[1].from))
	from_up=reshape(res[:,2], size(stl[1].from))
	from_mean=reshape(res_mean[:,1], size(stl[1].from))

	# to
	bb = hcat(map((i)->stl[i].to, 1:reps)...)
	res = zeros(size(bb)[1],2)
	res_mean = zeros(size(bb)[1],1)
	for i=1:size(bb)[1]
	    res[i, :] = sort(bb[i,:])[[ind_down, ind_up]]
			res_mean[i, 1] = mean(bb[i,:])
	end
	to_down=reshape(res[:,1], size(stl[1].to))
	to_up=reshape(res[:,2], size(stl[1].to))
	to_mean=reshape(res_mean[:,1], size(stl[1].to))

	# net
	bb = hcat(map((i)->stl[i].net, 1:reps)...)
	res = zeros(size(bb)[1],2)
	res_mean = zeros(size(bb)[1],1)
	for i=1:size(bb)[1]
	    res[i, :] = sort(bb[i,:])[[ind_down, ind_up]]
			res_mean[i, 1] = mean(bb[i,:])
	end
	net_down=reshape(res[:,1], size(stl[1].net))
	net_up=reshape(res[:,2], size(stl[1].net))
	net_mean=reshape(res_mean[:,1], size(stl[1].net))

	# overall
	bb = hcat(map((i)->stl[i].overall, 1:reps)...)
	res = sort(bb[1,:])[[ind_down, ind_up]]
	res_mean = mean(bb[1,:])
	overall_down=res[1]
	overall_up=res[2]
	overall_mean=res_mean

	return (SpilloverTable(names, table_mean, from_mean, to_mean, net_mean, pairwise_mean, overall_mean), SpilloverTable(names, table_up, from_up, to_up, net_up, pairwise_up, overall_up), SpilloverTable(names, table_down, from_down, to_down, net_down, pairwise_down, overall_down))
end
