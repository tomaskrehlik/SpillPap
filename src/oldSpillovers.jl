function spilloverDiebold09(data, lags, H, typ)
	estimate = varEstimate(data, lags, typ)
	FEVD = fevd(estimate, H)
	spill = 1-sum(diag(FEVD))/estimate.vars
	return spill::Float64
end

function spilloverDiebold09Rolling(data, lags, H, typ, window)
	return [spilloverDiebold09(data[(1:window)+i,:], lags, H, typ) for i=0:(size(data)[1]-window)]
end

function spilloverDiebold12(data, lags, H, typ)
	estimate = varEstimate(data, lags, typ)
	GFEVD = genFEVD(estimate, H)
	spill = 1-sum(diag(GFEVD))/estimate.vars
	return spill::Float64
end

function spilloverDiebold12Rolling(data, lags, H, typ, window)
	return [spilloverDiebold12(data[(1:window)+i,:], lags, H, typ) for i=0:(size(data)[1]-window)]
end

function spilloverDiebold12NoCorr(data, lags, H, typ)
	estimate = varEstimate(data, lags, typ)
	GFEVD = genFEVD(estimate, H, true)
	spill = 1-sum(diag(GFEVD))/estimate.vars
	return spill::Float64
end