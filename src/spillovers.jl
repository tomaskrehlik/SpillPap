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

