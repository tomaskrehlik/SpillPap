module SpillPap

	using Formatting
	using VARmodels
	using Distributions

	# Exports for the simulations of the Bi-variate two component model
	export q,
		s,
		simBivar

	# Exports of the standard spillovers from Diebold & Yilmaz
	export spilloverDiebold09,
		spilloverDiebold12
		
	# Exports of the Fourier decomposed spillovers.
	export fftSpillover09,
		fftSpillover12,
		fftSpilloverDec09,
		fftSpilloverDec12,
		fftSpilloverTableDec09,
		fftSpilloverTableDec12

	export spilloverTable,
		spilloverMeasure,
		bootSpillovers

	# Exports for the nullification of correlation
	export 	spilloverDiebold12NoCorr

	# Output spillovers.
	export latexSpilloverTable,
		latexSpilloverTablePresent

	include("biVariateGarch.jl")
	include("oldSpillovers.jl")
	include("fftSpillover.jl")
	include("spilloverTable.jl")
	include("spillovers.jl")

end # module
