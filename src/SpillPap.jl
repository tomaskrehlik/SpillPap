module SpillPap

	using Formatting
	using Wavelets
	using VARmodels

	# VAR exports that will be removed once I move the estimation infrastructure to the VAR
	# models
	# export 	VAREST,
	# 	fevd,
	# 	genFEVD,
	# 	getEnergies,
	# 	Psi,
	# 	restrictVAR,
	# 	restrictVAR2,
	# 	Phi,
	# 	testNonOverlapping

	# Exports for the simulations of the Bi-variate two component model
	export q,
		s,
		simBivar

	# Exports of wavelet spillovers, probably wrong estimation!
	export generalisedWavSpillover,
		generalisedWavDecompSpillover,
		spilloverOurs,
		spilloverOursRolling

	# Exports of the standatd spillovers from Diebold & Yilmaz
	export spilloverDiebold09,
		spilloverDiebold09Rolling,
		spilloverDiebold12,
		spilloverDiebold12Rolling
		
	# Exports of the Fourier decomposed spillovers.
	export fftFEVD, 
		fftGenFEVD,
		fftSpillover09,
		fftSpillover12,
		fftSpilloverDec09,
		fftSpilloverDec12,
		fftSpilloverTableDec09,
		fftSpilloverTableDec12

	# Exports for the nullification of correlation
	export 	spilloverDiebold12NoCorr

	# Output spillovers.
	export latexSpilloverTable


	# include("VAR.jl")
	include("biVariateGarch.jl")
	include("spillovers.jl")
	include("oldSpillovers.jl")
	include("newSpillovers.jl")
	include("fftSpillover.jl")
	include("spilloverTable.jl")

end # module
