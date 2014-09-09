module SpillPap

export 
	q,
	s,
	simBivar,
	spillover,
	generalisedSpillover,
	generalisedWavSpillover,
	testNonOverlapping,
	generalisedWavDecompSpillover,
	VAREST,
	fevd,
	getEnergies,
	Psi,
	restrictVAR,
	restrictVAR2,
	spilloverFinal,
	Phi

	include("VAR.jl")
	include("biVariateGarch.jl")
	include("spillovers.jl")


end # module
