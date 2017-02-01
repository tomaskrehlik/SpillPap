"""

The `SpilloverTableFrequency` function creates an object of the same name. It takes in the fevd tables
computes all the corresponding measures, both within and absolute variants as described in the Barunik, Krehlik (2017)
The arguments are:

	- `table` a fevd frequency based table produced by the function `spilloverTable`
	- `bands` a vector of frequency bounds.
	- `names` a vector of names that should have the same number of entries as there is series. Can be left empty.

The object `SpilloverTableFrequency` has the following components:

	- `bands` a vector of frequency bounds.
	- `tables_within` an object of type `SpilloverTable` that has the within measures inside
	- `tables_absolute` and object of type `SpilloverTable` that has the absolute measures inside.
	
The access to the individual spillovers is intuitive.

"""

type SpilloverTableFrequency
	bands::Vector{Float64}
	tables_within::Vector{SpilloverTable}
	tables_absolute::Vector{SpilloverTable}
end

function show(io::IO, tab::SpilloverTableFrequency)
	n = size(tab.bands)[1]
	print("WITHIN SPILLOVERS: \n")
	for i=1:(n-1)
		print("Band: " * string(tab.bands[i]) * " to " * string(tab.bands[i+1]) * "\n")
		print(tab.tables_within[i]) # This is a printing function for SpilloverTable, see SpilloverTable.jl for source.
	end
	print("\n\n")
	print("ABSOLUTE SPILLOVERS: \n")
	for i=1:(n-1)
		print("Band: " * string(tab.bands[i]) * " to " * string(tab.bands[i+1]) * "\n")
		print(tab.tables_absolute[i]) # This is a printing function for SpilloverTable, see SpilloverTable.jl for source.
	end
end

function SpilloverTableFrequency(table, bands)
	@assert typeof(table) <: Array{Float64, 3} "The argument should be a spillover table produced by function spilloverTable."
	SpilloverTableFrequency(table, bands, ["v" * j for j=[sprintf1("%d", i) for i=1:size(table)[1]]])
end

function SpilloverTableFrequency(table, bands::Vector{Float64}, names::Vector{String})
	@assert typeof(table) <: Array{Float64, 3} "The argument should be a spillover table produced by function spilloverTable."
	@assert size(table)[3] > 1 "There should only be one frequency band, for frequency based measures, use SpilloverTableFrequency"

	output_within = Vector{SpilloverTable}(size(table)[3])

	for i = 1:size(table)[3]
		output_within[i] = SpilloverTable(table[:,:,i], names)
	end

	output_absolute = deepcopy(output_within)

	for i = 1:size(table)[3]
		variance_share = sum(output_within[i].table)/(100*size(table)[2])
		output_absolute[i].overall *= variance_share
		output_absolute[i].to *= variance_share
		output_absolute[i].from *= variance_share
		output_absolute[i].pairwise *= variance_share
		output_absolute[i].net *= variance_share
	end
	SpilloverTableFrequency(bands, output_within, output_absolute)
end


