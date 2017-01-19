"""

The `SpilloverTable` function creates an object that holds the whole spillover table. That means 
it computes all measures such as from, to, net, and pairwise spillovers. The access to the individual
spillovers are intuitive.

	- `table` a fevd table produced by function spilloverTable
	- `names` a vector of names of the series, can be left empty, then prints names as "V_",
	where _ is a generic number

The object `SpilloverTable` has the followind attributes:

	- `table`
	- `to`
	- `from`
	- `net`
	- `pairwise`
	- `overall`

These correspond to individual spillovers.

"""
type SpilloverTable
	names::Vector{String}
	table::Matrix{Float64}
	from::Vector{Float64}
	to::Vector{Float64}
	net::Vector{Float64}
	pairwise::Matrix{Float64}
	overall::Float64
end

function SpilloverTable(table)
	if ((typeof(table) <: Array{Float64, 3}) && (size(table)[3] == 1))
		table = table[:,:,1]
	end
	SpilloverTable(table, ["V" * j for j=[sprintf1("%d", i) for i=1:size(table)[1]]])
end

function SpilloverTable(table, names)
	if ((typeof(table) <: Array{Float64, 3}) && (size(table)[3] == 1))
		table = table[:,:,1]
	end

	@assert typeof(table) <: Array{Float64, 2} "The argument should be a spillover table produced by function spilloverTable."
	@assert typeof(names) <: Vector{String} "The names argument should be a vector of names of the series."
	@assert length(names)==size(table)[1] "The names vector should have the same number of entries as there is of series in the system."

	overall = 1 - sum(diag(table))/sum(table)
	from = (vec(sum(table, 1)) - diag(table))/sum(table)
	to = (vec(sum(table, 2)) - diag(table))/sum(table)
	net = from-to
	pairwise = (triu(table) - tril(table)')
	SpilloverTable(names, table, from, to, net, pairwise, overall)
end


function show(io::IO, tab::SpilloverTable)
	printTable(tab, tab.names)
end


function printTable(tab::SpilloverTable, names::Vector{String})
	# Set the variables
	fmt = "%02.1f"
	n = size(tab.table)[1]

	function printShorter(s::String, max_l::Int)
		if length(s)<max_l
			print(repeat(" ", max_l-length(s)))
		end
		print(s[1:min(max_l,length(s))])
	end

	function padLeft(s::String, max_l::Int)
		if length(s)==max_l
			print(s)
		else 
			print(repeat(" ", max_l-length(s)))
			print(s)
		end
	end

	# Spillover table
	# Print the first row with names
	print(repeat(" ", 6))
	print("| ")
	for i=names
		printShorter(i, 6)
		print(" ")
	end
	print("| FROM")
	print("\n")

	# Print the division row
	print(repeat("-", 7*n + 6*2 + 3))
	print("\n")

	# Print the body of the table
	for i=1:n
		printShorter(names[i], 6)
		print("| ")
		for j=1:n
			s = sprintf1( fmt, tab.table[i,j] )
			padLeft(s, 6)
			print(" ")
			
		end
		print("|")
		s = sprintf1( fmt, 100*tab.from[i] )
		padLeft(s, 6)
		print("\n")
	end

	# Print the division row
	print(repeat("-", (6+1)*n + 6*2 + 3))
	print("\n")

	# Print the TO row + overall
	print("  TO  | ")
	for j=1:n
		s = sprintf1( fmt, 100*tab.to[j] )
		padLeft(s, 6)
		print(" ")
	end
	print("|")
	s = sprintf1( fmt, 100*tab.overall )
		padLeft(s, 6)
	print("\n")

	# Print the NET row
	print("  NET | ")
	for j=1:n
		s = sprintf1( fmt, 100*tab.net[j] )
		padLeft(s, 6)
		print(" ")
	end
	print("| \n\n")


	# The table with pairwise spillovers
	print("PAIRWISE\n")
	print(repeat(" ", 6))
	print("| ")
	for i=names
		printShorter(i, 6)
		print(" ")
	end
	print("\n")
	print(repeat("-", (6+1)*n + 6 + 3))
	print("\n")
	for i=1:n
		printShorter(names[i], 6)
		print("| ")
		for j=1:size(tab.pairwise)[2]
			s = sprintf1( fmt, tab.pairwise[i,j] )
			padLeft(s, 6)
			print(" ")
			
		end
		print("\n")
	end
end


