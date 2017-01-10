function latexSpilloverTable(data, lags, stepsAhead, typ, bounds, nams, caption = "", label = "")

	function pad_string(left, str, right)
		return map((x) -> join([left, x, right], ""), str)
	end

	total = mapslices(sum, fftSpilloverTableDec12(data, lags, stepsAhead, typ, bounds, false), [3])[:,:,1]
	total_ll = vec(sum(total,[1])) - diag(total)
	total = [total sum(total,[2])-diag(total)]
	abs = fftSpilloverTableDec12(data, lags, stepsAhead, typ, bounds, false)
	abs_ll = mapslices(x ->  vec(sum(x,[1]))-diag(x), abs, [1,2])
	abs = mapslices(x -> [x sum(x,[2])-diag(x)], abs, [1,2])
#==#
	prop = fftSpilloverTableDec12(data, lags, stepsAhead, typ, bounds, true)
	prop_ll = mapslices(x ->  vec(sum(x,[1]))-diag(x), prop, [1,2])
	prop = mapslices(x -> [x sum(x,[2])-diag(x)], prop, [1,2])

	total=map((x)->format(100*x, precision=1), total)
	abs=map((x)->format(100*x, precision=1), abs)
	prop=map((x)->format(100*x, precision=1), prop)
	total_ll=map((x)->format(100*x, precision=1), total_ll)
	abs_ll=map((x)->format(100*x, precision=1), abs_ll)
	prop_ll=map((x)->format(100*x, precision=1), prop_ll)


	total_l = mapslices((x) -> join(x, " & "), total, [2])[:,1]
	abs_l = mapslices((x) -> join(x, " & "), pad_string("\\tiny{(",mapslices((elem)->join(elem, " -- "), abs, [3]), ")}"), [2])[:,1,1]
	prop_l = mapslices((x) -> join(x, " & "), pad_string("\\tiny{(",mapslices((elem)->join(elem, " -- "), prop, [3]), ")}"), [2])[:,1,1]

	first = apply(vcat,[vec([j, " ", " "]) for j=[join(["\\multirow{3}{*}{",i,"}"]) for i=nams]])

	tab = vec([total_l abs_l prop_l]')
	tab = join(mapslices((x)->join(x, " & "), [first tab], [2]), "\\\\ \n")

	tt_t = format(100*sum(fftSpilloverDec12(data, lags, stepsAhead, typ, bounds, false)), precision=1)
	tt_d = map((x)->format(100*x, precision=1), fftSpilloverDec12(data, lags, stepsAhead, typ, bounds, false))
	tt_n = map((x)->format(100*x, precision=1), fftSpilloverDec12(data, lags, stepsAhead, typ, bounds, false, true))

	tt_t_l = join(["\\multirow{3}{*}{To}", total_ll, tt_t], " & ")
	tt_d_l = join(["", vec(mapslices(x->join(["\\tiny{(",join(x, "--"),")}"],""), abs_ll, [3])), join(["\\tiny{(",join(tt_d, " -- "),")}"], "")], " & ")
	tt_n_l = join(["", vec(mapslices(x->join(["\\tiny{(",join(x, "--"),")}"],""), prop_ll, [3])), join(["\\tiny{(",join(tt_n, " -- "),")}"], "")], " & ")

	header = join(["\\begin{table}[th]
		\\scriptsize
		\\centering
		\\begin{tabular}{c|",join(["c" for i=1:size(nams)[1]], ""),"|c}
		\\toprule
		  & ", join(nams, " & "), " & From \\\\
		  \\midrule\n"],"")

	footer = join(["\\\\ \n	\\bottomrule
		\\end{tabular}
		\\caption{", caption, "}
		\\label{", label, "}
	\\end{table}
	"],"")
	return join([header, join([tab, "\\midrule"], " \\\\ \n"), join([tt_t_l, tt_d_l, tt_n_l], " \\\\ \n"), footer], " \n")
end

function latexSpilloverTablePresent(data, lags, stepsAhead, typ, bounds, nams, caption = "", label = "")

	function pad_string(left, str, right)
		return map((x) -> join([left, x, right], ""), str)
	end

	total = mapslices(sum, fftSpilloverTableDec12(data, lags, stepsAhead, typ, bounds, false), [3])[:,:,1]
	total_ll = vec(sum(total,[1])) - diag(total)
	total = [total sum(total,[2])-diag(total)]
	abs = fftSpilloverTableDec12(data, lags, stepsAhead, typ, bounds, false)
	abs_ll = mapslices(x ->  vec(sum(x,[1]))-diag(x), abs, [1,2])
	abs = mapslices(x -> [x sum(x,[2])-diag(x)], abs, [1,2])
#==#
	prop = fftSpilloverTableDec12(data, lags, stepsAhead, typ, bounds, true)
	prop_ll = mapslices(x ->  vec(sum(x,[1]))-diag(x), prop, [1,2])
	prop = mapslices(x -> [x sum(x,[2])-diag(x)], prop, [1,2])

	total=map((x)->format(100*x, precision=0), total)
	abs=map((x)->format(100*x, precision=0), abs)
	prop=map((x)->format(100*x, precision=0), prop)
	total_ll=map((x)->format(100*x, precision=0), total_ll)
	abs_ll=map((x)->format(100*x, precision=0), abs_ll)
	prop_ll=map((x)->format(100*x, precision=0), prop_ll)


	total_l = mapslices((x) -> join(x, " & "), total, [2])[:,1]
	abs_l = mapslices((x) -> join(x, " & "), pad_string("\\tiny{(",mapslices((elem)->join(elem, "|"), abs, [3]), ")}"), [2])[:,1,1]
	prop_l = mapslices((x) -> join(x, " & "), pad_string("\\tiny{(",mapslices((elem)->join(elem, "|"), prop, [3]), ")}"), [2])[:,1,1]

	first = apply(vcat,[vec([j, " ", " "]) for j=[join(["\\multirow{3}{*}{",i,"}"]) for i=nams]])

	tab = vec([total_l abs_l prop_l]')
	tab = join(mapslices((x)->join(x, " & "), [first tab], [2]), "\\\\ \n")

	tt_t = format(100*sum(fftSpilloverDec12(data, lags, stepsAhead, typ, bounds, false)), precision=0)
	tt_d = map((x)->format(100*x, precision=0), fftSpilloverDec12(data, lags, stepsAhead, typ, bounds, false))
	tt_n = map((x)->format(100*x, precision=0), fftSpilloverDec12(data, lags, stepsAhead, typ, bounds, false, true))

	tt_t_l = join(["\\multirow{3}{*}{To}", total_ll, tt_t], " & ")
	tt_d_l = join(["", vec(mapslices(x->join(["\\tiny{(",join(x, "|"),")}"],""), abs_ll, [3])), join(["\\tiny{(",join(tt_d, "|"),")}"], "")], " & ")
	tt_n_l = join(["", vec(mapslices(x->join(["\\tiny{(",join(x, "|"),")}"],""), prop_ll, [3])), join(["\\tiny{(",join(tt_n, "|"),")}"], "")], " & ")

	header = join(["\\begin{table}[th]
		\\tiny
		\\centering
		\\begin{tabular}{c|",join(["c" for i=1:size(nams)[1]], ""),"|c}
		\\toprule
		  & ", join(nams, " & "), " & From \\\\
		  \\midrule\n"],"")

	footer = join(["\\\\ \n	\\bottomrule
		\\end{tabular}
		\\caption{", caption, "}
		\\label{", label, "}
	\\end{table}
	"],"")
	return join([header, join([tab, "\\midrule"], " \\\\ \n"), join([tt_t_l, tt_d_l, tt_n_l], " \\\\ \n"), footer], " \n")
end