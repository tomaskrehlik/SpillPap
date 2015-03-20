function latexSpilloverTable(data, lags, stepsAhead, typ, bounds, nams, caption, label)

	function pad_string(left, str, right)
		return map((x) -> join([left, x, right], ""), str)
	end


	total=map((x)->format(100*x, precision=1), mapslices(sum, fftSpilloverTableDec12(data, lags, stepsAhead, "Const and trend", bounds, false), [3])[:,:,1])
	abs=map((x)->format(100*x, precision=1), fftSpilloverTableDec12(data, lags, stepsAhead, "Const and trend", bounds, false))
	prop=map((x)->format(100*x, precision=1), fftSpilloverTableDec12(data, lags, stepsAhead, "Const and trend", bounds, true))

	total_l = mapslices((x) -> join(x, " & "), total, [2])[:,1]
	abs_l = mapslices((x) -> join(x, " & "), pad_string("\\tiny{(",mapslices((elem)->join(elem, " -- "), abs, [3]), ")}"), [2])[:,1,1]
	prop_l = mapslices((x) -> join(x, " & "), pad_string("\\tiny{(",mapslices((elem)->join(elem, " -- "), prop, [3]), ")}"), [2])[:,1,1]

	first = apply(vcat,[vec([j, " ", " "]) for j=[join(["\\multirow{3}{*}{",i,"}"]) for i=nams]])

	tab = vec([total_l abs_l prop_l]')
	tab = join(mapslices((x)->join(x, " & "), [first tab], [2]), " & \\\\ \n")

	tt_t = format(100*spilloverDiebold12(data, lags, stepsAhead, "Const and trend"), precision=1)
	tt_d = map((x)->format(100*x, precision=1), fftSpilloverDec12(data, lags, stepsAhead, "Const and trend", bounds, false))
	tt_n = map((x)->format(100*x, precision=1), fftSpilloverDec12(data, lags, stepsAhead, "Const and trend", bounds, false, true))

	tt_t_l = join(["   ", ["   " for i=1:size(nams)[1]], tt_t], " & ")
	tt_d_l = join(["   ", ["   " for i=1:size(nams)[1]], join(["\\tiny{(",join(tt_d, " -- "),")}"], "")], " & ")
	tt_n_l = join(["   ", ["   " for i=1:size(nams)[1]], join(["\\tiny{(",join(tt_n, " -- "),")}"], "")], " & ")

	header = join(["\\begin{table}[th]
		\\scriptsize
		\\centering
		\\begin{tabular}{c|",join(["c" for i=1:size(nams)[1]], ""),"|c}
		\\toprule
		  & ", join(nams, " & "), " & Total \\\\
		  \\midrule\n"],"")

	footer = join(["\\\\ \n	\\bottomrule
		\\end{tabular}
		\\caption{", caption, "}
		\\label{", label, "}
	\\end{table}
	"],"")
	return join([header, join([tab, "\\midrule"], " \\\\ \n"), join([tt_t_l, tt_d_l, tt_n_l], " \\\\ \n"), footer], " \n")
end