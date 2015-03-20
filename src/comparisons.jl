using Gadfly
using SpillPap
using Wavelets
using DataFrames

cd("/Users/tomaskrehlik/Documents/PHD/Spillovers_FIVAR/simulation study/images/")

set_default_plot_size(20cm, 20cm)

burnout = 100
sims=100

A = [0.01, 0.01]

B = [0.75 0.00;
	 0.00 0.75]

C = [1.00 0.00;
	 0.00 1.00]

R = [0.80 0.19;
	 0.19 0.80]

P = [0.01, 0.01]

O = [0.01, 0.03]

levels = 5
filt = ["haar", "la8", "d8"]
out = [begin
	dif = [begin
		(qdata, sdata) = simBivar(A, B, C, O, P, R);
		(qdata, sdata) = (qdata[:,burnout:end], sdata[:,burnout:end]);

		dat = modwt((qdata+sdata)', "la8", levels, "periodic");
		
		(D, S) = (dat.W, dat.V);
		decomMODWT = [apply(hcat, [D[:,j,i] for i=1:2, j=1:levels]) apply(hcat, [S[:,levels,i] for i=1:2])];
		(D, S) = mra(dat);
		decomMRA = [apply(hcat, [D[:,j,i] for i=1:2, j=1:levels]) apply(hcat, [S[:,levels,i] for i=1:2])];

		spilDecMRA = generalisedWavDecompSpillover((qdata+sdata)', levels, 1, 11, "la8", "mra");
		spilDecMODWT = generalisedWavDecompSpillover((qdata+sdata)', levels, 1, 11, "la8", "modwt");


		weights = [[var(D[:,i,1]) for i=1:levels], var(S[:,levels,1])]./var((qdata+sdata)[1,:]);
		weightedspilMRA = weights.*[spilDecMRA[i][2] for i=1:(levels+1)];
		weightedspilMODWT = weights.*[spilDecMODWT[i][2] for i=1:(levels+1)];
		[sum(weightedspilMODWT), sum(weightedspilMRA), generalisedSpillover(VAREST(decomMODWT, 1, "Const"), 11)[2], generalisedSpillover(VAREST(decomMRA, 1, "Const"), 11)[2], generalisedSpillover(VAREST((qdata+sdata)', 1, "Const"),11)[2]]
	end for i=1:100]
	mapslices(mean, apply(hcat, dif), 2)
end for j = filt]


level = 2
a=[begin


(qdata, sdata) = simBivar(A, B, C, O, P, R, 1000);
(qdata, sdata) = (qdata[:,burnout:end], sdata[:,burnout:end]);

dat = modwt((qdata+sdata)', "la8", level, "periodic");
		
(D, S) = (dat.W, dat.V);
decomMODWT = [apply(hcat, [D[:,j,i] for i=1:2, j=1:level]) apply(hcat, [S[:,level,i] for i=1:2])];
(D, S) = mra(dat);
decomMRA = [apply(hcat, [D[:,j,i] for i=1:2, j=1:level]) apply(hcat, [S[:,level,i] for i=1:2])];
weights = [[var(D[:,i,1]) for i=1:level], var(S[:,level,1])]./var((qdata+sdata)[1,:]);

spilDecMODWT = generalisedWavDecompSpillover((qdata+sdata)', level, 1, 11, "la8", "mra")

generalisedWavSpillover(VAREST(decomMRA, 1, "Const"), 11, level)[2] - sum(weights.*[spilDecMODWT[i][2] for i=1:(level+1)])
end for i=1:400]

plot(y=a, x=[1:400], Geom.point)