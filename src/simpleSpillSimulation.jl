# Non-decomposed simulation of the bi-variate GARCH

burnout = 100
sims=1000

# Simulation without any interconnections
withoutSpills = [begin
					(qdata, sdata) = simBivar(A, B, C, O, P, R)

					(qdata, sdata) = (qdata[:,burnout:end], sdata[:,burnout:end])

					spilloverFinal(VAREST((qdata+sdata)', 3, "Const"), 10)[2]
				end for i=1:sims]

(mean(withoutSpills), std(withoutSpills))

# Simulation with covariance in the residuals
covarianceSpills = [begin

					C = [1.00 j;
	 					j 1.00]

	 				[begin
						(qdata, sdata) = simBivar(A, B, C, O, P, R)

						(qdata, sdata) = (qdata[:,burnout:end], sdata[:,burnout:end])

						generalisedSpillover(VAREST((qdata+sdata)', 3, "Const"), 10)[2]
					end for i=1:sims]	
	
					end for j=[0, 0.3, 0.6, 0.9]]

[map(mean, covarianceSpills) map(std, covarianceSpills)]

C = [1.00 0.00;
	 0.00 1.00]

# Simulation with short term interconnections
shortTermSpills = [begin

					B = [0.75-j j;
	 					j 0.75-j]

	 				[begin
						(qdata, sdata) = simBivar(A, B, C, O, P, R)

						(qdata, sdata) = (qdata[:,burnout:end], sdata[:,burnout:end])

						generalisedSpillover(VAREST((qdata+sdata)', 3, "Const"), 10)[2]
					end for i=1:sims]	
	
					end for j=[0, 0.3, 0.6]]

[map(mean, shortTermSpills) map(std, shortTermSpills)]

B = [0.75 0.00;
	 0.00 0.75]

# Simulation with long term interconnections
longTermSpills = [begin

					R = [0.98-j j;
	 					j 0.98-j]

	 				[begin
						(qdata, sdata) = simBivar(A, B, C, O, P, R)

						(qdata, sdata) = (qdata[:,burnout:end], sdata[:,burnout:end])

						generalisedSpillover(VAREST((qdata+sdata)', 3, "Const"), 10)[2]
					end for i=1:sims]	
	
					end for j=[0, 0.3, 0.6, 0.9]]

[map(mean, longTermSpills) map(std, longTermSpills)]