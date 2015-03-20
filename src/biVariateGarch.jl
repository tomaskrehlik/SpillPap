using Distributions

function q(ω, ρ, φ, ql, sl, ε)
	return convert(Array{Float64,1}, ω + ρ*ql + φ.*(ε.^2 - ql - sl))
end

function s(α, β, ql, sl, ε)
	A = [α[1] 0.00;
		 0.00 α[2]]
	return convert(Array{Float64,1}, (A + β)*sl + α.*(ε.^2 - ql - sl))
end


function simBivar(A, B, C, O, P, R)
    return simBivar(A, B, C, O, P, R, 1000)
end

function simBivar(A, B, C, O, P, R, len)

    E = rand(MvNormal(C), len+1)

	qout = q(O, R, P, [0.0, 0.0], [0.0, 0.0], E[:,1])
	sout = s(A, B, [0.0, 0.0], [0.0, 0.0], E[:,1])

	qout = [qout q(O, R, P, qout, sout, E[:,1])]
	sout = [sout s(A, B, qout[:,1], sout, E[:,1])]

	for i = 2:len
		qout = [qout q(O, R, P, qout[:,i], sout[:,i], E[:,i+1])]
		sout = [sout s(A, B, qout[:,i], sout[:,i], E[:,i+1])]
	end

    return (qout, sout)
end