using DynamicalSystems, LinearAlgebra
function lyapunov_exponents(ds::DynamicalSystem, N, k, Δt)
	integ = tangent_integrators(ds, k)
	t0, λ = integ.t, zeros(k)
	for i in 2:N
		step!(integ, Δt)
		Y = get_deviations(integ)
		qrdec = qr(Y)
		for j in 1:k
			λ[j] += log(abs(qrdec.R[j, j]))
		end
		set_deviations!(integ, qrdec.Q)
	end
	λ ./= (integ.t - t0)
	return λ
end
