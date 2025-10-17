using LinearAlgebra, DifferentialEquations

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

function lyapunov(ds, u0, p, Δt, N, δ0)
	xrefr_prob = ODEProblem(ds, u0,                          (0.0, N*Δt), p)
	xtest_prob = ODEProblem(ds, u0 .+ (δ0/sqrt(length(u0))), (0.0, N*Δt), p)

	xrefr_int = init(xrefr_prob, Tsit5())
	xtest_int = init(xtest_prob, Tsit5())

	δ = zeros(N)
	for i in 1:N
		step!(xrefr_int, Δt, true)
		step!(xtest_int, Δt, true)
		δ[i] += norm(xrefr_int.u .- xtest_int.u)
		reinit!(xtest_int, xrefr_int.u .+ (xtest_int.u .- xrefr_int.u) .* (δ0/δ[i]))
	end

	λ = sum(log.(δ ./ δ0))
	λ /= (N*Δt)
end

u0 = [20.0, 20.0, 20.0]
p  = [10.0, 28.0, 8/3]
λ1 = lyapunov(lorenz63!, u0, p, 1.0, 1_000, 0.1)
