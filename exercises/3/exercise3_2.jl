using LinearAlgebra, DifferentialEquations, Plots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

function lyapunov(ds, u0, p, Δt, N, δ0)
	xrefr_prob = ODEProblem(ds, u0,            (0.0, Δt*N), p)
	xtest_prob = ODEProblem(ds, u0 .+ (δ0/√3), (0.0, Δt*N), p)

	xrefr_int = init(xrefr_prob, Tsit5())
	xtest_int = init(xtest_prob, Tsit5())

	δi = []

	for _ in 1:N
		step!(xrefr_int, Δt, true)
		step!(xtest_int, Δt, true)
		push!(δi, norm(xrefr_int.u .- xtest_int.u))
		reinit!(xtest_int, xrefr_int.u .+ (xtest_int.u .- xrefr_int.u) .* (δ0/δi[end]))
	end

	λ1 = (1/(N*Δt)) * sum(log.(δi ./ δ0))
end

u0 = [20.0, 20.0, 20.0]
p  = [10.0, 28.0, 8/3]
λ1 = lyapunov(lorenz63!, u0, p, 5.0, 1_000, 0.01)

prob1 = ODEProblem(lorenz63!, u0, (0.0,45.0), p)
sol1 = solve(prob1)

prob2 = ODEProblem(lorenz63!, u0 .+ 0.000001, (0.0,45.0), p)
sol2 = solve(prob2)

l = minimum([length(sol1.t), length(sol2.t)])

plot(
	sol1.t[1:l],
	log.(norm.(sol1.u[1:l] .- sol2.u[1:l])),
	title="Largest Lyapunov Exponent",
	xlabel="t", label="ln(δ(t))",
	linecolor=:deepskyblue
)
plot!(
	4:18, λ1 .* collect(0:14) .- 13,
	label="λ1=$(round(λ1, digits=2))", 
	linecolor=:mediumpurple
)
savefig("plots/exercise3_2.png")
