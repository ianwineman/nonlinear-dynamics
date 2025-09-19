using LinearAlgebra, DifferentialEquations, Plots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

function lyapunov(ds, u0, p, δ0, Δt, N)
	tis = Δt:Δt:(N*Δt)

	xref_prob = ODEProblem(ds, u0, (0.0,N*Δt), p)
	xref_sol = solve(xref_prob, tstops=tis)

	xref_tstops_i = [findfirst(x->x==ti, xref_sol.t) for ti in Δt:Δt:N*Δt]
	xrefs = getindex(xref_sol.u, xref_tstops_i)

	xtest_prob = ODEProblem(ds, u0 .+ δ0, (0.0,Δt), p)
	integrator = init(xtest_prob, Tsit5())

	for i=1:N
		solve!(integrator)
		reinit!(integrator, xrefs[i] .+ δ0; erase_sol=false, t0=i*Δt, tf=(i+1.0)*Δt)
	end

	xtest_sol = integrator.sol
	xtest_tstops_i = [findfirst(x->x==ti, xtest_sol.t) for ti in Δt:Δt:N*Δt]
	xtests = getindex(xtest_sol.u, xtest_tstops_i)

	δis = norm.(xrefs .- xtests)

	return sum(log.(δis ./ δ0)) / (N * Δt)
end

u0 = [20.0, 20.0, 20.0]
p = [10.0, 28.0, 8/3]
λ1 = lyapunov(lorenz63!, u0, p, 0.000001, 1.0, 45)

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
	4:19, 0.91 .* collect(0:15) .- 13,
	label="λ1=$(round(λ1, digits=2))", 
	linecolor=:mediumpurple
)

