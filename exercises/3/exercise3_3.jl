using DynamicalSystems, DifferentialEquations
import PredefinedDynamicalSystems as PDS
using LinearAlgebra, Plots

function lyapunov_discrete(f, x0, δ0, Δt, N)
	xref  = [x0]
	xtest = [x0 .+ δ0]
	δi    = []

	for i=1:(N*Δt)
		push!(xref, f(xref[end]))

		if i % Δt == 0
			push!(xtest, f(xtest[end]))
			push!(δi, norm(xref[end] .- xtest[end]))
			xtest[end] = xref[end] .+ δ0
		else
			push!(xtest, f(xtest[end]))
		end
	end
	sum(log.(δi ./ norm(δ0))) / (N * Δt)
end

function lyapunov_continuous(ds, u0, p, Δt, N, δ0)
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

function standard_map(x)
	θ, v= x
	k  = 0.6
	(θ + v + k*sin(θ), v + k*sin(θ)) .% 2π
end

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

my_λ1_sm = lyapunov_discrete(standard_map, (0.1, 0.11), (0.0001, 0.0001), 50, 10_000)
sm = PDS.standardmap([0.1, 0.11]; k=0.6)
# BUG in DynamicalSystems.jl prevents using lyapunov(), but lyapunovspectrum()[1] works...
λ1_sm = lyapunovspectrum(sm, 10_000)[1]

u0 = [20.0, 20.0, 20.0]
p  = [10.0, 28.0, 8/3]
my_λ1_lo = lyapunov_continuous(lorenz63!, u0, p, 5.0, 1_000, 0.01)
lo = PDS.lorenz([20.0, 20.0, 20.0]; σ=10.0, ρ=28.0, β=8/3)
λ1_lo = lyapunov(lo, 1_000)

println("-----------------+-------+-----------------------")
println("Dynamical System | My λ1 | DynamicalSystems.jl λ1")
println("-----------------+-------+-----------------------")
println("Standard map     | $(round(my_λ1_sm; digits=2))  | $(round(λ1_sm; digits=2))")
println("-----------------+-------+-----------------------")
println("Lorenz 63        | $(round(my_λ1_lo; digits=2))  | $(round(λ1_lo; digits=2))")
println("-----------------+-------+-----------------------")


#=
julia> include("exercise3_3.jl")
-----------------+-------+-----------------------
Dynamical System | My λ1 | DynamicalSystems.jl λ1
-----------------+-------+-----------------------
Standard map     | 0.13  | 0.09
-----------------+-------+-----------------------
Lorenz 63        | 1.03  | 0.9
-----------------+-------+-----------------------
=#
