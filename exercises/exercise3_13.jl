using LinearAlgebra, DifferentialEquations, Plots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

function lyapunov_convergance(ds, u0, p, Δt, N, δ0)
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

	λN = [sum(log.(δ[1:i] ./ δ0))/(i*Δt) for i in 1:length(δ)]
end

let u0 = [20.0, 20.0, 20.0], p  = [10.0, 28.0, 8/3]
	λN = lyapunov_convergance(lorenz63!, u0, p, 1.0, 1_000, 0.1)

	plot(
		1:length(λN), λN,
		xlabel="Iterations, \$N\$", ylabel="\$λ_1\$",
		label="\$λ_1(N)\$",
		linecolor=:orange,
		title="Converagance of \$λ_1\$"
	)
	scatter!(
		[length(λN)], [λN[end]], 
		label="\$λ_1 = $(round(λN[end], digits=4))\$",
		markercolor=:red
	)
	savefig("plots/exercise3_13_1.png")
end

let Δts = [0.1, 0.5, 1.0, 5.0, 10.0]
	λNs = []
	for Δt in Δts
		u0 = [20.0, 20.0, 20.0]
		p  = [10.0, 28.0, 8/3]
		λN = lyapunov_convergance(lorenz63!, u0, p, Δt, 1_000, 0.1)
		push!(λNs, λN)
	end
	plot(
		1:1_000,
		hcat(λNs...),
		label=hcat(["\$ Δt = $Δt\$" for Δt in Δts]...),
		xlabel="Iterations, \$N\$", ylabel="\$λ_1\$", 
		title="Converagance of \$λ_1\$"
	)
	scatter!(
		[length(λNs[3])], [λNs[3][end]], 
		label="\$λ_1 = $(round(λNs[3][end], digits=4))\$",
		markercolor=:red
	)
	savefig("plots/exercise3_13_2.png")
end

let δ0s = [0.01, 0.05, 0.1, 1.0, 5.0]
	λNs = []
	for δ0 in δ0s
		u0 = [20.0, 20.0, 20.0]
		p  = [10.0, 28.0, 8/3]
		λN = lyapunov_convergance(lorenz63!, u0, p, 1.0, 1_000, δ0)
		push!(λNs, λN)
	end
	plot(
		1:1_000,
		hcat(λNs...),
		label=hcat(["\$ δ_0 = $δ0\$" for δ0 in δ0s]...),
		xlabel="Iterations, \$N\$", ylabel="\$λ_1\$", 
		title="Converagance of \$λ_1\$"
	)
	scatter!(
		[length(λNs[3])], [λNs[3][end]], 
		label="\$λ_1 = $(round(λNs[3][end], digits=4))\$",
		markercolor=:red
	)
	savefig("plots/exercise3_13_3.png")
end
