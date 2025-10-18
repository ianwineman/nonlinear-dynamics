### Doing 3.16 with Hénon–Heiles system instead of H(q,p), thus we reproduce Fig. 3.7 here.
###
using DifferentialEquations, LinearAlgebra, PredefinedDynamicalSystems, Plots

function lyapunov(ds, u0, p, Δt, N, δ0)
	xrefr_prob = ODEProblem(ds, u0,                          (0.0, Δt*N), p)
	xtest_prob = ODEProblem(ds, u0 .+ (δ0/sqrt(length(u0))), (0.0, Δt*N), p)

	xrefr_int = init(xrefr_prob, Tsit5())
	xtest_int = init(xtest_prob, Tsit5())

	δ = zeros(N)
	for i in 1:N
		step!(xrefr_int, Δt, true)
		step!(xtest_int, Δt, true)
		δ[i] = norm(xrefr_int.u .- xtest_int.u)
		reinit!(xtest_int, xrefr_int.u .+ (xtest_int.u .- xrefr_int.u) .* (δ0/δ[i]))
	end

	λ = (1/(N*Δt)) * sum(log.(δ ./ δ0))
end

let u0s = PredefinedDynamicalSystems.henonheiles_ics(0.13, 24), T = 1_250.0, Δt = 1.0, pc_ϵ=0.1
	poincare_section = Vector{Tuple{Vector{Float64}, Float64}}()

	function henon_heiles_rule(u)
		x, y, vx, vy = u
		return [
			vx
			vy
			-x - 2*x*y
			-y - (x^2 - y^2)
		]
	end
	function henon_heiles_rule_de!(du, u, p, t)
		x, y, vx, vy = u
		du[1] = vx
		du[2] = vy
		du[3] = -x - 2*x*y
		du[4] = -y - (x^2 - y^2)
	end

	for (i, u0) in enumerate(u0s)
		println(vcat([" " for _ in 1:(length(string(length(u0s)))-length(string(i)))], ["$i", " of $(length(u0s))"])...)
		λ1 = lyapunov(henon_heiles_rule_de!, u0, [], 1.0, Int(T), 0.01)

		prob = ODEProblem(henon_heiles_rule_de!, u0, (0.0,T))
		sol  = solve(prob)
		
		fa = findall(x->((abs(x[1])<=pc_ϵ) && (x[3]>=0.0)), sol.u)
		points = sol.u[fa]

		push!(poincare_section, [(p,λ1) for p in points]...)
	end

	marker_colors = [
		:spring,
	    :imola,
	    :Purples,
	    :Greens,
	    :inferno,
	    :winter,
	    :nuuk,
	    :lipari
	]
	for color in marker_colors
		scatter(
			[p[1][4] for p in poincare_section], 
			[p[1][2] for p in poincare_section],
			marker_z=[p[2] for p in poincare_section],
			markerstrokewidth=0.0,
			markersize=2,
			mc=color,
			label=false,
			backgroundinside=:black,
			backgroundoutside=:black,
			grid=false,
			xtick=false,
			ytick=false,
			cbar=false,
			ratio=:equal,
			size=(1179,2556)
		)
		savefig("plots/fun/exercise3_16_fun_$(color).png")
	end
end
