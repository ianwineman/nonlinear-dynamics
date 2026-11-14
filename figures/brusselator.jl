using DifferentialEquations
using LinearSolve, SparseArrays
using SparseConnectivityTracer, ADTypes
using Random

function brusselator!(du, u, p, t)
	a, b, d, L, M  = p
	h = L / (M - 1)

	for m_i in 1:M
		for m_j in 1:M
			fps_1 = (u[1 + mod(m_i, M - 1), m_j, 1] + u[1 + mod(m_i - 2, M - 1), m_j, 1] + u[m_i, 1 + mod(m_j, M - 1), 1] + u[m_i, 1 + mod(m_j - 2, M - 1), 1] - 4*u[m_i, m_j, 1]) / h^2
			fps_2 = (u[1 + mod(m_i, M - 1), m_j, 2] + u[1 + mod(m_i - 2, M - 1), m_j, 2] + u[m_i, 1 + mod(m_j, M - 1), 2] + u[m_i, 1 + mod(m_j - 2, M - 1), 2] - 4*u[m_i, m_j, 2]) / h^2
			du[m_i, m_j, 1] = 1 - ((b + 1) * u[m_i, m_j, 1]) + (a * (u[m_i, m_j, 1])^2 * u[m_i, m_j, 2]) + fps_1
			du[m_i, m_j, 2] = b * u[m_i, m_j, 1] - (a * (u[m_i, m_j, 1])^2 * u[m_i, m_j, 2]) + d*fps_2
		end
	end
end

a = 9.0
b = 9.8
d = 3.0
L = 50.0
M = 100
p = (a, b, d, L, M)

Random.seed!(4321)
u0 = rand(M, M, 2)
tspan = (0.0, 3000.0)

detector = TracerSparsityDetector()
du0 = similar(u0)
jac_sparsity = jacobian_sparsity(
	(du, u) -> brusselator!(du, u, p, 0.0), du0, u0, detector
)

f = ODEFunction(brusselator!; jac_prototype=float.(jac_sparsity))
prob_sparse = ODEProblem(f, u0, tspan, p)

alg = FBDF(linsolve=KLUFactorization())
stats = @timed solve(prob_sparse, alg)
solv = stats.value

es = [extrema(solv.u[i][:, :, 1] .- solv.u[i][:, :, 2]) for i in 1:length(solv.u)]
cmin = minimum(first.([extrema(solv.u[i][:, :, 1] .- solv.u[i][:, :, 2]) for i in 1:length(solv.u)]))
cmax = maximum(last.([extrema(solv.u[i][:, :, 1] .- solv.u[i][:, :, 2]) for i in 1:length(solv.u)]))
clim = maximum(abs.([cmin, cmax]))

anim = @animate for (i, t) in enumerate(solv.t)
	u1 = solv.u[i][:, :, 1]
	u2 = solv.u[i][:, :, 2]
	u = u1 .- u2
   heatmap(
   		u, 
   		title="\$ t=$(round(t; digits=2))\$", 
   		xlabel="\$ x_1\$", 
   		ylabel="\$ x_2\$", 
   		colorbar_title="\$ u\$",
   		c=cgrad([:blue, :red]),
   		clims=(-clim, clim)
   	)
end

println(stats.time)
gif(anim, "figures/brusselator.gif", fps=30)
