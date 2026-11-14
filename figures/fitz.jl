using DifferentialEquations
using LinearSolve, SparseArrays
using SparseConnectivityTracer, ADTypes
using Random
using Plots

function fitzhughnagumo!(du, u, p, t)
	a, b, d, ϵ, L, M  = p
	h = L / (M - 1)

	for m_i in 1:M
		for m_j in 1:M
			m_i_p = ifelse(m_i + 1 > M, m_i, m_i + 1)
			m_i_m = ifelse(m_i - 1 < 1, 1, m_i - 1)
			m_j_p = ifelse(m_j + 1 > M, m_j, m_j + 1)
			m_j_m = ifelse(m_j - 1 < 1, 1, m_j - 1)

			fps = (u[m_i_p, m_j, 1] + u[m_i_m, m_j, 1] + u[m_i, m_j_p, 1] + u[m_i, m_j_m, 1] - 4*u[m_i, m_j, 1]) / h^2
			
			du[m_i, m_j, 1] = a * u[m_i, m_j, 1] * (u[m_i, m_j, 1] - b) * (1 - u[m_i, m_j, 1]) - u[m_i, m_j, 2] + d * fps
			du[m_i, m_j, 2] = ϵ * (u[m_i, m_j, 1] - u[m_i, m_j, 2])
		end
	end
end


a = 3.0
b = 0.2
d = 1.0
ϵ = 0.01
L = 300.0
M = 150
p = (a, b, d, ϵ, L, M)

# Figure A
Random.seed!(4321)
u0 = zeros(M, M, 2)
u0[:, :, 1] .+= (rand(M, M) * 0.01)
u0[40:41, 40:41, 1] .= 0.999
u0[75:76, 75:76, 1] .= 0.999

tspan = (0.0, 250.0)

detector = TracerSparsityDetector()
du0 = similar(u0)
jac_sparsity = jacobian_sparsity(
	(du, u) -> fitzhughnagumo!(du, u, p, 0.0), du0, u0, detector
)

f = ODEFunction(fitzhughnagumo!; jac_prototype=float.(jac_sparsity))
prob_sparse = ODEProblem(f, u0, tspan, p)

alg = FBDF(linsolve=KLUFactorization())
stats = @timed solve(prob_sparse, alg, saveat=range(tspan...))
solv = stats.value

es = [extrema(solv.u[i][:, :, 1]) for i in 1:length(solv.u)]
cmin = minimum(first.(es))
cmax = maximum(last.(es))

x1 = [((m1 - 1) * L)/(M - 1) for m1 in 1:M]
x2 = [((m2 - 1) * L)/(M - 1) for m2 in 1:M]
anim = @animate for (i, t) in enumerate(solv.t)
    u = solv.u[i][:, :, 1]
	heatmap(
		x1,
		x2,
   		u, 
   		title="\$ t=$(round(t; digits=2))\$", 
   		xlabel="\$ x_1\$", 
   		ylabel="\$ x_2\$", 
   		colorbar_title="\$ u\$",
   		c=:haline,
   		clims=(cmin, cmax),
   		aspect_ratio=:equal,
   		xlim=[0.0, L]
   	)
end

println("solved in $(stats.time) s")
mp4(anim, "figures/fitzA.mp4", fps=60)

# Figure B
Random.seed!(4321)
u0 = zeros(M, M, 2)
u0[:, :, 1] .+= (rand(M, M) * 0.01)
u0[:, 10:11, 1] .= 0.999

tspan = (0.0, 300.0)

prob_sparse = ODEProblem(f, u0, tspan, p)

alg = FBDF(linsolve=KLUFactorization())
stats1 = @timed solve(prob_sparse, alg, saveat=range(tspan...))
solv1 = stats1.value

u0 = solv1.u[end]
u0[55:73, 1:70, 1] .= 0.999

tspan = (300.0, 1000.0)

prob_sparse = ODEProblem(f, u0, tspan, p)

alg = FBDF(linsolve=KLUFactorization())
stats2 = @timed solve(prob_sparse, alg, saveat=range(tspan...))
solv2 = stats2.value

solvu = vcat(solv1.u, solv2.u)
solvt = vcat(solv1.t, solv2.t)

es = [extrema(solvu[i][:, :, 1]) for i in 1:length(solvu)]
cmin = minimum(first.(es))
cmax = maximum(last.(es))

x1 = [((m1 - 1) * L)/(M - 1) for m1 in 1:M]
x2 = [((m2 - 1) * L)/(M - 1) for m2 in 1:M]
anim = @animate for (i, t) in enumerate(solvt)
    u = solvu[i][:, :, 1]
	heatmap(
		x1,
		x2,
   		u, 
   		title="\$ t=$(round(t; digits=2))\$", 
   		xlabel="\$ x_1\$", 
   		ylabel="\$ x_2\$", 
   		colorbar_title="\$ u\$",
   		c=:haline,
   		clims=(cmin, cmax),
   		aspect_ratio=:equal,
   		xlim=[0.0, L]
   	)
end

println("solved in $(stats1.time + stats2.time) s")
mp4(anim, "figures/fitzB.mp4", fps=60)
