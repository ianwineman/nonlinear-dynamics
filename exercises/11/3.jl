using DifferentialEquations
using LinearSolve, SparseArrays
using SparseConnectivityTracer, ADTypes
using Plots
using Statistics

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

u0 = zeros(M, M, 2)
u0[1:100,  :  , 1] .= 0.5
u0[ :   , 1:86, 2]  .= 0.5

tspan = (0.0, 1000.0)

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
println("solved in $(stats.time) s")

es = [extrema(solv.u[i][:, :, 1]) for i in 1:length(solv.u)]
cmin = minimum(first.(es))
cmax = maximum(last.(es))

x1 = [((m1 - 1) * L)/(M - 1) for m1 in 1:M]
x2 = [((m2 - 1) * L)/(M - 1) for m2 in 1:M]

p1 = heatmap(
	x1, 
	x2, 
	solv.u[end][:, :, 1], 
	aspect_ratio=:equal, 
	xlim=[0.0, L],
	ylim=[0.0, L],
   	xlabel="\$ x_1\$", 
   	ylabel="\$ x_2\$", 
   	colorbar_title="\$ u\$",
   	c=:haline,
   	clims=(cmin, cmax)
)

es = [extrema(solv.u[i][:, :, 2]) for i in 1:length(solv.u)]
cmin = minimum(first.(es))
cmax = maximum(last.(es))

p2 = heatmap(
	x1, 
	x2, 
	solv.u[end][:, :, 2], 
	aspect_ratio=:equal, 
	xlim=[0.0, L],
	ylim=[0.0, L],
   	xlabel="\$ x_1\$", 
   	ylabel="\$ x_2\$", 
   	colorbar_title="\$ w\$",
   	c=:turku,
   	clims=(cmin, cmax),
)

u = solv.u[end][:, :, 1]
w = solv.u[end][:, :, 2]
uμ = mean(u)
wμ = mean(w)
uσ =  std(u)
wσ =  std(w)

α = atan.((w .- wμ) ./ wσ, (u .- uμ) ./ uσ)

p3 = heatmap(
	x1, 
	x2, 
	α, 
	c=:lightrainbow,
	xlabel="\$ x_1\$", 
   	ylabel="\$ x_2\$", 
   	colorbar_title="\$ α\$",
)

p4 = heatmap(
	x2[75:100], 
	x1[120:145], 
	α[120:145, 75:100], 
	c=:lightrainbow,
	xlabel="\$ x_1\$", 
   	ylabel="\$ x_2\$", 
   	colorbar_title="\$ α\$",
)

plot(p1, p2, p3, p4)
savefig("plots/11.3.png")

#=
ENV["GKSwstype"] = "100"

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

mp4(anim, "plots/11.3.mp4", fps=60)
=#


