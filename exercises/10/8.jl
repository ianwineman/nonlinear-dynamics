using DifferentialEquations
using Plots
using Plots.PlotMeasures
using Statistics

function grid!(du, u, p, t)
	N, P, K, γ = p

	# u[1:N] is x1:xN
	# u[N+1:2N] is y1:yN

    for i in 1:N
    	du[i]     = u[N + i] # ẋ_i
    	du[i + N] = P[i] - (γ * u[i + N]) + sum([K[i, j] * sin(u[j] - u[i]) for j in 1:N]) # ẏ_i
    end
end

# original
K = [0 1 0 1 1 0 0 0;
 	 1 0 1 0 0 0 0 0;
 	 0 1 0 1 0 0 1 0;
	 1 0 1 0 0 0 0 0;
	 1 0 0 0 0 1 0 1;
	 0 0 0 0 1 0 1 0;
	 0 0 1 0 0 1 0 1;
	 0 0 0 0 1 0 1 0] .|> Float64
K .*= 1.03

P = [1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0]

N = first(size(K))

γ = 1.0

p = (N, P, K, γ)

prob = ODEProblem(grid!, zeros(2 * N), (0.0, 20.0), p)
solv = solve(prob)

p1 = plot(
   solv,
   idxs=1:8,
   lc=hcat([P[i] > 0.0 ? :blue : :red for i in 1:N]...),
   label=false,
   ylim=[-2π, 2π],
   xlabel="\$ t\$",
   ylabel="\$ ϕ_i\$",
   ytick=(-2π:π:2π, ["\$ -2π\$", "\$ -π\$", "\$ 0\$", "\$ π\$", "\$ 2π\$"]),
   title="original"
)
hspan!([π, 2π], color=:gray, alpha=0.1, label=false)
hspan!([-π, -2π], color=:gray, alpha=0.1, label=false)

# added capacity 1 → 5
K = [0 1 0 1 1 0 0 0;
 	 1 0 1 0 0 0 0 0;
 	 0 1 0 1 0 0 1 0;
	 1 0 1 0 0 0 0 0;
	 1 0 0 0 0 1 0 1;
	 0 0 0 0 1 0 1 0;
	 0 0 1 0 0 1 0 1;
	 0 0 0 0 1 0 1 0] .|> Float64
K .*= 1.03
K[1, 5] *= 2
K[5, 1] *= 2

P = [1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0]

N = first(size(K))

γ = 1.0

p = (N, P, K, γ)

prob = ODEProblem(grid!, zeros(2 * N), (0.0, 20.0), p)
solv = solve(prob)

p2 = plot(
   solv,
   idxs=1:8,
   lc=hcat([P[i] > 0.0 ? :blue : :red for i in 1:N]...),
   label=false,
   ylim=[-2π, 2π],
   xlabel="\$ t\$",
   ylabel="\$ ϕ_i\$",
   ytick=(-2π:π:2π, ["\$ -2π\$", "\$ -π\$", "\$ 0\$", "\$ π\$", "\$ 2π\$"]),
   title="added capacity"
)
hspan!([ π,  2π], color=:gray, alpha=0.1, label=false)
hspan!([-π, -2π], color=:gray, alpha=0.1, label=false)

# added link 1 → 6
K = [0 1 0 1 1 1 0 0;
 	 1 0 1 0 0 0 0 0;
 	 0 1 0 1 0 0 1 0;
	 1 0 1 0 0 0 0 0;
	 1 0 0 0 0 1 0 1;
	 1 0 0 0 1 0 1 0;
	 0 0 1 0 0 1 0 1;
	 0 0 0 0 1 0 1 0] .|> Float64
K .*= 1.03

P = [1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0]

N = first(size(K))

γ = 1.0

p = (N, P, K, γ)

prob = ODEProblem(grid!, zeros(2 * N), (0.0, 20.0), p)
solv = solve(prob)

p3 = plot(
   solv,
   idxs=1:8,
   lc=hcat([P[i] > 0.0 ? :blue : :red for i in 1:N]...),
   label=false,
   ylim=[-2π, 2π],
   xlabel="\$ t\$",
   ylabel="\$ ϕ_i\$",
   ytick=(-2π:π:2π, ["\$ -2π\$", "\$ -π\$", "\$ 0\$", "\$ π\$", "\$ 2π\$"]),
   title="added link"
)
hspan!([π, 2π], color=:gray, alpha=0.1, label=false)
hspan!([-π, -2π], color=:gray, alpha=0.1, label=false)

# added capacity heatmap
S = zeros(N, N)

for i in 1:N
	for j in 1:N
		if i != j
			K = [0 1 0 1 1 1 0 0;
			 	 1 0 1 0 0 0 0 0;
			 	 0 1 0 1 0 0 1 0;
				 1 0 1 0 0 0 0 0;
				 1 0 0 0 0 1 0 1;
				 1 0 0 0 1 0 1 0;
				 0 0 1 0 0 1 0 1;
				 0 0 0 0 1 0 1 0] .|> Float64
			K[i, j] += 1
			K[j, i] += 1

			P = [1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0]

			N = first(size(K))

			γ = 1.0

			p = (N, P, K, γ)

			prob = ODEProblem(grid!, zeros(2 * N), (0.0, 20.0), p)
			solv = solve(prob)

			sync = maximum([mean(abs.(u[1:8])) for u in solv.u]) < π

			if sync
				S[i, j] = 1
				S[j, i] = 1
			else
				S[i, j] = -1
				S[j, i] = -1
			end
		end
	end
end

p4 = heatmap(
	S, 
	c=cgrad([:red, :black, :blue], 3, categorical=true),
	cbar_title="red = no sync, blue = sync",
	colorbar_ticks=false,
	xtick=1:8,
	ytick=1:8,
	title="added capacity"
)

# reduced capacity heatmap
S = zeros(N, N)

for i in 1:N
	for j in 1:N
		if i != j
			K = [0 1 0 1 1 1 0 0;
			 	 1 0 1 0 0 0 0 0;
			 	 0 1 0 1 0 0 1 0;
				 1 0 1 0 0 0 0 0;
				 1 0 0 0 0 1 0 1;
				 1 0 0 0 1 0 1 0;
				 0 0 1 0 0 1 0 1;
				 0 0 0 0 1 0 1 0] .|> Float64
			K[i, j] -= ifelse(K[i, j] > 0, 1, 0)
			K[j, i] -= ifelse(K[j, i] > 0, 1, 0)

			P = [1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0]

			N = first(size(K))

			γ = 1.0

			p = (N, P, K, γ)

			prob = ODEProblem(grid!, zeros(2 * N), (0.0, 20.0), p)
			solv = solve(prob)

			sync = maximum([mean(abs.(u[1:8])) for u in solv.u]) < π

			if sync
				S[i, j] = 1
				S[j, i] = 1
			else
				S[i, j] = -1
				S[j, i] = -1
			end
		end
	end
end

p5 = heatmap(
	S, 
	c=cgrad([:red, :black, :blue], 3, categorical=true),
	cbar_title="red = no sync, blue = sync",
	colorbar_ticks=false,
	xtick=1:8,
	ytick=1:8,
	title="reduced capacity"
)

# plot all plots
layout = @layout [a b c; d e]
plot(p1, p2, p3, p4, p5, layout=layout, size=(900, 600), margin=3mm)
savefig("plots/10.8.png")
