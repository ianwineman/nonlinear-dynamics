using DifferentialEquations
using Plots
using LinearAlgebra

function rössler!(du, u, p, t)
	a, b, c = p
	x, y, z = u

    du[1] = -y - z
    du[2] = x + a*y
    du[3] = b + z*(x - c)
end

function embed(X, τ, d)
	return [[X[i + n*τ] for n in 0:(d-1)] for i in 1:(length(X) - d - 1)]
end

function recurrencematrix(X, ϵ; n=length(X))
	R = zeros(Int, n, n)
	for i in 1:n
		for j in 1:n
			if norm(X[i] - X[j]) < ϵ
				R[i, j] = 1
			end
		end
	end
	return R
end

chaotic_prob = ODEProblem(rössler!, [1, -2, 0.1], (0.0, 1_000.0), [0.2, 0.2, 5.7])
chaotic_solv = solve(prob)

periodic_prob = ODEProblem(rössler!, [1, -2, 0.1], (0.0, 1_000.0), [0.2, 0.2, 2.5])
periodic_solv = solve(prob)

rand_solv = [rand(3) for _ in 1:4_000]

ns = collect(100:100:1_000)
ϵs = collect(range(0.01, 1.0; length=10))
ds = collect(1:10)

anim = @animate for i in 1:10
	p1 = heatmap(
			recurrencematrix(embed(first.(chaotic_solv.u),  1, 3), 0.25; n=ns[i]), 
			aspect_ratio=:equal,
			title="length=$(ns[i])",
			titlefontsize=10,
			cbar=false,
			xlim=[0, ns[i]],
			ylim=[0, ns[i]],
			xtick=[0, ns[i]],
			ytick=[0, ns[i]],
			c=:devon
		)
	p2 = heatmap(
			recurrencematrix(embed(first.(periodic_solv.u), 1, 3), 0.25; n=ns[i]), 
			aspect_ratio=:equal,
			title="length=$(ns[i])",
			titlefontsize=10,
			cbar=false,
			xlim=[0, ns[i]],
			ylim=[0, ns[i]],
			xtick=[0, ns[i]],
			ytick=[0, ns[i]],
			c=:devon
		)
	p3 = heatmap(
			recurrencematrix(embed(first.(rand_solv),       1, 3), 0.25; n=ns[i]), 
			aspect_ratio=:equal,
			title="length=$(ns[i])",
			titlefontsize=10,
			cbar=false,
			xlim=[0, ns[i]],
			ylim=[0, ns[i]],
			xtick=[0, ns[i]],
			ytick=[0, ns[i]],
			c=:devon
		)

	p4 = heatmap(
			recurrencematrix(embed(first.(chaotic_solv.u),  1, 3), ϵs[i]; n=1_000), 
			aspect_ratio=:equal,
			title="ϵ=$(ϵs[i])",
			titlefontsize=10,
			cbar=false,
			xlim=[0, 1_000],
			ylim=[0, 1_000],
			xtick=[0, 1_000],
			ytick=[0, 1_000],
			c=:devon
		)
	p5 = heatmap(
			recurrencematrix(embed(first.(periodic_solv.u), 1, 3), ϵs[i]; n=1_000), 
			aspect_ratio=:equal,
			title="ϵ=$(ϵs[i])",
			titlefontsize=10,
			cbar=false,
			xlim=[0, 1_000],
			ylim=[0, 1_000],
			xtick=[0, 1_000],
			ytick=[0, 1_000],
			c=:devon
		)
	p6 = heatmap(
			recurrencematrix(embed(first.(rand_solv),       1, 3), ϵs[i]; n=1_000), 
			aspect_ratio=:equal,
			title="ϵ=$(ϵs[i])",
			titlefontsize=10,
			cbar=false,
			xlim=[0, 1_000],
			ylim=[0, 1_000],
			xtick=[0, 1_000],
			ytick=[0, 1_000],
			c=:devon
		)

	p7 = heatmap(
			recurrencematrix(embed(first.(chaotic_solv.u),  1, ds[i]), 0.25; n=1_000), 
			aspect_ratio=:equal,
			title="d=$(ds[i])",
			xlabel="Chaotic Rössler",
			titlefontsize=10,
			cbar=false,
			xlim=[0, 1_000],
			ylim=[0, 1_000],
			xtick=[0, 1_000],
			ytick=[0, 1_000],
			c=:devon
		)
	p8 = heatmap(
			recurrencematrix(embed(first.(periodic_solv.u), 1, ds[i]), 0.25; n=1_000), 
			aspect_ratio=:equal,
			title="d=$(ds[i])",
			xlabel="Periodic Rössler",
			titlefontsize=10,
			cbar=false,
			xlim=[0, 1_000],
			ylim=[0, 1_000],
			xtick=[0, 1_000],
			ytick=[0, 1_000],
			c=:devon
		)
	p9 = heatmap(
			recurrencematrix(embed(first.(rand_solv),       1, ds[i]), 0.25; n=1_000), 
			aspect_ratio=:equal,
			title="d=$(ds[i])",
			xlabel="Random numbers",
			titlefontsize=10,
			cbar=false,
			xlim=[0, 1_000],
			ylim=[0, 1_000],
			xtick=[0, 1_000],
			ytick=[0, 1_000],
			c=:devon
		)

	plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(3, 3), size=(900, 600))
end
gif(anim, "plots/exercise8_11.gif", fps=1)
