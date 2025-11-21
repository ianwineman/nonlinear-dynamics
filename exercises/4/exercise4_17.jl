using Plots, Statistics

function intermittent_map(x,ϵ)
	((1+ϵ)*x + (1-ϵ)*x^2) % 1.0
end

function intermittent_map_orbit(x,ϵ,N)
	orbit    = zeros(N)
	orbit[1] = x
	for i in 2:N
		orbit[i] = intermittent_map(orbit[i-1],ϵ)
	end
	return orbit
end

function laminar(orbit)
	laminars = []
	for (i,x) in enumerate(orbit)
		if i < length(orbit) && abs(x - orbit[i+1]) < 0.01
			push!(laminars, (i,-0.01))
		end
	end
	return laminars
end

function mean_laminar_length(laminar_indices::Vector{Int64})
	lengths = []
	current_length = 1
	for (j,li) in enumerate(laminar_indices)
		if j < length(laminar_indices) && laminar_indices[j+1] == li + 1
			current_length += 1
		else
			if current_length != 1
				push!(lengths, current_length)
				current_length = 1
			end
		end
	end
	return mean(lengths)
end

N = 10_000
o = intermittent_map_orbit(0.1, 0.0001, N)
scatter(
	1:N, o,
	msw=0,
	ms=2,
	label=false,
	xlim=[0,1000],
	xlabel="\$n\$",
	ylabel="\$x_n\$",
	title="\$x_{n+1}=[(1+ϵ)x_n + (1-ϵ)x_n^2] \$ \$mod\$ \$1\$"
)
l = laminar(o)
scatter!(first.(l), last.(l), msw=0, ms=2, label="Laminar, \$⟨ℓ⟩=$(round(mean_laminar_length(first.(l)), digits=2))\$")
savefig("plots/exercise4_17.png")

let ϵs = 0.0001:0.0001:0.05
	N = 10_000
	mlls = []
	for ϵ in ϵs
		o = intermittent_map_orbit(0.1, ϵ, N)
		l = laminar(o)
		mll = mean_laminar_length(first.(l))
		push!(mlls, mll)
	end
	plot(
		ϵs,mlls,
		lw=2,
		xlim=[0.0001,0.05],
		ylim=[0,226],
		yticks=[Int(round(mlls[end])),50,100,150,200,226],
		label=false,
		xlabel="\$ϵ\$",
		ylabel="\$⟨ℓ⟩\$",
		title="\$x_{n+1}=[(1+ϵ)x_n + (1-ϵ)x_n^2] \$ \$mod\$ \$1\$"
	)
	savefig("plots/exercise4_17_mlls.png")
end
