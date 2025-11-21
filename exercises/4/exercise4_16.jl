using Plots, Statistics

function logistic_map(x,r)
	r * x * (1 - x)
end

function logsitic_map_orbit(x,r,N)
	orbit    = zeros(N)
	orbit[1] = x
	for i in 2:N
		orbit[i] = logistic_map(orbit[i-1],r)
	end
	return orbit
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

rc = 1 + √8

r = rc - 0.0002
N = 10000
o = logsitic_map_orbit(0.4,r,N)
scatter(
	1:N, o, 
	label=false, 
	xlim=[0,500],
	xlabel="\$n\$", 
	ylabel="\$x_n\$",
	title="Intermittency in Logistic map \$(r=1+√8-0.0002)\$",
	msw=0,
	ms=2,
	mc=:slateblue,
	legendposition=:topleft
)
plot!(1:N, o, lα=0.2, lc=:black, label="\$x_n\$")

δ = ones(N)
for (n,xn) in enumerate(o)
	if n ≥ 4
		δ[n] = abs(xn-o[n-3])
	end
end
plot!(1:N, δ, lc=:darkgreen, lα=0.5, label="\$δ\$")

laminar = []
for i in 2:length(δ)-1
	ϵ = 0.005
	if (δ[i] < ϵ) && (δ[i-1] < ϵ) && (δ[i+1] < ϵ)
		push!(laminar, (i,0.0))
	end
end
scatter!(
	first.(laminar), last.(laminar), 
	label="Laminar, \$⟨ℓ⟩=$(round(mean_laminar_length(first.(laminar)), digits=2))\$",
	mc=:lightblue, 
	msw=0
)
savefig("plots/exercise4_16.png")

#---------------------------------------------------------------------------------

mean_laminar_lenghts = []
δrs = 0.00001:0.00001:0.001
for δr in δrs
	local r = 1 + √8 - δr
	local N = 10000
	local o = logsitic_map_orbit(0.4,r,N)

	local δ = ones(N)
	for (n,xn) in enumerate(o)
		if n ≥ 4
			δ[n] = abs(xn-o[n-3])
		end
	end

	local laminar = []
	for i in 2:length(δ)-1
		ϵ = 0.005
		if (δ[i] < ϵ) && (δ[i-1] < ϵ) && (δ[i+1] < ϵ)
			push!(laminar, (i,0.0))
		end
	end

	mll = mean_laminar_length(first.(laminar))
	push!(mean_laminar_lenghts, mll)
end
plot(
	δrs, mean_laminar_lenghts,
	lw=2,
	ylim=[0,225],
	yticks=[Int(round(mean_laminar_lenghts[end])),50,100,150,200,225],
	label=false,
	xlabel="\$δr\$",
	ylabel="\$⟨ℓ⟩\$",
	title="Logistic map \$⟨ℓ⟩\$ for \$r=1 + √8 - δr\$"
)
savefig("plots/exercise4_16_mlls.png")
