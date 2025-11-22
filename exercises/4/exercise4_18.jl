using Plots, FFTW, Statistics, LsqFit

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

function laminar_lengths(laminar_indices::Vector{Int64})
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
	return lengths
end

let
	N = 1_000
	o = intermittent_map_orbit(0.1, 0.0001, N)
	P = abs2.(rfft(o .- mean(o)))
	P ./= maximum(P)
	ν = rfftfreq(length(o))
	plot(
		ν,P,
		label="Power spectrum",
		title="Power spectrum",
		lc=:red
	)

	@. model(x,p) = p[1]*x^(-p[2])-1
	fit = curve_fit(model, ν[2:end], P[2:end], [1.0,0.5])
	c, α = coef(fit)
	y = zeros(length(ν))
	@. y = c*ν^-α-1
	plot!(
		ν,y,
		label="Power law, \$α=$(round(α, digits=3))\$",
		lc=:blue
	)
	savefig("plots/exercise4_18_pl.png")
end

let ϵs = 0.0001:0.0001:0.05
	N = 100_000

	αs = zeros(length(ϵs))
	for (i,ϵ) in enumerate(ϵs)
		o = intermittent_map_orbit(0.1, ϵ, N)
		P = abs2.(rfft(o .- mean(o)))
		P ./= maximum(P)
		ν = rfftfreq(length(o))

		@. model(x,p) = p[1]*x^(-p[2])-1
		fit = curve_fit(model, ν[2:end], P[2:end], [1.0,0.5])
		_, α = coef(fit)
		αs[i] = α
	end

	plot(
		ϵs,αs,
		title="Power spectrum power law exponent \$α\$",
		label=false,
		xlabel="\$ϵ\$",
		ylabel="\$α\$",
		lc=:orange
	)
	savefig("plots/exercise4_18_ae.png")
end

let
	N = 100_000
	o = intermittent_map_orbit(0.1, 0.0001, N)
	l = laminar(o)
	lls = sort(laminar_lengths(first.(l)))
	P = abs2.(rfft(lls .- mean(lls)))
	P ./= maximum(P)
	ν = rfftfreq(length(lls))

	plot(
		ν,P,
		label="Power spectrum",
		title="Laminar length distribution power law",
		lc=:red
	)

	@. model(x,p) = p[1]*x^(-p[2])-1
	fit = curve_fit(model, ν[2:end], P[2:end], [1.0,0.5])
	c, α = coef(fit)
	y = zeros(length(ν))
	@. y = c*ν^-α-1
	plot!(
		ν,y,
		label="Power law, \$α=$(round(α, digits=3))\$",
		lc=:blue
	)
	savefig("plots/exercise4_18_llpl.png")
end

let ϵs = 0.0001:0.0001:0.05
	N = 100_000

	αs = zeros(length(ϵs))
	for (i,ϵ) in enumerate(ϵs)
		o = intermittent_map_orbit(0.1, ϵ, N)
		l = laminar(o)
		lls = sort(laminar_lengths(first.(l)))
		P = abs2.(rfft(lls .- mean(lls)))
		P ./= maximum(P)
		ν = rfftfreq(length(lls))

		@. model(x,p) = p[1]*x^(-p[2])-1
		fit = curve_fit(model, ν[2:end], P[2:end], [1.0,0.5])
		_, α = coef(fit)
		αs[i] = α
	end

	plot(
		ϵs,αs,
		title="Laminar length distribution power law exponent \$α\$",
		label=false,
		xlabel="\$ϵ\$",
		ylabel="\$α\$",
		lc=:orange
	)
	savefig("plots/exercise4_18_llae.png")
end
