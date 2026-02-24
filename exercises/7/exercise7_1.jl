using Plots
using Statistics, StatsBase, Random

"Amplitude bins and their associated probabilities for 1d timeseries `X` with bin width `ϵ`"
function amplitudebinning(X::Vector{Float64}, ϵ)
	min, max = extrema(X)
	bins = ceil(Int, (max - min)/ϵ)
	r = range(min, max; length=bins)
	Ps = zeros(Int, bins - 1)
	for x in X
		i = findfirst(y -> y > x, r)
		if !isnothing(i) && i != 1
			Ps[i - 1] += 1
		end
	end
	return Ps ./ length(X)
end


"Amplitude bins and their associated probabilities for 2d timeseries `X` with bin widths `ϵ`"
function amplitudebinning(X::Vector{Tuple{Float64, Float64}}, ϵ)
	Y = last.(X)
	X = first.(X)

	xmin, xmax = extrema(X)
	ymin, ymax = extrema(Y)

	xbins = ceil(Int, (xmax - xmin)/ϵ)
	ybins = ceil(Int, (ymax - ymin)/ϵ)

	xr = range(xmin, xmax; length=xbins)
	yr = range(ymin, ymax; length=ybins)

	Ps = zeros((xbins - 1, ybins - 1))

	for (x, y) in zip(X, Y)
		i = findfirst(z -> z > x, xr)
		j = findfirst(z -> z > y, yr)
		if !isnothing(i) && !isnothing(j) && i != 1 && j != 1
			Ps[i - 1, j - 1] += 1
		end
	end 
	return vec(Ps ./ length(X))
end

"Shannon entropy of `X` in bits"
function entropy(X; ϵ=0.1)
	bins = amplitudebinning(X, ϵ)
	bins = bins[findall(!=(0), bins)]
	return -sum(bins .* log2.(bins))
end

"Mutual information of `X` and `Y` in bits"
function mutualinformation(X::Vector{Float64}, Y::Vector{Float64}; ϵ=0.1)
	return entropy(X; ϵ=ϵ) + entropy(Y; ϵ=ϵ) - entropy(collect(zip(X, Y)); ϵ=ϵ)
end


#==#

ϕs = range(0, π; length=100)
cors = []
for ϕ in ϕs
	t = range(0, 2π * 10; length=1_000)
	x = cos.(t)
	y = cos.(t .+ ϕ)
	pcor = cor(x, y)
	scor = corspearman(x, y)
	mi = mutualinformation(x, y)

	null = zeros(10_000)
	Threads.@threads for i in 1:10_000
		sx = shuffle(x)
		sy = shuffle(y)
		null[i] = mutualinformation(sx, sy)
	end
	μ = mean(null)
	σ = std(null)
	conf = abs(mi - μ) > 2*σ

	push!(cors, [pcor, scor, mi, Int(conf)])
end
p1 = plot(
	ϕs, 
	[getindex.(cors, 1) getindex.(cors, 2) getindex.(cors, 3) getindex.(cors, 4)],
	label=["Pearson" "Spearman" "Mutual information" "\$1 \$ if \$ |mi - μ| > 2σ \$, \$ 0\$ otherwise"],
	xlabel="\$ ϕ\$",
	xlim=[0, π],
	ls=[:solid :solid :solid :dash],
)

t = range(0, 2π * 10; length=1_000)
x = cos.(t)
y = cos.(t .+ π/3)
p2 = plot(
	t, 
	[x y],
	xlabel="\$ t\$",
	ylabel="\$ f(t)\$",
	label=["\$ cos(t)\$" "\$ cos(t + π/3)\$"]
)

m = mutualinformation(x, y)
null = zeros(10_000)
Threads.@threads for i in 1:10_000
	sx = shuffle(x)
	sy = shuffle(y)
	null[i] = mutualinformation(sx, sy)
end

μ = mean(null)
σ = std(null)
p3 = histogram(
	null,
	label="Null PDF",
	xlabel="Mutual information",
	linecolor=false,
	legendtitle="\$ ϕ=π/3\$"
)
vline!([μ], label="\$ μ\$", lw=2)
vline!([μ + 2*σ, μ - 2*σ], label="\$ μ ± 2σ\$", lw=1.5, ls=:dash)
vline!([m], label="\$ m\$", lw=1.5)

plot(
	p1, p2, p3,
	layout=@layout([
		a
		b c
	]),
	size=(900, 600)
)
savefig("plots/exercise7_1.png")
