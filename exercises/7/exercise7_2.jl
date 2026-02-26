using CSV, DataFrames
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

"""Test whether `X` and `Y` are correlated using Mutual information. \n
`σs` determines required threshold to reject null hypothesis (`X` and `Y` are uncorrelated)"""
function micorrelation(X, Y; N::Int64=10_000, σs=2)
	mi = mutualinformation(X, Y)

	null = zeros(N)
	Threads.@threads for i in 1:N
		SX = shuffle(X)
		SY = shuffle(Y)
		null[i] = mutualinformation(SX, SY)
	end
	μ = mean(null)
	σ = std(null)
	verdict = abs(mi - μ) > σs*σ
	return (null, μ, σ, mi, σs, verdict)
end

df = CSV.read("book/exercise_data/11.csv", DataFrame; header=false)
X, Y = df[:, 1], df[:, 2]
null, μ, σ, mi, σs, verdict = micorrelation(X, Y)

histogram(
	null,
	label="Null PDF",
	xlabel="Mutual information",
	title="Dataset 11",
	color=:lavender,
	background_color=RGB(255/255, 248/255, 231/255)
)
vline!(
	[μ],
	lw=1.5,
	ls=:dash,
	lc=:purple4,
	label="\$ μ\$"
)
vline!(
	[μ - σs*σ, μ + σs*σ],
	lw=1.5,
	ls=:dashdot,
	lc=:purple4,
	label="\$ μ \$ \$± \$ \$$(σs)σ\$"
)
vline!(
	[mi],
	lw=2,
	lc=:purple1,
	label="\$ mi\$"
)
savefig("plots/exercise7_2.png")
