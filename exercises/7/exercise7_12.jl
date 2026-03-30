using CSV, DataFrames
using DynamicalSystems
using FFTW
using Plots
using Statistics

df = CSV.read("book/exercise_data/11.csv", DataFrame; header=false)
x, y = df[:, 1], df[:, 2]

function detrend(x)
	y = zeros(length(x) - 1)
	for i in 1:(length(x) - 1)
		y[i] = x[i + 1] - x[i]
	end
	return y
end

function surrogate(X; alg=:FT)
	if alg == :FT
		irfft(rfft(X) .* [exp(im * 2π * rand(Float64)) for _ in 1:length(rfft(X))], length(X))
	elseif alg == :AAFT
		S = irfft(rfft(X) .* [exp(im * 2π * rand(Float64)) for _ in 1:length(rfft(X))], length(X))
		S[sortperm(S)] .= sort(X)
		return S
	else
		@error "Unknown algorithm: $alg"
	end
end

# x
stat = sum(selfmutualinfo(x, 1:6)) / 6
surrogate_count = 1_000
null = zeros(surrogate_count)

for i in 1:surrogate_count
	s = surrogate(detrend(x); alg=:AAFT)
	null[i] = sum(selfmutualinfo(s, 1:6)) / 6
end

μ, σ = mean(null), std(null)
p1 = histogram(
	null,
	label="\$ null\$",
	xlabel="\$ SMI\$",
	c=:lightgray
)
vline!(
	[μ],
	label="μ",
	lw=2,
	lc=:blue
)
vline!(
	[μ - 2*σ, μ + 2*σ],
	label="μ ± 2σ",
	lw=2,
	ls=:dot,
	lc=:black
)
vline!(
	[stat],
	label="\$ SMI(x)\$",
	lw=2,
	lc=:red
)

# y
stat = sum(selfmutualinfo(y, 1:6)) / 6
surrogate_count = 1_000
null = zeros(surrogate_count)

for i in 1:surrogate_count
	s = surrogate(detrend(y); alg=:AAFT)
	null[i] = sum(selfmutualinfo(s, 1:6)) / 6
end

μ, σ = mean(null), std(null)
p2 = histogram(
	null,
	label="\$ null\$",
	xlabel="\$ SMI\$",
	c=:lightgray
)
vline!(
	[μ],
	label="μ",
	lw=2,
	lc=:blue
)
vline!(
	[μ - 3*σ, μ + 3*σ],
	label="μ ± 3σ",
	lw=2,
	ls=:dot,
	lc=:black
)
vline!(
	[stat],
	label="\$ SMI(y)\$",
	lw=2,
	lc=:red
)

plot(p1, p2)
savefig("plots/exercise7_12.png")
