using Plots
using Statistics, StatsBase, Random
using DifferentialEquations
using ProgressMeter

"Permutation for `X`"
function permra(X::Vector{Float64})
	return [findfirst(==(y), sort(X)) for y in X]
end

"Permutation for `X`"
function permra(X::Vector{Tuple{Float64, Float64}})
	A, B = first.(X), last.(X)
	return [(findfirst(==(A[i]), sort(A)), findfirst(==(B[i]), sort(B))) for i in 1:length(X)]
end

"""Shannon entropy of `X` in bits 
using permutation of relative amplitude. 
Permutations of length `d`
"""
function entropy(X, d)
	pras = [permra(X[i:(i+d-1)]) for i in 1:(length(X)-d+1)]
	upras = unique(pras)
	Pis = [count(==(up), upras) |> float for up in upras]
	Pis ./= length(Pis)
	return -sum(Pis .* log2.(Pis))
end


"""Mutual information of `X` and `Y` in bits
using permutations of relative amplitude.
Permutations of length `d`
"""
function mutualinformation(X::Vector{Float64}, Y::Vector{Float64}, d)
	return entropy(X, d) + entropy(Y, d) - entropy(collect(zip(X, Y)), d)
end

"""Test whether `X` and `Y` are correlated using Mutual information. \n
`σs` determines required threshold to reject null hypothesis (`X` and `Y` are uncorrelated)"""
function micorrelation(X, Y, d; N::Int64=10_000, σs=2)
	mi = mutualinformation(X, Y, d)

	null = zeros(N)
	@showprogress Threads.@threads for i in 1:N
		SX = shuffle(X)
		SY = shuffle(Y)
		null[i] = mutualinformation(SX, SY, d)
	end
	μ = mean(null)
	σ = std(null)
	verdict = abs(mi - μ) > σs*σ
	return (null, μ, σ, mi, σs, verdict)
end

function rössler!(du, u, p, t)
	a, b, c = p
	x, y, z = u

	du[1] = -y - z
	du[2] = x + a * y
	du[3] = b + z * (x - c)
end

prob = ODEProblem(rössler!, [1, -2, 0.1], (0.0, 200.0), [0.2, 0.2, 5.7])
solv = solve(prob; dtmax=0.1)
X, Y = getindex.(solv.u, 1), getindex.(solv.u, 2)

null, μ, σ, mi, σs, verdict = micorrelation(X, Y, 5)

p1 = histogram(
	null,
	label="Null PDF",
	xlabel="Mutual information",
	color=:lavender,
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

p2 = plot(
	solv, 
	idxs=(1,2,3),
	label="Rössler system",
	xlabel="\$ x\$",
	ylabel="\$ y\$",
	zlabel="\$ z\$",
	lα=0.5
)

p3 = plot(
	solv.t, 
	[X, Y],
	label=["\$ x\$" "\$ y\$"],
	xlabel="\$ t\$"
)

plot(
	p1, p2, p3, 
	size=(900, 600), 
	layout=@layout([
		a
		b c
	])
)
savefig("plots/exercise7_3.png")
