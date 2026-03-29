using DifferentialEquations
using Statistics, Random
using Plots, StatsPlots
using FFTW
using LinearAlgebra
using DynamicalSystems
using ProgressMeter

function rössler!(du, u, p, t)
	a, b, c = p
	x, y, z = u

	du[1] = -y - z
	du[2] = x + a * y
	du[3] = b + z * (x - c)
end

prob = ODEProblem(rössler!, [1, 1, 1], (0.0, 500.0), [0.165, 0.2, 10.0])
solv = solve(prob; dtmax=0.1)

X = first.(solv.u)
X ./= std(X)
X .+= randn(length(X)) * std(X) * 0.1
X .+= 2 # for plotting (per Fig. 7.5 code)

p1 = plot(1:1_000, X[1:1_000], xlim=[0, 1_000], label="Rössler", lc=RGB(159/255, 178/255, 97/255), xlabel="\$ t\$")

Random.seed!(77163)
rng = Random.MersenneTwister(77163)

η = randn(rng, 5000)
S = ones(5000)
for n in 4:5000
    S[n] = 1.625 * S[n - 1] - 0.284 * S[n - 2] - 0.355 * S[n - 3] + η[n] - 0.96 * η[n - 1]
end
S ./= std(S)
S .-= 2 # for plotting (per Fig. 7.5 code)

plot!(1:1_000, S[1:1_000], label="ARMA", lc=RGB(96/255, 77/255, 158/255))

function surrogate(X; alg=:FT)
	if alg == :FT
		S = rfft(X)
		shuffle!(S)
		return irfft(S, length(X))
	elseif alg == :AAFT
		S = rfft(X)
		shuffle!(S)
		S = irfft(S, length(X))
		S[sortperm(S)] .= sort(X)
		return S
	else
		@error "Unknown algorithm: $alg"
	end
end

function takensbestestimator(X, ϵmax)
	η  = 0
	ij = [] # TODO ∀ i,j st. i < j && norm(X[i] - X[j]) < ϵmax

	for j in 1:length(X)
		for i in 1:(j-1)
			if norm(X[i] - X[j]) < ϵmax
				push!(ij, (i, j))
			end
		end
	end

	N = length(ij)

	for (i, j) in ij
		η += log(norm(X[i] - X[j]) / ϵmax)
	end

	η /= N
	return -1 / η
end

function delayembed(w::Vector{Float64}, τ::Int, d::Int)
	L = length(w) - (d-1)*τ
	masks = [[i + j*τ for j in 0:(d-1)] for i in 1:L]
	return [w[m] for m in masks]
end

#ee = 10 .^ (-5.0:0.1:1.0)
#CC = correlationsum(X, ee)
#plot(log.(ee), log.(CC))
X_ϵmax = 3.1622776601683795

#ee = 10 .^ (-5.0:0.1:1.0)
#CC = correlationsum(S, ee)
#plot(log.(ee), log.(CC))
S_ϵmax = 3.1622776601683795

CX = takensbestestimator(delayembed(X, 3, 15), X_ϵmax)
CS = takensbestestimator(delayembed(S, 4, 17), S_ϵmax)

box_Xft, box_Xaaft, box_Sft, box_Saaft = [], [], [], []
@showprogress for i in 1:100
	Xft = delayembed(surrogate(X[1:1_000]), 3, 15)
	Xaaft = delayembed(surrogate(X[1:1_000]; alg=:AAFT), 3, 15)
	Sft = delayembed(surrogate(S[1:1_000]), 4, 17)
	Saaft = delayembed(surrogate(S[1:1_000]; alg=:AAFT), 4, 17)

	push!(box_Xft, takensbestestimator(Xft, X_ϵmax))
	push!(box_Xaaft, takensbestestimator(Xaaft, X_ϵmax))
	push!(box_Sft, takensbestestimator(Sft, 0.25))
	push!(box_Saaft, takensbestestimator(Saaft, 0.25))
end

p2 = boxplot(
	box_Xft,
	label="FT",
	title="Rössler",
	xtick=false,
	ylabel="\$ Δ^{(T)}\$",
	c=:red;
	outliers=false
)
boxplot!(box_Xaaft, label="AAFT", c=:orange; outliers=false)
hline!([CX], lw=1, lc=:blue, ls=:dash, label=false)

p3 = boxplot(
	box_Sft,
	label="FT",
	title="ARMA",
	xtick=false,
	ylabel="\$ Δ^{(T)}\$",
	c=:red;
	outliers=false
)
boxplot!(box_Saaft, label="AAFT", c=:orange; outliers=false)
hline!([CS], lw=1, lc=:blue, ls=:dash, label=false)

plot(p1, p2, p3, layout=@layout[a{0.5w} b{0.25w} c{0.25w}], size=(800, 400))


