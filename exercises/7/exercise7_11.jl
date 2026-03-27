using DifferentialEquations
using Statistics, Random
using Plots
using FFTW
using LinearAlgebra
using DynamicalSystems

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

p1 = plot(1:1_000, X[1:1_000], xlim=[0, 1_000], label="Rössler")

Random.seed!(77163)
rng = Random.MersenneTwister(77163)

η = randn(rng, 5000)
S = ones(5000)
for n in 4:5000
    S[n] = 1.625 * S[n - 1] - 0.284 * S[n - 2] - 0.355 * S[n - 3] + η[n] - 0.96 * η[n - 1]
end
S ./= std(S)
S .-= 2 # for plotting (per Fig. 7.5 code)

plot!(1:1_000, S[1:1_000], label="ARMA")

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

CX, CS
