using DifferentialEquations
using Plots

function SIR!(du, u, p, t)
	β, γ, N = p
	S, I, R = u

	du[1] = (-β * S * I)/N
	du[2] = ( β * S * I)/N - (γ * I)
	du[3] = γ * I
end

N = 10_000
β = 0.3
γ = 0.1
S0 = 9999.0
u0 = [S0, N - S0, 0.0]

prob = ODEProblem(SIR!, u0, (0.0, 150.0), [β, γ, N])
solv = solve(prob)

plot(solv, label=["S" "I" "R"], xlabel="\$ t\$")

C = N - (((γ*N)/β) * log(S0)) # integration constant
Smax = (γ*N)/β
Imax = Smax * log(Smax) - Smax + C

hline!([Imax], label="\$ I_{max} = $(round(Imax; digits=3))\$", ls=:dash, lc=:black)
savefig("plots/10.11.png")
