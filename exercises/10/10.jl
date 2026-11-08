using DifferentialEquations
using Plots

function I(t, N, β, γ, I0)
	c1 = N - ((N * γ) / β)
	c2 = β - γ
	c3 = (c1 / I0) - 1

	return (c1 * exp(c2 * t)) / (c3 + exp(c2 * t))
end

function SI!(du, u, p, t)
	β, γ, N = p
	S, I = u

	du[1] = (-β * S * I)/N + (γ * I)
	du[2] = ( β * S * I)/N - (γ * I)
end


I0 = 1
N = 10_000
β = 0.3
γ = 0.1

ts = range(0.0, 150.0; length=250)
Is = zeros(length(ts))
Is[1] = I0

for i in 2:length(ts)
	Is[i] = I(ts[i], N, β, γ, I0)
end

u0 = [9999.0, 1.0]
prob = ODEProblem(SI!, u0, (0.0, 150.0), [β, γ, N])
solv = solve(prob)

plot(
	ts, 
	Is, 
	label="\$ I\$ (analytic)", 
	lc=:red,
	xlabel="\$ t\$"
)

scatter!(solv.t, last.(solv.u), label="\$ I\$ (numerical)", mc=:red, ms=2, msw=0)

plot!(ts, N .- Is, label="\$ S\$ (analytic)", lc=:blue)

scatter!(solv.t, first.(solv.u), label="\$ S\$ (numerical)", mc=:blue, ms=2, msw=0)

hline!([N * (1 - (γ/β))], label="\$ I\$ asymptote", ls=:dash, lc=:black)
savefig("plots/10.10.png")

# if γ=0, I(t) → N, as t → ∞
# İ = 0 if I = 0 or I = N(1 - γ/β)
# if γ > β, I(t) → 0, as t → ∞ (since I > 0)
# if γ < β, I(t) → N(1 - γ/β), as t → ∞
