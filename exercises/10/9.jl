using DifferentialEquations, Plots

function SIR!(du, u, p, t)
	S, I, R = u
	β, γ = p

	N = S + I + R

	if 15.0 < t < 25.0
		β /= 2
	end

	du[1] = -(β * S * I) / N
	du[2] =  (β * S * I) / N - γ*I
	du[3] = γ*I
end

prob = ODEProblem(SIR!, [99.0, 1.0, 0.0], (0.0, 100.0), [0.3, 0.1])
solv = solve(prob)

plot(solv, xtick=vcat(0:20:100, 15, 25), label=["S" "I" "R"])
vspan!([15.0, 25.0], color=:orange, alpha=0.2, label="quarantine \$ β → β/2\$") 
savefig("plots/10.9.png")
