using DifferentialEquations, LinearAlgebra
using Plots

function correlation_sum(X, ϵ, w)
	N = length(X)
	C = 0

	for i in 1:(N - w - 1)
		for j in (1 + w + i):N
			C += ifelse(norm(X[i] .- X[j]) < ϵ, 1, 0)
		end
	end

	return (2*C)/((N - w)*(N - w - 1))
end


Ds = 4:10
Δs = zeros(length(Ds))
for D in Ds
	println("D=$D")
	function lorenz96!(du, u, p, t)
		F = p

	    for i in 1:D
	    	du[i] = ((u[mod1(i+1, D)] - u[mod1(i-2, D)]) * u[mod1(i-1, D)]) - u[i] + F
	    end
	end

	prob = ODEProblem(lorenz96!, rand(-10:10, D), (0.0,100.0), 24.0)
	solv = solve(prob)

	ϵs = [exp(i) for i in -10.0:0.1:5.0]
	Cs = [correlation_sum(solv.u, ϵ, 3) for ϵ in ϵs]

	i = findfirst(x -> !isinf(x), log.(Cs)) + round(Int, length(Cs)/10)
	j = findfirst(x -> x == 0.0, log.(Cs))  - round(Int, length(Cs)/10)

	xstart = ϵs[i]
	xend = ϵs[j]
	Δ = (log.(Cs)[j] - log.(Cs)[i]) / (log.(ϵs)[j] - log.(ϵs)[i])
	Δs[D - first(Ds) + 1] = Δ
end

scatter(
   Ds, Δs,
   xlabel="Lorenz-96 dimension, \$D\$",
   ylabel="Fractal dimension, \$Δ\$",
   label=false,
   ms=4, mc=:purple1, msw=0,
   background_color=RGB(255/255, 248/255, 231/255),
   gridalpha=0.2
)
plot!(
   Ds, (Δs[end] - Δs[1])/(Ds[end] - Ds[1]) .* Ds,
   lw=2,
   lc=:lightblue1,
   label="\$Δ/D\$ \$≈\$ \$$(round((Δs[end] - Δs[1])/(Ds[end] - Ds[1]); digits=2))\$"
)
savefig("plots/exercise5_10.png")
