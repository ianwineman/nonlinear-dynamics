using Plots

function f(θ; Ω=0.3, K=1.0)
	θ + Ω - (K/2π)*sin(2π*θ)
end

function nonlinearcirclemap(θ; Ω=0.3, K=1.0)
	return mod(θ + Ω - (K/2π)*sin(2π*θ), 1)
end

function trajectory(θ1, n, f; Ω=0.3, K=1.0)
	θ = zeros(n)
	θ[1] = θ1
	for i in 2:length(θ)
		θ[i] = f(θ[i-1]; Ω=Ω, K=K)
	end
	return θ
end

function W(θ1, Ω, K; n=10_000)
	θ = trajectory(θ1, n, f; Ω=Ω, K=K)
	return (θ[end] - θ[1]) / n
end

function orbitdiagram(θ1, Ωr, K; n=1000, t=100)
	points = []
	for Ω in Ωr
		orbitp = trajectory(θ1, n, nonlinearcirclemap; Ω=Ω, K=K)[end-t:end]
		for point in orbitp
			push!(points, (Ω, point))
		end
	end
	return points
end

Ws = []
Ωs = 0.0:0.001:2.0
for Ω in Ωs
	w1 = W(0.7, Ω, 0.1)
	w2 = W(0.7, Ω, 0.5)
	w3 = W(0.7, Ω, 1.0)
	push!(Ws, (w1, w2, w3))
end

p1 = plot(
	Ωs, 
	getindex.(Ws, 1), 
	label="\$ K=0.1\$",
	xlabel="\$ Ω\$",
	ylabel="\$ W\$",
	lc=:blue
)
plot!(Ωs, getindex.(Ws, 2), label="\$ K=0.5\$", lc=:red)
plot!(Ωs, getindex.(Ws, 3), label="\$ K=1.0\$", lc=:green)

od = orbitdiagram(0.7, Ωs, 0.1; n=10_000)
p2 = scatter(first.(od), last.(od), msw=0, ms=0.1, mc=:blue, label=false, xlabel="\$ Ω\$", ylabel="\$ θ\$")

od = orbitdiagram(0.7, Ωs, 0.5; n=10_000)
p3 = scatter(first.(od), last.(od), msw=0, ms=0.1, mc=:red, label=false, xlabel="\$ Ω\$")

od = orbitdiagram(0.7, Ωs, 1.0; n=10_000)
p4 = scatter(first.(od), last.(od), msw=0, ms=0.1, mc=:green, label=false, xlabel="\$ Ω\$")

plot(p1, p2, p3, p4, layout=@layout [a; b c d])
savefig("plots/9.10.png")
