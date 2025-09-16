using LinearAlgebra, Plots

function standard_map(x)
	θ, v= x
	k  = 0.6
	(θ + v + k*sin(θ), v + k*sin(θ)) .% 2π
end

function largest_lyapunov_exponent(f, x0, δ0, Δt, N)
	xref  = [x0]
	xtest = [x0 .+ δ0]
	δi    = []

	for i=1:(N*Δt)
		push!(xref, f(xref[end]))

		if i % Δt == 0
			push!(xtest, f(xtest[end]))
			push!(δi, norm(xref[end] .- xtest[end]))
			xtest[end] = xref[end] .+ δ0
		else
			push!(xtest, f(xtest[end]))
		end
	end
	sum(log.(δi ./ norm(δ0))) / (N * Δt)
end

λ1 = largest_lyapunov_exponent(standard_map, (0.1, 0.11), (0.0001, 0.0001), 50, 10000)

xrefs = [(0.1, 0.11)]
xtests = [(0.1, 0.11) .+ (0.0001, 0.0001)]
δis = [norm((0.0001, 0.0001))]

for _=1:200
	push!(xrefs, standard_map(xrefs[end]))
	push!(xtests, standard_map(xtests[end]))
	push!(δis, norm(xrefs[end] .- xtests[end]))
end

plot(
	1:201, log.(δis),
	title="Largest Lyapunov Exponent",
	xlabel="t",
	label="ln(δ(t))",
	linecolor=:deepskyblue
)
plot!(
	1:50, λ1 .* collect(1:50) .+ -8.86376678169621,
	label="λ1 = $(round(λ1, digits=2))",
	linecolor=:mediumpurple
)
savefig("plots/exercise3_1.png")
