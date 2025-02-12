using Plots
using LinearAlgebra

function fermi_ulam(u, θ, ϵ)
	return (abs(u + ϵ * sin(θ)), mod(θ + 2π / abs(u + ϵ * sin(θ)),2π))
end

function n_iterations(u, θ, ϵ, n)
	last = (u, θ)
	for step=1:n
		next = fermi_ulam(last..., ϵ)
		last = next
	end
	return last
end

function iterate_points(points, ϵ, n)
	iterated_points = zeros(Float64,length(points))
	Threads.@threads for i=1:length(points)
		iterated_points[i] = norm(n_iterations(points[i]...,ϵ,n)...)
	end
	return iterated_points
end

us = 0.0:0.01:2π
θs = 0.0:0.01:2π

epsilon::Float64 = 1.0
iters::Int64 = 100

points = [(u,θ) for u=us for θ=θs]
iterated_points = iterate_points(points, epsilon, iters)

heatmap(
	us, θs, reshape(iterated_points, length(us), length(θs))', 
	title="Fermi-Ulam map, \$ ϵ=$epsilon\$, \$||(u_{$iters},θ_{$iters})||\$",
	c=cgrad(:managua),
	xlim=[0,2π],
	ylim=[0,2π],
	xguide="\$u_0\$",
	yguide="\$θ_0\$"
)

savefig("plots/exercise1_11.png")