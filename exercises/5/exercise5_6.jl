using LinearAlgebra, Plots

function next_line_iterate(v11, v12, center)
	dir = v12 .- v11                       # direction from v11 to v12
	len = norm(dir)                        # length of dir
	unit_dir = dir ./ len                  # unit vector in direction from v11 to v12
	unit_pr = (unit_dir[2], -unit_dir[1])  # unit vector perpindicular to unit_dir

	v21 = v11 .+ unit_dir .* (len/3)       # new vertex colinear with v11 & v12
	v22 = v11 .+ unit_dir .* (2*len/3)     # new vertex colinear with v11 & v12
	v2h = v11 .+ unit_dir .* (len/2)       # mid-point between v11 & v12

	v23 = v2h .+ unit_pr .* (len/(2*√3))   # new vertex non-colinear with v11 & v12

	return v21, v23, v22
end

function correlation_dimension(X::Vector{Tuple{Float64, Float64}}, ϵ::Real, w::Int)
	N = length(X)
	C = 0

	for i in 1:(N - w - 1)
		for j in (1 + w + i):N
			C += ifelse(norm(X[i] .- X[j]) < ϵ,1,0)
		end
	end

	CS = (2*C)/((N - w)*(N - w - 1)) # correlation sum
	return log(CS)/log(ϵ)
end

function next_koch_snowflake_step(vertices)
	center = (sum(first.(vertices))/length(vertices), sum(last.(vertices))/length(vertices))
	new_vertices = []
	for (i,v) in enumerate(vertices)
		push!(new_vertices, vertices[i])
		j = ifelse(i<length(vertices),i+1,1)
		nv = collect(next_line_iterate(vertices[i], vertices[j], center))
		push!(new_vertices, nv...)
	end
	return new_vertices
end

steps = [[1/√3 .* (cos(π/2 + θ), sin(π/2 + θ)) for θ in [0, 2*π/3, 4*π/3]]]
max_step = 8

println("0/$max_step Koch step")
for i in 1:max_step
	println("$i/$max_step Koch step")
	push!(steps, next_koch_snowflake_step(steps[end]))
end

dims = []
for (i,step) in enumerate(steps)
	println("$(i-1)/$(length(steps)-1) Dimension step")
	push!(dims, correlation_dimension(step, 0.1, 3))
end

dims = [ifelse(isnan(x) || isinf(x), 1, x) for x in dims]
println("Dim ≈ $(round(dims[end]; digits=2)) ($(max_step+1) steps)")

display_step = 3
ks = steps[display_step+1]
ks = vcat(ks, ks[1])
plot(
	first.(ks), 
	last.(ks), 
	xlim=[-0.55,0.55],
	lc=:blue, 
	aspect_ratio=:equal,
	size=(500,500),
	grid=false,
	xaxis=false,
	yaxis=false,
	xtick=false,
	ytick=false,
	label="Step $(display_step), \$Δ^{(C)} ≈ $(round(dims[display_step+1]; digits=2))\$",
	title="Koch snowflake",
	legendposition=:inside,
	foregroundcolorlegend = nothing
)
savefig("plots/exercise5_6.png")
