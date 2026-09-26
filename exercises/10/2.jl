using LinearAlgebra, StatsBase, Plots

function nodedegree(A, i)
	A[:, i] |> sum
end

function nodedegreedistribution(A)
	[nodedegree(A, i) for i in 1:first(size(A))]
end

function barabasialbert(N, n0, n)
	@assert N > n0 >= n
	A = zeros(Int, N, N)
	A[1:n0, 1:n0] = ones(Int, n0, n0)
	A .-= Diagonal(A)

	for i in (n0+1):N
		newly_connected_node_indexes = []
		for _ in 1:n
			potential_new_connected_nodes = setdiff(1:N, i, newly_connected_node_indexes)
			Pjs = [nodedegree(A, j) / sum([nodedegree(A, l) for l in potential_new_connected_nodes if l != j]) for j in potential_new_connected_nodes]
			Pjs = ifelse.(isinf.(Pjs), 1.0, Pjs)
			newly_connected_node_index = sample(potential_new_connected_nodes, Weights(Pjs), 1)[1]
			push!(newly_connected_node_indexes, newly_connected_node_index)
			A[i, newly_connected_node_index] = 1
			A[newly_connected_node_index, i] = 1
		end
	end
	A
end

coeffs = []
for n in 1:5
	local A = barabasialbert(100, 5, n)
	local ndd = nodedegreedistribution(A)

	local x = sort(unique(ndd))
	local y = [count(==(i), ndd) for i in x]
	local lx, ly = log.(x), log.(y)
	local X = [ones(length(lx)) lx]
	local β = X \ ly
	local a, b = exp(β[1]), β[2]

	push!(coeffs, abs(b))
end
p1 = plot(1:5, coeffs, marker=:circ, xlabel="\$ n\$", ylabel="Power-law coefficient")

A = barabasialbert(100, 5, 5)
ndd = nodedegreedistribution(A)
p2 = histogram(ndd, bins=minimum(ndd):maximum(ndd), label=false, title="\$ N=1000, n0=5, n=5\$")

x = sort(unique(ndd))
y = [count(==(i), ndd) for i in x]
lx, ly = log.(x), log.(y)
X = [ones(length(lx)) lx]
β = X \ ly
a, b = exp(β[1]), β[2]

xs = range(minimum(x), maximum(x), length=200)
plot!(xs, a .* xs .^ b, lw=2)

plot(p1, p2)
savefig("plots/10.2.png")
