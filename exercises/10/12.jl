using LinearAlgebra, Plots, Statistics, StatsBase

function erdosrenyi(n, p)
	A = rand(n, n) .<= p
	A = (A |> LowerTriangular) .+ transpose(A |> LowerTriangular)
	A .-= Diagonal(A)
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

function wattsstrogatz(n, K, β)
	A = zeros(Int, n, n)

	for i in 1:n
		for j in (i-K+1):(i+K-1)
			if i != mod(j, n)+1
				A[i, mod(j, n)+1] = 1
				A[mod(j, n)+1, i] = 1
			end
		end
	end

	edges = findall(==(1), A) .|> Tuple
	uniqueedges = sort(edges, by=z -> z[1]^2 + z[2]^2)[1:2:end]
	for ue in uniqueedges
		i, j = ue
		if rand() < β
			A[i, j] = 0
			A[j, i] = 0

			k = rand(setdiff(1:n, [i, j]))
			A[i, k] = 1
			A[k, i] = 1
		end
	end
	A
end

function nodedegree(A, i)
	A[:, i] |> sum
end

@enum State susceptible infectious recovered

mutable struct Agent
	state::State
	days_infected::Int
end

struct Graph
	nodes::Vector{Agent}
	adjacency::Matrix{Int}
end

function Graph(A)
	G = Graph([Agent(susceptible, 0) for _ in 1:(size(A) |> first)], A)
	i = rand(1:(size(A) |> first))
	G.nodes[i].state = infectious
	G
end

function step!(G::Graph, τ, p)
	I = findall(a -> a.state == infectious, G.nodes)
	for i in I
		a = G.nodes[i]
		if a.days_infected >= p
			a.state = recovered
			a.days_infected = 0
		else
			neighbors = findall(==(1), G.adjacency[i, :])
			for i in neighbors
				n = G.nodes[i]
				if (n.state == susceptible) && (τ > rand())
					n.state = infectious
				end
			end
		end
	end
	I = findall(a -> a.state == infectious, G.nodes)
	for i in I
		a = G.nodes[i]
		a.days_infected += 1
	end
	return nothing
end

function infectiouscount(G)
	count(a -> a.state == infectious, G.nodes)
end

τ = 0.1
p = 14
N = 100
trials = 100

Is = []
for i in 1:trials
	i |> println

	A = erdosrenyi(100, 0.1)
	G = Graph(A)
	I = zeros(Int, N)
	I[1] = 1
	for j in 2:N
		step!(G, τ, p)
		I[j] = G |> infectiouscount
	end
	push!(Is, I)
end
Is = Is |> stack
means = [mean(Is[i, :]) for i in 1:size(Is)[1]]
σs = [std(Is[i, :]) for i in 1:size(Is)[1]]

plot(
	1:N, 
	means, 
	ribbon=σs, 
	fillalpha=0.3, 
	lw=2,
	xlim=[0, 100],
	label="ER",
	xlabel="\$ t\$",
	ylabel="infectious agents (count)"
)

Is = []
for i in 1:trials
	i |> println

	A = barabasialbert(100, 3, 2)
	G = Graph(A)
	I = zeros(Int, N)
	I[1] = 1
	for j in 2:N
		step!(G, τ, p)
		I[j] = G |> infectiouscount
	end
	push!(Is, I)
end
Is = Is |> stack
means = [mean(Is[i, :]) for i in 1:size(Is)[1]]
σs = [std(Is[i, :]) for i in 1:size(Is)[1]]

plot!(
	1:N, 
	means, 
	ribbon=σs, 
	fillalpha=0.3, 
	lw=2,
	xlim=[0, 100],
	label="BA"
)

Is = []
for i in 1:trials
	i |> println

	A = wattsstrogatz(100, 2, 0.1)
	G = Graph(A)
	I = zeros(Int, N)
	I[1] = 1
	for j in 2:N
		step!(G, τ, p)
		I[j] = G |> infectiouscount
	end
	push!(Is, I)
end
Is = Is |> stack
means = [mean(Is[i, :]) for i in 1:size(Is)[1]]
σs = [std(Is[i, :]) for i in 1:size(Is)[1]]

plot!(
	1:N, 
	means, 
	ribbon=σs, 
	fillalpha=0.3, 
	lw=2,
	xlim=[0, 100],
	label="WS"
)
savefig("plots/10.12.png")
