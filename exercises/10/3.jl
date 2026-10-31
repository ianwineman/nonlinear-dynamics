using Plots

function wattsstrogatz(n, K, β)	
	edges = Tuple{Int, Int}[]
	for i in 1:n
		for j in (i-K):(i+K)
			if i != j
				j = mod(j, 1:n)
				if j < i
					push!(edges, (i, j) |> sort)
					# sort means all edges are stored (small, large)
				end
			end
		end
	end

	for (e, edge) in enumerate(edges)
		if β > rand()
			new_edge = nothing
			while true # rejection sampling
				a = rand(1:n)
				b = rand(1:n)
				a == b && continue
				sort((a, b)) in edges && continue
				new_edge = (a, b) |> sort # ensure (small, large)
				break
			end
			edges[e] = new_edge # replace old edge with new edge
		end
	end
	edges

	A = zeros(Int, n, n)
	for edge in edges
		i, j = edge
		A[i, j] = 1
		A[j, i] = 1
	end
	A
end

function clusteringcoefficient(A)
	n = first(size(A))
	C = 0
	for i in 1:n
		Ci = 0
		for j in 1:n
			for ξ in 1:n
				Ci += A[i, j] * A[j, ξ] * A[ξ, i]
			end
		end
		k = sum([A[i, l] for l in 1:n])
		if k > 1
			Ci /= (k * (k - 1)) #k < 2 ? 0.0 : Ci / (k * (k - 1))
			C += Ci
		end
	end
	return C / n
end

function averagepathlength(A)
	n = first(size(A))
	ds = []
	for i in 1:n
		for j in 1:n
			if i != j
				# i → j
				paths = map(x -> [x], findall(==(1), A[i, :]))
				while true
					complete_paths = findall(x -> x[end] == j, paths)
					if isempty(complete_paths)
						new_paths = []
						for path in paths
							path_end = path[end]
							neighbors = findall(==(1), A[path_end, :])
							for neighbor in neighbors
								if neighbor ∉ path
									new_path = vcat(path, neighbor)
									push!(new_paths, new_path)
								end
							end
						end

						if isempty(new_paths)
							break
						else
							paths = new_paths
						end
					else
						d = first(length.(paths[complete_paths]))
						push!(ds, d)
						break
					end
				end
			end
		end
	end
	return sum(ds) / (n * (n - 1))
end


βs = [10^x for x in range(-2.0, 0.0; length=200)]
Cs = [clusteringcoefficient(wattsstrogatz(100, 5, β)) for β in βs]

p1 = plot(
	βs, 
	Cs, 
	label="\$ K=5\$",
	xlabel="\$ β\$",
	ylabel="\$ C\$",
	xscale=:log10
)

Cs = [clusteringcoefficient(wattsstrogatz(100, 10, β)) for β in βs]
plot!(βs, Cs, label="\$ K=10\$", xscale=:log10)

Cs = [clusteringcoefficient(wattsstrogatz(100, 20, β)) for β in βs]
plot!(βs, Cs, label="\$ K=20\$", xscale=:log10)

ℓs = [averagepathlength(wattsstrogatz(10, 1, β)) for β in βs]
p2 = plot(
	βs, 
	ℓs, 
	label="\$ K=1\$",
	xlabel="\$ β\$",
	ylabel="\$ ℓ\$",
	xscale=:log10
)

ℓs = [averagepathlength(wattsstrogatz(10, 2, β)) for β in βs]
plot!(βs, ℓs, label="\$ K=2\$")

ℓs = [averagepathlength(wattsstrogatz(10, 3, β)) for β in βs]
plot!(βs, ℓs, label="\$ K=3\$")

plot(p1, p2)
savefig("plots/10.3.png")
