using LinearAlgebra, Plots

function erdosrenyi(n, p)
	A = rand(n, n) .>= p
	A = (A |> LowerTriangular) .+ transpose(A |> LowerTriangular)
	A .-= Diagonal(A)
end

function barabasialbert(n, n0)
	@assert n > n0
	A = zeros(Int, n, n)
	A[1:n0, 1:n0] = ones(Int, n0, n0)
	A .-= Diagonal(A)

	for i in (n0+1):n
		for j in 1:n
			if i != j
				Pj = nodedegree(A, j) / sum([nodedegree(A, l) for l in 1:n])
				if Pj > rand()
					A[i, j] = 1
					A[j, i] = 1
				end
			end
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

function nodedegreedistribution(A)
	[nodedegree(A, i) for i in 1:first(size(A))]
end

A = erdosrenyi(1000, 0.5)
ndd = nodedegreedistribution(A)
p1 = histogram(
	ndd, 
	bins=minimum(ndd):maximum(ndd),
	lw=0,
	label=false,
	title="Erdös-Rényi",
	normalize=true,
	xlabel="\$ k\$",
	ylabel="\$ P(k)\$"
)

A = barabasialbert(100, 4)
ndd = nodedegreedistribution(A)
p2 = histogram(
	ndd, 
	bins=minimum(ndd):maximum(ndd),
	lw=0,
	label=false,
	title="Barabási-Albert",
	normalize=true,
	xlabel="\$ k\$",
	ylabel="\$ P(k)\$"
)

A = wattsstrogatz(1000, 2, 0.1)
ndd = nodedegreedistribution(A)
p3 = histogram(
	ndd, 
	bins=minimum(ndd):maximum(ndd),
	lw=0,
	label=false,
	title="Watts-Strogatz",
	normalize=true,
	xlabel="\$ k\$",
	ylabel="\$ P(k)\$"
)

plot(p1, p2, p3, layout=(1, 3))
savefig("plots/10.1.png")
