using Plots; ENV["GKSwstype"] = "100"


function initialize(M, N)
	@assert 2*N < M^2
	
	A = rand(M, M)
	r1 = sort(vec(A))[N]
	r2 = sort(vec(A))[2*N]
	A = map(x -> ifelse(       x <= r1,  1.0, x), A)
	A = map(x -> ifelse(       x <= r2, -1.0, x), A)
	A = map(x -> ifelse(-1.0 < x <  1.0, 0.0, x), A)
	A
end

function update!(A, k)
	n = size(A, 1)
	for i in 1:n
		for j in 1:n
			if A[i, j] != 0.0
				neighbors = [get(A, (i + x, j + y), 0.0) for x in -1:1 for y in -1:1 if x != 0.0 || y != 0.0]
				same_count = count(==(A[i, j]), neighbors)
				if same_count < k
					empty_cells = findall(iszero, A)
					u, v = Tuple(empty_cells[rand(1:length(empty_cells))])
					A[u, v] = A[i, j]
					A[i, j] = 0.0
				end
			end
		end
	end
	return nothing
end

A = initialize(20, 100)

anim = @animate for i in 1:25
	println(i)
	update!(A, 3)
	heatmap(A, c = cgrad([:orange, :black, :green], 3, categorical=true), title="$i")
end

gif(anim, "plots/10.15.gif", fps=5)
