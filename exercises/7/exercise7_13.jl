using DifferentialEquations
using Plots, Measures
using LinearAlgebra

function lorenz96!(du, u, p, t)
	D, F = Int(p[1]), p[2]

    for i in 1:D
    	du[i] = ((u[mod1(i+1, D)] - u[mod1(i-2, D)]) * u[mod1(i-1, D)]) - u[i] + F
    end
end

function delay_embed(w::Vector{Float64}, τ::Int, d::Int)
	L = length(w) - (d-1)*τ
	masks = [[i + j*τ for j in 0:(d-1)] for i in 1:L]
	return [w[m] for m in masks]
end

function nearest_neighbors(x, X; k=1, w=5)
	i = findfirst(==(x), X)
	if i == nothing
		norms = [(j, norm(x - y)) for (j, y) in enumerate(X)]
		return [i for i in first.(sort(norms; by=x->x[2])[1:k]) if i < length(X)]
	else
		norms = [(j, norm(x - y)) for (j, y) in enumerate(X)]
		window = maximum([1, i-w]):minimum([length(X), i+w])
		norms[window] .= [(typemax(Int), Inf) for _ in window]
		return [i for i in first.(sort(norms; by=x->x[2])[1:k]) if i < length(X)]
	end
end

function nearest_neighbors_forecasting(x, R; k=1, w=5)
	return mean(R[nearest_neighbors(x, R; k=k, w=w) .+ 1])
end

function MSE(x, y)
	return sum((x .- y) .^ 2) / length(x)
end

function NRMSE(x, y, x̄)
	return sqrt(MSE(x, y) / MSE(x, x̄))
end

u0 = [0.12013623147677799, 0.274707399757551, 0.4738961112067477, 0.14933883808016557, 0.4274098963989116] # rand(5)
prob = ODEProblem(lorenz96!, u0, (0.0, 100.0), [5.0, 8.0])
solv = solve(prob)

#plot(solv, idxs=(1,2,3))

function ragwitz_criterion(x, i)
	mean_NRMSEs = zeros(length(4:10), length(1:20))
	for d in 4:10
		for τ in 1:20
			h = getindex.(x, i) 
			R = delay_embed(h, τ, d)

			ground_truth = R[(end-99):end]
			prediction = similar(R, 100)
			for j in 1:100
				q = j == 1 ? R[end - 100] : prediction[j - 1]
				q̃ = nearest_neighbors_forecasting(q, R; k=5) 
				push!(h, first(q̃))
				R = delay_embed(h, τ, d)
				prediction[j] = q̃
			end

			NRMSEs = zeros(100)
			for j in 1:100
				NRMSEs[j] = NRMSE(ground_truth[j], prediction[j], mean(ground_truth))
			end
			mean_NRMSEs[d - 3, τ] = mean(NRMSEs)
		end
	end
	return mean_NRMSEs
end

mn = ragwitz_criterion(solv.u, 1)
p1 = heatmap(
	mn, 
	xlabel="\$ τ\$", 
	ylabel="\$ d\$", 
	colorbar_title="\nMean NRMSE", 
	ytick=(1:7, 4:10), margin=5mm,
	title="\nLorenz96 \$ x_1\$"
)

mn = ragwitz_criterion(solv.u, 2)
p2 = heatmap(
	mn, 
	xlabel="\$ τ\$", 
	ylabel="\$ d\$", 
	colorbar_title="\nMean NRMSE", 
	ytick=(1:7, 4:10), margin=5mm,
	title="\nLorenz96 \$ x_2\$"
)

mn = ragwitz_criterion(solv.u, 3)
p3 = heatmap(
	mn, 
	xlabel="\$ τ\$", 
	ylabel="\$ d\$", 
	colorbar_title="\nMean NRMSE", 
	ytick=(1:7, 4:10), margin=5mm,
	title="\nLorenz96 \$ x_3\$"
)

plot(p1, p2, p3, layout=(3, 1), size=(400, 800))
savefig("plots/exercise7_13.png")

# TODO use println("n_c = $(findfirst(>=(0.5), NRMSEs))") to use n_c not mean(NRMSEs) to rate pair of params
