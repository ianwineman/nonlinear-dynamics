using CSV, DataFrames
using LinearAlgebra
using Plots
using Base.Threads
using RecurrenceAnalysis

function recurrencematrix(X, ϵ; n=length(X))
	R = zeros(Int, n, n)
	@threads for i in 1:n
		for j in 1:n
			if norm(X[i] - X[j]) < ϵ
				R[i, j] = 1
			end
		end
	end
	return R
end

df = CSV.read("book/exercise_data/16.csv", DataFrame; header=false, delim=" ", ignorerepeated=true)

ϵs = range(0.1, 0.2; length=20)
n = 2_000

dlxs = []
dlys = []
dlzs = []

vlxs = []
vlys = []
vlzs = []

anim = @animate for ϵ in ϵs
	"ϵ=$ϵ" |> println

	x, y, z = df[:, 2], df[:, 3], df[:, 5]
	x ./= maximum(x)
	y ./= maximum(y)
	z ./= maximum(z)
	
	Rx = recurrencematrix(x, ϵ; n=n)
	p1 = heatmap(
		Rx,
		cbar=false,
		c=cgrad(:devon, 3, categorical=true)[1:2],
		aspect_ratio=:equal,
		xlim=[0, n],
		ylim=[0, n],
		xtick=false,
		ytick=false
	)
	dlx = dl_entropy(RecurrenceMatrix(x, ϵ))
	push!(dlxs, (ϵ, dlx))
	vlx = vl_entropy(RecurrenceMatrix(x, ϵ))
	push!(vlxs, (ϵ, vlx))
	p4 = plot(
		first.(dlxs), 
		last.(dlxs), 
		xlim=[0.1, 0.2],
		ylim=[0, 10],
		label="\$ H_ℓ\$",
		legendposition=:topleft,
		xlabel="\$ ϵ\$",
		lc=cgrad(:devon, 3, categorical=true)[1],
		xtick=[0.1, 0.2]
	)
	plot!(
		first.(vlxs), 
		last.(vlxs),
		label="\$ H_r\$",
		lc=cgrad(:devon, 3, categorical=true)[2]
	)

	Ry = recurrencematrix(y, ϵ; n=n)
	p2 = heatmap(
		Ry,
		cbar=false,
		c=cgrad(:devon, 3, categorical=true)[1:2],
		aspect_ratio=:equal,
		xlim=[0, n],
		ylim=[0, n],
		xtick=false,
		ytick=false
	)
	dly = dl_entropy(RecurrenceMatrix(y, ϵ))
	push!(dlys, (ϵ, dly))
	vly = vl_entropy(RecurrenceMatrix(y, ϵ))
	push!(vlys, (ϵ, vly))
	p5 = plot(
		first.(dlys), 
		last.(dlys), 
		xlim=[0.1, 0.2],
		ylim=[0, 10],
		label="\$ H_ℓ\$",
		legendposition=:topleft,
		xlabel="\$ ϵ\$",
		lc=cgrad(:devon, 3, categorical=true)[1],
		xtick=[0.1, 0.2]
	)
	plot!(
		first.(vlys), 
		last.(vlys),
		label="\$ H_r\$",
		lc=cgrad(:devon, 3, categorical=true)[2]
	)

	Rz = recurrencematrix(z, ϵ; n=n)
	p3 = heatmap(
		Rz,
		cbar=false,
		c=cgrad(:devon, 3, categorical=true)[1:2],
		aspect_ratio=:equal,
		xlim=[0, n],
		ylim=[0, n],
		xtick=false,
		ytick=false
	)
	dlz = dl_entropy(RecurrenceMatrix(z, ϵ))
	push!(dlzs, (ϵ, dlz))
	vlz = vl_entropy(RecurrenceMatrix(z, ϵ))
	push!(vlzs, (ϵ, vlz))
	p6 = plot(
		first.(dlzs), 
		last.(dlzs), 
		xlim=[0.1, 0.2],
		ylim=[0, 10],
		label="\$ H_ℓ\$",
		legendposition=:topleft,
		xlabel="\$ ϵ\$",
		lc=cgrad(:devon, 3, categorical=true)[1],
		xtick=[0.1, 0.2]
	)
	plot!(
		first.(vlzs), 
		last.(vlzs),
		label="\$ H_r\$",
		lc=cgrad(:devon, 3, categorical=true)[2]
	)

	plot(p1, p2, p3, p4, p5, p6, size=(900, 600))
end
gif(anim, "plots/exercise8_12.gif", fps=5)
