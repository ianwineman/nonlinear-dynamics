using CSV, DataFrames
using LinearAlgebra
using Plots
using Base.Threads
using DynamicalSystems
using Images
using Statistics
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

function recurrencematrix(X; q=0.1, n=length(X))
	@assert 0.0 < q < 1.0
	R = zeros(n, n)
	@threads for i in 1:n
		for j in 1:n
			R[i, j] = norm(X[i] - X[j])
		end
	end

	qu = quantile(vec(R), q)
	Rqu = zeros(Int, n, n)

	@threads for i in 1:n
		for j in 1:n
			if R[i, j] <= qu
               Rqu[i, j] = 1
           end
		end
	end
	return Rqu
end

df = CSV.read("book/exercise_data/16.csv", DataFrame; header=false, delim=" ", ignorerepeated=true)

x, y, z = df[:, 2], df[:, 3], df[:, 5]
x ./= maximum(x)
y ./= maximum(y)
z ./= maximum(z)

τ = 50 #plot(1:100, [selfmutualinfo(x, 1:100), selfmutualinfo(y, 1:100), selfmutualinfo(z, 1:100)], xtick=1:5:100, label=["x" "y" "z"])
d = 4  #plot(1:10, [delay_afnn(x, τ, 1:10; w=5), delay_afnn(y, τ, 1:10; w=5), delay_afnn(z, τ, 1:10; w=5)], label=["x" "y" "z"])

dx = embed(x, d, τ)
dy = embed(y, d, τ)
dz = embed(z, d, τ)

p_x_1 = recurrencematrix(x[1:length(dx)])
p_x_d = recurrencematrix(dx)

p_y_1 = recurrencematrix(y[1:length(dy)])
p_y_d = recurrencematrix(dy)

p_z_1 = recurrencematrix(z[1:length(dz)])
p_z_d = recurrencematrix(dz)

#save("plots/p_x_1.png", p_x_1 |> rotr90 .|> Float64)
#save("plots/p_x_d.png", p_x_d |> rotr90 .|> Float64)

#save("plots/p_y_1.png", p_y_1 |> rotr90 .|> Float64)
#save("plots/p_y_d.png", p_y_d |> rotr90 .|> Float64)

#save("plots/p_z_1.png", p_z_1 |> rotr90 .|> Float64)
#save("plots/p_z_d.png", p_z_d |> rotr90 .|> Float64)

####

function plot_R(R; ℓ=(size(R) |> first), title="")
	heatmap(
		R[1:ℓ, 1:ℓ],
		cbar=false,
		showaxis=false,
		aspect_ratio=:equal,
		c=cgrad(:devon, 3, categorical=true)[1:2],
		grid=false,
		title=title
	)
end

ell = 3_000
p1 = plot_R(p_x_1; ℓ=ell, title="\$ x\$ (original) \n\$ H_{diag}=$(round(RecurrenceMatrix(x, GlobalRecurrenceRate(0.1)) |> dl_entropy; digits=2))\$ \$ H_{vert}=$(round(RecurrenceMatrix(x, GlobalRecurrenceRate(0.1)) |> vl_entropy; digits=2))\$")
p2 = plot_R(p_x_d; ℓ=ell, title="\$ x\$ (embedded) \n\$ H_{diag}=$(round(RecurrenceMatrix(dx, GlobalRecurrenceRate(0.1)) |> dl_entropy; digits=2))\$ \$ H_{vert}=$(round(RecurrenceMatrix(dx, GlobalRecurrenceRate(0.1)) |> vl_entropy; digits=2))\$")

p3 = plot_R(p_y_1; ℓ=ell, title="\$ y\$ (original) \n\$ H_{diag}=$(round(RecurrenceMatrix(y, GlobalRecurrenceRate(0.1)) |> dl_entropy; digits=2))\$ \$ H_{vert}=$(round(RecurrenceMatrix(y, GlobalRecurrenceRate(0.1)) |> vl_entropy; digits=2))\$")
p4 = plot_R(p_y_d; ℓ=ell, title="\$ y\$ (embedded) \n\$ H_{diag}=$(round(RecurrenceMatrix(dy, GlobalRecurrenceRate(0.1)) |> dl_entropy; digits=2))\$ \$ H_{vert}=$(round(RecurrenceMatrix(dy, GlobalRecurrenceRate(0.1)) |> vl_entropy; digits=2))\$")

p5 = plot_R(p_z_1; ℓ=ell, title="\$ z\$ (original) \n\$ H_{diag}=$(round(RecurrenceMatrix(z, GlobalRecurrenceRate(0.1)) |> dl_entropy; digits=2))\$ \$ H_{vert}=$(round(RecurrenceMatrix(z, GlobalRecurrenceRate(0.1)) |> vl_entropy; digits=2))\$")
p6 = plot_R(p_z_d; ℓ=ell, title="\$ z\$ (embedded) \n\$ H_{diag}=$(round(RecurrenceMatrix(dz, GlobalRecurrenceRate(0.1)) |> dl_entropy; digits=2))\$ \$ H_{vert}=$(round(RecurrenceMatrix(dz, GlobalRecurrenceRate(0.1)) |> vl_entropy; digits=2))\$")

plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(600, 1200))
savefig("plots/exercise8_10.png")
