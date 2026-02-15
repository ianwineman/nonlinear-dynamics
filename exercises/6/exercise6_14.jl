using CSV, DataFrames
using Plots, Measures

"Permutation of relative amplitude"
function permra(x)
	return [findfirst(==(y), sort(x)) for y in x]
end

"Permutation entropy"
function perment(X, d)
	pras = [permra(X[i:(i+d-1)]) for i in 1:(length(X)-d+1)]
	upras = unique(pras)
	Pis = [count(==(up), upras) |> float for up in upras]
	Pis ./= length(Pis)
	return -sum(Pis .* log2.(Pis))
end

"""Permutation entropy \n
for each `W` length window of `X` where each window has an overlap of `w`
with its pre-/suc-ceeding windows. \n\n
Last window may be shorter than `W`.
"""
function perment(X, d, w, W)
	window_ranges = []
	i = 1
	j = i + W - 1
	while true
		push!(window_ranges, i:j)
		i += (W - w)
		j += (W - w)
		if j > length(X)
			push!(window_ranges, i:length(X))
			break
		end
	end
	
	perments = []
	for wr in window_ranges
		Y = X[wr]
		push!(perments, (wr, perment(Y, d)))
	end
	return perments
end


df = CSV.read("book/exercise_data/7.csv", DataFrame; header=false)
X = df[:, 1]

anim = @animate for d in 2:10
	w, W = 100, 1_000
	pes = perment(X, d, w, W)

	plt = plot(xlabel="\$T\$", ylabel="Permutation entropy", legendtitle="\$ d=$d\$", legendposition=:outerright)
	for (i, perment) in enumerate(pes)
		xs, pe = perment
		plot!(plt, xs, collect(range(pe, pe; length=length(xs))), label=ifelse(i==1,"\$ (w, W) = ($w, $W)\$", ""), lc=:blue)
	end

	w, W = 10, 100
	pes = perment(X, d, w, W)

	for (i, perment) in enumerate(pes)
		xs, pe = perment
		plot!(plt, xs, collect(range(pe, pe; length=length(xs))), label=ifelse(i==1,"\$ (w, W) = ($w, $W)\$", ""), lc=:red)
	end

	w, W = 25, 100
	pes = perment(X, d, w, W)

	for (i, perment) in enumerate(pes)
		xs, pe = perment
		plot!(plt, xs, collect(range(pe, pe; length=length(xs))), label=ifelse(i==1,"\$ (w, W) = ($w, $W)\$", ""), lc=:orange)
	end
	plot(
		plt, 
		size=(800, 400), 
		left_margin=5mm, 
		bottom_margin=5mm, 
		xlim=[0, length(X)], 
		ylim=[0, 10]
	)
end
gif(anim, "plots/exercise6_14.gif", fps=0.5)
