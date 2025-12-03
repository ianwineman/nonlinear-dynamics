using Plots

function logistic_map(x)
	return 4*x*(1-x)
end

function tent_map(x)
	if x ≤ 0.5
		return 2*x
	else
		return 2 - 2*x
	end
end

function trajectory(f, x0, N)
	t = zeros(typeof(x0), N)
	t[1] = x0
	for n in 2:N
		t[n] = f(t[n-1])
	end
	return t
end



let 
	setprecision(BigFloat, 10_000)
	x0 = BigFloat(π)/4
	N = 10_000
	t = trajectory(tent_map, x0, N)
	p1 = plot(t, 1:N, 
		xlabel="\$x_n\$", ylabel="\$n\$", 
		lα=0.1, lw=0.5, label=false,
		xlim=[0,1], ylim=[0,N],
		grid=false,
		title="Tent map: \$μ=2\$"
	)

	ϵ = 0.01
	bins = zeros(Int(1.0/ϵ))
	for point in t
		i = Int(point ÷ ϵ + 1)
		bins[i] += 1
	end
	bins ./= N

	bin_centers = (ϵ/2):ϵ:(1.0-(ϵ/2))
	p2 = bar(
		bin_centers, bins,
		xlim=[0,1],
		label=false,
		lw=0,
		ylim=[minimum(bins), maximum(bins)],
		title="Histogram: \$ϵ=$ϵ\$",
		xlabel="\$x\$"
	)

	plot(p1,p2)
	savefig("plots/exercise5_3_tent_hist.png")

	H = -sum(bins .* log.(bins))
	println("Tent map     Shannon entropy: $H") # ~4.601
end

let
	x0 = 0.4
	N = 10_000
	t = trajectory(logistic_map, x0, N)
	p1 = plot(t, 1:N, 
		xlabel="\$x_n\$", ylabel="\$n\$", 
		lα=0.1, lw=0.5, label=false,
		xlim=[0,1], ylim=[0,N],
		grid=false,
		title="Logistic map: \$r=4\$"
	)

	ϵ = 0.01
	bins = zeros(Int(1.0/ϵ))
	for point in t
		i = Int(point ÷ ϵ + 1)
		bins[i] += 1
	end
	bins ./= N

	bin_centers = (ϵ/2):ϵ:(1.0-(ϵ/2))
	p2 = bar(
		bin_centers, bins,
		xlim=[0,1],
		label=false,
		lw=0,
		ylim=[minimum(bins), maximum(bins)],
		title="Histogram: \$ϵ=$ϵ\$",
		xlabel="\$x\$"
	)

	plot(p1,p2)
	savefig("plots/exercise5_3_logistic_hist.png")

	H = -sum(bins .* log.(bins))
	println("Logistic map Shannon entropy: $H") # ~4.391
end
