using Plots

function henon(u; p=[1.4, 0.3])
	x, y = u
	a, b = p
	return [1 - a*x^2 + y, b*x]
end

tr_len = 15
trs = []
for i in 1:100_000
	u0 = [rand() * 2 - 1, rand() * 2 - 1]
	tr = [u0]

	for j in 1:tr_len
		push!(tr, henon(tr[end]))
	end
	push!(trs, tr)
end


anim = @animate for frame in 0:tr_len
	println(frame + 1)
	bin_dims = 100
	bins = zeros(bin_dims, bin_dims)
	xbins = range(-1.5, 1.5; length=(bin_dims + 1)) # each number in xbins is the left side of bin, except the last num which is the right bound of the last box
	ybins = range(-1.5, 1.5; length=(bin_dims + 1))

	#ends = last.(trs)
	ends = getindex.(trs, frame + 1)

	for e in ends
		x, y = e
		if (-1.5 <= x <= 1.5) && (-1.5 <= x <= 1.5)
			xbin = findfirst(>=(x), xbins) - 1
			ybin = findfirst(>=(y), ybins) - 1
			bins[ybin, xbin] += 1
		end
	end

	bins ./= maximum(bins)
	cgrad_break = bins[findfirst(>(0.0), bins)]

	heatmap(
		xbins,
		ybins,
		bins, 
		c=cgrad([:white, :blue, :orange], [0.0, cgrad_break, cgrad_break + std(bins)]), 
		xlim=[-1.5, 1.5], 
		ylim=[-1.5, 1.5],
		aspect_ratio=:equal,
		xlabel="\$ x\$",
		ylabel="\$ y\$",
		colorbar_title="\$ ρ\$",
		title="\$ n = $frame\$"
	)
end
gif(anim, "plots/exercise8_7.gif", fps=1)
