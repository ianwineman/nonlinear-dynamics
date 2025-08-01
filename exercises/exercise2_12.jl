using Plots, FFTW, Statistics, CSV, DataFrames

df = CSV.read("exercises/exercise2_12.csv", DataFrame; header=false);

plots = []
Δt = 0.05
for col in names(df)
	r = df[:, col]
	P = abs2.(rfft(r .- mean(r)))
	ν = rfftfreq(length(r), 1/Δt)

	p = plot(ν, P ./ maximum(P), label="$col")
	push!(plots, p)
end

plot(
	plots..., 
	layout=(2,5), 
	ylim=[10.0^(-5), 1.0], 
	yscale=:log10, 
	size=(1200,600),
	linewidth=1,
	linecolor=:cornflowerblue,
	plot_title="Normalized power spectra"
)
savefig("plots/figure2_12.png")
