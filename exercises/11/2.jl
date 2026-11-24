using Plots

q = range(-2, 2; length=100)
σ = -1 .+ (2 .* (q .^ 2)) .- (q .^ 4)  # u^* = 0

plot(
	q, 
	σ, 
	xlabel="\$ q\$",
	ylabel="\$ σ(q)\$",
	label="\$ u^* = 0, r = r_c = 0\$"
)
scatter!([-1, 1], [0, 0], label="\$ q_c\$")
savefig("plots/11.2.png")
