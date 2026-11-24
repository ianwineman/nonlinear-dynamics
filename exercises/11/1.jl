using Plots

q = range(-5, 5; length=100)
y1 = 0 .- q.^2  # u^* = 0, r = r_c = 0
y2 = -2 .- q.^2 # u^* = ±√r, r = 1 (no r_c)

plot(
	q, 
	y1, 
	xlabel="\$ q\$",
	ylabel="\$ σ(q)\$",
	label="\$ u^* = 0, r = r_c = 0\$"
)
plot!(q, y2, label="\$ u^* = ±√r, r = 1 \$ (no \$ r_c\$)")
savefig("plots/11.1.png")

# u^* = 0 stable when r < r_c = 0, unstable otherwise
# u^* = ±√r always stable
