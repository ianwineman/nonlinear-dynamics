using LinearAlgebra, Plots
re_λ1s, re_λ2s, im_λ1s, im_λ2s = [], [], [], []
μ_range = -10.0:0.01:10.0
for μ in μ_range
   J = [0 1; -1 -μ]
   λ1, λ2 = eigen(J).values
   push!(re_λ1s, real(λ1))
   push!(re_λ2s, real(λ2))
   push!(im_λ1s, imag(λ1))
   push!(im_λ2s, imag(λ2))
end
plot(
	μ_range, [re_λ1s, re_λ2s, im_λ1s, im_λ2s], 
	xtick=-10:1:10, 
	xlabel="\$μ\$",
	label=["\$Re(λ_1)\$" "\$Re(λ_2)\$" "\$Im(λ_1)\$" "\$Im(λ_2)\$"], 
	lc=[:red :orange :blue :purple]
)
savefig("plots/exercise4_4.png")
