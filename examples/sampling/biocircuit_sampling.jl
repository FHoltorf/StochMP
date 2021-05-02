cd(@__DIR__)
include("..\\util\\sampling.jl")

#=
reaction system
1:  DNA → DNA + mRNA
2:  mRNA → ∅ (degradation)
3:  mRNA → mRNA + P
4:  P → ∅ (degradation)
5:  P + DNA → P:DNA (binding of repressor to promoter) (negative feedback)
6:  P:DNA → P + DNA (unbinding of repressor and promoter) (negative feedback)
=#

# number of species
N = 4
# initial state
D_T1 = 20
x0 = [10,0,D_T1,0]
# stoichiometry
S = [[1,0,0,0], [-1,0,0,0], [0,1,0,0], [0,-1,0,0], [0,-1,-1,1], [0,1,1,-1]]
# kinetic parameters
c = [0.2, log(2)/5, 0.5, log(2)/20, 5.0, 1.0] # [1/min] # paper
# options
Tf_range = [50.0]
trange = range(0, stop=Tf_range[end], length=500)
x_MC, var_MC = KMC(x0, c, [[3=>1], [1=>1], [1=>1], [2=>1], [2=>1, 3=>1], [4=>1]],
                          [[i => s[i] for i in 1:N] for s in S], Tf_range, 1e5, trange)


save("results\\biocircuit_sampling.jld", "x_MC", x_MC,
                                         "var_MC", var_MC,
                                         "trange", trange)
