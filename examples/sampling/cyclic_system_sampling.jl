cd(@__DIR__)
include("..\\util\\sampling.jl")

# number of species
N = 4
# initial state
x0 = [20.0,10.0,10.0,0.0]
# stoichiometry
S = [[-1,-1,1,0],[0,0,-1,1],[1,1,0,-1]]
# kinetic parameters
c = [1, 2.1, 0.3]
# options
Tf_range = [2.5]
trange = range(0, stop=Tf_range[end], length=50)
x_MC, var_MC = KMC(x0, c, [[i=> -s[i] for i in 1:N if s[i] < 0] for s in S],
                          [[i => s[i] for i in 1:N] for s in S], Tf_range, 1e5, trange)

save("results\\cyclic_system_sampling.jld", "x_MC", x_MC,
                                           "var_MC", var_MC,
                                           "trange", trange)
