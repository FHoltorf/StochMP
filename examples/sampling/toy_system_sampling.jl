cd(@__DIR__)
include("..\\util\\sampling.jl")

# number of species
N = 4
# initial state
x0 = [40.0,41.0,0.0,0.0]
# stoichiometry
S = [[-1,-1,1,0],[0,0,-1,1],[0,0,1,-1]]
# kinetic parameters
c = [1, 2.1, 0.3]
# Propensities
Tf_range = [5.0]
trange = range(0, stop=Tf_range[end], length=1000)
x_MC, var_MC = KMC(x0, c, [[i=> -s[i] for i in 1:N if s[i] < 0] for s in S],
                          [[i => s[i] for i in 1:N] for s in S], Tf_range, 1e5, trange)

save("results\\toy_system_sampling.jld", "x_MC", x_MC,
                                        "var_MC", var_MC,
                                        "trange", trange)
