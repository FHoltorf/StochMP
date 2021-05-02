cd(@__DIR__)
include("..\\util\\sampling.jl")

# number of species
N = 1
# initial condition
x0 = [2.0]
# stoichiometry
S = [[1], [-2]]
# kinetic data
c = [1, 0.01]
# final time
Tf_range = [50.0]
# time points
trange = range(0, stop=Tf_range[end], length=50)
x_MC, var_MC = KMC(x0, c, [[i=> -s[i] for i in 1:N if s[i] < 0] for s in S],
                          [[i => s[i] for i in 1:N] for s in S], Tf_range, 1e5, trange)

save("results\\birth_death_sampling.jld", "x_MC", x_MC,
                                         "var_MC", var_MC,
                                         "trange", trange)
