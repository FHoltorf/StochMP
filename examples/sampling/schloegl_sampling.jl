cd(@__DIR__)
include("..\\util\\sampling.jl")

# number of species
N = 1
# initial state
x0 = [20.0]
# stoichiometry
S = [[1], [-1], [1], [-1]]
# kinetic parameters
c = [0.15, 1.5e-3, 20.0, 2.0]
# options
Tf_range = [5.0]
trange = range(0, stop=Tf_range[end], length=50)
x_MC, var_MC = KMC(x0, c, [[1=>2], [1=>3], [], [1=>1]],
                          [[i => s[i] for i in 1:N] for s in S], Tf_range, 1e5, trange)

save("results\\schloegl_sampling.jld", "x_MC", x_MC,
                                      "var_MC", var_MC,
                                      "trange", trange)
