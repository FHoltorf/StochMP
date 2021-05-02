cd(@__DIR__)
include("..\\util\\sampling.jl")

# number of species
N = 3
# initial state
x0 = [1.0,0,0.0]
# stoichiometry
S = [[1, -1, 0], [-1, 0, 0], [0, 1, 0], [0, -1, -1], [0, 0, 1], [0, 0, -1]]
# kinetic parameters
c = [0.025, 0.25, 1.0, 7.5e-6, 1000.0, 1.99]
Tf_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40, 50.0]
#=
reaction system
1:  genome → template
2:  template → ∅ (degradation)
3:  template → genome (but template is not lost)
4:  genome + struct → ∅ (virus formation)
5:  template → struct (but template is not lost)
6:  struct → ∅ (degradation)
=#
trange = range(0, stop=Tf_range[end], length=1000)
x_MC, var_MC = KMC(x0, c, [[2=>1], [1=>1], [1=>1], [2=>1,3=>1], [1=>1], [3=>1]],
                   [[i => s[i] for i in 1:N] for s in S], Tf_range, 20000, trange)

save("results\\viral_infection_sampling.jld", "x_MC", x_MC,
                                             "var_MC", var_MC,
                                             "trange", trange)
