cd(@__DIR__)
include("..\\util\\sampling.jl")

# initial point
x0 = [53.0, 0, 0, 0, 0, 53.0, 0, 0, 0, 0]
# number of species
N = 10
# stoichiometry
S = Array{Array}(undef,14)
S[1] = [-1, 1, 1, zeros(Int64,7)...]
S[2] = -S[1]
S[3] = [0, -1, 0, 2, zeros(Int64,6)...]
S[4] = -S[3]
S[5] = [0, 0, -1, 0, 1, zeros(Int64,5)...]
S[6] = -S[5]
S[7] = [0, 0, 0, 0, -1, -1, 1, zeros(Int64,3)...]
S[8] = -S[7]
S[9] = [-1, zeros(Int64,6)..., 1, 0, 0]
S[10] = -S[9]
S[11] = [zeros(Int64,7)..., -1, 1, 0]
S[12] = [zeros(Int64,8)..., -1, 1]
S[13] = -S[12]
S[14] = [1, zeros(Int64,8)..., -1]
c = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1e4, 1e4, 1.0, 1.0, 1.0, 1e5, 1e5, 1.0]
Tf_range = [0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]

trange = range(0, stop=Tf_range[end], length=50)

x_MC, var_MC = KMC(x0, c, [[i => -s[i] for i in 1:N if s[i] < 0] for s in S],
                   [[i => s[i] for i in 1:N] for s in S], Tf_range, 5000, trange)

fig, ax = subplots()
ax.plot(trange, [x[1] for x in x_MC], color="k")
display(fig)

fig, ax = subplots()
ax.plot(trange, [x[1] for x in var_MC], color="k")
display(fig)

save("results\\large_system_sampling.jld", "x_MC", x_MC,
                                          "var_MC", var_MC,
                                          "trange", trange)
