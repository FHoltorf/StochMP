cd(@__DIR__)
include("..\\util\\compare_poly.jl")
include("..\\util\\compare_delta.jl")
include("..\\util\\figformat.jl")

## Process Model
#=
Species:
    Substrate => S => 1
    Enzyme => E => 2
    Substrate/Enzyme complex => S:E => 3
    Product => P => 4

Reaction Network:
    1/2:    S + E ⇌ S:E
    3:      S:E → P + E
    4:      P → S
=#
# Stoichiometry Data
x0 = [100.0,100.0,0.0,0.0]
S = [-1 -1 1 0;
     1 1 -1 0;
     0 1 -1 1;
     1 0 0 -1]

# Reachable Set (eliminate reaction invariants)
ds = [2,3]
x, x_full, is = jump_χ(S,x0,ds)
x0, S = x0[is], S[:,is]
χ = ReachableSet(g=x_full)

# Kinetic Data
c = [1.0, 1.0, 1.0, 1]
a = [c[1]*x_full[1]*x_full[2],
     c[2]*x_full[3],
     c[3]*x_full[3],
     c[4]*x_full[4]]

# Assemble Stochastic Process
jump = JumpProcess(x, a, S, χ)

## Moment Problem Setup
# USER CHOICES!!!
m = 4
n = 1
r = 6
nR = 3
nT = 20
Tf_range = [0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3,
            0.4, 0.5, 0.625, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=[50,50])

## Bound trajectories
# Holtorf, Barton
objs_weak = [(1,x[1],nT),(1,x[2],nT)]
bounds_weak, status_weak, time_weak = trajectory_bounds_δ(mp, x0, objs_weak, n, nT, nR, Tf_range)

# Dowdy, Barton
objs_poly = [(1,x[1],nT),(1,x[2],nT)]
bounds_poly, status_poly, time_poly = trajectory_bounds_poly(mp, x0, objs_poly, r, nT, Tf_range)

## Save/Load results
save("results\\michaelis_menten.jld", "user_choices", (m,n,nT,nR,Tf_range),
                                 "objs_weak", objs_weak,
                                 "objs_poly", objs_poly,
                                 "bounds_weak", bounds_weak,
                                 "bounds_poly", bounds_poly,
                                 "time_weak", time_weak,
                                 "time_poly", time_poly,
                                 "status_weak", status_weak,
                                 "status_poly", status_poly)

sampling_results = load("..\\sampling\\results\\michaelis_menten_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[1]], marker="o", color="magenta", label=L"⟨S⟩^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[1]], marker="o", color="red", label=L"⟨S⟩^{L/ U}_{poly}"))
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[2]], marker="o", color="cyan", label=L"⟨P⟩^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[2]], marker="o", color="dodgerblue", label=L"⟨P⟩^{L/ U}_{poly}"))
push!(lns, ax.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨S⟩_{SSA}"))
push!(lns, ax.plot(trange, [x[4] for x in x_MC], linestyle="dashed", color = "k", label=L"⟨P⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns], loc="upper right")
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\michaelis_menten_mean.pdf")
display(fig)

#=
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[3]], marker="o", color="magenta", label=L"V(S)^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[3]], marker="o", color="red", label=L"V(S)^{L/ U}_{poly}"))
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[4]], marker="o", color="cyan", label=L"V(P)^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[4]], marker="o", color="dodgerblue", label=L"V(P)^{L/ U}_{poly}"))
push!(lns, ax.plot(trange, [x[1] for x in var_MC], color = "k", label=L"V(S)_{SSA}"))
push!(lns, ax.plot(trange, [x[4] for x in var_MC], linestyle="dashed", color = "k", label=L"V(P)_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\michaelis_menten_variance_useless.pdf")
display(fig)

lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[3]], marker="o", color="magenta", label=L"V(S)^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[4]], marker="o", color="cyan", label=L"V(P)^{L/ U}_{weak}"))
push!(lns, ax.plot(trange, [x[1] for x in var_MC], color = "k", label=L"V(S)_{SSA}"))
push!(lns, ax.plot(trange, [x[4] for x in var_MC], linestyle="dashed", color = "k", label=L"V(P)_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns], loc="right")
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\michaelis_menten_variance_useful.pdf")
display(fig)
=#
