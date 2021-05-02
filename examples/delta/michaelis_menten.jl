cd(@__DIR__)
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
n = 2
nR = 3
nT = 20
Tf_range = [0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3,
            0.4, 0.5, 0.625, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=[50,50])

## Bound trajectories
# Holtorf, Barton
objs_HB = [(1,x[1],nT),(1,x[2],nT),(2,x[1],nT),(2,x[2],nT)]
bounds_HB, status_HB, time_HB = trajectory_bounds_δ(mp, x0, objs_HB, n, nT, nR, Tf_range)

# Dowdy, Barton
objs_DB = [(1,x[1],1),(1,x[2],1),(2,x[1],1),(2,x[2],1)]
bounds_DB, status_DB, time_DB = trajectory_bounds_δ(mp, x0, objs_DB, 1, 1, nR, Tf_range)

## Save/Load results
save("results/michaelis_menten.jld", "user_choices", (m,n,nT,nR,Tf_range),
                                 "objs_HB", objs_HB,
                                 "objs_DB", objs_DB,
                                 "bounds_HB", bounds_HB,
                                 "bounds_DB", bounds_DB,
                                 "time_HB", time_HB,
                                 "time_DB", time_DB,
                                 "status_HB", status_HB,
                                 "status_DB", status_DB)

sampling_results = load("..\\sampling\\results\\michaelis_menten_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
red = PyPlot.cm.Reds(0.7)
blue = PyPlot.cm.Blues(0.7)
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[1]], marker="o", color=blue, label=L"⟨S⟩^{L/ U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[1]], marker="o", color=red,  label=L"⟨S⟩^{L/ U}_{HB}"))
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[2]], marker="o", color=blue, linestyle="dashed", label=L"⟨P⟩^{L/ U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[2]], marker="o", color=red, linestyle="dashed", label=L"⟨P⟩^{L/ U}_{HB}"))
push!(lns, ax.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨S⟩_{SSA}"))
push!(lns, ax.plot(trange, [x[4] for x in x_MC], linestyle="dashed", color = "k", label=L"⟨P⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns], loc="upper right")
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\michaelis_menten_mean.pdf")
display(fig)

lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[3]], marker="o", color=blue, label=L"Var(S)^{U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[3]], marker="o", color=red, label=L"Var(S)^{U}_{HB}"))
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[4]], marker="o", color=blue, linestyle="dashed", label=L"Var(P)^{U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[4]], marker="o", color=red, linestyle="dashed", label=L"Var(P)^{U}_{HB}"))
push!(lns, ax.plot(trange, [x[1] for x in var_MC], color = "k", label=L"Var(S)_{SSA}"))
push!(lns, ax.plot(trange, [x[4] for x in var_MC], linestyle="dashed", color = "k", label=L"Var(P)_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\michaelis_menten_variance_useless.pdf")
display(fig)

lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[3]], marker="o", color=red, label=L"Var(S)^{U}_{HB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[4]], marker="o", color=red, linestyle="dashed", label=L"Var(P)^{U}_{HB}"))
push!(lns, ax.plot(trange, [x[1] for x in var_MC], color = "k", label=L"Var(S)_{SSA}"))
push!(lns, ax.plot(trange, [x[4] for x in var_MC], linestyle="dashed", color = "k", label=L"Var(P)_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns], loc="right")
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\michaelis_menten_variance_useful.pdf")
display(fig)
