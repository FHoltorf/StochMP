cd(@__DIR__)
include("..\\util\\compare_delta.jl")
include("..\\util\\figformat.jl")

## Process Model
#=
Species:
    in alphabetic order

Reaction Network:
    1:  A + B → C
    2:  C → D
    3:  D → A + B
=#
# Stoichiometry Data
x0 = [20.0,10.0,10.0,0.0]
S = [-1 -1 1 0;
     0 0 -1 1;
     1 1 0 -1]

# Reachable Set (eliminate reaction invariants)
ds = [2,4] # dependent species
x, x_full, is = jump_χ(S,x0,ds)
x0, S = x0[is], S[:,is]
χ = ReachableSet(g=x_full)

# Kinetic Data
c = [1, 2.1, 0.3]
a = [c[1]*x_full[1]*x_full[2],
     c[2]*x_full[3],
     c[3]*x_full[4]]

# Assemble Stochastic Process
jump = JumpProcess(x, a, S, χ)

## Moment Problem Setup
# USER CHOICES!!!
m = 2
n = 2
nR = 3
nT = 10
Tf_range = [0.0,0.05, 0.075, 0.1,0.15,0.2,0.25,0.3,0.35,
            0.4,0.45,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5]

## Bound trajectories
# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=20*ones(2))

# Holtorf, Barton
objs_HB = [(1,x[1],nT),(1,x[2],nT),(2,x[1],nT),(2,x[2],nT)]
bounds_HB, status_HB, time_HB = trajectory_bounds_δ(mp, x0, objs_HB, n, nT, nR, Tf_range)

# Dowdy, Barton
objs_DB = [(1,x[1],1),(1,x[2],1),(2,x[1],1),(2,x[2],1)]
bounds_DB, status_DB, time_DB = trajectory_bounds_δ(mp, x0, objs_DB, 1, 1, nR, Tf_range)

## Save/Load results
save("results/cylic_system.jld", "user_choices", (m,n,nT,nR,Tf_range),
                         "objs_HB", objs_HB,
                         "objs_DB", objs_DB,
                         "bounds_HB", bounds_HB,
                         "bounds_DB", bounds_DB,
                         "time_HB", time_HB,
                         "time_DB", time_DB,
                         "status_HB", status_HB,
                         "status_DB", status_DB)

sampling_results = load("..\\sampling\\results\\cyclic_system_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
red = PyPlot.cm.Reds(0.7)
blue = PyPlot.cm.Blues(0.7)
lns = []
fig, axA = subplots()
axC = axA.twinx()
push!(lns, axA.plot(Tf_range, bounds_DB[objs_DB[1]], marker="o", color=blue, label=L"⟨A⟩^{L/ U}_{DB}"))
push!(lns, axA.plot(Tf_range, bounds_HB[objs_HB[1]], marker="o", color=red, label=L"⟨A⟩^{L/ U}_{HB}"))
push!(lns, axC.plot(Tf_range, bounds_DB[objs_DB[2]], marker="o", color=blue, linestyle="dashed", label=L"⟨C⟩^{L/ U}_{DB}"))
push!(lns, axC.plot(Tf_range, bounds_HB[objs_HB[2]], marker="o", color=red, linestyle="dashed", label=L"⟨C⟩^{L/ U}_{HB}"))
push!(lns, axA.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨A⟩_{SSA}"))
push!(lns, axC.plot(trange, [x[3] for x in x_MC], linestyle="dashed", color = "k", label=L"⟨C⟩_{SSA}"))
axA.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
axA.set_xlabel("time [s]")
axA.set_ylabel("molecular count of A [-]")
axC.set_ylabel("molecular count of C [-]")
fig.savefig("figures\\cyclic_system_mean.pdf")
display(fig)

lns = []
fig, axA = subplots()
axC = axA.twinx()
push!(lns, axA.plot(Tf_range, bounds_DB[objs_DB[3]], marker="o", color=blue, label=L"Var(A)^{U}_{DB}"))
push!(lns, axA.plot(Tf_range, bounds_HB[objs_HB[3]], marker="o", color=red, label=L"Var(A)^{U}_{HB}"))
push!(lns, axC.plot(Tf_range, bounds_DB[objs_DB[4]], marker="o", color=blue, linestyle="dashed", label=L"Var(C)^{U}_{DB}"))
push!(lns, axC.plot(Tf_range, bounds_HB[objs_HB[4]], marker="o", color=red, linestyle="dashed", label=L"Var(C)^{U}_{HB}"))
push!(lns, axA.plot(trange, [x[1] for x in var_MC], color = "k", label=L"Var(A)_{SSA}"))
push!(lns, axC.plot(trange, [x[3] for x in var_MC], linestyle="dashed", color = "k", label=L"Var(C)_{SSA}"))
axA.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
axA.set_xlabel("time [s]")
axA.set_ylabel("molecular count of A [-]")
axC.set_ylabel("molecular count of C [-]")
fig.savefig("figures\\cyclic_system_variance.pdf")
display(fig)
