cd(@__DIR__)
include("..\\util\\compare.jl")
include("..\\util\\figformat.jl")

## Process Model
#=
Species:
    in alphabetic order

Reaction Network:
    1:  ∅ → A
    2:  2A → ∅
=#
# Stoichiometry Data
x0 = [2.0]
S = [1, -2]

# ReachableSet (no reaction invariants)
@polyvar(x)
χ = ReachableSet(g=[x])

# Kinetic Data
c = [1, 0.01]
a = [c[1]*x^0, c[2]*x*(x-1)]

# Assemble Stochastic Process
jump = JumpProcess(x, a, S, χ)

## Moment Problem Setup
# USER CHOICES!!!
m = 8
nI = 2
nR = 3
nT = 20
Tf_range = [0.0, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 5.0, 7.5,
            10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=[7.0])

## Bound trajectories
# Holtorf, Barton
objs_HB = [(1,x,nT),(2,x,nT)]
bounds_HB, status_HB, time_HB = trajectory_bounds(mp, x0, objs_HB, nI, nT, nR, Tf_range)

# Dowdy, Barton
objs_DB = [(1,x,1),(2,x,1)]
bounds_DB, status_DB, time_DB = trajectory_bounds(mp, x0, objs_DB, 1, 1, nR, Tf_range)

## Save/Load results
save("results/birth_death.jld", "user_choices", (m,nI,nT,nR,Tf_range),
                                 "objs_HB", objs_HB,
                                 "objs_DB", objs_DB,
                                 "bounds_HB", bounds_HB,
                                 "bounds_DB", bounds_DB,
                                 "time_HB", time_HB,
                                 "time_DB", time_DB,
                                 "status_HB", status_HB,
                                 "status_DB", status_DB)

sampling_results = load("..\\sampling\\results\\birth_death_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
red = PyPlot.cm.Reds(0.7)
blue = PyPlot.cm.Blues(0.7)
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[1]], marker="o", color=blue, label=L"⟨A⟩^{L/ U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[1]], marker="o", color=red, label=L"⟨A⟩^{L/ U}_{HB}"))
push!(lns, ax.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨A⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\birth_death_mean.pdf")
display(fig)

lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[2]], marker="o", color=blue, label=L"Var(A)^{U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[2]], marker="o", color=red, label=L"Var(A)^{U}_{HB}"))
push!(lns, ax.plot(trange, [x[1] for x in var_MC], color = "k", label=L"Var(A)_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\birth_death_variance.pdf")
display(fig)
