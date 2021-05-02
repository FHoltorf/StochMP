cd(@__DIR__)
include("..\\util\\compare.jl")
include("..\\util\\figformat.jl")

## Process Model
#=
Species:
    Template => T  => 1,
    Genome => G => 2,
    Struct => S => 3

reaction system
    1:  G → T
    2:  T → ∅ (degradation)
    3:  T → G (but template is not lost)
    4:  G + S → ∅ (virus formation)
    5:  T → S (but template is not lost)
    6:  S → ∅ (degradation/secretion)
=#

# Stoichiometry Data
x0 = [1.0,0,0.0]
S = [1 -1 0;
     -1 0 0;
     0 1 0;
     0 -1 -1;
     0 0 1;
     0 0 -1]


# Reachable Set (eliminate reaction invariants)
@polyvar(x[1:3])
χ = ReachableSet(g=x)

# Kinetic Data
c = [0.025, 0.25, 1.0, 7.5e-6, 1000.0, 1.99]
a = [c[1]*x[2],
     c[2]*x[1],
     c[3]*x[1],
     c[4]*x[2]*x[3],
     c[5]*x[1],
     c[6]*x[3]]

# Assemble Stochastic Process
jump = JumpProcess(x, a, S, χ)

## Moment Problem Setup
# USER CHOICES!!!
m = 4
nI = 2
nR = 3
nT = 10 # 20
Tf_range = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40, 50.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=[20,200.0,20000.0])

## Bound trajectories
# Holtorf, Barton
objs_HB = [(1,x[1],nT),(1,x[2],nT),(1,x[3],nT)]
bounds_HB, status_HB, time_HB = trajectory_bounds(mp, x0, objs_HB, nI, nT, nR, Tf_range)

# Dowdy, Barton
objs_DB = [(1,x[1],1),(1,x[2],1),(1,x[3],1)]
bounds_DB, status_DB, time_DB = trajectory_bounds(mp, x0, objs_DB, 1, 1, nR, Tf_range)

## Save/Load results
save("results\\viral_infection.jld", "user_choices", (m,nI,nT,nR,Tf_range),
                                    "objs_HB", objs_HB,
                                    "objs_DB", objs_DB,
                                    "bounds_HB", bounds_HB,
                                    "bounds_DB", bounds_DB,
                                    "time_HB", time_HB,
                                    "time_DB", time_DB,
                                    "status_HB", status_HB,
                                    "status_DB", status_DB)

sampling_results = load("..\\sampling\\results\\viral_infection_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
red = PyPlot.cm.Reds(0.7)
blue = PyPlot.cm.Blues(0.7)
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[1]], marker="o", color=blue, label=L"⟨T⟩^{L/ U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[1]], marker="o", color=red, label=L"⟨T⟩^{L/ U}_{HB}"))
push!(lns, ax.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨T⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\viral_infection_mean_template.pdf")
display(fig)

lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[2]], marker="o", color=red, label=L"⟨G⟩^{L/ U}_{HB}"))
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[2]], marker="o", color=blue, label=L"⟨G⟩^{L/ U}_{DB}"))
push!(lns, ax.plot(trange, [x[2] for x in x_MC], color = "k", label=L"⟨G⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\viral_infection_mean_genome.pdf")
display(fig)

lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[3]], marker="o", color=red, label=L"⟨S⟩^{L/ U}_{HB}"))
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[3]], marker="o", color=blue, label=L"⟨S⟩^{L/ U}_{DB}"))
push!(lns, ax.plot(trange, [x[3] for x in x_MC], color = "k", label=L"⟨S⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\viral_infection_mean_struct.pdf")
display(fig)
