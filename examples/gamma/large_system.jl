cd(@__DIR__)
include("..\\util\\compare_gamma.jl")
include("..\\util\\figformat.jl")

## Process Model
#=
Species:
    in alphabetic order

Reaction Network:
    1/2:        A ⇋ B + C
    3/4:        B ⇋ 2D
    5/6:        C ⇋ E
    7/8:        E + F ⇋ G
    9/10:       A ⇋ H
    11:         H → I
    12/13:      I ⇋ J
    14:         J → A
=#
# Stoichiometry Data
x0 = [53.0, 0, 0, 0, 0, 53.0, 0, 0, 0, 0]
S = [  -1   1   1   0   0   0   0   0   0   0;
        1  -1  -1   0   0   0   0   0   0   0;
        0  -1   0   2   0   0   0   0   0   0;
        0   1   0  -2   0   0   0   0   0   0;
        0   0  -1   0   1   0   0   0   0   0;
        0   0   1   0  -1   0   0   0   0   0;
        0   0   0   0  -1  -1   1   0   0   0;
        0   0   0   0   1   1  -1   0   0   0;
       -1   0   0   0   0   0   0   1   0   0;
        1   0   0   0   0   0   0  -1   0   0;
        0   0   0   0   0   0   0  -1   1   0;
        0   0   0   0   0   0   0   0  -1   1;
        0   0   0   0   0   0   0   0   1  -1;
        1   0   0   0   0   0   0   0   0  -1]

# Reachable Set (eliminate reaction invariants)
ds = [5,7,9] # dependent species
x, x_full, is = jump_χ(S,x0,ds)
x0, S = x0[is], S[:,is]
χ = ReachableSet(g=x_full)

# Kinetic Data
c = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1e4, 1e4, 1.0, 1.0, 1.0, 1e5, 1e5, 1.0]
a = [c[1]*x_full[1],
     c[2]*x_full[2]*x_full[3],
     c[3]*x_full[2],
     c[4]*x_full[4]*(x_full[4]-1.0),
     c[5]*x_full[3],
     c[6]*x_full[5],
     c[7]*x_full[5]*x_full[6],
     c[8]*x_full[7],
     c[9]*x_full[1],
     c[10]*x_full[8],
     c[11]*x_full[8],
     c[12]*x_full[9],
     c[13]*x_full[10],
     c[14]*x_full[10]]

# Assemble Stochastic Process
jump = JumpProcess(x, a, S, χ)

## Moment Problem Setup
# USER CHOICES!!!
m = 2
nI = 2
nR = 5
nT = 5
Tf_range = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75,
            1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=53*ones(7))

## Bound trajectories
# Holtorf, Barton
objs_HB = [(1,x[1],nT),(1,x[6],nT)]
bounds_HB, status_HB, time_HB = trajectory_bounds_γ(mp, x0, objs_HB, nI, nT, nR, Tf_range)

# Dowdy, Barton
objs_DB = [(1,x[1],1),(1,x[6],1)]
bounds_DB, status_DB, time_DB = trajectory_bounds_γ(mp, x0, objs_DB, 0, 1, nR, Tf_range)

## Save/Load results
save("results/large_system.jld", "user_choices", (m,nI,nT,nR,Tf_range),
                                 "objs_HB", objs_HB,
                                 "objs_DB", objs_DB,
                                 "bounds_HB", bounds_HB,
                                 "bounds_DB", bounds_DB,
                                 "time_HB", time_HB,
                                 "time_DB", time_DB,
                                 "status_HB", status_HB,
                                 "status_DB", status_DB)

sampling_results = load("..\\sampling\\results\\large_system_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
red = PyPlot.cm.Reds(0.7)
blue = PyPlot.cm.Blues(0.7)
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[1]], marker="o", color=blue, label=L"⟨A⟩^{L/ U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[1]], marker="o", color=red, label=L"⟨A⟩^{L/ U}_{HB}"))
push!(lns, ax.plot(Tf_range, bounds_DB[objs_DB[2]], marker="o", color=blue, linestyle="dashed", label=L"⟨H⟩^{L/ U}_{DB}"))
push!(lns, ax.plot(Tf_range, bounds_HB[objs_HB[2]], marker="o", color=red, linestyle="dashed", label=L"⟨H⟩^{L/ U}_{HB}"))
push!(lns, ax.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨A⟩_{SSA}"))
push!(lns, ax.plot(trange, [x[8] for x in x_MC], linestyle="dashed", color = "k", label=L"⟨H⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\large_system_mean.pdf")
display(fig)
