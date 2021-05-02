cd(@__DIR__)
include("..\\util\\compare_poly.jl")
include("..\\util\\compare_delta.jl")
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
r = 4
n = 1
nR = 5
nT = 5
Tf_range = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75,
            1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=53*ones(7))

## Bound trajectories
# Holtorf, Barton
objs_weak = [(1,x[1],nT),(1,x[6],nT)]
bounds_weak, status_weak, time_weak = trajectory_bounds_δ(mp, x0, objs_weak, n, nT, nR, Tf_range)

# Dowdy, Barton
objs_poly = [(1,x[1],nT),(1,x[6],nT)]
bounds_poly, status_poly, time_poly = trajectory_bounds_poly(mp, x0, objs_poly, r, nT, Tf_range)

## Save/Load results
save("results/large_system.jld", "user_choices", (m,n,nT,nR,Tf_range),
                                 "objs_weak", objs_weak,
                                 "objs_poly", objs_poly,
                                 "bounds_weak", bounds_weak,
                                 "bounds_poly", bounds_poly,
                                 "time_weak", time_weak,
                                 "time_poly", time_poly,
                                 "status_weak", status_weak,
                                 "status_poly", status_poly)

sampling_results = load("..\\sampling\\results\\large_system_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[1]], marker="o", color="magenta", label=L"⟨A⟩^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[1]], marker="o", color="red", label=L"⟨A⟩^{L/ U}_{poly}"))
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[2]], marker="o", color="cyan", label=L"⟨H⟩^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[2]], marker="o", color="dodgerblue", label=L"⟨H⟩^{L/ U}_{poly}"))
push!(lns, ax.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨A⟩_{SSA}"))
push!(lns, ax.plot(trange, [x[8] for x in x_MC], linestyle="dashed", color = "k", label=L"⟨H⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\large_system_mean.pdf")
display(fig)
