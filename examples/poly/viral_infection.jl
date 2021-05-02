cd(@__DIR__)
include("..\\util\\compare_poly.jl")
include("..\\util\\compare_delta.jl")
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
r = 6
n = 1
nR = 3
nT = 10 # 20
Tf_range = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40, 50.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=[20,200.0,20000.0])

## Bound trajectories
# Holtorf, Barton
objs_weak = [(1,x[1],nT),(1,x[2],nT),(1,x[3],nT)]
bounds_weak, status_weak, time_weak = trajectory_bounds_δ(mp, x0, objs_weak, n, nT, nR, Tf_range)

# Dowdy, Barton
objs_poly = [(1,x[1],nT),(1,x[2],nT),(1,x[3],nT)]
bounds_poly, status_poly, time_poly = trajectory_bounds_poly(mp, x0, objs_poly, r, nT, Tf_range)

## Save/Load results
save("results\\viral_infection.jld", "user_choices", (m,n,nT,nR,Tf_range),
                                 "objs_weak", objs_weak,
                                 "objs_poly", objs_poly,
                                 "bounds_weak", bounds_weak,
                                 "bounds_poly", bounds_poly,
                                 "time_weak", time_weak,
                                 "time_poly", time_poly,
                                 "status_weak", status_weak,
                                 "status_poly", status_poly)

sampling_results = load("..\\sampling\\results\\viral_infection_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[1]], marker="o", color="magenta", label=L"⟨T⟩^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[1]], marker="o", color="red", label=L"⟨T⟩^{L/ U}_{poly}"))
push!(lns, ax.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨T⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\toy_system_mean_T.pdf")
display(fig)

lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[2]], marker="o", color="magenta", label=L"⟨G⟩^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[2]], marker="o", color="red", label=L"⟨G⟩^{L/ U}_{poly}"))
push!(lns, ax.plot(trange, [x[2] for x in x_MC], color = "k", label=L"⟨G⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\toy_system_mean_G.pdf")
display(fig)

lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[3]], marker="o", color="magenta", label=L"⟨S⟩^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[3]], marker="o", color="red", label=L"⟨S⟩^{L/ U}_{poly}"))
push!(lns, ax.plot(trange, [x[3] for x in x_MC], color = "k", label=L"⟨S⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\toy_system_mean_S.pdf")
display(fig)
