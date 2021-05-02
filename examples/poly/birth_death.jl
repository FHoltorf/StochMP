cd(@__DIR__)
include("..\\util\\compare_poly.jl")
include("..\\util\\compare_delta.jl")
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
m = 4
r = 6
n = 1
nR = 3
nT = 20
Tf_range = [0.0, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 5.0, 7.5,
            10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=[7.0])

## Bound trajectories
# Holtorf, Barton
objs_weak = [(1,x,nT)]
bounds_weak, status_weak, time_weak = trajectory_bounds_δ(mp, x0, objs_weak, n, nT, nR, Tf_range)

# Dowdy, Barton
objs_poly = [(1,x,nT)]
bounds_poly, status_poly, time_poly = trajectory_bounds_poly(mp, x0, objs_poly, r, nT, Tf_range)


## Save/Load results
save("results/birth_death.jld", "user_choices", (m,n,nT,nR,Tf_range),
                                 "objs_weak", objs_weak,
                                 "objs_poly", objs_poly,
                                 "bounds_weak", bounds_weak,
                                 "bounds_poly", bounds_poly,
                                 "time_weak", time_weak,
                                 "time_poly", time_poly,
                                 "status_weak", status_weak,
                                 "status_poly", status_poly)

sampling_results = load("..\\sampling\\results\\birth_death_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[1]], marker="o", color="magenta", label=L"⟨A⟩^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[1]], marker="o", color="red", label=L"⟨A⟩^{L/ U}_{poly}"))
push!(lns, ax.plot(trange, [x[1] for x in x_MC], color = "k", label=L"⟨A⟩_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\birth_death_mean.pdf")
display(fig)

#=
lns = []
fig, ax = subplots()
push!(lns, ax.plot(Tf_range, bounds_weak[objs_weak[2]], marker="o", color="magenta", label=L"V(A)^{L/ U}_{weak}"))
push!(lns, ax.plot(Tf_range, bounds_poly[objs_poly[2]], marker="o", color="red", label=L"V(A)^{L/ U}_{poly}"))
push!(lns, ax.plot(trange, [x[1] for x in var_MC], color = "k", label=L"V(A)_{SSA}"))
ax.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns])
ax.set_xlabel("time [s]")
ax.set_ylabel("molecular count [-]")
fig.savefig("figures\\birth_death_variance.pdf")
display(fig)
=#
