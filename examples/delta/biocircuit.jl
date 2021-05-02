cd(@__DIR__)
include("..\\util\\compare_delta.jl")
include("..\\util\\figformat.jl")

## Process Model
#=
Species:
    mRNA => 1
    Protein => P => 2
    DNA => 3
    Protein/DNA complex => P:DNA => 4

Reaction Network:
    1:  DNA → DNA + mRNA
    2:  mRNA → ∅ (degradation)
    3:  mRNA → mRNA + P
    4:  P → ∅ (degradation)
    5:  P + DNA → P:DNA (binding of repressor to promoter) (negative feedback)
    6:  P:DNA → P + DNA (unbinding of repressor and promoter) (negative feedback)
=#
# Stoichiometry Data
D_T1 = 20
x0 = [10,0,D_T1,0]
S = [1 0 0 0;
     -1 0 0 0;
     0 1 0 0;
     0 -1 0 0;
     0 -1 -1 1;
     0 1 1 -1]

# Reachable Set (eliminate reaction invariants)
ds = [4] # dependent species
x, x_full, is = jump_χ(S,x0,ds)
x0, S = x0[is], S[:,is]
χ = ReachableSet(g=x_full)

# Kinetic Data
c = [0.2, log(2)/5, 0.5, log(2)/20, 5, 1.0]
a = [c[1]*x_full[3],
     c[2]*x_full[1],
     c[3]*x_full[1],
     c[4]*x_full[2],
     c[5]*x_full[2]*x_full[3],
     c[6]*x_full[4]]

# Assemble Stochastic Process
jump = JumpProcess(x, a, S, χ)

## Moment Problem Setup
# USER CHOICES!!!
m = 4
n = 2
nR = 3
nT = 10
Tf_range = [0.0, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.5,
            10.0, 12.5, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=D_T1*ones(3))

## Bound trajectories
# Holtorf, Barton
objs_HB = [(1,x[1],nT),(1,x[2],nT),(2,x[1],nT),(2,x[2],nT)]
bounds_HB, status_HB, time_HB = trajectory_bounds_δ(mp, x0, objs_HB, n, nT, nR, Tf_range)

# Dowdy, Barton
objs_DB = [(1,x[1],1),(1,x[2],1),(2,x[1],1),(2,x[2],1)]
bounds_DB, status_DB, time_DB = trajectory_bounds_δ(mp, x0, objs_DB, 1, 1, nR, Tf_range)
## Save/Load results
save("results\\bio_circuit.jld", "user_choices", (m,n,nT,nR,Tf_range),
                                "objs_HB", objs_HB,
                                "objs_DB", objs_DB,
                                "bounds_HB", bounds_HB,
                                "bounds_DB", bounds_DB,
                                "time_HB", time_HB,
                                "time_DB", time_DB,
                                "status_HB", status_HB,
                                "status_DB", status_DB)

sampling_results = load("..\\sampling\\results\\biocircuit_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

## Plot results
red = PyPlot.cm.Reds(0.7)
blue = PyPlot.cm.Blues(0.7)
lns = []
fig, axP = subplots()
axmRNA = axP.twinx()
push!(lns, axP.plot(Tf_range, bounds_DB[objs_DB[2]], color = blue, marker="o", label=L"⟨P⟩^{L/U}_{DB}"))
push!(lns, axP.plot(Tf_range, bounds_HB[objs_HB[2]], color = red, marker="o", label=L"⟨P⟩^{L/U}_{HB}"))
push!(lns, axmRNA.plot(Tf_range, bounds_DB[objs_DB[1]], color = blue, linestyle="dashed",  marker="o", label=L"⟨mRNA⟩^{L/U}_{DB}"))
push!(lns, axmRNA.plot(Tf_range, bounds_HB[objs_HB[1]], color = red, linestyle="dashed", marker="o", label=L"⟨mRNA⟩^{L/U}_{HB}"))
push!(lns, axP.plot(trange, [x[2] for x in x_MC], color = "k", label = L"⟨P⟩_{SSA}"))
push!(lns, axmRNA.plot(trange, [x[1] for x in x_MC], linestyle = "dashed", color = "k", label = L"⟨mRNA⟩_{SSA}"))
axP.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns], bbox_to_anchor=(0.5, -0.35), ncol = 3, loc="center")

axP.set_xlabel("time [min]")
axP.set_ylabel("molecular count of P [-]")
axmRNA.set_ylabel("molecular count of mRNA [-]")
fig.tight_layout()
fig.savefig("figures\\biocircuit_mean.pdf")
display(fig)

# generate plot with two axes
lns = []
fig, axP = subplots()
axmRNA = axP.twinx()
push!(lns, axP.plot(Tf_range, bounds_DB[objs_DB[4]], color = blue, marker="o", label=L"Var(P)^{U}_{DB}"))
push!(lns, axP.plot(Tf_range, bounds_HB[objs_HB[4]], color = red, marker="o", label=L"Var(P)^{U}_{HB}"))
push!(lns, axmRNA.plot(Tf_range, bounds_DB[objs_DB[3]], color = blue, linestyle="dashed", marker="o", label=L"Var(mRNA)^{U}_{DB}"))
push!(lns, axmRNA.plot(Tf_range, bounds_HB[objs_HB[3]], color = red, linestyle="dashed", marker="o", label=L"Var(mRNA)^{U}_{HB}"))
push!(lns, axP.plot(trange, [x[2] for x in var_MC], color = "k", label = L"Var(P)_{SSA}"))
push!(lns, axmRNA.plot(trange, [x[1] for x in var_MC], linestyle = "dashed", color = "k", label = L"Var(mRNA)_{SSA}"))
axP.legend([ln[1] for ln in lns], [ln[1].get_label() for ln in lns], bbox_to_anchor=(0.5, -0.35), ncol = 3, loc="center")
axP.set_xlabel("time [min]")
axP.set_ylabel("molecular count of P [-]")
axmRNA.set_ylabel("molecular count of mRNA [-]")
fig.tight_layout()
fig.savefig("figures\\biocircuit_variance.pdf")
display(fig)
