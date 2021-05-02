using StochMP, PyPlot, LaTeXStrings, LinearAlgebra, DynamicPolynomials, JLD, MosekTools
cd(@__DIR__)
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

n = 2
nR = 1
Tf = 20.0
ϕ_range = range(0, stop = 2π, length=250)

x_scale = 7.0

## "True" solution
sampling_results = load("..\\sampling\\results\\birth_death_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]

true_y1 = x_MC[20]+ (x_MC[21]-x_MC[20])/(trange[21]-trange[20])*(Tf-trange[20])
true_var = var_MC[20]+ (var_MC[21]-var_MC[20])/(trange[21]-trange[20])*(Tf-trange[20])
true_y2 = true_var .+ true_y1.^2
## Varying Truncation Order/Fixed Time Discretization
m_range = [2,4,6,8]
nT = 5
contours = Dict()
sol_stats = Dict()
for m in m_range
    mp = MomentProblem(jump,m;x_scale=[x_scale])
    R = -unique(sort(round.(svd(mp.A).S, digits=3)))[1:min(nR,end)]
    T = range(Tf/nT, stop=Tf, length=nT)
    contours[m] = []
    sol_stats[m] = []
    obj = [(3,Dict([x=>0, x^2=>1]),nT)]
    d = []
    for Tf_test in range(0.5, stop = Tf, length=20)
        T = range(Tf_test/nT, stop=Tf_test, length=nT)
        sol, d, model = transient_bounds(mp, x0, R, T, n, obj, Mosek.Optimizer; d=d)
    end
    for ϕ in ϕ_range
        obj = [(3,Dict([x=>sin(ϕ)/x_scale, x^2=>cos(ϕ)/x_scale^2]),nT)]
        sol, d, model = transient_bounds(mp, x0, R, T, n, obj, Mosek.Optimizer; d=d)
        push!(contours[m], [sol[obj[1]][2][2], sol[obj[1]][2][3]])
        push!(sol_stats[m], sol[obj[1]][3])
    end
end
fig, ax = subplots()
for m in m_range
    ax.plot([c[1] for c in contours[m]], [c[2] for c in contours[m]],
             label=latexstring("n_{T} = ",nT,", m = ", m))
end
ax.plot(true_y1, true_y2, marker="*", markersize=10, color="k", linestyle="", label="true moments")
ax.legend()
ax.set_xlabel(L"y_{1}(t_f)")
ax.set_ylabel(L"y_{2}(t_f)")
fig.savefig("figures\\set_projection_truncation_order.pdf")
display(fig)

## Fixed Truncation Order/Varying Time Discretization
m = 2
nT_range = [5,10,15,20]

mp = MomentProblem(jump,m; x_scale=[x_scale])
R = -unique(sort(round.(svd(mp.A).S, digits=3)))[1:min(nR,end)]
contours = Dict()
sol_stats = Dict()
for nT in nT_range
    T = range(Tf/nT, stop=Tf, length=nT)
    contours[nT] = []
    sol_stats[nT] = []
    obj = [(3,Dict([x=>1, x^2=>0]),nT)]
    sol, d, model = transient_bounds(mp, x0, R, T, n, obj, Mosek.Optimizer)
    for ϕ in ϕ_range
        obj = [(3,Dict([x=>sin(ϕ)/x_scale, x^2=>cos(ϕ)/x_scale^2]),nT)]
        sol, d, model = transient_bounds(mp, x0, R, T, n, obj, Mosek.Optimizer; d=d)
        push!(contours[nT], [sol[obj[1]][2][2], sol[obj[1]][2][3]])
        push!(sol_stats[nT], sol[obj[1]][3])
    end
end
fig, ax = subplots()
for nT in nT_range
    ax.plot([c[1] for c in contours[nT]], [c[2] for c in contours[nT]],
             label=latexstring("n_{T} = ", nT,", m = ", m))
end
ax.plot(true_y1, true_y2, marker="*", markersize=10, color="k", linestyle="", label="true moments")
ax.legend()
ax.set_xlabel(L"y_{1}(t_f)")
ax.set_ylabel(L"y_{2}(t_f)")
fig.savefig("figures\\set_projection_time_discretization.pdf")
display(fig)

## Varying Truncation Order/Varying Time Discretization
pars = [(2,5), (4,5), (4,10), (6,10)]

contours = Dict()
sol_stats = Dict()
for p in pars
    m, nT = p[1], p[2]
    mp = MomentProblem(jump, m; x_scale=[x_scale])
    R = -unique(sort(round.(svd(mp.A).S, digits=3)))[1:min(nR,end)]
    T = range(Tf/nT, stop=Tf, length=nT)
    contours[p] = []
    sol_stats[p] = []
    obj = [(3,Dict([x=>1, x^2=>0]),nT)]
    sol, d, model = transient_bounds(mp, x0, R, T, n, obj, Mosek.Optimizer)
    for ϕ in ϕ_range
        obj = [(3,Dict([x=>sin(ϕ)/x_scale, x^2=>cos(ϕ)/x_scale^2]),nT)]
        sol, d, model = transient_bounds(mp, x0, R, T, n, obj, Mosek.Optimizer; d=d)
        push!(contours[p], [sol[obj[1]][2][2], sol[obj[1]][2][3]])
        push!(sol_stats[p], sol[obj[1]][3])
    end
end
fig, ax = subplots()
for p in pars
    ax.plot([c[1] for c in contours[p]], [c[2] for c in contours[p]],
             label=latexstring("n_{T} = ", p[2],", m = ",p[1]))
end
ax.plot(true_y1, true_y2, marker="*", markersize=10, color="k", linestyle="", label="true moments")
ax.legend()
ax.set_xlabel(L"y_{1}(t_f)")
ax.set_ylabel(L"y_{2}(t_f)")
fig.savefig("figures\\set_projection_joint.pdf")
display(fig)
