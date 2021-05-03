cd(@__DIR__)
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


parameter_sweep = true
if parameter_sweep
## Moment Problem Setup
m_min, Δm, m_max = 2, 2, 10
n_min, Δn, n_max= 2, 1, 10
nR_min, ΔnR, nR_max = 1, 1, 1
nT_min, ΔnT, nT_max = 2, 2, 20

m_range = m_min:Δm:m_max
n_range = n_min:Δn:n_max
nR_range = nR_min:ΔnR:nR_max
nT_range = nT_min:ΔnT:nT_max

R = []
Tf_range = [0.0, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 5.0, 7.5,
            10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]

objs = [(1,x),(2,x)]


dir = string("figures\\birth_death_parsweep")
try
    readdir(dir)
catch err
    mkdir(dir)
end

global results = Dict()
sampling_results = load("..\\sampling\\results\\birth_death_sampling.jld")
trange, x_MC, var_MC = sampling_results["trange"], sampling_results["x_MC"], sampling_results["var_MC"]


## traces for m, n fixed, nT increasing
m = m_min
n = n_min
nR = nR_max
bounds, status, time = parameter_sweep_δ(m -> MomentProblem(jump, m; x_scale=[7.0]),
                                              x0, objs, [m], [n], nT_range, [nR], Tf_range; R = R)

results[("traces", "nT")] = ([bounds, status, time], [[m], [n], nT_range, [nR], Tf_range])
fig, ax = subplots()
Z = [[0.01,0.01],[0.01,0.01]]
levels = (nT_min:ΔnT:nT_max+ΔnT)
CS3 = ax.contourf(Z, levels, cmap= PyPlot.cm.Reds)

fig, ax = subplots()
lns = []
rM = length(nT_min:ΔnT:nT_max)
r = 1
for nT in nT_min:ΔnT:nT_max
    obj = (1,x,nT)
    bnds = copy(bounds[m,n,nT,nR][obj])
    push!(lns, ax.plot(Tf_range, bnds, color = PyPlot.cm.Reds(0.2 + 0.6*r/rM), marker="o"))
    global r += 1
end
push!(lns, ax.plot(trange, x_MC, color="black", linestyle="dashed"))
ax.legend([lns[end][1]], ["sample average"])
ax.set_xlabel(L"t \ [s]")
ax.set_ylabel(L"⟨A(t)⟩")
cbar = fig.colorbar(CS3, label=L"n_{\mathsf{T}}", ticks= 1 .+ (nT_min:ΔnT:nT_max), format="%i")
cbar.ax.set_yticklabels(nT_min:ΔnT:nT_max)
fig.tight_layout()
fig.savefig(string(dir,"\\mean_nT.pdf"))
display(fig)

fig, ax = subplots()
lns = []
rM = length(nT_min:ΔnT:nT_max)
r = 1
for nT in nT_min:ΔnT:nT_max
    obj = (2,x,nT)
    bnds = bounds[m,n,nT,nR][obj]
    push!(lns, ax.plot(Tf_range, bnds, color = PyPlot.cm.Reds(0.2 + 0.6*r/rM), marker="o"))
    global r += 1
end
push!(lns, ax.plot(trange, var_MC, color="black", linestyle="dashed"))
ax.legend([lns[end][1]], ["sample variance"])
ax.set_xlabel(L"t \ [s]")
ax.set_ylabel(L"Var(A(t))")
cbar = fig.colorbar(CS3, label=L"n_{\mathsf{T}}", ticks= 1 .+ (nT_min:ΔnT:nT_max), format="%i")
cbar.ax.set_yticklabels(nT_min:ΔnT:nT_max)
fig.tight_layout()
fig.savefig(string(dir,"\\var_nT.pdf"))
display(fig)

## traces for nT,n fixed, m increasing
nT = nT_min
n = n_min
nR = nR_max
bounds, status, time = parameter_sweep_δ(m -> MomentProblem(jump, m; x_scale=[7.0]),
                            x0, objs, m_range, [n], [nT], [nR], Tf_range; R = R)


results[("traces", "m")] = ([bounds, status, time], [m_range, [n], [nT], [nR], Tf_range])
fig, ax = subplots()
Z = [[0.01,0.01],[0.01,0.01]]
levels = (m_min:Δm:m_max+Δm)
CS3 = ax.contourf(Z, levels, cmap= PyPlot.cm.Reds)

fig, ax = subplots()
lns = []
rM = length(m_min:Δm:m_max)
r = 1
for m in m_min:Δm:m_max
    obj = (1,x,nT)
    bnds = copy(bounds[m,n,nT,nR][obj])
    push!(lns, ax.plot(Tf_range, bnds, color = PyPlot.cm.Reds(0.2 + 0.6*r/rM), marker="o"))
    global r += 1
end
push!(lns, ax.plot(trange, x_MC, color="black", linestyle="dashed"))
ax.legend([lns[end][1]], ["sample average"])
ax.set_xlabel(L"t \ [s]")
ax.set_ylabel(L"⟨A(t)⟩")
cbar = fig.colorbar(CS3, label=L"m", ticks = 1.0 .+ (m_min:Δm:m_max), format="%i")
cbar.ax.set_yticklabels(m_min:Δm:m_max)
fig.tight_layout()
fig.savefig(string(dir,"\\mean_m.pdf"))
display(fig)

fig, ax = subplots()
lns = []
rM = length(m_min:Δm:m_max)
r = 1
for m in m_min:Δm:m_max
    obj = (2,x,nT)
    bnds = copy(bounds[m,n,nT,nR][obj])
    push!(lns, ax.plot(Tf_range, bnds, color = PyPlot.cm.Reds(0.2 + 0.6*r/rM), marker="o"))
    global r += 1
end
push!(lns, ax.plot(trange, var_MC, color="black", linestyle="dashed"))
ax.legend([lns[end][1]], ["sample variance"])
ax.set_xlabel(L"t \ [s]")
ax.set_ylabel(L"Var(A(t))")
cbar = fig.colorbar(CS3, label=L"m", ticks= 1.0 .+ (m_min:Δm:m_max), format="%i")
cbar.ax.set_yticklabels(m_min:Δm:m_max)
fig.tight_layout()
fig.savefig(string(dir,"\\var_m.pdf"))
display(fig)

## traces for m,n_T fixed, n increasing
m = m_min
nT = nT_min
nR = nR_max

bounds, status, time = parameter_sweep_δ(m -> MomentProblem(jump, m; x_scale=[7.0]),
                            x0, objs, [m], n_range, [nT], [nR], Tf_range; R = R)

results[("traces", "n")] = ([bounds, status, time], [ [m], n_range, [nT], [nR], Tf_range])

fig, ax = subplots()
Z = [[0.01,0.01],[0.01,0.01]]
levels = (n_min:Δn:n_max+Δn)
CS3 = ax.contourf(Z, levels, cmap= PyPlot.cm.Reds)

fig, ax = subplots()
lns = []
rM = length(n_min:Δn:n_max)
r = 1
for n in n_min:Δn:n_max
    obj = (1,x,nT)
    bnds = copy(bounds[m,n,nT,nR][obj])
    push!(lns, ax.plot(Tf_range, bnds, color = PyPlot.cm.Reds(0.2 + 0.6*r/rM), marker="o"))
    global r += 1
end
push!(lns, ax.plot(trange, x_MC, color="black", linestyle="dashed"))
ax.legend([lns[end][1]], ["sample average"])
ax.set_xlabel(L"t \ [s]")
ax.set_ylabel(L"⟨A(t)⟩")
cbar = fig.colorbar(CS3, label=L"n_I", ticks= 0.5 .+ (n_min:Δn:n_max), format="%i")
cbar.ax.set_yticklabels((n_min:Δn:n_max))
fig.tight_layout()
fig.savefig(string(dir,"\\mean_n.pdf"))
display(fig)

fig, ax = subplots()
lns = []
rM = length(n_min:Δn:n_max)
r = 1
for n in n_min:Δn:n_max
    obj = (2,x,nT)
    bnds = copy(bounds[m,n,nT,nR][obj])
    push!(lns, ax.plot(Tf_range, bnds, color = PyPlot.cm.Reds(0.2 + 0.6*r/rM), marker="o"))
    global r += 1
end
push!(lns, ax.plot(trange, var_MC, color="black", linestyle="dashed"))
ax.legend([lns[end][1]], ["sample variance"])
ax.set_xlabel(L"t \ [s]")
ax.set_ylabel(L"Var(A(t))")
cbar = fig.colorbar(CS3, label=L"n_I", ticks=0.5 .+ (n_min:Δn:n_max), format="%i")
cbar.ax.set_yticklabels((n_min:Δn:n_max))
fig.tight_layout()
fig.savefig(string(dir,"\\var_n.pdf"))
display(fig)

## HEATMAPS
## max_{t_i} (ub-lb) heat map for and m and nT, n fixed
n = n_min
nR = nR_max

bounds, status, time = parameter_sweep_δ(m -> MomentProblem(jump, m; x_scale=[7.0]),
                                         x0, objs, m_range, [n], nT_range, [nR], Tf_range; R = R)

results[("heatmap", "m", "nT")] = ([bounds, status, time], [m_range, [n], nT_range, [nR], Tf_range])

t = length(Tf_range)
i, j = 1, 1
δ = zeros(length(m_min:Δm:m_max), length(nT_min:ΔnT:nT_max))
for m in m_min:Δm:m_max
    global j = 1
    for nT in nT_min:ΔnT:nT_max
        obj = (1,x,nT)
        bnds = copy(bounds[m,n,nT,nR][obj])
        δ[i,j] = maximum([b[2] for b in bnds] .- [b[1] for b in bnds])
        global j += 1
    end
    global i += 1
end
fig, ax = subplots()
c = ax.imshow(δ, cmap="RdBu", vmin=0, origin="lower", aspect="auto")
ax.set_yticks(0:length(m_range)-1)
ax.set_xticks(0:length(nT_range)-1)
ax.set_yticklabels(m_range)
ax.set_xticklabels(nT_range)
ax.set_ylabel(L"m")
ax.set_xlabel(L"n_{\mathsf{T}}")
fig.colorbar(c, ax = ax, label=L"\max_{i} \quad ⟨A(t_i)⟩^U - ⟨A(t_i)⟩^L")
fig.tight_layout()
fig.savefig(string(dir,"\\heatmap_m_nT.pdf"))
display(fig)

## max_{t_i} (ub-lb) heat map for and m and n, nT fixed
nT = nT_min
nR = nR_max

bounds, status, time = parameter_sweep_δ(m -> MomentProblem(jump, m; x_scale=[7.0]),
                                         x0, objs, m_range, n_range, [nT], [nR], Tf_range; R = R)

results[("heatmap", "m", "n")] = ([bounds, status, time], [m_range, n_range, [nT], [nR], Tf_range])
t = length(Tf_range)
i, j = 1, 1
δ = zeros(length(m_min:Δm:m_max), length(n_min:Δn:n_max))
for m in m_min:Δm:m_max
    global j = 1
    for n in n_min:Δn:n_max
        obj = (1,x,nT)
        bnds = copy(bounds[m,n,nT,nR][obj])
        δ[i,j] = maximum([b[2] for b in bnds] .- [b[1] for b in bnds])
        global j += 1
    end
    global i += 1
end
fig, ax = subplots()
c = ax.imshow(δ, cmap="RdBu", vmin=0, origin="lower", aspect="auto")
ax.set_yticks(0:length(m_range)-1)
ax.set_xticks(0:length(n_range)-1)
ax.set_yticklabels(m_range)
ax.set_xticklabels(n_range)
ax.set_ylabel(L"m")
ax.set_xlabel(L"n_I")
fig.colorbar(c, ax = ax, label=L"\max_{i} \quad ⟨A(t_i)⟩^U - ⟨A(t_i)⟩^L")
fig.tight_layout()
fig.savefig(string(dir,"\\heatmap_m_n.pdf"))
display(fig)

## max_{t_i} (ub-lb) heat map for and n and nT, m fixed
m = m_min
nR = nR_max

bounds, status, time = parameter_sweep_δ(m -> MomentProblem(jump, m; x_scale=[7.0]),
                                         x0, objs, [m], n_range, nT_range, [nR], Tf_range; R = R)

results[("heatmap", "n", "nT")] = ([bounds, status, time], [[m], n_range, nT_range, [nR], Tf_range])

t = length(Tf_range)
i, j = 1, 1

δ = zeros(length(n_min:Δn:n_max), length(nT_min:ΔnT:nT_max))
for n in n_min:Δn:n_max
    global j = 1
    for nT in nT_min:ΔnT:nT_max
        obj = (1,x,nT)
        bnds = copy(bounds[m,n,nT,nR][obj])
        δ[i,j] = maximum([b[2] for b in bnds] .- [b[1] for b in bnds])
        global j += 1
    end
    global i += 1
end
fig, ax = subplots()
c = ax.imshow(δ, cmap="RdBu", vmin=0, origin="lower", aspect="auto")
ax.set_xticks(0:length(nT_range)-1)
ax.set_yticks(0:length(n_range)-1)
ax.set_xticklabels(nT_range)
ax.set_yticklabels(n_range)
ax.set_xlabel(L"n_{\mathsf{T}}")
ax.set_ylabel(L"n_{I}")
fig.colorbar(c, ax = ax, label=L"\max_{i} \quad ⟨A(t_i)⟩^U - ⟨A(t_i)⟩^L")
fig.tight_layout()
fig.savefig(string(dir,"\\heatmap_n_nT.pdf"))
display(fig)

io = open(string(dir,"\\readme.txt"), "w")
println(io, "m_min = ", m_min)
println(io, "delta_m = ", Δm)
println(io, "m_max = ", m_max)
println(io, "nT_min = ", nT_min)
println(io, "delta_nT = ", ΔnT)
println(io, "nT_max = ", nT_max)
println(io, "n_min = ", n_min)
println(io, "delta_n = ", Δn)
println(io, "n_max = ", n_max)
println(io, "m gradient figures are at n = ", n_min,", nT = ", nT_min, ", nR = ", nR_max)
println(io, "n gradient figures are at nT = ", nT_min, ", m = ", m_min, ", nR = ", nR_max)
println(io, "nT gradient figures are at n = ", n_min, ", m = ", m_min, ", nR = ", nR_max)
println(io, "heatmap m, nT is at n = ", n_min, ", nR = ", nR_max)
println(io, "heatmap n, nT is at n = ", m_min, ", nR = ", nR_max)
println(io, "heatmap m, n is at nT = ", nT_min, ", nR = ", nR_max)
close(io)

save("results/birth_death_parsweep.jld", "results", results,
                                         "n_range", n_range,
                                         "nT_range", nT_range,
                                         "m_range", m_range,
                                         "Tf_range", Tf_range,
                                         "nR_range", nR_range,
                                         "x0", x0,
                                         "R", R)

#
else
## Moment Problem Setup
# USER CHOICES!!!
m = 8
n = 2
nR = 3
nT = 20
Tf_range = [0.0, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 5.0, 7.5,
            10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]

# Assemble Moment Problem
mp = MomentProblem(jump, m; x_scale=[7.0])

## Bound trajectories
# Holtorf, Barton
objs_HB = [(1,x,nT),(2,x,nT)]
bounds_HB, status_HB, time_HB = trajectory_bounds_δ(mp, x0, objs_HB, n, nT, nR, Tf_range)

# Dowdy, Barton
objs_DB = [(1,x,1),(2,x,1)]
bounds_DB, status_DB, time_DB = trajectory_bounds_δ(mp, x0, objs_DB, 1, 1, nR, Tf_range)

## Save/Load results
save("results/birth_death.jld", "user_choices", (m,n,nT,nR,Tf_range),
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
end
