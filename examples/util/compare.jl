using StochMP, LaTeXStrings, LinearAlgebra, DynamicPolynomials, PyPlot, MathOptInterface, JLD, MosekTools

try readdir("results")
catch error
    @warn "created 'results' directory for output files"
    mkdir("results")
end
try readdir("figures")
catch error
    @warn "created 'figures' directory for output figures"
    mkdir("figures")
end

function trajectory_bounds(mp, x0, objs, n, nT, nR, Tf_range; R=[])
    # basis function selection
    if isempty(R)
        R = -unique(sort(round.(svd(mp.A).S, digits=3)))[1:min(nR,end)]
    end
    bounds, status, time = Dict(), Dict(), Dict()
    for obj in objs
        d = []
        if typeof(obj[2]) <: Dict
            spec = [key for key in keys(obj[2])][1]
            idx = length(spec.name) > 1 ? parse(Int64, spec.name[3:end-1]) : 1
        else
            idx = length(obj[2].name) > 1 ? parse(Int64,obj[2].name[3:end-1]) : 1
        end
        for Tf in Tf_range
            if Tf == 0
                bounds[obj] = [obj[1] == 1 ? x0[idx]*ones(2) : 0.0]
                status[obj] = [obj[1] == 1 ? [1,1] :
                                              1]
                time[obj] = [obj[1] == 1 ? 0.0*ones(2) : 0.0]
            else
                T = range(Tf/nT, stop=Tf, length=nT)
                sol, d, m = transient_bounds(mp, x0, R, T, n, [obj], Mosek.Optimizer; d=d)
                push!(bounds[obj], sol[obj][1])
                push!(status[obj], typeof(obj[2]) <: Dict ? sol[obj][3] : sol[obj][2])
                push!(time[obj], typeof(obj[2]) <: Dict ? sol[obj][4] : sol[obj][3])
            end
        end
    end
    return bounds, status, time
end
