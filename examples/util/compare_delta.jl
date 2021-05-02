using StochMP, LaTeXStrings, LinearAlgebra, DynamicPolynomials, PyPlot, MathOptInterface, JLD, MosekTools

function trajectory_bounds_δ(mp, x0, objs, n, nT, nR, Tf_range; R=[])
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
            idx = length(obj[2].name) > 1 ? parse(Int64, obj[2].name[3:end-1]) : 1
        end
        for Tf in Tf_range
            if Tf == 0
                bounds[obj] = [obj[1] == 1 ? x0[idx]*ones(2) : 0.0]
                status[obj] = [obj[1] == 1 ? [1, 1] :
                                              1]
                time[obj] = [obj[1] == 1 ? 0.0*ones(2) : 0.0]
            else
                T = range(Tf/nT, stop=Tf, length=nT)
                sol, d, m = transient_bounds_δ(mp, x0, R, T, n, [obj], Mosek.Optimizer; d=d)
                push!(bounds[obj], sol[obj][1])
                push!(status[obj], typeof(obj[2]) <: Dict ? sol[obj][3] : sol[obj][2])
                push!(time[obj], typeof(obj[2]) <: Dict ? sol[obj][4] : sol[obj][3])
            end
        end
    end
    return bounds, status, time
end

function parameter_sweep_δ(mp_fxn, x0, objs_raw, m_range, n_range, nT_range, nR_range, Tf_range; R=[])
    # basis function selection
    bounds, status, time = Dict(), Dict(), Dict()
    for m in m_range, n in n_range, nT in nT_range, nR in nR_range
        objs = [(obj...,nT) for obj in objs_raw]
        mp = mp_fxn(m)
        b, s, t = trajectory_bounds_δ(mp, x0, objs, n, nT, nR, Tf_range; R=R)
        bounds[m,n,nT,nR] = b
        status[m,n,nT,nR] = s
        time[m,n,nT,nR] = t
    end
    return bounds, status, time
end
