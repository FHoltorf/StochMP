using StochMP, LaTeXStrings, LinearAlgebra, DynamicPolynomials, PyPlot, MathOptInterface, JLD, MosekTools

function trajectory_bounds_poly(mp, x0, objs, r, nT, Tf_range; R=[])
    # basis function selection
    bounds, status, time = Dict(), Dict(), Dict()
    for obj in objs
        if typeof(obj[2]) <: Dict
            spec = [key for key in keys(obj[2])][1]
            idx = length(spec.name) > 1 ? parse(Int64, spec.name[3:end-1]) : 1
        else
            idx = length(obj[2].name) > 1 ? parse(Int64,obj[2].name[3:end-1]) : 1
        end
        for Tf in Tf_range
            if Tf == 0
                bounds[obj] = [obj[1] == 1 ? x0[idx]*ones(2) : 0.0]
                status[obj] = [obj[1] == 1 ? [MathOptInterface.FEASIBLE_POINT, MathOptInterface.FEASIBLE_POINT] :
                                              MathOptInterface.FEASIBLE_POINT]
                time[obj] = [obj[1] == 1 ? 0.0*ones(2) : 0.0]
            else
                T = range(Tf/nT, stop=Tf, length=nT)
                sol, m = poly_relaxation(mp, x0, T, r, objs, Mosek.Optimizer)
                push!(bounds[obj], sol[obj][1])
                push!(status[obj], typeof(obj[2]) <: Dict ? sol[obj][3] : sol[obj][2])
                push!(time[obj], typeof(obj[2]) <: Dict ? sol[obj][4] : sol[obj][3])
            end
        end
    end
    return bounds, status, time
end
