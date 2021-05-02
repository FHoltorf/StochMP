module StochMP

using LinearAlgebra, JuMP, SumOfSquares, DynamicPolynomials, MultivariatePolynomials

export ReachableSet, StochasticProcess, StochasticOCP, MomentProblem,
       ItoProcess, JumpProcess, JumpDiffProcess, jump_χ,
       transient_bounds, transient_bounds_γ, transient_bounds_δ, stationary_bounds,
       poly_restriction, poly_relaxation

abstract type AbstractStochasticProcess end
abstract type AbstractReachableSet end
abstract type AbstractMomentProblem end
include("util.jl")

mutable struct ReachableSet{T1,T2<:AbstractPolynomialLike} <: AbstractReachableSet
      g::Array{T1,1}
      h::Array{T2,1}
end

mutable struct StochasticProcess{F<:Function} <: AbstractStochasticProcess
      nx::Int64
      ny::Int64
      x::Array{<:AbstractPolynomialLike,1}
      y::Array{<:AbstractPolynomialLike,1}
      q::Int64
      G::F
      χ::ReachableSet
      class::String
end

mutable struct MomentProblem{T1,T2,T3 <: Number} <: AbstractMomentProblem
      nL::Int64
      nH::Int64
      b2s::Dict # basis
      s2b::Dict #
      D::Union{Diagonal,UniformScaling{Bool}}
      A::Array{T1,2} # dynamics
      B::Array{Array{T2,2},1} # algebraic constraints
      M::Array{Array{Array{T3,1},2},1} # LMIs
      Bdim::Array{Tuple{Int64,Int64}}
      Mdim::Array{Tuple{Int64,Int64}}
end

function ReachableSet(; g::Array{<:AbstractPolynomialLike,1}=Array{Polynomial,1}(),
                        h::Array{<:AbstractPolynomialLike,1}=Array{Polynomial,1}())
      return ReachableSet(g,h)
end

function MomentProblem(sp::StochasticProcess, m::Int64;
                       x_scale=ones(sp.nx), y_scale=ones(sp.ny))
      z = isempty(sp.y) ? sp.x : cat(sp.x,sp.y,dims=1)
      z_scale = cat(x_scale,y_scale,dims=1)

      #basis
      bx, bxL = sort(monomials(sp.x,0:m+sp.q)), sort(monomials(sp.x,0:m))
      by, byL = sort(monomials(z,0:m+sp.q,p -> !all(degree.(p,sp.y) .== 0))),
                sort(monomials(z,0:m,p -> !all(degree.(p,sp.y) .== 0)))
      bz, bzL = cat(bx,by,dims=1), cat(bxL,byL,dims=1)
      nL = length(bxL)
      nH = length(bz) - nL

      # scaling factor
      D = all(z_scale .== 1) ? I : Diagonal([i <= nL ? 1.0*bz[i](x_scale...) : 1.0*bz[i](z_scale...)
                                             for i in 1:nL+nH])

      # dynamics
      A = zeros(nL, nL+nH)
      for i in 1:nL
            A[i,:] = coefficients(sp.G(exponents(bxL[i])),bz)
      end
      A = D==I ? A : inv(Diagonal(D[1:nL,1:nL]))*A*D

      # algebraic constraints
      B = Array{Array{Float64,2},1}()
      Bdim = Array{Tuple{Int64,Int64},1}()
      for h in sp.χ.h
            push!(B, assemble_alg_con(h,m+sp.q,z,bz; D=D))
            push!(Bdim, (binomial(m-maxdegree(h)+sp.nx, sp.nx), size(B[end],1)))
      end

      # support
      M = [assemble_lmi_con(sp.x[1]^0,m+sp.q,sp.x,sp.y,bz; D=D)]
      Mdim = [(binomial(floor(Int64,m/2)+sp.nx,sp.nx), size(M[1],1))]
      for g in sp.χ.g
            push!(M, assemble_lmi_con(g,m+sp.q,sp.x,sp.y,bz; D=D))
            push!(Mdim, (binomial(floor(Int64,(m-maxdegree(g))/2)+sp.nx,sp.nx), size(M[end],1)))
      end
      return MomentProblem(nL,nH,Dict([bz[i] => i for i in 1:length(bz)]),
                                 Dict([i => bz[i] for i in 1:length(bz)]),
                                 D,A,B,M,Bdim,Mdim)
end

function transient_bounds(mp::MomentProblem, x0, R, T, n::Int64, objs, solver; d=[])
      # convenient shorthands
      nT = length(T)
      Tf = T[end]
      # initial moments
      μ0 = [mp.s2b[i](x0...)/mp.D[i,i] for i in 1:mp.nL]

      # setting up model
      m = Model(solver)
      if isempty(d)
            μ = Dict([t => [@variable(m, base_name=var_name(t,i)) for i in 1:mp.nL]
                      for t in 1:nT])
            z = Dict([(ρ,l,t) => (l > 0 ? [T[t]^(l-1)/factorial(l-1)*@variable(m, base_name=var_name(ρ,l,t,i)) for i in 1:mp.nL+mp.nH] :
                                           exp(ρ*(Tf-T[t]))*μ[t])
                      for ρ in R, l in 0:n, t in 1:nT])

      else
            μ = Dict([t => [d[t][i]*@variable(m, base_name=var_name(t,i)) for i in 1:mp.nL]
                      for t in 1:nT])
            z = Dict([(ρ,l,t) => (l > 0 ? [d[ρ,l,t][i]*T[t]^(l-1)/factorial(l-1)*@variable(m, base_name=var_name(ρ,l,t,i)) for i in 1:mp.nL+mp.nH] :
                                          exp(ρ*(Tf-T[t]))*μ[t])
                      for ρ in R, l in 0:n, t in 1:nT])
      end

      # dynamics
      for ρ in R, l in 1:n, t in 1:nT, i in 1:mp.nL
            if l == 1 && ρ == 0 && i == 1
                  @constraint(m, z[ρ,n,t][1] == T[t]^n/factorial(n))
                  @constraint(m, μ[t][1] == 1)
            else
                  if t == 1
                        Δt = T[t]
                        scale = 1#C(0,T[t],l-1,0)
                        @constraint(m, (mp.A[i,:]'*z[ρ,l,t] - ρ*z[ρ,l,t][i]) ./ scale  == (z[ρ,l-1,t][i] - T[t]^(l-1)/factorial(l-1)*exp(ρ*Tf)*μ0[i]) ./ scale )
                  else
                        Δt = T[t] - T[t-1]
                        scale = 1 #C(T[t-1],T[t],l,0)
                        @constraint(m, (mp.A[i,:]'*(z[ρ,l,t]-(T[t]/T[t-1])^(l-1)*z[ρ,l,t-1]) - ρ*(z[ρ,l,t][i] - (T[t]/T[t-1])^(l-1)*z[ρ,l,t-1][i])) ./ scale ==
                                       (z[ρ,l-1,t][i] - (T[t]/T[t-1])^(l-1)*z[ρ,l-1,t-1][i]) ./ scale)
                  end
            end
      end

      # closure

      ## Support
      # affine constraints
      if !isempty(mp.B)
            @constraint(m, [l in 1:n, ρ in R, t in 1:nT, i in 1:length(mp.B)],
                           mp.B[i]*z[ρ,l,t] .== 0)
            @constraint(m, [t in 1:nT, i in 1:length(mp.B)],
                           mp.B[i][1:mp.Bdim[i][1],1:mp.nL]*μ[t] .== 0)
      end

      # LMIs
      for p in 1:length(mp.M)
            if mp.Mdim[p][1] > 1
                  @constraint(m, [t in 1:nT], set_LMI(mp.M[p],μ[t];dim=mp.Mdim[p][1]) in PSDCone())
            else
                  @constraint(m, [t in 1:nT], set_LMI(mp.M[p],μ[t];dim=mp.Mdim[p][1]) .>= 0)
            end
            for l in 1:n-1
                  for k in 0:l
                        @constraint(m, [ρ in R, t in 1:nT],
                                       (t > 1 ? set_LMI(mp.M[p], Ω(k,l-k,T[t-1],T[t],
                                                        [z[ρ,i,t-1] for i in 1:l+1],
                                                        [z[ρ,i,t] for i in 1:l+1])) :
                                                set_LMI(mp.M[p], Ω(k,l-k,0.0,T[t],
                                                        [zeros(mp.nL+mp.nH) for i in 1:l+1],
                                                        [z[ρ,i,t] for i in 1:l+1])))
                                                in PSDCone())
                  end
            end
            if n == 1 # recover Dowdy & Barton
                  @constraint(m, [ρ in R, t in 1:nT],
                                 (t > 1 ? set_LMI(mp.M[p],z[ρ,1,t] - z[ρ,1,t-1]) :
                                          set_LMI(mp.M[p],z[ρ,1,t])) in PSDCone())
            end
      end

      # objectives
      sol = Dict()
      for obj in objs
            if obj[1] == 1 # means
                  s = mp.b2s[obj[2]]
                  t = obj[3]
                  @objective(m, Min, μ[t][s])
                  optimize!(m)
                  lb, stat_lb, t_lb = get_results(m; scale=mp.D[s,s])

                  @objective(m, Max, μ[t][s])
                  optimize!(m)
                  ub, stat_ub, t_ub = get_results(m; scale=mp.D[s,s])

                  sol[obj] = ([lb,ub], [stat_lb,stat_ub], [t_lb,t_ub])
            elseif obj[1] == 2 # variance
                  s1 = mp.b2s[obj[2]]
                  s2 = mp.b2s[obj[2]^2]
                  t = obj[3]
                  slack = @variable(m)
                  scale = isempty(d) ? 1 : d[t][s2]
                  @constraint(m, [μ[t][s2]-slack μ[t][s1]
                                  μ[t][s1] 1] ./ mp.D[s2,s2] in PSDCone())
                  @objective(m, Max, slack)
                  optimize!(m)

                  sol[obj] = get_results(m; scale= mp.D[s2,s2])
            elseif obj[1] == 3 # tailored lin. combination
                  @objective(m, Min,
                                sum(mp.D[mp.b2s[i],mp.b2s[i]]*obj[2][i]*μ[obj[3]][mp.b2s[i]]
                                    for i in keys(obj[2])) )
                  optimize!(m)
                  sol[obj] = get_results(m,mp.D[1:mp.nL,1:mp.nL]*μ[obj[3]])
            else
                  @error "objective option not implemented"
            end
      end
      # get scaling factors here
      d = get_scaling_factors(μ,z,d)
      return sol, d, m
end

# Averaging z by polynomial test function, implied by the formulation beneath
function transient_bounds_γ(mp::MomentProblem, x0, R, T, n::Int64, objs, solver; d=[])
      # convenient shorthands
      nT = length(T)
      Tf = T[end]
      n -= 2

      # initial moments
      μ0 = [mp.s2b[i](x0...)/mp.D[i,i] for i in 1:mp.nL]

      # setting up model
      m = Model(solver)
      if isempty(d)
            μ = Dict([t => [@variable(m, base_name=var_name(t,i)) for i in 1:mp.nL]
                      for t in 1:nT])
            z = Dict([(ρ,t) => [@variable(m, base_name=var_name(ρ,t,i)) for i in 1:mp.nL+mp.nH] for ρ in R, t in 1:nT])
            γ = Dict([(ρ,k,l,t) => [@variable(m, base_name=var_name(ρ,k,l,t,i)) for i in 1:mp.nL+mp.nH] for ρ in R, t in 1:nT, l in 0:n, k in 0:n])
      else
            μ = Dict([t => [d[t][i]*@variable(m, base_name=var_name(t,i)) for i in 1:mp.nL]
                      for t in 1:nT])
            z = Dict([(ρ,t) => [d[ρ,t][i]*@variable(m, base_name=var_name(ρ,t,i)) for i in 1:mp.nL+mp.nH] for ρ in R, t in 1:nT])
            γ = Dict([(ρ,k,l,t) => [d[ρ,k,l,t][i]*@variable(m, base_name=var_name(ρ,k,l,t,i)) for i in 1:mp.nL+mp.nH] for ρ in R, t in 1:nT, l in 0:n, k in 0:n])
      end

      # dynamics
      for ρ in R, t in 1:nT, i in 1:mp.nL
            if ρ == 0 && i == 1
                  @constraint(m, μ[t][1] == 1)
                  @constraint(m, z[ρ,t][1] == T[t])
            else
                  @constraint(m, mp.A[i,:]'*z[ρ,t] - ρ*z[ρ,t][i] == exp(ρ*(Tf-T[t]))*μ[t][i] - exp(ρ*Tf)*μ0[i])
            end
            for l in 0:n, k in 0:n
                  t1, t2 = (t > 1 ? T[t-1] : 0), T[t]
                  if ρ == 0 && i == 1
                        @constraint(m, γ[ρ,k,l,t][i] == ((1+k)*t1+(1+l)*t2)/(2+k+l))
                  else
                        rhs = AffExpr(0.0)
                        if k == 0
                              rhs += (t2-t1)^l/C(t1,t2,k,l)*z[ρ,t][i]
                        end
                        if l == 0 && t > 1
                              rhs -= (t2-t1)^k/C(t1,t2,k,l)*z[ρ,t-1][i]
                        end
                        if k > 0
                              rhs += k*C(t1,t2,k-1,l)/C(t1,t2,k,l)*γ[ρ,k-1,l,t][i]
                        end
                        if l > 0
                              rhs -= l*C(t1,t2,k,l-1)/C(t1,t2,k,l)*γ[ρ,k,l-1,t][i]
                        end
                        @constraint(m, mp.A[i,:]'*γ[ρ,k,l,t] - ρ*γ[ρ,k,l,t][i] == rhs - exp(ρ*Tf)*μ0[i] )
                  end
            end
      end

      # closure

      ## Support
      # affine constraints
      if !isempty(mp.B)
            @constraint(m, [ρ in R, t in 1:nT, i in 1:length(mp.B)],
                           mp.B[i]*z[ρ,t] .== 0)
            @constraint(m, [k in 0:n, l in 0:n, ρ in R, t in 1:nT, i in 1:length(mp.B)],
                           mp.B[i]*γ[ρ,k,l,t] .== 0)
            @constraint(m, [t in 1:nT, i in 1:length(mp.B)],
                           mp.B[i][1:mp.Bdim[i][1],1:mp.nL]*μ[t] .== 0)
      end

      # LMIs
      for p in 1:length(mp.M)
            if mp.Mdim[p][1] > 1
                  @constraint(m, [t in 1:nT], set_LMI(mp.M[p],μ[t];dim=mp.Mdim[p][1]) in PSDCone())
            else
                  @constraint(m, [t in 1:nT], set_LMI(mp.M[p],μ[t];dim=mp.Mdim[p][1]) .>= 0)
            end
            for ρ in R, t in 1:nT
                  z1, z2 = (t > 1 ? z[ρ,t-1] : zeros(mp.nL+mp.nH)), z[ρ,t]
                  for k in 0:n, l in 0:n
                        @constraint(m, (t > 1 ? set_LMI(mp.M[p], γ[ρ,k,l,t] - z[ρ,t-1]) : set_LMI(mp.M[p], γ[ρ,k,l,t])) in PSDCone())
                        @constraint(m, set_LMI(mp.M[p], z[ρ,t] - γ[ρ,k,l,t]) in PSDCone())
                  end
                  if n < 0
                        @constraint(m, (t > 1 ? set_LMI(mp.M[p], z[ρ,t] - z[ρ,t-1]) : set_LMI(mp.M[p], z[ρ,t])) in PSDCone())
                  end
            end
      end

      # objectives
      sol = Dict()
      for obj in objs
            if obj[1] == 1 # means
                  s = mp.b2s[obj[2]]
                  t = obj[3]
                  @objective(m, Min, μ[t][s])
                  optimize!(m)
                  lb, stat_lb, t_lb = get_results(m; scale=mp.D[s,s])

                  @objective(m, Max, μ[t][s])
                  optimize!(m)
                  ub, stat_ub, t_ub = get_results(m; scale=mp.D[s,s])

                  sol[obj] = ([lb,ub], [stat_lb,stat_ub], [t_lb,t_ub])
            elseif obj[1] == 2 # variance
                  s1 = mp.b2s[obj[2]]
                  s2 = mp.b2s[obj[2]^2]
                  t = obj[3]
                  slack = @variable(m)
                  @constraint(m, [μ[t][s2]-slack μ[t][s1]
                                  μ[t][s1] 1] in PSDCone())
                  @objective(m, Max, slack)
                  optimize!(m)

                  sol[obj] = get_results(m; scale=mp.D[s2,s2])
            elseif obj[1] == 3 # tailored lin. combination
                  @objective(m, Min,
                                sum(mp.D[mp.b2s[i],mp.b2s[i]]*obj[2][i]*μ[obj[3]][mp.b2s[i]]
                                    for i in keys(obj[2])) )
                  optimize!(m)
                  sol[obj] = get_results(m,mp.D[1:mp.nL,1:mp.nL]*μ[obj[3]])
            else
                  @error "objective option not implemented"
            end
      end
      # get scaling factors here
      d = get_scaling_factors(μ,z,γ,d)
      return sol, d, m
end

function transient_bounds_δ(mp::MomentProblem, x0, R, T, n::Int64, objs, solver; d=[])
      # convenient shorthands
      nT = length(T)
      Tf = T[end]
      n -= 1

      # initial moments
      μ0 = [mp.s2b[i](x0...)/mp.D[i,i] for i in 1:mp.nL]

      # setting up model
      m = Model(solver)
      if isempty(d)
            μ = Dict([t => [@variable(m, base_name=var_name(t,i)) for i in 1:mp.nL]
                      for t in 1:nT])
            z = Dict([(ρ,k,l,t) => [@variable(m, base_name=var_name(ρ,k,l,t,i)) for i in 1:mp.nL+mp.nH]
                      for ρ in R, l in 0:n, k in 0:n, t in 1:nT if min(max(1-l,0),n) <= k <= n-l])
      else
            μ = Dict([t => [d[t][i]*@variable(m, base_name=var_name(t,i)) for i in 1:mp.nL]
                      for t in 1:nT])
            z = Dict([(ρ,k,l,t) => [d[ρ,k,l,t][i]*@variable(m, base_name=var_name(ρ,k,l,t,i)) for i in 1:mp.nL+mp.nH]
                      for ρ in R, l in 0:n, k in 0:n, t in 1:nT if min(max(1-l,0),n) <= k <= n-l])
      end

      # dynamics
      for t in 1:nT
            t1, t2 = (t > 1 ? T[t-1] : 0), T[t]
            @constraint(m, μ[t][1] == 1)
            for ρ in R, i in 1:mp.nL, k in 0:n, l in min(max(1-k,0),n):n-k
                  if ρ == 0 && i == 1
                        @constraint(m, z[ρ,k,l,t][1] == 1)
                  else
                        rhs = AffExpr(0.0)
                        if k == 0
                              rhs += (t2-t1)^l*exp(ρ*(Tf-t2))/C(t1,t2,l,k)*μ[t][i]
                        elseif k == 1 && l == 0
                              rhs += k*C(t1,t2,l,k-1)/C(t1,t2,l,k)*(z[ρ,0,1,t][i] + z[ρ,1,0,t][i])/2
                        else
                              rhs += k*C(t1,t2,l,k-1)/C(t1,t2,l,k)*z[ρ,k-1,l,t][i]
                        end
                        if l == 0
                              rhs -= (t2-t1)^k*exp(ρ*(Tf-t1))/C(t1,t2,l,k)*(t > 1 ? μ[t-1][i] : μ0[i])
                        elseif l == 1 && k == 0
                              rhs -= l*C(t1,t2,l-1,k)/C(t1,t2,l,k)*(z[ρ,0,1,t][i] + z[ρ,1,0,t][i])/2
                        else
                              rhs -= l*C(t1,t2,l-1,k)/C(t1,t2,l,k)*z[ρ,k,l-1,t][i]
                        end
                        rhs += ρ*z[ρ,k,l,t][i]
                        @constraint(m, mp.A[i,:]'*z[ρ,k,l,t] == rhs)
                  end
            end
      end

      ## Support
      # affine constraints
      if !isempty(mp.B)
            @constraint(m, [ρ in R, k in 0:n, l in min(max(1-k,0),n):n-k, t in 1:nT, i in 1:length(mp.B)],
                           mp.B[i]*z[ρ,k,l,t] .== 0)
            @constraint(m, [t in 1:nT, i in 1:length(mp.B)],
                           mp.B[i][1:mp.Bdim[i][1],1:mp.nL]*μ[t] .== 0)
      end

      # LMIs
      for p in 1:length(mp.M)
            if mp.Mdim[p][1] > 1
                  @constraint(m, [t in 1:nT], set_LMI(mp.M[p],μ[t];dim=mp.Mdim[p][1]) in PSDCone())
            else
                  @constraint(m, [t in 1:nT], set_LMI(mp.M[p],μ[t];dim=mp.Mdim[p][1]) .>= 0)
            end
            @constraint(m, [ρ in R, k in 0:n, l in min(max(1-k,0),n):n-k, t in 1:nT],
                           set_LMI(mp.M[p],z[ρ,k,l,t]) in PSDCone())
      end

      # objectives
      sol = Dict()
      for obj in objs
            if obj[1] == 1 # means
                  s = mp.b2s[obj[2]]
                  t = obj[3]
                  @objective(m, Min, μ[t][s])
                  optimize!(m)
                  lb, stat_lb, t_lb = get_results(m; scale=mp.D[s,s])

                  @objective(m, Max, μ[t][s])
                  optimize!(m)
                  ub, stat_ub, t_ub = get_results(m; scale=mp.D[s,s])

                  sol[obj] = ([lb,ub], [stat_lb,stat_ub], [t_lb,t_ub])
            elseif obj[1] == 2 # variance
                  s1 = mp.b2s[obj[2]]
                  s2 = mp.b2s[obj[2]^2]
                  t = obj[3]
                  slack = @variable(m)
                  @constraint(m, [μ[t][s2]-slack μ[t][s1]
                                  μ[t][s1] 1] in PSDCone())
                  @objective(m, Max, slack)
                  optimize!(m)

                  sol[obj] = get_results(m; scale=mp.D[s2,s2])
            elseif obj[1] == 3 # tailored lin. combination
                  @objective(m, Min,
                                sum(mp.D[mp.b2s[i],mp.b2s[i]]*obj[2][i]*μ[obj[3]][mp.b2s[i]]
                                    for i in keys(obj[2])) )
                  optimize!(m)
                  sol[obj] = get_results(m, mp.D[1:mp.nL,1:mp.nL]*μ[obj[3]])
            else
                  @error "objective option not implemented"
            end
      end
      # get scaling factors here
      d = get_scaling_factors(μ,z,d)
      return sol, d, m
end

function poly_restriction(mp::MomentProblem, x0, T, r, objs, solver)
      # convenient shorthands
      nT = length(T)
      Δt = [k > 1 ? T[k]-T[k-1] : T[1] for k in 1:nT]
      println(Δt)

      # initial condition
      μ0 = [mp.s2b[i](x0...)/mp.D[i,i] for i in 1:mp.nL+mp.nH]
      println(μ0)

      @polyvar(t)
      @polyvar(y[1:maximum([d[2] for d in mp.Mdim])])
      bt = sort(monomials(t, 0:r))
      bt_red = isodd(r) ? (bt[1:r], bt[1:r]) : (bt[1:r+1], bt[1:r-1])
      bt_sos = isodd(r) ? (bt[1:div(r+1,2)], bt[1:div(r+1,2)]) : (bt[1:div(r,2)+1], bt[1:div(r,2)])

      m = Model(solver)
      μ = Dict([k => @variable(m, [1:mp.nL+mp.nH], Poly(bt)) for k in 1:nT])
      Y, Q = Dict(), Dict()
      for i in 1:length(mp.M)
            dim = mp.Mdim[i][2]
            bty = (cat([y[k].*bt_sos[1] for k in 1:dim]..., dims=1),
                   cat([y[k].*bt_sos[2] for k in 1:dim]..., dims=1))
            for k in 1:nT
                  Y[i,k] = (@variable(m, [1:dim,1:dim], Poly(bt_red[1])),
                            @variable(m, [1:dim,1:dim], Poly(bt_red[2])))
                  Q[i,k] = (@variable(m, [1], SOSPoly(bty[1]))[1],
                            @variable(m, [1], SOSPoly(bty[2]))[1])
                  @constraint(m, y[1:dim]'*Y[i,k][1]*y[1:dim] == Q[i,k][1])
                  @constraint(m, y[1:dim]'*Y[i,k][2]*y[1:dim] == Q[i,k][2])
                  @constraint(m, set_LMI(mp.M[i],μ[k]) .==  Λ(Y[i,k], t, r, Δt[k]) )
            end
      end

      # dynamics
      @constraint(m, ode[k in 1:nT, i in 1:mp.nL],
                     ∇(μ[k][i],t) == mp.A[i,:]'*μ[k])
      @constraint(m, continuity[k in 1:nT, i in 1:mp.nL],
                     μ[k][i](0.0) == (k > 1 ? μ[k-1][i](Δt[k-1]) : μ0[i]))

      # affine constraints
      if !isempty(mp.B)
            @constraint(m, [k in 1:nT, i in 1:length(mp.B)], mp.B[i]*μ[k] .== 0)
      end

      sol = Dict()
      for obj in objs
            if obj[1] == 1 # means
                  s = mp.b2s[obj[2]]
                  k = findfirst(x -> x >= obj[3], T)
                  t_obj = obj[3] - (k > 1 ? T[k-1] : 0)
                  println(t_obj)
                  println(k)
                  @objective(m, Min, μ[k][s](t_obj))
                  optimize!(m)
                  lb, stat_lb, t_lb = get_results(m; scale=mp.D[s,s])

                  @objective(m, Max, μ[k][s](t_obj))
                  optimize!(m)
                  ub, stat_ub, t_ub = get_results(m; scale=mp.D[s,s])

                  sol[obj] = ([lb,ub], [stat_lb,stat_ub], [t_lb,t_ub])
            elseif obj[1] == 2 # variance
                  s1, s2 = mp.b2s[obj[2]], mp.b2s[obj[2]^2]
                  k = findfirst(x -> x >= obj[3], T)
                  t_obj = obj[3] - (k > 1 ? T[k-1] : 0)
                  println(t_obj)
                  println(k)
                  slack = @variable(m)
                  @constraint(m, [μ[k][s2](t_obj)-slack μ[k][s1](t_obj)
                                  μ[k][s1](t_obj) 1] in PSDCone())
                  @objective(m, Max, slack)
                  optimize!(m)

                  sol[obj] = get_results(m; scale=mp.D[s2,s2])
            elseif obj[1] == 3 # tailored lin. combination
                  k = findfirst(x -> x >= obj[3], T)
                  t_obj = obj[3] - (k > 1 ? T[k-1] : 0)
                  @objective(m, Min,
                                sum(mp.D[mp.b2s[i],mp.b2s[i]]*obj[2][i]*μ[k][mp.b2s[i]](t_obj)
                                    for i in keys(obj[2])) )
                  optimize!(m)

                  sol[obj] = get_results(m)
            else
                  @error "objective option not supported yet"
            end
      end
      return sol, m
end

function poly_relaxation(mp::MomentProblem, x0, T, r, objs, solver)
      # convenient shorthands
      nT = length(T)
      Δt = [k > 1 ? T[k]-T[k-1] : T[1] for k in 1:nT]

      # initial condition
      μ0 = [mp.s2b[i](x0...)/mp.D[i,i] for i in 1:mp.nL]

      @polyvar(t)
      @polyvar(y[1:maximum([d[2] for d in mp.Mdim])])
      bt = sort(monomials(t, 0:r))
      bt_red = isodd(r) ? (bt[1:r], bt[1:r]) : (bt[1:r+1], bt[1:r-1])
      bt_sos = isodd(r) ? (bt[1:div(r+1,2)], bt[1:div(r+1,2)]) : (bt[1:div(r,2)+1], bt[1:div(r,2)])

      m = Model(solver)
      λ = Dict([k => @variable(m, [1:mp.nL], Poly(bt)) for k in 1:nT])
      Y, Q, S = Dict(), Dict(), Dict()
      for i in 1:length(mp.M)
            dim = mp.Mdim[i][2]
            bty = (cat([y[k].*bt_sos[1] for k in 1:dim]..., dims=1),
                   cat([y[k].*bt_sos[2] for k in 1:dim]..., dims=1))
            for k in 1:nT
                  Y[i,k] = (@variable(m, [1:dim,1:dim], Poly(bt_red[1])),
                            @variable(m, [1:dim,1:dim], Poly(bt_red[2])))
                  Q[i,k] = (@variable(m, [1], SOSPoly(bty[1]))[1],
                            @variable(m, [1], SOSPoly(bty[2]))[1])
                  S[i,k] = @variable(m, [1:dim,1:dim], Poly(bt))
                  @constraint(m, y[1:dim]'*Y[i,k][1]*y[1:dim] == Q[i,k][1])
                  @constraint(m, y[1:dim]'*Y[i,k][2]*y[1:dim] == Q[i,k][2])
                  @constraint(m, S[i,k] .== Λ(Y[i,k], t, r, Δt[k]))
            end
      end

      K = cat(Matrix{Float64}(I,mp.nL,mp.nL), zeros(mp.nH,mp.nL), dims=1)
      # dynamics
      @constraint(m, dae[k in 1:nT],
                     K*∇.(λ[k],t) .== - mp.A'*λ[k] .+ sum(set_dual_LMI(mp.M[j], S[j,k]) for j in 1:length(mp.M)))
      @constraint(m, continuity[k in 1:nT-1, i in 1:mp.nL],
                     λ[k][i](Δt[k]) == λ[k+1][i](0))

      if !isempty(mp.B)
            @error "Equality constraints are not supported"
      end

      terminal_condition = @constraint(m, [i in 1:mp.nL], λ[nT][i](Δt[end]) == 0)

      @objective(m, Max, sum(μ0[i]*λ[1][i](0) for i in 1:mp.nL))
      slack = @variable(m, [1:2,1:2], PSD)
      sol = Dict()

      for obj in objs
            if obj[1] == 1 # means
                  s = mp.b2s[obj[2]]
                  c = zeros(mp.nL)
                  c[s] = 1.0
                  delete(m, terminal_condition)
                  terminal_condition = @constraint(m, [i in 1:mp.nL], λ[nT][i](Δt[nT]) == c[i])
                  optimize!(m)
                  lb, stat_lb, t_lb = get_results(m; scale=mp.D[s,s])

                  c[s] = -1.0
                  delete(m,terminal_condition)
                  terminal_condition = @constraint(m, [i in 1:mp.nL], λ[nT][i](Δt[nT]) == c[i])
                  optimize!(m)
                  ub, stat_ub, t_ub = get_results(m; scale=mp.D[s,s])
                  sol[obj] = ([lb,-ub], [stat_lb,stat_ub], [t_lb,t_ub])
            else
                  s1 = mp.b2s[obj[2]]
                  s2 = mp.b2s[obj[2]^2]
                  delete(m, terminal_condition)
                  terminal_condition = @constraint(m, [i in 1:mp.nL], λ[nT][i](Δt[nT]) == -(i == s1 ? tr([0.0 1.0; 1.0 0.0]*slack) : (i == s2 ? tr([1.0 0.0; 0.0 0.0]*slack) : 0.0)))
                  slack_con = @constraint(m, tr([-1 0; 0 0]*slack) == -1) # ⟺ slack[1,1] == 1
                  @objective(m, Max, sum(μ0[i]*λ[1][i](0) for i in 1:mp.nL) - tr(slack*[0 0; 0 1])) #⟺ set_objective_coefficient(m, slack[2,2], -1.0)
                  optimize!(m)
                  ub, stat_ub, t_ub = get_results(m; scale=mp.D[s2,s2])
                  sol[obj] = (-ub, stat_ub, t_ub)
                  # clean up
                  @objective(m, Max, sum(μ0[i]*λ[1][i](0) for i in 1:mp.nL))#set_objective_coefficient(m, slack[2,2], 0.0)
                  delete(m, slack_con)
            end
      end
      return sol, m
end

function stationary_bounds(mp::MomentProblem, objs, solver; d=[])
      m = Model(solver)
      μ = []

      m = Model(solver)
      if isempty(d)
            μ = [@variable(m, base_name=var_name(i)) for i in 1:mp.nL+mp.nH]
      else
            μ = [d[i]*@variable(m, base_name=var_name(i)) for i in 1:mp.nL+mp.nH]
      end

      # dynamics
      @constraint(m, μ[1] == 1)
      @constraint(m, dynamics[i in 2:mp.nL], mp.A[i,:]'*μ == 0)

      # support
      # affine constraints
      if !isempty(mp.B)
            @constraint(m, [i in 1:length(mp.B)], mp.B[i]*μ .== 0)
      end
      # LMIs
      for p in 1:length(mp.M)
            if mp.Mdim[p][1] > 1
                  @constraint(m, set_LMI(mp.M[p],μ) in PSDCone())
            else
                  @constraint(m, set_LMI(mp.M[p],μ) .>= 0)
            end
      end

      sol = Dict()
      for obj in objs
            if obj[1] == 1 # means
                  s = mp.b2s[obj[2]]
                  @objective(m, Min, μ[s])
                  optimize!(m)
                  lb, stat_lb, t_lb = get_results(m; scale=mp.D[s,s])

                  @objective(m, Max, μ[s])
                  optimize!(m)
                  ub, stat_ub, t_ub = get_results(m; scale=mp.D[s,s])

                  sol[obj] = ([lb,ub], [stat_lb,stat_ub], [t_lb,t_ub])
            elseif obj[1] == 2 # variance
                  s1 = mp.b2s[obj[2]]
                  s2 = mp.b2s[obj[2]^2]
                  slack = @variable(m)
                  @constraint(m, [μ[s2]-slack μ[s1]
                                  μ[s1] 1] in PSDCone())
                  @objective(m, Max, slack)
                  optimize!(m)

                  sol[obj] = get_results(m; scale=mp.D[s2,s2])
            elseif obj[1] == 3 # tailored lin. combination
                  @objective(m, Min,
                                sum(mp.D[mp.b2s[i],mp.b2s[i]]*obj[2][i]*μ[mp.b2s[i]]
                                    for i in keys(obj[2])))
                  optimize!(m)

                  sol[obj] = get_results(m)
            else
                  @error "objective option not supported"
            end
      end
      d = get_scaling_factors(μ,d)
      return sol, d, m
end

function jump_G(midx::Union{Array{Int64,1}, Tuple}, a::Array{<:T1,1},
                x::Array{<:T2,1}, S::Array{<:Number,2}) where {T1,T2<:AbstractPolynomialLike}
      return sum(a[r]*(polyPow(midx, x .+ S[r,:]) - polyPow(midx, x)) for r in 1:length(a))
end

function jump_G(midx::Union{Array{Int64,1}, Tuple, Int64}, a::Array{<:T1,1},
                x::T2, S::Array{<:Number,1}) where {T1,T2 <: AbstractPolynomialLike}
      return sum(a[r]*((x + S[r])^midx[1] - x^midx[1]) for r in 1:length(a))
end

function ito_G(midx::Union{Array{Int64,1}, Tuple}, μ::Array{T1,1},
               σ::Array{T2}, x::Array{T3,1}) where {T1,T2,T3 <:AbstractPolynomialLike}
      f = polyPow(midx,x)
      return dot(μ,∇(f,x)) + 1/2 * tr(∆(f,x)*σ)
end

function ito_G(midx::Union{Array{Int64,1}, Int64, Tuple},
               μ::T1, σ::T2, x::T3) where {T1,T2,T3 <:AbstractPolynomialLike}
      f = x^midx[1]
      return μ*∇(f,x) + 1/2*∆(f,x)*σ
end

function jumpdiff_G()
      # tbd
end

function JumpProcess(x::Array{T1,1}, a::Array{T2,1}, S::Array{<:Number,2},
                     χ::ReachableSet = ReachableSet()) where {T1,T2<:AbstractPolynomialLike}
      G(midx) = jump_G(midx, a, x, S)
      q = maximum(maxdegree.(a)) - 1
      return StochasticProcess(length(x),0,x,Array{AbstractPolynomialLike,1}(),q,G,χ,"Jump")
end

function JumpProcess(x::T1, a::Array{T2,1}, S::Array{<:Number,1},
                     χ::ReachableSet = ReachableSet()) where {T1,T2<:AbstractPolynomialLike}
      G(midx) = jump_G(midx, a, x, S)
      q = maximum(maxdegree.(a)) - 1
      return StochasticProcess(1,0,[x],Array{AbstractPolynomialLike,1}(),q,G,χ,"Jump")
end

function ItoProcess(x::Array{T1,1}, y::Array{T2,1}, μ::Array{T3,1}, σ::Array{T4},
                    χ::ReachableSet = ReachableSet()) where {T1,T2,T3,T4 <: AbstractPolynomialLike}
      G(midx) = ito_G(midx, μ, σ, x)
      μ_deg = maximum(maxdegree.(μ))
      σ_deg = maximum(maxdegree.(σ))
      q = max(μ_deg-1, σ_deg-2, 0)
      return StochasticProcess(length(x),length(y),x,y,q,G,χ,"Ito")
end

function ItoProcess(x::T1, y::Array{T2,1}, μ::T3, σ::T4,
                    χ::ReachableSet = ReachableSet()) where {T1,T2,T3,T4 <: AbstractPolynomialLike}
      G(midx) = ito_G(midx, μ, σ, x)
      μ_deg = maximum(maxdegree.(μ))
      σ_deg = maximum(maxdegree.(σ))
      q = max(μ_deg-1, σ_deg-2, 0)
      return StochasticProcess(1,length(y),[x],y,q,G,χ,"Ito")
end

function ItoProcess(x::Array{T1,1}, y::T2, μ::Array{T3,1}, σ::Array{T4,2},
                    χ::ReachableSet = ReachableSet()) where {T1,T2,T3,T4 <: AbstractPolynomialLike}
      G(midx) = ito_G(midx, μ, σ, x)
      μ_deg = maximum(maxdegree.(μ))
      σ_deg = maximum(maxdegree.(σ))
      q = max(μ_deg-1, σ_deg-2, 0)
      return StochasticProcess(length(x),1,x,[y],q,G,χ,"Ito")
end

function ItoProcess(x::T1, y::T2, μ::T3, σ::T4,
                    χ::ReachableSet = ReachableSet()) where {T1,T2,T3,T4 <: AbstractPolynomialLike}
      G(midx) = ito_G(midx, μ, σ, x)
      μ_deg = maximum(maxdegree.(μ))
      σ_deg = maximum(maxdegree.(σ))
      q = max(μ_deg-1, σ_deg-2, 0)
      return StochasticProcess(1,1,[x],[y],q,G,χ,"Ito")
end

function ItoProcess(x::T1, μ::T2, σ::T3,
                    χ::ReachableSet = ReachableSet()) where {T1,T2,T3 <: AbstractPolynomialLike}
      G(midx) = ito_G(midx, μ, σ, x)
      μ_deg = maximum(maxdegree.(μ))
      σ_deg = maximum(maxdegree.(σ))
      q = max(μ_deg-1, σ_deg-2, 0)
      return StochasticProcess(1,0,[x],Array{AbstractPolynomialLike,1}(),q,G,χ,"Ito")
end

function ItoProcess(x::Array{T1,1}, μ::Array{T2,1}, σ::Array{T3,2},
                    χ::ReachableSet = ReachableSet()) where {T1,T2,T3 <: AbstractPolynomialLike}
      G(midx) = ito_G(midx, μ, σ, x)
      μ_deg = maximum(maxdegree.(μ))
      σ_deg = maximum(maxdegree.(σ))
      q = max(μ_deg-1, σ_deg-2, 0)
      return StochasticProcess(length(x),0,x,Array{AbstractPolynomialLike,1}(),q,G,χ,"Ito")
end
end
