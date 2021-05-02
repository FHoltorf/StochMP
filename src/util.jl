var_name(ρ,l,t,i) = string("z^{",l,"}_{",i,"}(",ρ,";t_",t,")")
var_name(t,i) = string("μ_{",i,"}(t_",t,")")
var_name(i) = string("μ_{",i,"}")
#
var_name(ρ,t,i) = string("z_{",i,"}(",ρ,";t_",t,")")
var_name(ρ,k,l,t,i) = string("γ^{",l,",",k,"}_{",i,"}(",ρ,";t_",t,")")
C(t1,t2,k,l) = (t2-t1)^(1+k+l)*factorial(l)*factorial(k)/factorial(1+k+l) # coefficient

function get_results(m; scale::Number=1.0)
      if dual_status(m) == MOI.FEASIBLE_POINT
            status = 1
      else
            status = 0
      end
      return scale*dual_objective_value(m), status, MOI.get(m, MOI.SolveTime())
end

function get_results(m, vars; scale::Number=1.0)
      if dual_status(m) == MOI.FEASIBLE_POINT
            status = 1
      else
            status = 0
      end
      return scale*dual_objective_value(m), value.(vars), status, MOI.get(m, MOI.SolveTime())
end

function get_scaling_factors(μ,z,d)
      try
            if isempty(d)
                  return merge(Dict([key => max.(0.5,abs.(value.(μ[key]))) for key in keys(μ)]),
                         Dict([key => max.(0.5,abs.(value.(z[key]))) for key in keys(z)]))
            else
                  return merge(Dict([key => max.(0.5,abs.(value.(μ[key]./d[key]))) .* d[key] for key in keys(μ)]),
                         Dict([key => max.(0.5,abs.(value.(z[key]./d[key]))) .* d[key] for key in keys(z)]))
            end
      catch
            @error "scaling factor extraction failed"
            return []
      end
end

function get_scaling_factors(μ,z,γ,d)
      if isempty(d)
            return merge(Dict([key => max.(0.5,abs.(value.(μ[key]))) for key in keys(μ)]),
                   Dict([key => max.(0.5,abs.(value.(z[key]))) for key in keys(z)]),
                   Dict([key => max.(0.5,abs.(value.(γ[key]))) for key in keys(γ)]))
      else
            return merge(Dict([key => max.(0.5,abs.(value.(μ[key]./d[key]))) .* d[key] for key in keys(μ)]),
                         Dict([key => max.(0.5,abs.(value.(z[key]./d[key]))) .* d[key] for key in keys(z)]),
                         Dict([key => max.(0.5,abs.(value.(γ[key]./d[key]))) .* d[key] for key in keys(γ)]))
      end
end

function get_scaling_factors(μ,d)
      if isempty(d)
            return max.(0.5, abs.(value.(μ)))
      else
            return max.(0.5, abs.(value.(μ./d))) .* d
      end
end

function Ω(m::Int64,n::Int64,t1::Float64,t2::Float64,z1,z2)
      expr = AffExpr.(zeros(size(z1[1])))
      for k in 0:m
            expr .+= (-1)^k * binomial(n+k,n) * (t2-t1)^(m-k)/factorial(m-k)*z2[n+1+k]
      end
      for k in 0:n
            expr .+= (-1)^(m+1) * binomial(m+k,m) * (t2-t1)^(n-k)/factorial(n-k)*z1[m+1+k]
      end
      #scale = maximum([maximum(abs.(e.terms.vals)) for e in expr])
      return expr ./ (binomial(m+n, m) * (t2-t1)^(max(n,m)))
end


function set_LMI(M::Array{Array{N,1},2},v::Array{T,1}; dim=size(M,1)) where
                                    {T <: Union{VariableRef, GenericAffExpr}, N <: Number}
      expr = AffExpr.(zeros(dim,dim))
      for i in 1:dim
            for j in i:dim
                  nz = findall(x->x!=0, M[i,j])
                  for k in nz
                        add_to_expression!(expr[i,j],M[i,j][k]*v[k])
                        if i!=j
                              add_to_expression!(expr[j,i],M[i,j][k]*v[k])
                        end
                  end
            end
      end
      return expr
end

function set_LMI(M::Array{Array{N,1},2},v::Array{Polynomial{true,T},1}; dim=size(M,1)) where
                                          {T <: Union{VariableRef,GenericAffExpr}, N <: Number}
      expr = zeros(Polynomial{true,GenericAffExpr{Float64,VariableRef}},dim,dim)
      for i in 1:dim
            for j in i:dim
                  nz = findall(x->x!=0, M[i,j])
                  for k in nz
                        expr[i,j] += 1.0*M[i,j][k]*v[k]
                        if i!=j
                              expr[j,i] += 1.0*M[i,j][k]*v[k]
                        end
                  end
            end
      end
      return expr
end

function set_LMI(M::Array{Array{N1,1},2},v::Array{N2,1}; dim=size(M,1)) where
                                                             {N1,N2 <: Number}
      expr = zeros(typeof(M[1,1][1]*v[1]),dim,dim)
      for i in 1:dim
            for j in i:dim
                  nz = findall(x->x!=0, M[i,j])
                  for k in nz
                        expr[i,j] += M[i,j][k]*v[k]
                        if i!=j
                              expr[j,i] += M[i,j][k]*v[k]
                        end
                  end
            end
      end
      return expr
end

function set_dual_LMI(M::Array{Array{N,1},2},S::Array{Polynomial{true,T},2}; dim=size(M,1)) where
                                             {T <: Union{VariableRef,GenericAffExpr}, N <: Number}
      expr = zeros(Polynomial{true,GenericAffExpr{Float64,VariableRef}},length(M[1,1]))
      for i in 1:length(M[1,1])
            Mi = [M[j,k][i] for j in 1:dim, k in 1:dim]
            expr[i] = 1.0*tr(Mi*S)
      end
      return expr
end

function set_dual_LMI(M::Array{Array{N,1},2}, S::Array{T,2}; dim=size(M,1)) where
                                              {T <: Union{VariableRef,GenericAffExpr}, N <: Number}
      expr = zeros(GenericAffExpr{Float64,VariableRef},length(M[1,1]))
      for i in 1:length(M[1,1])
            Mi = [M[j,k][i] for j in 1:dim, k in 1:dim]
            expr[i] = 1.0*tr(Mi*S)
      end
      return expr
end

function set_dual_LMI(M::Array{Array{N,1},2}, S::Symmetric{T1,Array{T2,2}}; dim=size(M,1)) where
                                              {T1,T2 <: Union{VariableRef,GenericAffExpr}, N <: Number}
      expr = zeros(GenericAffExpr{Float64,VariableRef},length(M[1,1]))
      for i in 1:length(M[1,1])
            Mi = [M[j,k][i] for j in 1:dim, k in 1:dim]
            expr[i] = 1.0*tr(Mi*S)
      end
      return expr
end


function assemble_alg_con(h::T1, order::Int64, x::Array{T2,1}, b::Array{T3,1}; D=I) where
                          {T1,T2,T3<:AbstractPolynomialLike}
    dh = order - maxdegree(h)
    if dh < 1
        @error "Truncation order incompatible with LMI generating function.
                Increase truncation order or remove LMI generating function!"
    end
    bh = sort(monomials(x, 0:dh))
    Mh = h*bh

    return row_scaling(cat(map(x -> coefficients(x,b), Mh)...,dims=2)*D)
end

function assemble_lmi_con(g::T1, order::Int64, x::Array{T2,1}, y::Array{T3,1}, b::Array{T4,1}; D=I) where
                         {T1,T2,T3,T4<:AbstractPolynomialLike}
    dg = floor(Int64, (order-maxdegree(g))/2)
    if dg < 1
        @error "Truncation order incompatible with LMI generating function.
                Increase truncation order or remove LMI generating function!"
    end
    if isempty(y)
          bg = sort(monomials(x, 0:dg))
    else
          z = cat(x,y,dims=1)
          bg = cat(sort(monomials(x, 0:dg)),
                   sort(monomials(z, 0:dg, p -> !all(degree.(p,y) .== 0))), dims = 1)
    end
    g_scaled = coefficients(g,b)'*D*b
    Mg = g_scaled*bg*bg'
    return lmi_row_scaling(map(x -> coefficients(x,b), Mg))
end

function ∇(p::Union{AbstractPolynomialLike, Array{<:AbstractPolynomialLike}},
           x::Union{PolyVar, Array{<:PolyVar}})
    return differentiate(p,x)
end

function ∆(p::AbstractPolynomialLike, x::Union{PolyVar, Array{<:PolyVar}}) where {C}
    return ∇(∇(p,x),x)
end

function Λ(Y, t, r, Δt)
    return isodd(r) ? t*Y[1] + (Δt-t)*Y[2] : Y[1] + t*(Δt-t)*Y[2]
end

function polyPow(midx::Union{Array{Int64,1}, Tuple}, p::Array{<: AbstractPolynomialLike,1})
    return prod(p[i]^midx[i] for i in 1:length(midx))
end

function jump_χ(S,x0)
  N = nullspace(S)
  if isempty(N)
        @error "Stoichometry matrix has trivial nullspace.
                There are no reaction invariants."
  end
  R = round.(qr(N').R, digits=5)
  ds = [1]; j = 2
  for i in 2:size(R,2)
        if R[i,j] == 0
              continue
        else
              push!(ds, i)
              j += 1
              if j > size(R,1)
                    break
              end
        end
  end

  # independent species
  is = setdiff(1:size(S,2),ds)
  Nds, Nis = N[ds,:]', N[is,:]'

  @polyvar(x[1:length(is)])

  xds = x0[ds] .+ round.(inv(Nds)*Nis, digits=5)*(x0[is] .- x)
  x_full = Array{Polynomial,1}(undef,size(S,2))
  x_full[is] = x
  x_full[ds] = xds
  return x, xds, is
end

function jump_χ(S,x0,ds)
  N = nullspace(S)
  R = round.(qr(N').R, digits=5)
  is = setdiff(1:size(S,2),ds)
  Nds, Nis = N[ds,:]', N[is,:]'
  @polyvar(x[1:length(is)])

  xds = x0[ds] .+ round.(inv(Nds)*Nis, digits=5)*(x0[is] .- x)
  x_full = Array{Polynomial{true,Float64},1}(undef,size(S,2))
  x_full[is] = x
  x_full[ds] = xds
  return x, x_full, is, ds
end

function row_scaling(A::Array{T,2}) where {T<:Number}
      for i in 1:size(A,1)
            s = maximum(abs.(A[i,:])) < 1e-6 ? 1.0 : maximum(abs.(A[i,:]))
            A[i,:] = A[i,:] ./ s
      end
      return A
end

function lmi_row_scaling(M)
      s = maximum(maximum.(map(x->abs.(x),M)))
      return map(x -> x/s, M)
end
