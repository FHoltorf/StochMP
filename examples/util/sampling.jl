using DifferentialEquations, ProgressMeter, JLD

function KMC(x0, c, reac_S, net_S, Tf_range, n_samples, trange; moments=[])
      kinetics = MassActionJump(c, reac_S, net_S, scale_rates=false)
      prob = DiscreteProblem(x0,(0,Tf_range[end]))
      jump_prob = JumpProblem(prob, Direct(), kinetics)
      nt = length(trange)
      X = [zeros(size(x0)) for i = 1:nt]
      X_sq = [zeros(size(x0)) for i = 1:nt]
      X_m = Dict(m=>[zeros(size(x0)) for i = 1:nt] for m in moments)
      @showprogress for i = 1:n_samples
            sol = DifferentialEquations.solve(jump_prob,SSAStepper())
            for k = 1:nt
                  X[k] += sol(trange[k])/n_samples
                  X_sq[k] += sol(trange[k]).^2/n_samples
                  for m in moments
                        X_m[m][k] += sol(trange[k]).^m/n_samples
                  end
            end
      end
      X_var = [X_sq[i] - X[i].^2 for i in eachindex(X)]
      if isempty(moments)
            return X, X_var
      else
            return X, X_var, X_m
      end
end
