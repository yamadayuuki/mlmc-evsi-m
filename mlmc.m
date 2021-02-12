% function [P, Nl] = mlmc(Lmin,Lmax,N0,eps,mlmc_l, alpha,beta,gamma)
%
% multi-level Monte Carlo estimation
%
% P     = value
% Nl    = number of samples at each level
%
% Lmin  = minimum level of refinement       >= 2
% Lmax  = maximum level of refinement       >= Lmin
% N0    = initial number of samples         > 0
% eps   = desired accuracy (rms error)      > 0 
%
% alpha -> weak error is  O(2^{-alpha*l})
% beta  -> variance is    O(2^{-beta*l})
% gamma -> sample cost is O(2^{gamma*l})    > 0
%
% if alpha, beta are not positive then they will be estimated
%
% mlmc_l = function for level l estimator 
%
% sums = mlmc_fn(l,N)     low-level routine
%
% inputs:  l = level
%          N = number of paths
%
% output: sums(1) = sum(Y)
%         sums(2) = sum(Y.^2)
%         where Y are iid samples with expected value:
%         E[P_0]           on level 0
%         E[P_l - P_{l-1}] on level l>0

function [P, Nl] = mlmc(Lmin,Lmax,N0,eps,mlmc_l, alpha_0,beta_0,gamma)

%
% check input parameters
%
  if (Lmin<2)
    error('error: needs Lmin >= 2');
  end

  if (Lmax<Lmin)
    error('error: needs Lmax >= Lmin');
  end

  if (N0<=0 || eps<=0 || gamma <= 0)
    error('error: needs N>0, eps>0, gamma>0 \n');
  end

%
% initialisation
%
  alpha = max(0, alpha_0);
  beta  = max(0, beta_0);

  theta = 0.25;

  L = Lmin;

  Nl(1:L+1)       = 0;
  suml(1:2,1:L+1) = 0;
  dNl(1:L+1)      = N0;

  while sum(dNl) > 0

%
% update sample sums
%   

      if dNl(0+1) > 0
        sums        = feval(mlmc_l,0,dNl(0+1));
        Nl(0+1)     = Nl(0+1) + dNl(0+1);
        suml(1,0+1) = suml(1,0+1) + sums(5);
        suml(2,0+1) = suml(2,0+1) + sums(6);
      end
    
    for l=1:L
      if dNl(l+1) > 0
        sums        = feval(mlmc_l,l,dNl(l+1));
        Nl(l+1)     = Nl(l+1) + dNl(l+1);
        suml(1,l+1) = suml(1,l+1) + sums(1);
        suml(2,l+1) = suml(2,l+1) + sums(2);
      end
    end


%
% compute absolute average and variance
%
    ml = abs(   suml(1,:)./Nl);
    Vl = max(0, suml(2,:)./Nl - ml.^2);

%
% fix to cope with possible zero values for ml and Vl
% (can happen in some applications when there are few samples)
%
    for l = 3:L+1
      ml(l) = max(ml(l), 0.5*ml(l-1)/2^alpha);
      Vl(l) = max(Vl(l), 0.5*Vl(l-1)/2^beta);
    end

%
% use linear regression to estimate alpha, beta if not given
%
    if alpha_0 <= 0
      A     = repmat((1:L)',1,2).^repmat(1:-1:0,L,1);
      x     = A \ log2(ml(2:end))';
      alpha = max(0.5,-x(1));
    end

    if beta_0 <= 0
      A    = repmat((1:L)',1,2).^repmat(1:-1:0,L,1);
      x    = A \ log2(Vl(2:end))';
      beta = max(0.5,-x(1));
    end
%
% set optimal number of additional samples
%
    Cl  = 16 * 2.^(gamma*(0:L));
    Ns  = ceil( sqrt(Vl./Cl) * sum(sqrt(Vl.*Cl)) ...
                             / ((1-theta)*eps^2) );
    dNl = max(0, Ns-Nl);
%
% if (almost) converged, estimate remaining error and decide 
% whether a new level is required
%
    if sum( dNl > 0.01*Nl ) == 0
      rem = ml(L+1) / (2^alpha - 1);

      if rem > sqrt(theta)*eps
        if (L==Lmax)
          fprintf(1,'*** failed to achieve weak convergence *** \n');
        else
          L       = L+1;
          Vl(L+1) = Vl(L) / 2^beta;
          Nl(L+1) = 0;
          suml(1:4,L+1) = 0;

          Cl  = 2.^(gamma*(0:L));
          Ns  = ceil( sqrt(Vl./Cl) * sum(sqrt(Vl.*Cl)) ...
                                   / ((1-theta)*eps^2) );
          dNl = max(0, Ns-Nl);
        end
      end
    end
  end

%
% finally, evaluate multilevel estimator
%
  P = sum(suml(1,:)./Nl);
end
