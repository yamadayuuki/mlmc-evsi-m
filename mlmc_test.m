%
% function mlmc_test(mlmc_fn,M, N,L, N0,Eps,Lmin,Lmax, fp)
%
% multilevel Monte Carlo test routine
%
% sums = mlmc_fn(l,N)     low-level routine
%
% inputs:  l = level
%          N = number of paths
%
% output: sums(1) = sum(Pf-Pc)
%         sums(2) = sum((Pf-Pc).^2)
%         sums(3) = sum((Pf-Pc).^3)
%         sums(4) = sum((Pf-Pc).^4)
%         sums(5) = sum(Pf)
%         sums(6) = sum(Pf.^2)
%
% M      = refinement cost factor (2^gamma in general MLMC Thm)
%
% N      = number of samples for convergence tests
% L      = number of levels for convergence tests
%
% N0     = initial number of samples for MLMC calcs
% Eps    = desired accuracy array for MLMC calcs
%

function mlmc_test(mlmc_fn,M, N,L, N0,Eps,Lmin,Lmax, fp)

%
% first, convergence tests
%

PRINTF2(fp,'\n');
PRINTF2(fp,'**********************************************************\n');
PRINTF2(fp,'*** Convergence tests, kurtosis, telescoping sum check ***\n');
PRINTF2(fp,'**********************************************************\n');
PRINTF2(fp,'\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)');
PRINTF2(fp,'    kurtosis     check \n-------------------------');
PRINTF2(fp,'--------------------------------------------------\n');

rng('default');    % reset random number generator

del1 = [];
del2 = [];
var1 = [];
var2 = [];
kur1 = [];
chk1 = [];
cost = [];

for l = 0:L
%  disp(sprintf('l = %d',l))
  tic;
  sums = feval(mlmc_fn,l,N);
  cost = [ cost toc ];
  sums = sums/N;
  if (l==0)
    kurt = 0.0;
  else
    kurt = (     sums(4)             ...
             - 4*sums(3)*sums(1)     ...
             + 6*sums(2)*sums(1)^2   ...
             - 3*sums(1)*sums(1)^3 ) ...
           / (sums(2)-sums(1)^2)^2;
  end

  del1 = [del1 sums(1)];
  del2 = [del2 sums(5)];
  var1 = [var1 sums(2)-sums(1)^2 ];
  var2 = [var2 sums(6)-sums(5)^2 ];
  var2 = max(var2, 1e-10);  % fix for cases with var=0
  kur1 = [kur1 kurt ];

  if l==0
    check = 0;
  else
    check = abs(       del1(l+1)  +      del2(l)  -      del2(l+1)) ...
      /    ( 3.0*(sqrt(var1(l+1)) + sqrt(var2(l)) + sqrt(var2(l+1)) )/sqrt(N));
  end
  chk1 = [chk1 check];

  PRINTF2(fp,'%2d   %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e \n', ...
          l,del1(l+1),del2(l+1),var1(l+1),var2(l+1),kur1(l+1),chk1(l+1));
end

%
% print out a warning if kurtosis or consistency check looks bad
%

if ( kur1(end) > 100.0 )
  PRINTF2(fp,'\n WARNING: kurtosis on finest level = %f \n',kur1(end));
  PRINTF2(fp,' indicates MLMC correction dominated by a few rare paths; \n');
  PRINTF2(fp,' for information on the connection to variance of sample variances,\n');
  PRINTF2(fp,' see http://mathworld.wolfram.com/SampleVarianceDistribution.html\n\n');
end

if ( max(chk1) > 1.0 )
  PRINTF2(fp,'\n WARNING: maximum consistency error = %f \n',max(chk1));
  PRINTF2(fp,' indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n\n');
end

%
% use linear regression to estimate alpha, beta and gamma
%

L1 = 2;
L2 = L+1;
pa    = polyfit(L1:L2,log2(abs(del1(L1:L2))),1);  alpha = -pa(1);
pb    = polyfit(L1:L2,log2(abs(var1(L1:L2))),1);  beta  = -pb(1);

% just last two points to minimise effect of MATLAB overhead
gamma = log2(cost(end)/cost(end-1));

PRINTF2(fp,'\n******************************************************\n');
PRINTF2(fp,'*** Linear regression estimates of MLMC parameters ***\n');
PRINTF2(fp,'******************************************************\n');
PRINTF2(fp,'\n alpha = %f  (exponent for MLMC weak convergence)\n',alpha);
PRINTF2(fp,' beta  = %f  (exponent for MLMC variance) \n',beta);
PRINTF2(fp,' gamma = %f  (exponent for MLMC cost) \n',gamma);

%
% second, mlmc complexity tests
%

PRINTF2(fp,'\n');
PRINTF2(fp,'***************************** \n');
PRINTF2(fp,'*** MLMC complexity tests *** \n');
PRINTF2(fp,'***************************** \n\n');
PRINTF2(fp,'   eps       value    mlmc_cost   std_cost  savings     N_l \n');
PRINTF2(fp,'----------------------------------------------------------- \n');
 
rng('default');    % reset random number generator

alpha = max(alpha,0.5);
beta  = max(beta,0.5);
gamma = log2(M);
theta = 0.25;

for i = 1:length(Eps)
  eps = Eps(i);
  [P, Nl] = mlmc(Lmin,Lmax,N0,eps,mlmc_fn,alpha,beta,gamma);
  l = length(Nl)-1;
%  mlmc_cost = (1+1/M)*sum(Nl.*M.^(0:l));
%  std_cost  = sum(var2(end)*M.^(0:l))/((1.0-theta)*eps^2);
  mlmc_cost = sum(16 * Nl.*M.^(0:l));
  l2 = min(l+1,length(var2));
  std_cost  = var2(l2)*16* M^l / ((1.0-theta)*eps^2);
%  std_cost  = var2(l2)/eps^4;

  PRINTF2(fp,'%.3e  %.3e  %.3e  %.3e  %7.2f ', ...
	  eps, P, mlmc_cost, std_cost, std_cost/mlmc_cost);
  PRINTF2(fp,'%12d',Nl);
  PRINTF2(fp,'\n');
end

PRINTF2(fp,'\n');

end

%
% function to print to both a file and stdout 
%

function PRINTF2(fp,varargin)
  fprintf(fp,varargin{:});
  fprintf( 1,varargin{:});
end
