function [alpha beta] = betapar(par_est)
    alpha = (1 - par_est)/0.01 - par_est;
    beta = alpha*(1/par_est - 1);
end