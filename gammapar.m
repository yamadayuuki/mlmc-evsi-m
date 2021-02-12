function [alpha beta scale] = gammapar(par_est)
    alpha = 100;
    beta = 100/par_est;
    scale = par_est/100;
end