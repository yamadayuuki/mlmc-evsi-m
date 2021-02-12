function [m,s] = in_beta_par(a,b)
m = a/(a + b);
s = sqrt(a*b/(((a + b)^2)*(a + b + 1)));
end