function [a,b] = beta_par(m,s)
a = m*( (m*(1-m)/s^2) - 1 );
b = (1-m)*( (m*(1-m)/s^2) -1 );
end