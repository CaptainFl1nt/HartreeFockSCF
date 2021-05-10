function [f0] = F0(t)
%F0 See function defined in equation A.32 (pp. 415) of Szabo and Ostlund.
if abs(t) < 10^(-7)
    f0 = 1;
else
    f0 = 0.5*sqrt(pi/t)*erf(sqrt(t));
end
end

