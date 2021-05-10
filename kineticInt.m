function integral = kineticInt(a,Ra,b,Rb)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    A = (4*a*b/pi^2) ^ 0.75;
    fact = a*b/(a+b);
    integral = A * fact * (3 - 2*fact*difference(Ra,Rb))*(pi/(a+b))^1.5*exp(-fact*difference(Ra,Rb));
end

