function integral = overlapInt(a,Ra,b,Rb)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    A = (4*a*b/pi^2) ^ 0.75;
    integral = A * (pi/(a+b))^1.5 * exp(-a*b/(a+b) * difference(Ra,Rb));
end

