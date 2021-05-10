function integral = nucAttInt(a,Ra,b,Rb,Rc)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    A = (4*a*b/pi^2) ^ 0.75;
    Rp = (a*Ra + b*Rb)/(a + b);
    integral = -2*pi/(a+b)*exp(-difference(Ra,Rb)*a*b/(a+b))*F0((a+b)*difference(Rp,Rc))*A;
end

