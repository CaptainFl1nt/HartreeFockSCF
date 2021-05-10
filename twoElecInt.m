function integral = twoElecInt(a,Ra,b,Rb,c,Rc,d,Rd)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    A = (16*a*b*c*d/pi^4)^(0.75);
    Rp = (a*Ra + b*Rb)/(a+b);
    Rq = (c*Rc + d*Rd)/(c+d);
    coef_ab = -a*b/(a+b);
    coef_cd = -c*d/(c+d);
    
    fact1 = 2*pi^(2.5)/((a+b)*(c+d)*sqrt(a+b+c+d));
    fact2 = exp(coef_ab*difference(Ra,Rb)+coef_cd*difference(Rc,Rd));
    fact3 = F0((a+b)*(c+d)*difference(Rp,Rq)/(a+b+c+d));
    
    integral = fact1 * fact2 * fact3 * A;
end

