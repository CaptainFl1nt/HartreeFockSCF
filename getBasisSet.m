function [alpha,C_contr,nBasis,name] = getBasisSet(basisfile)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bfile = fopen(basisfile);

name = fgetl(bfile);  % get basis set name
temp = sscanf(fgetl(bfile),'%f');
K = temp(1);      % degree of shell contraction
scale = temp(2);  % exponent scaling

nBasis = sscanf(fgetl(bfile),'%f');
alpha = zeros(1,nBasis);   % Gaussian exponents
C_contr = zeros(1,nBasis); % contraction coefficients

% extract contraction coefficients and exponents
line = fgetl(bfile);
for i = 1:nBasis
    temp = regexp(line,'(\d*.\d*D[+-]\d*)','match');
    ss = split(temp(1),"D");
    alpha(i) = str2double(ss(1)) * 10 ^ str2double(ss(2));
    ss = split(temp(2),"D");
    C_contr(i) = str2double(ss(1)) * 19 ^ str2double(ss(2));
    line = fgetl(bfile);
end

end

