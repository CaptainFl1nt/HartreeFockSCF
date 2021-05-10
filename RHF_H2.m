function [E0,Etot,P] = RHF_H2(R,alpha,C_contr,nBasis,delta,max_steps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = 2;              % number of orbitals

% compute overlap, kinetic, nuclear attraction, and two-electron integrals
S = zeros(N,N);
T = zeros(N,N);
V = zeros(N,N);
tETensor = zeros(N,N,N,N); % two electron integral tensro

% evaluate overlap, kinetic, and nuclear attraction integrals
for mu = 1:N
    Rmu = R(mu,:);
    for nu = 1:N
        Rnu = R(nu,:);
        for i = 1:nBasis
            c_mu_i = C_contr(i);
            a_mu = alpha(i);
            for j = 1:nBasis
                c_nu_j = C_contr(j);
                a_nu = alpha(j);
                fact = c_mu_i*c_nu_j;
                
                S(mu,nu) = S(mu,nu) + fact*overlapInt(a_mu,Rmu,a_nu,Rnu);
                T(mu,nu) = T(mu,nu) + fact*kineticInt(a_mu,Rmu,a_nu,Rnu);
                
                for c = 1:2
                    V(mu,nu) = V(mu,nu) + fact*nucAttInt(a_mu,Rmu,a_nu,Rnu,R(c,:));
                end
            end
        end
        
        for lambda = 1:N
            Rlam = R(lambda,:);
            for sigma = 1:N
                Rsig = R(sigma,:);
                for p = 1:nBasis
                    c_mu_p = C_contr(p);
                    a_mu = alpha(p);
                    for q = 1:nBasis
                        c_nu_q = C_contr(q);
                        a_nu = alpha(q);
                        for k = 1:nBasis
                            c_lam_k = C_contr(k);
                            a_lam = alpha(k);
                            for m = 1:nBasis
                                c_sig_m = C_contr(m);
                                a_sig = alpha(m);
                                fact = c_mu_p*c_nu_q*c_lam_k*c_sig_m;
                                integral = twoElecInt(a_mu,Rmu,a_nu,Rnu,a_lam,Rlam,a_sig,Rsig);
                                tETensor(mu,nu,lambda,sigma) = tETensor(mu,nu,lambda,sigma) + fact*integral;
                            end
                        end
                    end
                end
            end
        end      
    end
end
    
H_core = T + V;

% SCF Procedure

% diagonalize matrix S
[U,Seig] = eig(S);
s = transpose(U)*S*U;
X = U*inv(sqrt(s))*transpose(U);

% guess of density matrix P
P = zeros(N,N);
P_new = zeros(N,N);

n_steps = 0;
while true
    P = P_new;
    % calculate G matrix
    G = zeros(N,N);
    for mu = 1:N
        for nu = 1:N
            for lambda = 1:N
                for sigma = 1:N
                    mnsl = tETensor(mu,nu,sigma,lambda);
                    mlsn = tETensor(mu,lambda,sigma,nu);
                    G(mu,nu) = G(mu,nu) + P(lambda,sigma)*(mnsl - 0.5*mlsn);
                end
            end
        end
    end

    % Fock matrix
    F = H_core + G;

    % transformed Fock matrix
    F_prime = transpose(X) * F * X;
    
    % calculate C' and epsilon
    [C_prime, eps] = eig(F_prime);
    [eps,index] = sort(diag(eps));
    C_prime = C_prime(:,index);

    % calculate C
    C = X*C_prime;

    % calculate new P
    P_new = zeros(N,N);
    for mu = 1:N
        for nu = 1:N
            for a = 1:(N/2)
                P_new(mu,nu) = P_new(mu,nu) + 2*C(mu,a)*C(nu,a);
            end
        end
    end
    
    % check for convergence
    
    total = 0;
    for mu = 1:N
        for nu = 1:N
            total = total + (P_new(mu,nu) - P(mu,nu))^2;
        end
    end
    total = (total / N^2)^(1/2);
    if total < delta
        break
    end
    
    n_steps = n_steps + 1;
    if (n_steps >= max_steps)
        disp("Maximum Steps Exceeded.");
        break;
    end
end

E0 = 0;
for mu = 1:N
    for nu = 1:N
        E0 = E0 + P(nu,mu)*(H_core(mu,nu)+F(mu,nu));
    end
end

E0 = E0 * 0.5; 

Enuc = 0;
for a = 1:N
    for b = (a+1):N
        Enuc = Enuc + 1/sqrt(difference(R(a,:),R(b,:)));
    end
end

Etot = E0 + Enuc;
end

