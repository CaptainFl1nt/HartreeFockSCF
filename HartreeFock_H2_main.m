% Hartree-Fock Implementation for H2
% All units are in atomic units.
%% Initial Setup
delta = 10^(-7);
max_steps = 10;
N = 50;
a0 = linspace(0.4,4,N);
Etot = zeros(1,N);

a0_exact = [0.4,0.6,1,1.2,1.3,1.35,1.375,1.4,1.425,1.45,1.5,1.6,2.0,2.4,3.2,4.0];
E_exact = [-0.0749825,-0.7247206,-1.084500,-1.124541,-1.131293,-1.135239,-1.133103,-1.137474,-1.134544,-1.136994,-1.130748,-1.125689,-1.091085,-1.044736,-0.978582,-0.916097];
%% Basis Set Comparison
hold on
basisfiles = ["STO1G.txt", "STO2G.txt", "STO3G.txt", "STO6G.txt"];
for k = 1:4
    [alpha,C_contr,nBasis,name] = getBasisSet(basisfiles(k));
    for i = 1:N
        R = [[0 0 0]; [a0(i) 0 0]];
        [E0, Etot(i), P] = RHF_H2(R,alpha,C_contr,nBasis,delta,max_steps);
    end
    plot(a0,Etot,'DisplayName',name)
end

plot(a0_exact,E_exact,'o','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName',"QM MC");
legend
xlabel("Internuclear Separation (a_0)");
ylabel("Energy (Hartree)");
title("Energy vs. Separation for STO-nG Basis Sets");
xlim([0.4,4]);
hold off

%% Density Plot
N = 2;
a = 1.4;
R = [[0 0 0]; [1.4 0 0]];
[alpha,C_contr,nBasis,name] = getBasisSet("STO6G.txt");
[E0, E_total, P] = RHF_H2(R,alpha,C_contr,nBasis,delta,max_steps);

xlim = [-2 + a/2, 2 + a/2];
ylim = [-2 2];
zlim = [0 0];

nx = 200;
ny = 200;
nz = 1;

[X, Y, Z] = meshgrid(linspace(xlim(1),xlim(2),nx),linspace(ylim(1),ylim(2),ny),linspace(zlim(1),zlim(2),nz));
rho = zeros(nx,ny,nz);

for mu = 1:N
    Rmu = R(mu,:);
    for nu = 1:N
        Rnu = R(nu,:);
        sum1 = 0;
        sum2 = 0;
        for j = 1:nBasis
            a_j = alpha(j);
            c_j = C_contr(j);
            sum1 = sum1 + (2*a_j/pi)^(0.75) * c_j*exp(-a_j*((X-Rmu(1)).^2 + (Y-Rmu(2)).^2 + (Z-Rmu(3)).^2));
        end
        for k = 1:nBasis
            a_k = alpha(k);
            c_k = C_contr(k);
            sum2 = sum2 + (2*a_k/pi)^(0.75) * c_k*exp(-a_k*((X-Rmu(1)).^2 + (Y-Rmu(2)).^2 + (Z-Rmu(3)).^2));
        end
        rho = rho + P(mu,nu).*sum1.*sum2;
    end
end
%% Density Plot
colormap turbo
axis square
imagesc(X(:),Y(:),rho)
cb = colorbar;
cb.Ticks = [0.01,0.45];
cb.TickLabels = {'Low','High'};
cb.Title.String = "Density";
xlabel('$x$ ($a_0$)','Interpreter','latex');
title("Electron Density for r = 1.4 a_0");
ylabel('$y$ ($a_0$)','Interpreter','latex');