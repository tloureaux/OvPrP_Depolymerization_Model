clear all;
clc;


alpha=1e-2;
rho=0.7;
omega=0.1;
sigma=1e-4;
beta=150;
mu=0.01;
lambda=230;
gamma=8;
tau=0.8;
delta=0.8;
k1=250;
C0 = 0.367;

% Sizes
nT=2; n1=6; n2=12;
N=100;
sizes=[1,nT,n1,n2+nT*(0:N-1)];

t_early = linspace(0,30,400);

%% STAGE 1: POLYMERIZATION
M = 100;
Z0 = zeros(N+2,1);

F1 = @(t,Z) F_ode(Z,M,n1,n2,N,alpha,sigma,rho,omega,beta,mu,lambda, ...
                  5,0.1,0.1,10);

[~,Y] = ode23s(F1,[0 300],Z0);
Y = Y';

z0 = Y(:,end)/M;
z0(1)=0; z0(2)=0;

M_agreg = sizes(2:end)*z0;
z0_adim = z0 / M_agreg;

%% SWEEP
Z0 = z0_adim * C0;

F2 = @(t,Z) F_ode(Z,C0,n1,n2,N,alpha,sigma,rho,omega,beta,mu,lambda, ...
                              gamma,tau,delta,k1);

[t,Y_late] = ode15s(F2,[0 500],Z0);
 Y_late = Y_late';

T  = Y_late(1,:);
W  = Y_late(2,:);
S  = Y_late(3,:);
Ei = Y_late(4:end,:);

sizes_red = [2,n1,n2+2*(0:N-1)];
V = C0 - sizes_red * Y_late;

C_late = zeros(N+3,length(t));
C_late(2:end,:) = Y_late / C0;
C_late(1,:)     = 1 - sizes(2:end)*Y_late / C0;

Wbar = mean(W(floor(end/2):end));
Sbar = mean(S(floor(end/2):end));
Tbar = mean(T(floor(end/2):end));
Vbar = mean(V(floor(end/2):end));

Eibar = [];
for i=1:N-1
    Eibar = [Eibar mean(Ei(i,floor(end/2):end))];
end

Sstar = @(W) (gamma*W + lambda*W.^2) ./ (2*(mu + k1*W));

X = @(W) (4*beta^3*Sstar(W).^4)./(27*sigma*(mu + k1*W).^2);
u_plus  = @(W) ((mu + k1*W)./(2*sigma)) .* (1 + sqrt(1 + X(W))) .* Sstar(W);
u_minus = @(W) ((mu + k1*W)./(2*sigma)) .* (1 - sqrt(1 + X(W))) .* Sstar(W);
Tstar = @(W) nthroot(u_plus(W),3) + nthroot(u_minus(W),3);

A = @(W) delta*(alpha + tau*Sstar(W) + rho*Tstar(W)) + tau*(omega*Tstar(W) + 3*lambda*W.^2);
B = @(W) 4*delta*tau*(alpha + rho*Tstar(W)).*(omega*Tstar(W) + 3*lambda.*W.^2);
Vstar_sq_minus = @(W) (A(W)./(2*tau*(alpha + rho*Tstar(W)))) .* (1 - sqrt(1 - B(W)./A(W).^2));
Vstarminus = @(W) sqrt(Vstar_sq_minus(W));

qminus = @(W) (tau*Vstarminus(W).^2)./delta;
Ei_starm = @(i,W) (qminus(W)).^i .* Sstar(W);

fWminus = @(W) Vstarminus(W) + 2*Tstar(W) + 6*W ...
              + 12*Sstar(W) ...
              + ((14*qminus(W) - 12*qminus(W).^2)./((1-qminus(W)).^2)).*Sstar(W) ...
              - C0;

Wmin = 0.0;
Wmax = C0/6;
nScan = 1000000;

Wgrid = linspace(Wmin, Wmax, nScan);

Fgridm = fWminus(Wgrid);
idx = find(Fgridm(1:end-1).*Fgridm(2:end) < 0);
brackets = [Wgrid(idx).'  Wgrid(idx+1).'];
Wstarm = brackets(:,1);
if Wstarm(1) == 0
    opts = optimoptions('fsolve','Display','iter','TolFun',1e-12,'TolX',1e-12);
    f_trans = @(x) fWminus(x.^2);   % W = x^2 ensures W >= 0
    [xroot, Fval, exitflag] = fsolve(f_trans, 0, opts);
    Wrootm = xroot^2;
    Wstarm(1) = Wrootm;
end

disp(Wstarm);

Sstarm = Sstar(Wstarm(end));
Tstarm = Tstar(Wstarm(end));
Vstarm = Vstarminus(Wstarm(end));
Eistarm = @(i) (qminus(Wstarm(end))).^i .* Sstar(Wstarm(end));

M2 = C_late(1,:) + sizes(2:end).^2 * Y_late / C0;

Mwstarm = [];
    Ceq = [ Tstarm ; Wstarm(end) ; Sstarm];
    for j=1:N-1
        x = Eistarm(j);
        Ceq = [Ceq ; x];
    end
    Mwstarm = [Mwstarm 1- (sizes(2:end)*Ceq)/C0 ...
            + (sizes(2:end).^2 * Ceq) / C0;];

%% TRAJECTORY NEAR EQUILIBRIUM

eps = 1e-4;   % small perturbation
Z0_eq = Ceq + eps*randn(size(Ceq));   % perturbed initial condition

[t_eq, Y_eq] = ode15s(F2, [0 200], Z0_eq);
Y_eq = Y_eq';

W_eq  = Y_eq(2,:);
S_eq  = Y_eq(3,:);
Ei_eq = Y_eq(4:end,:);

Mw_eq = 1 - (sizes(2:end)*Y_eq)/C0 + (sizes(2:end).^2 * Y_eq) / C0;
Mw_eq = Mw_eq(1,:);

t_start = 0;
t_max   = 1000;
idx = find(t_eq >= t_start & t_eq <= t_max);
W_eq_r  = W_eq(idx);
S_eq_r  = S_eq(idx);
Mw_eq_r = Mw_eq(idx);
t_eq_r  = t_eq(idx);

t_norm = (t_eq_r - min(t_eq_r)) / (max(t_eq_r) - min(t_eq_r));
cmap = jet(length(t_norm));

%% ASYMPTOTICS

%Assuming kW << mu
Wapp1 = (2*sigma/gamma)*((C0/2)-(1/8)*(omega/alpha)*(sqrt(1+(8*alpha*C0/omega))-1))^3;
Sapp1 = gamma*Wapp1/(2*mu);
Tapp1 = (gamma*Wapp1/(2*sigma))^(1/3);
Vapp1 = sqrt((omega/alpha)*Tapp1);
Eiapp1 = zeros(1,N);
for i=1:N
    Eiapp1(i) = Sapp1*(tau*Vapp1^2/delta)^i;
end

Capp1 = [ Tapp1 ; Wapp1 ; Sapp1; Eiapp1(1:N-1)'];
Mwapp1 = 1- (sizes(2:end)*Capp1)/C0 + (sizes(2:end).^2 * Capp1) / C0;

%Assuming kW >> mu
K = 18/((24*sqrt(3))^(2/3))*(beta*lambda/(rho*k1^2));
Cconst = sqrt(omega/rho)+6*(gamma/k1)+7*(gamma/k1)*(tau*omega/(delta*rho));
Cw = (K*gamma/(2*sqrt(omega/rho)))+((24*sqrt(3))^(2/3)/3)*(k1^2/(beta*gamma))+6+6*(lambda/k1)+7*(lambda/k1)*(tau*omega/(delta*rho));
Wapp3 = (C0-Cconst)/Cw;
Sapp3 = (lambda*Wapp3+gamma)/(2*k1);
Tapp3 = ((24*sqrt(3))^(2/3)/6)*(k1^2/beta)*(Wapp3/(lambda*Wapp3+gamma));
Vapp3 = sqrt((omega/rho)+(18/(24*sqrt(3))^(2/3))*(beta*lambda/(rho*k1^2))*(lambda*Wapp3^2+gamma*Wapp3));
Eiapp3 = zeros(1,N);
for i=1:N
    Eiapp3(i) = Sapp3*(tau*Vapp3^2/delta)^i;
end

Capp3 = [ Tapp3 ; Wapp3 ; Sapp3; Eiapp3(1:N-1)'];
Mwapp3 = 1- (sizes(2:end)*Capp3)/C0 + (sizes(2:end).^2 * Capp3) / C0;


disp('W');
disp([Wstarm' Wapp1 Wapp3]);
disp('S');
disp([Sstarm' Sapp1 Sapp3]);
disp('T');
disp([Tstarm' Tapp1 Tapp3]);
disp('V');
disp([Vstarm' Vapp1 Vapp3]);
disp('E1');
disp([Eistarm(1)' Eiapp1(1) Eiapp3(1)]);
disp('E2');
disp([Eistarm(2)' Eiapp1(2) Eiapp3(2)]);
disp('Mw');
disp([Mwstarm Mwapp1 Mwapp3]);

figure;
plot(t, W, 'black', 'LineWidth', 1.5); hold on;
colors = lines(length(Wstarm));
k = 1;
for i = 1:length(Wstarm)
    yline(Wstarm(i), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    k = k + 1;
end
xlabel('Time t');
ylabel('W');
title('W(t)');

figure;
plot(t, S, 'black', 'LineWidth', 1.5); hold on;
colors = lines(length(Wstarm));
k = 1;
for i = 1:length(Sstarm)
    yline(Sstarm(i), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    k = k + 1;
end
xlabel('Time t');
ylabel('S');
title('S(t)');

figure;
plot(t, T, 'black', 'LineWidth', 1.5); hold on;
colors = lines(length(Wstarm));
k = 1;
for i = 1:length(Tstarm)
    yline(Tstarm(i), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    k = k + 1;
end
xlabel('Time t');
ylabel('T');
title('T(t)');

figure;
plot(t, V, 'black', 'LineWidth', 1.5); hold on;
colors = lines(length(Wstarm));
k = 1;
for i = 1:length(Vstarm)
    yline(Vstarm(i), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    k = k + 1;
end
xlabel('Time t');
ylabel('V');
title('V(t)');

figure;
plot(t,M2,'black'); hold on;
colors = lines(length(Wstarm));
k = 1;
for i = 1:length(Vstarm)
    yline(Mwstarm(i), '-', 'Color', colors(k,:), 'LineWidth', 1.5);
    k = k + 1;
end
xlabel('Time t');
ylabel('Mw');
title(sprintf('Mw(t)   |   C0 = %.2f,   \\delta = %.2f,   \\tau = %.2f', ...
              C0, delta, tau));

tail = floor(0.8*length(t)) : length(t);

W_tail = W(tail);
S_tail = S(tail);
T_tail = T(tail);
V_tail = V(tail);
Mw_tail = M2(tail);

amp = max(M2(tail)) - min(M2(tail));
fprintf('C0 = %.2f, tail amplitude in W = %.6e\n', C0, amp);

Ceq = [ Tstarm ; Wstarm(end) ; Sstarm ];
for j = 1:N-1
    x = Eistarm(j);
    Ceq = [Ceq ; x(1)];
end

diff_full = Y_late - Ceq;                 % broadcast over columns
dist_full = vecnorm(diff_full,2,1);       % Euclidean norm at each time

[min_dist_full, idx_min_full] = min(dist_full);

fprintf('Min full-state distance to equilibrium: %.3e\n', min_dist_full);
fprintf('At time t = %.3f\n', t(idx_min_full));


%fprintf('Max distance from equilibrium: W = %.3e, S = %.3e\n', max(dist_W), max(dist_S));

%figure;
%plot3(W_tail, S_tail, T_tail, 'k', 'LineWidth', 1.5);
%grid on;
%xlabel('W');
%ylabel('S');
%zlabel('T');
%title('3D limit cycle in (W, S, T)');

figure; hold on;
plot(W_tail, Mw_tail, 'k', 'LineWidth', 2.0);
plot(Wstarm, Mwstarm, 'r.', 'MarkerSize', 18);
plot(W_eq_r, Mw_eq_r, 'Color', [0 0.4 1 0.25], 'LineWidth', 0.75);
grid on;
xlabel('W'); ylabel('Mw');
title('(W, Mw) projection with restrained trajectory');

figure; hold on;
plot(W_tail, S_tail, 'k', 'LineWidth', 2.0);
plot(Wstarm, Sstarm, 'r.', 'MarkerSize', 18);
plot(W_eq_r, S_eq_r, 'Color', [0 0.4 1 0.25], 'LineWidth', 0.75);
grid on;
xlabel('W'); ylabel('S');
title('(W, S) projection with restrained trajectory');

figure; hold on;
plot3(W_tail, S_tail, Mw_tail, 'k', 'LineWidth', 1.5);
plot3(Wstarm, Sstarm, Mwstarm, 'r.', 'MarkerSize', 14);
plot3(W_eq_r, S_eq_r, Mw_eq_r, ...
      'Color', [0 0.4 1 0.25], ...
      'LineWidth', 0.75);
grid on;
xlabel('W');
ylabel('S');
zlabel('Mw');
title('3D limit cycle in (W, S, Mw)');
view(45,25);

figure;
plot(Wgrid,fWminus(Wgrid)); hold on;
yline(0,'-r'); hold on;
plot(Wstarm,fWminus(Wstarm),'r.','MarkerSize',15);
xlabel('W')
ylabel('Mass conservation function');
title('EQUILIBRIA W^{*}');
xlim([-1*1e-5,7*1e-5]);
ylim([-0.04,0.38]);


%DfWminus = @(W) (fWminus(W + 1e-6) - fWminus(W - 1e-6)) / (2e-6);
%Df_vals = DfWminus(Wgrid);
%figure;
%plot(Wgrid, Df_vals, 'LineWidth', 1.5); hold on;
%plot(Wstarm, DfWminus(Wstarm), 'r.', 'MarkerSize', 10);
%xlabel('W');
%ylabel('d/dW fWminus(W)');
%title('Derivative of fWminus(W)');
%grid on;


%figure;
%plot(Wgrid,Tstar(Wgrid));
%title('T');

%F = @(W) (2*lambda)./(2*lambda.*W + gamma) ...
%       + (3*k)./(k.*W + mu) ...
%       - (2*lambda.*W + gamma)./(lambda.*W.^2 + gamma.*W) ...
%       - (2*lambda.*k.*W + 2*lambda*mu)./(lambda.*k.*W.^2 + 2*lambda*mu.*W + gamma*mu);

%figure;
%plot(Wgrid,F(Wgrid));
%title('F');

%figure;
%i_plot = 1:10;  % indices of Ei to plot
%[t_grid, i_grid] = meshgrid(t, i_plot);
%Ei_plot = Ei(i_plot, :);                 % select first 10 Ei

figure; hold on;
i_plot = 1:20;           % indices of Ei to plot
colors = lines(length(i_plot));
for j = 1:length(i_plot)
    plot3(t, i_plot(j)*ones(size(t)), Ei(i_plot(j),:), ...
    'Color', colors(j,:), 'MarkerSize', 3, 'LineWidth', 1);
end
xlabel('Time t');
ylabel('Index i');
zlabel('Ei(t)');
title('Discrete 3D plot of Ei(t)');
grid on;
view(45,30);  % adjust viewing angle

%Vapp0 = @(W) sqrt((omega/alpha)*Tstar(W));
%Wgrid0 = linspace(Wmin, 1, nScan);
%figure;
%plot(Wgrid0,Vstarminus(Wgrid0)); hold on;
%title('V');

%figure;
%plot(Tstar(Wgrid),Vstarminus(Wgrid));
%title('T vs V');