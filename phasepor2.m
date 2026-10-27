clear all; clc;

%% ---------------------------------------------------------
%  PARAMETERS (fixed)
%% ---------------------------------------------------------
alpha  = 1e-2;
rho    = 0.7;
omega  = 0.1;
sigma  = 1e-4;
beta   = 150;
mu     = 0.01;
lambda = 230;
gamma  = 8;
tau    = 0.8;
delta  = 0.8;
k1     = 250;

% Sweep over these C0 values
%C0_list = [linspace(0.01,0.36,36),linspace(0.366,0.368,21),linspace(0.4,10,97)];
%C0_list = [linspace(0.01,0.36,36),linspace(0.366,0.368,21)];
%C0_list = linspace(0.366,0.368,21);
%C0_list = linspace(0.4,6,57);
C0_list = linspace(0.1,4.5,45);

% Sizes
nT  = 2;
n1  = 6;
n2  = 12;
N   = 100;
sizes = [1, nT, n1, n2 + nT*(0:N-1)];

%% ---------------------------------------------------------
%  STAGE 1: Polymerization (run ONCE)
%% ---------------------------------------------------------
M  = 100;
Z0 = zeros(N+2,1);

F1 = @(t,Z) F_ode(Z, M, n1, n2, N, alpha, sigma, rho, omega, ...
                  beta, mu, lambda, 5, 0.1, 0.1, 10);

[~,Y] = ode23s(F1, [0 300], Z0);
Y = Y';

z0 = Y(:,end)/M;
z0(1) = 0;
z0(2) = 0;

M_agreg  = sizes(2:end)*z0;
z0_adim  = z0 / M_agreg;

%% ---------------------------------------------------------
%  3D FIGURE
%% ---------------------------------------------------------
figure; hold on; grid on;
xlabel('t'); ylabel('Mw'); zlabel('C0');
title('3D Trajectories in (t, Mw, C0)');
view(45,25);

%% ---------------------------------------------------------
%  LOOP OVER C0 VALUES
%% ---------------------------------------------------------
for cidx = 1:length(C0_list)

    C0 = C0_list(cidx);
    fprintf('Running C0 = %.3f\n', C0);

    %% ---------------------------------------------------------
    %  LATE-STAGE DYNAMICS (base IC)
    %% ---------------------------------------------------------
    Z0_base = z0_adim * C0;

    F2 = @(t,Z) F_ode(Z, C0, n1, n2, N, alpha, sigma, rho, omega, ...
                      beta, mu, lambda, gamma, tau, delta, k1);

    [t_base, Y_late] = ode15s(F2, [0 500], Z0_base);
    Y_late = Y_late';

    % W_base = Y_late(2,:); % no longer needed
    Mw_base = 1 - (sizes(2:end)*Y_late)/C0 + (sizes(2:end).^2 * Y_late)/C0;

    %% ---------------------------------------------------------
    %  EQUILIBRIA FOR THIS C0 (optional, can remain the same)
    %% ---------------------------------------------------------
    Sstar = @(W) (gamma*W + lambda*W.^2) ./ (2*(mu + k1*W));
    X = @(W) (4*beta^3*Sstar(W).^4)./(27*sigma*(mu + k1*W).^2);
    u_plus  = @(W) ((mu + k1*W)./(2*sigma)) .* (1 + sqrt(1 + X(W))) .* Sstar(W);
    u_minus = @(W) ((mu + k1*W)./(2*sigma)) .* (1 - sqrt(1 + X(W))) .* Sstar(W);
    Tstar = @(W) nthroot(u_plus(W),3) + nthroot(u_minus(W),3);

    A = @(W) delta*(alpha + tau*Sstar(W) + rho*Tstar(W)) + ...
             tau*(omega*Tstar(W) + 3*lambda*W.^2);

    B = @(W) 4*delta*tau*(alpha + rho*Tstar(W)).*(omega*Tstar(W) + 3*lambda.*W.^2);

    Vstar_sq_minus = @(W) (A(W)./(2*tau*(alpha + rho*Tstar(W)))) .* ...
                          (1 - sqrt(1 - B(W)./A(W).^2));

    Vstarminus = @(W) sqrt(Vstar_sq_minus(W));
    qminus     = @(W) (tau*Vstarminus(W).^2)./delta;

    Ei_starm = @(i,W) (qminus(W)).^i .* Sstar(W);

    fWminus = @(W) Vstarminus(W) + 2*Tstar(W) + 6*W ...
                   + 12*Sstar(W) ...
                   + ((14*qminus(W) - 12*qminus(W).^2)./((1-qminus(W)).^2)).*Sstar(W) ...
                   - C0;

    % Root finding
    Wmin = 0.0; Wmax = C0/6;
    nScan = 20000;
    Wgrid = linspace(Wmin, Wmax, nScan);
    Fgrid = fWminus(Wgrid);

    idx = find(Fgrid(1:end-1).*Fgrid(2:end) < 0);
    Wstarm = Wgrid(idx);
    Sstarm = Sstar(Wstarm);

    % Mw* for equilibria
    Mwstarm = [];
    for i = 1:length(Wstarm)
        Ceq = [ Tstar(Wstarm(i)) ; Wstarm(i) ; Sstarm(i) ];
        for j = 1:N-1
            Ceq = [Ceq ; Ei_starm(j, Wstarm(i))];
        end
        Mw_eq = 1 - (sizes(2:end)*Ceq)/C0 + (sizes(2:end).^2 * Ceq)/C0;
        Mwstarm = [Mwstarm Mw_eq];
    end

    %% ---------------------------------------------------------
    %  MULTIPLE TRAJECTORIES FOR THIS C0
    %% ---------------------------------------------------------
    idxW = 2; idxS = 3;

    nIC  = 1;
    pert = linspace(0, 0, nIC);

    colors = lines(nIC);

    for k = 1:nIC
        Zk = Z0_base;
        Zk(idxW) = Zk(idxW) * (1 + pert(k));
        Zk(idxS) = Zk(idxS) * (1 - pert(k));

        [t_traj, Y_traj] = ode15s(F2, [0 500], Zk);
        Y_traj = Y_traj';

        % W_traj  = Y_traj(2,:); % not used
        Mw_traj = 1 - (sizes(2:end)*Y_traj)/C0 + (sizes(2:end).^2 * Y_traj)/C0;

        % PLOT IN (t, Mw, C0) SPACE
        plot3(t_traj, C0*ones(size(t_traj)), Mw_traj, ...
              'Color', colors(k,:), 'LineWidth', 1.2);
    end

end