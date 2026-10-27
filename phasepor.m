clear all; clc;

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

% Sweep over C0
C0_list = linspace(0.1, 5, 491);

nT  = 2;
n1  = 6;
n2  = 12;
N   = 100;
sizes = [1, nT, n1, n2 + nT*(0:N-1)];

%% STAGE 1
M  = 100;
Z0 = zeros(N+2,1);

F1 = @(t,Z) F_ode(Z, M, n1, n2, N, alpha, sigma, rho, omega, ...
                  beta, mu, lambda, 5, 0.1, 0.1, 10);

[~, Y] = ode23s(F1, [0 300], Z0);
Y = Y';

z0 = Y(:,end) / M;
z0(1) = 0;
z0(2) = 0;

M_agreg = sizes(2:end) * z0;
z0_adim = z0 / M_agreg;

%% FIGURES
figure; hold on; grid on;
xlabel('C0'); ylabel('W'); zlabel('Mw');
title('3D Trajectories in (C0, W, Mw)');
view(45,25);

CMW_store = [];

%% LOOP
for cidx = 1:length(C0_list)

    C0 = C0_list(cidx);
    fprintf('Running C0 = %.3f\n', C0);

    %% 2ND PHASE
    Z0_base = z0_adim * C0;

    F2 = @(t,Z) F_ode(Z, C0, n1, n2, N, alpha, sigma, rho, omega, ...
                      beta, mu, lambda, gamma, tau, delta, k1);

    [t_base, Y_late] = ode15s(F2, [0 500], Z0_base);
    Y_late = Y_late';

    W_base  = Y_late(2,:);
    Mw_base = 1 - (sizes(2:end) * Y_late)/C0 + (sizes(2:end).^2 * Y_late)/C0;

    %% EQUILIBRIA
    Sstar = @(W) (gamma*W + lambda*W.^2) ./ (2*(mu + k1*W));

    X = @(W) (4*beta^3 * Sstar(W).^4) ./ (27*sigma*(mu + k1*W).^2);

    u_plus  = @(W) ((mu + k1*W)./(2*sigma)) .* (1 + sqrt(1 + X(W))) .* Sstar(W);
    u_minus = @(W) ((mu + k1*W)./(2*sigma)) .* (1 - sqrt(1 + X(W))) .* Sstar(W);

    Tstar = @(W) nthroot(u_plus(W),3) + nthroot(u_minus(W),3);

    A = @(W) delta*(alpha + tau*Sstar(W) + rho*Tstar(W)) + ...
             tau*(omega*Tstar(W) + 3*lambda*W.^2);

    B = @(W) 4*delta*tau*(alpha + rho*Tstar(W)) .* ...
             (omega*Tstar(W) + 3*lambda*W.^2);

    Vstar_sq_minus = @(W) (A(W)./(2*tau*(alpha + rho*Tstar(W)))) .* ...
                          (1 - sqrt(1 - B(W)./A(W).^2));

    Vstarminus = @(W) sqrt(Vstar_sq_minus(W));
    qminus     = @(W) (tau * Vstarminus(W).^2) ./ delta;

    Ei_starm = @(i,W) (qminus(W)).^i .* Sstar(W);

    fWminus = @(W) Vstarminus(W) + 2*Tstar(W) + 6*W ...
                   + 12*Sstar(W) ...
                   + ((14*qminus(W) - 12*qminus(W).^2) ./ ((1 - qminus(W)).^2)) .* Sstar(W) ...
                   - C0;

    % Root finding
    Wmin = 0.0; 
    Wmax = C0/6;
    nScan = 20000;

    Wgrid = linspace(Wmin, Wmax, nScan);
    Fgrid = fWminus(Wgrid);

    idx = find(Fgrid(1:end-1).*Fgrid(2:end) < 0);
    Wstarm = Wgrid(idx);
    Sstarm = Sstar(Wstarm);

    % Mw* for equilibria
    Mwstarm = [];
    for i = 1:length(Wstarm)
        Ceq = [Tstar(Wstarm(i)); Wstarm(i); Sstarm(i)];
        for j = 1:N-1
            Ceq = [Ceq; Ei_starm(j, Wstarm(i))];
        end
        Mw_eq = 1 - (sizes(2:end) * Ceq)/C0 + (sizes(2:end).^2 * Ceq)/C0;
        Mwstarm = [Mwstarm Mw_eq];
    end

    %% MULTIPLE TRAJECTORIES
    idxW = 2; 
    idxS = 3;

    nIC  = 1;
    pert = linspace(0, 0, nIC);

    colors = lines(nIC);

    for k = 1:nIC
        Zk = Z0_base;
        Zk(idxW) = Zk(idxW) * (1 + pert(k));
        Zk(idxS) = Zk(idxS) * (1 - pert(k));

        [t_traj, Y_traj] = ode15s(F2, [0 500], Zk);
        Y_traj = Y_traj';

        W_traj  = Y_traj(2,:);
        Mw_traj = 1 - (sizes(2:end) * Y_traj)/C0 + (sizes(2:end).^2 * Y_traj)/C0;

        % NEW: Store all (C0, W, Mw)
        CMW_store = [CMW_store; ...
                     C0 * ones(length(W_traj),1), ...
                     W_traj(:), ...
                     Mw_traj(:)];

        plot3(C0*ones(size(W_traj)), W_traj, Mw_traj, ...
              'Color', colors(k,:), 'LineWidth', 1.2);

        plot3(C0*ones(size(Wstarm)), Wstarm, Mwstarm, 'r.','MarkerSize',14);
    end

end

%% PLOT TRAJECTORIES

figure; hold on; grid on;

C0_all = CMW_store(:,1);
W_all  = CMW_store(:,2);
Mw_all = CMW_store(:,3);

C0_unique = unique(C0_all);

colormap(turbo);
caxis([min(Mw_all) max(Mw_all)]);
colorbar;
ylabel(colorbar, 'Mw intensity');

for i = 1:length(C0_unique)

    C0_val = C0_unique(i);

    % Extract trajectory points
    idx = (C0_all == C0_val);
    W_traj  = W_all(idx);
    Mw_traj = Mw_all(idx);

    X = [C0_val * ones(1,length(W_traj)); 
         C0_val * ones(1,length(W_traj))];

    Y = [W_traj'; 
         W_traj'];

    Z = [Mw_traj'; 
         Mw_traj'];

    C = Z;   % color = Mw

    surface(X, Y, Z, C, ...
        'FaceColor', 'none', ...
        'EdgeColor', 'interp', ...
        'LineWidth', 2);
end

xlabel('C0');
ylabel('W');
zlabel('Mw');
title('Trajectories');
view(45,25);