clear all;
clc;

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

% Concentrations to compare
C0_list = [0.5, 0.4, 0.37];
colors  = {'g', 'b', 'r'};

% Sizes
nT = 2; n1 = 6; n2 = 12;
N = 100;
sizes = [1, nT, n1, n2 + nT*(0:N-1)];

% Initial conditions (adimensional)
IC_list = [0.1   0.1   0.02;    % C0 = 0.5 highest
            0.1   0.1   0.013;    % C0 = 0.4 lowest
            0.1   0.1   0.022 ];  % C0 = 0.37 mid

assert(size(IC_list,1) == numel(C0_list), ...
    'IC_list must have one row per entry of C0_list.');

figure;
hold on;
grid on;

legend_labels = cell(1, numel(C0_list));

for k = 1:length(C0_list)
    C0 = C0_list(k);
% Initial condition for this C0
    T0 = IC_list(k,1);
    W0 = IC_list(k,2);
    S0 = IC_list(k,3);
    z0_adim = zeros(N+2, 1);
    z0_adim(1) = T0;
    z0_adim(2) = W0;
    z0_adim(3) = S0;
% Scale to concentration C0
    Z0 = C0 * z0_adim;
    F2 = @(t,Z) F_ode(Z, C0, n1, n2, N, ...
        alpha, sigma, rho, omega, beta, mu, ...
        lambda, gamma, tau, delta, k1);
    [t, Y] = ode15s(F2, [0 2000], Z0);
    Y = Y';
% Compute Mw
    C = zeros(N+3, length(t));
    C(2:end,:) = Y / C0;
    C(1,:) = 1 - sizes(2:end) * Y / C0;
    Mw = C(1,:) + sizes(2:end).^2 * Y / C0;
    h(k) = plot(t, Mw, 'LineWidth', 2, 'Color', colors{k});
    legend_labels{k} = sprintf('C_0 = %.2f (T_0=%.2f, W_0=%.2f, S_0=%.3f)', ...
        C0, T0, W0, S0);
end
xlabel('Time');
ylabel('Mw');
title('Mw trajectories for different initial conditions');
legend(h(end:-1:1), legend_labels(end:-1:1), 'Location', 'best');