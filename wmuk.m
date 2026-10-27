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

%% BIFURCATION CURVE W*(C0)

C0_values = linspace(0.1,0.5,41);

Wstar_values = NaN(size(C0_values));

for ic = 1:length(C0_values)
    disp(C0_values(ic));

    C0 = C0_values(ic);

    %Equilibria

    Sstar = @(W) (gamma*W + lambda*W.^2) ./ (2*(mu + k1*W));

    X = @(W) (4*beta^3*Sstar(W).^4)./(27*sigma*(mu + k1*W).^2);

    u_plus  = @(W) ((mu + k1*W)./(2*sigma)) .* ...
        (1 + sqrt(1 + X(W))) .* Sstar(W);

    u_minus = @(W) ((mu + k1*W)./(2*sigma)) .* ...
        (1 - sqrt(1 + X(W))) .* Sstar(W);

    Tstar = @(W) nthroot(u_plus(W),3) + nthroot(u_minus(W),3);

    A = @(W) delta*(alpha + tau*Sstar(W) + rho*Tstar(W)) ...
           + tau*(omega*Tstar(W) + 3*lambda*W.^2);

    B = @(W) 4*delta*tau*(alpha + rho*Tstar(W)) ...
        .*(omega*Tstar(W) + 3*lambda.*W.^2);

    Vstar_sq_minus = @(W) ...
        (A(W)./(2*tau*(alpha + rho*Tstar(W)))) ...
        .* (1 - sqrt(1 - B(W)./A(W).^2));

    Vstarminus = @(W) sqrt(Vstar_sq_minus(W));

    qminus = @(W) (tau*Vstarminus(W).^2)./delta;

    fWminus = @(W) ...
          Vstarminus(W) ...
        + 2*Tstar(W) ...
        + 6*W ...
        + 12*Sstar(W) ...
        + ((14*qminus(W)-12*qminus(W).^2)./(1-qminus(W)).^2).*Sstar(W) ...
        - C0;

    %Root search

    Wgrid = linspace(0,C0/6,5000);

    Fgrid = fWminus(Wgrid);

    idx = find(Fgrid(1:end-1).*Fgrid(2:end) < 0);

    if ~isempty(idx)

        % Keep largest equilibrium (change "end" to "1" if desired)
        a = Wgrid(idx(end));
        b = Wgrid(idx(end)+1);

        Wstar_values(ic) = fzero(fWminus,[a b]);

    end

end

%% FIGURE

figure;
plot(C0_values,Wstar_values,'k','LineWidth',2);
hold on;

yline(mu/k1,'r--','LineWidth',2);

xlabel('C_0');
ylabel('W^*');
title('EVOLUTION OF W^{*}(C_0)')

legend('W^*','\mu/k','Location','northwest');

grid on;