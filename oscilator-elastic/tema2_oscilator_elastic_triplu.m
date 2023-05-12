clear; close all; clc;
% Pentru simularea in timp real.
TR = 0;

% Parametrii fizici ai sistemului:
% kg; Masele corpurilor.
m1 = 0.35;
m2 = 0.15;
m3 = 0.45;

% N/m; Constantele elastice ale resorturilor.
ka = 10;
kb = 5;
kc = 7;
kd = 6;

% m; Lungimile resorturilor nedeformate.
La = 0.5;
Lb = 0.5;
Lc = 0.5;
Ld = 0.5;

% Deplasari si viteze initiale:
% m
eta10 = 0.2;
eta20 = -0.2;
eta30 = -0.1;

% m/s
v10 = -0.13;
v20 = 0.13;
v30 = -0.2;

% Pulsatii:
omega11 = sqrt(ka/m1);
omega21 = sqrt(kb/m1);
omega31 = sqrt(kc/m1);
omega41 = sqrt(kd/m1);

omega12 = sqrt(ka/m2);
omega22 = sqrt(kb/m2);
omega32 = sqrt(kc/m2);
omega42 = sqrt(kd/m2);

omega13 = sqrt(ka/m3);
omega23 = sqrt(kb/m3);
omega33 = sqrt(kc/m3);
omega43 = sqrt(kd/m3);

% "timp" caracteristic sistemului
Tmax = pi * 1.5/2;  % ales arbitrar
ti = 0;
tf = 5*Tmax;
N = 100000;
t = linspace(ti, tf, N);
dt = t(2) - t(1);

% Prealocare deplasari
eta1 = zeros(1, N);
eta2 = zeros(1, N);
eta3 = zeros(1, N);

% prealocare viteze
v1 = zeros(1, N);
v2 = zeros(1, N);
v3 = zeros(1, N);

% Valori de start:
% Deplasari initiale - valori de start pas 1
eta1(1) = eta10;
eta2(1) = eta20;
eta3(1) = eta30;

% Deplasari la t=dt - valori de start pas 2
eta1(2) = eta10 + v10*dt;
eta2(2) = eta20 + v20*dt;
eta3(2) = eta30 + v30*dt;

% Viteze de start
v1(1) = v10;
v2(1) = v20; 
v3(1) = v30;

for i = 2:N-1

    % Recurente de ordinul II

    aux1 = -omega11^2 * eta1(i) - omega21^2 * (eta1(i) - eta2(i)) - omega41^2 * (eta1(i) - eta3(i));

    aux2 = -omega22^2 * (eta2(i) - eta1(i)) - omega32^2 * (eta2(i) - eta3(i));

    aux3 = -omega43^2 * (eta3(i) - eta1(i)) - omega33^2 * (eta3(i) - eta2(i));

    eta1(i+1) = 2 * eta1(i) - eta1(i-1) + dt^2 * aux1;
    eta2(i+1) = 2 * eta2(i) - eta2(i-1) + dt^2 * aux2;
    eta3(i+1) = 2 * eta3(i) - eta3(i-1) + dt^2 * aux3;

    v1(i) = (eta1(i+1) - eta1(i))/dt;
    v2(i) = (eta2(i+1) - eta2(i))/dt;
    v3(i) = (eta3(i+1) - eta3(i))/dt;
end

v1(N) = v1(N-1);
v2(N) = v2(N-1);
v3(N) = v3(N-1);

T1 = 1/2 * m1 * v1.^2;
T2 = 1/2 * m2 * v2.^2;
T3 = 1/2 * m3 * v3.^2;

% Energia cinetica.
T = T1+T2+T3;

Ua = 1/2 * ka * eta1.^2;
Ub = 1/2 * kb * (eta2 - eta1).^2;
Uc = 1/2 * kc * (eta3 - eta2).^2;
Ud = 1/2 * kd * eta3.^2;

% Energia elastica.
U = Ua+Ub+Uc+Ud;

% Hamiltoniana sistemului (energia totala)
H = T+U;

% Coordonate x.
xs = 0;
x1 = La+eta1;
x2 = La+Lb+eta2;
x3 = La+Lb+Lc+eta3;
xd = La+Lb+Lc+Ld;

% Grosimile resorturilor nedeformate (acelasi material):
ga = sqrt(ka*La);
gb = sqrt(kb*Lb);
gc = sqrt(kc*Lc);
gd = sqrt(kd*Ld);

% Grosimile resorturilor deformate (vectori linie):
ga = ga * La./(La + eta1);
gb = gb * La./(La + eta2 - eta1);
gc = gc * La./(La - eta2);
gd = gd * La./(La - eta3);

% Controleaza dimensiunea grafica a corpurilor.
coef = 80;

% Raze grafice.
rg1 = coef*m1^(1/3);
rg2 = coef*m2^(1/3);
rg3 = coef*m3^(1/3);

figure('units', 'normalized', 'outerposition', [0.2 0.2 0.5 0.7]);

% Grafic pentru legile de miscare.
subplot(2, 2, [1 2]);
plot(t, eta1,'-r', t, eta2, '-b'); hold on;
plot(t, eta3, '-', 'color', [0 0.5 0]); hold off;

eta1min = min(eta1); eta2min = min(eta2); eta3min = min(eta3);
etamin = min([eta1min eta2min eta3min]);

eta1max = max(eta1); eta2max = max(eta2); eta3max = max(eta3);
etamax = max([eta1max eta2max eta3max]);

xlabel('t / s'); ylabel('deplasari / m');
legend('Oscilator 1','Oscilator 2', 'Oscilator 3');
title('Legile de miscare ale componentelor');
axis([ti tf 1.1*etamin 1.5*etamax]); % [xmin xmax ymin ymax]


% Porneste cronometrul si initializeaza timpul simularii.
tic; simt = 0;

while simt <= tf

    hold off;

    % Graficul simularii dinamice.
    subplot(2, 2, [3 4]);

    % Cauta cel mai apropiat t de simt.
    j = abs(t-simt)==min(abs(t-simt));

    % Resort a.
    plot([xs x1(j)], [0 0], '-m', 'LineWidth', ga(j)); hold on;
    % Resort b.
    plot([x1(j) x2(j)], [0 0], '-c', 'LineWidth', gb(j));
    % Resort c.
    plot([x2(j) x3(j)], [0 0], '-o', 'LineWidth', gc(j));
    % Resort d.
    plot([x3(j) xd], [0 0], '-', 'color', [0.4940 0.1840 0.5560], 'LineWidth', gd(j));

    % Oscilator 1.
    plot(x1(j), 0, '.r', 'MarkerSize', rg1);
    % Oscilator 2.
    plot(x2(j), 0, '.b', 'MarkerSize', rg2);
    % Oscilator 3.
    plot(x3(j), 0, '.', 'color', [0 0.5 0], 'MarkerSize', rg3);
    
    % Pozitiile de echilibru.
    plot(La, 0, '+k', La+Lb, 0, '+k', La+Lb+Lc, 0, '+k');

    title('Simularea dinamica');
    xlabel('x / m');
    axis([xs xd -5 5]);

    if TR == 1 % Simulare in timp real.

        % Actualizeaza timpul simularii cu ceasul sistemului.
        simt = toc;
        text(0.8*xd, 4, ['t = ', num2str(round(t(j))), ' s']);

    else
        % Incrementeaza timpul simularii.
        simt = simt+1e-2;

        text(0.1*xd, 4.0, ['t = ', num2str(round(t(j)*10)), ' ds']);
        text(0.8*xd, 4.0, ['T = ', num2str(round(T(j)*1e3)), ' mJ']);
        text(0.8*xd, 3.5, ['U = ', num2str(round(U(j)*1e3)), ' mJ']);
        text(0.8*xd, 3.0, ['H = ', num2str(round(H(j)*1e3)), ' mJ']);

        text(x1(j), 1.0, ['T1 = ', num2str(round(T1(j)*1e3)), ' mJ'], 'color', 'r');
        text(x2(j), -1.0, ['T2 = ', num2str(round(T2(j)*1e3)), ' mJ'], 'color', 'b');
        text(x3(j), 1.0, ['T3 = ', num2str(round(T3(j)*1e3)), ' mJ'], 'color', [0 0.5 0]);

        text(0.1*xd, -2, ['Ua = ', num2str(round(Ua(j)*1e3)), ' mJ'], 'color', 'm');
        text(0.3*xd, -2, ['Ub = ', num2str(round(Ub(j)*1e3)), ' mJ'], 'color', 'c');
        text(0.55*xd, -2, ['Uc = ', num2str(round(Uc(j)*1e3)), ' mJ'], 'color', [0.9290 0.6940 0.1250]);
        text(0.75*xd, -2, ['Ud = ', num2str(round(Ud(j)*1e3)), ' mJ'],'color', [0.4940 0.1840 0.5560]);

    end

    pause(1e-6)

end
