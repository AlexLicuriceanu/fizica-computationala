clc;
clear;
close all;
% toate marimile sunt exprimate in unitati SI

% Rata la care sa se desfasoare simularea dinamica. Eg. simrate = 3,
% simuleaza de trei ori mai repede decat viteza reala. simrate 1 =
% simularea in timp real.
simrate = 4;

% Parametrii fizici:
g = 9.81;           % acceleratia gravitationala (N/kg)
m0 = 194;           % masa totala initiala (kg)
D = 0.18;           % diametrul rachetei (m)
mc0 = 0.72 * m0;    % masa de combustibil lichid (kg)
ro0 = 1.22;         % densitatea aerului (kg/m^3)

% Conditii initiale:
v0 = 16;        % viteza initiala (m/s)
alpha0 = 53;    % unghiul de lansare (deg)
tau = 57;       % timp de ardere (s)
u = 3880;       % viteza de evacuare a gazelor (m/s)

eta = 1.81 * 1e-5;      % coeficientul de vascozitate (Pa*s)
c1 = 6.54;
c2 = 0.64;
b1 = c1 * eta * D;      % coeficientul liniar
b2 = c2 * ro0 * D^2;    % coeficientul patratic

% Definirea intervalului de timp de inters
t0 = 0;                         % momentul initial
tf = 1000;                      % momentul final
N = 10000;                      % numarul momentelor de timp
t = linspace(t0, tf, N);        % vectorul momentelor
dt = t(2) - t(1);               % pasul

% Prealocare si valori de inceput:
vx = zeros(1,N);
vy = zeros(1, N);
x = zeros(1,N);
y = zeros(1, N);
mass = zeros(1, N);

vx(1) = v0 * cosd(alpha0);
vy(1) = v0 * sind(alpha0);


m = m0;         % masa rachetei care se modifica cand se arde carburant
mc = mc0;       % masa carburantului
engineon = 0;   % pentru a salva pana la ce moment este pornit motorul

for i = 1:N-1
    % Calculez componentele fortei de frecare.
    Frx = -b1 * vx(i) - b2 * sqrt(vx(i)^2 + vy(i)^2) * vx(i);
    Fry = -b1 * vy(i) - b2 * sqrt(vx(i)^2 + vy(i)^2) * vy(i);

    if t(i) <= tau && mc > 0
        % Calculez debitul de gaze de evacuare.
        q = mc / tau;

        % Actualizez masa rachetei si masa carburantului.
        m = m - q *dt;
        mc = mc - q * dt;

        % Calculez componentele fortei de tractiune.
        Fx = q * u * (vx(i)/sqrt(vx(i)^2 + vy(i)^2));
        Fy = q * u * (vy(i)/sqrt(vx(i)^2 + vy(i)^2));

        % Calculez componentele acceleratiei.
        ax = (Fx - Frx) / m;
        ay = (Fy - Fry - m * g) / m;

        aux = 1 + (q * (1 + u/sqrt(vx(i)^2 + vy(i)^2) ) - b1 - b2 * sqrt(vx(i)^2 + vy(i)^2)) * dt/(m - q*t(i));

        % Salvez ultimul moment in care este pornit motorul.
        engineon = i;
    else
        % Calculez componentele acceleratiei.
        ax = Frx / m;
        ay = (Fry - m * g) / m;

        aux = (1 + dt*(-b1/m - b2/m * sqrt(vx(i)^2 + vy(i)^2)));
    end

    % Aproximez urmatoarea pozitie a rachetei.
    vx(i+1) = vx(i)*aux + ax * dt;
    vy(i+1) = vy(i)*aux + ay * dt;
    x(i+1) = x(i) + vx(i+1) * dt;
    y(i+1) = y(i) + vy(i+1) * dt;

    % Salvez valorile masei rachetei.
    mass(i) = m;

    % Verifica daca racheta a ajuns la sol
    if y(i+1) < 0
        break;
    end
end

% Elimin valorile in surplus.
t = t(1:i);
vx = vx(1:i);
vy = vy(1:i);
x = x(1:i);
y = y(1:i);
mass = mass(1:i);

% Afisarea unor marimi de interes
tf = t(i);      % timpul de zbor
b = x(i);       % bataia
h = max(y);     % altitudinea maxima
tu = t(y==h);   % timpul de urcare
tc = tf - tu;   % timpul de coborare
Q = 1/2 * m * (v0^2 - vx(i)^2 - vy(i)^2);     % caldura produsa prin frecare

afis=['Timpul de zbor: ', num2str(tf),' s']; disp(afis);
afis=['Bataia rachetei: ', num2str(b/1e3),' km']; disp(afis);
afis=['Altitudinea maxima: ', num2str(h/1e3),' km']; disp(afis);
afis=['Timpul de urcare: ', num2str(tu),' s']; disp(afis);
afis=['Timpul de coborare: ', num2str(tc),' s']; disp(afis);
afis=['Caldura produsa: ', num2str(Q/1e6),' MJ']; disp(afis);

figure('units', 'normalized', 'outerposition', [0 0 0.8 0.8]);

% Grafic pentru componenta verticala a vitezei
subplot(3, 3, 1)
plot(t(1:engineon), vx(1:engineon), '-b', 'LineWidth', 2);
hold on;
plot(t(engineon:i), vx(engineon:i), '-r', 'LineWidth', 2);
grid;
xlabel({'Timp (s)'});
ylabel({'Vx (m/s)'});
title({'Componenta verticala a vitezei'});
hold off;

% Grafic pentru componenta orizontala a vitezei
subplot(3, 3, 2)
plot(t(1:engineon), vy(1:engineon), '-b', 'LineWidth', 1.5);
hold on;
plot(t(engineon:i), vy(engineon:i), '-r', 'LineWidth', 1.5);
grid;
xlabel({'Timp (s)'});
ylabel({'Vy (m/s)'});
title({'Componenta orizontala a vitezei'});
hold off;

% Grafic pentru variatia masei
subplot(3, 3, 3)
plot(t(1:engineon), mass(1:engineon), '-b', 'LineWidth', 1.5);
hold on;
plot(t(engineon:i), mass(engineon:i), '-r', 'LineWidth', 1.5);
grid;
ylim([0 m0]);
xlabel({'Timp (s)'});
ylabel({'Masa (kg)'});
title({'Variatia masei'});
hold off;

% Grafic pentru traiectorie
subplot(3, 3, [6 9])
plot(x(1:engineon)/1e3, y(1:engineon)/1e3, '-b', 'LineWidth', 1.5);
hold on;
plot(x(engineon:i)/1e3, y(engineon:i)/1e3, '-r', 'LineWidth', 1.5);
grid;
xlim([0 inf]);
ylim([0 inf]);
xlabel({'Bataie (km)'});
ylabel({'Altitudine (km)'});
title({'Traiectorie'});
hold off;

% Grafic pentru simularea dinamica

tic;        % porneste cronometrul
simt = 0;   % retine timpul initial

while simt < tf
    subplot(3,3, [4 5 7 8]);
    plot(x(1:engineon)/1e3, y(1:engineon)/1e3, '-b', 'LineWidth', 1);
    hold on;
    plot(x(engineon:i)/1e3, y(engineon:i)/1e3, '-r', 'LineWidth', 1);
    hold on;
    
    xlabel('Bataie (km)'); ylabel('Altitudine (km)');
    grid;
    title('Simulare dinamica');
    hold on;

    % Cauta cel mai apropiat t din discretizare.
    index = abs(t-simt) == min(abs(t-simt)); 

    if simt <= tau
        plot(x(index)/1e3, y(index)/1e3, '*', 'MarkerSize', 15, 'Color', 'k');
    else
        plot(x(index)/1e3, y(index)/1e3, '.', 'MarkerSize', 25, 'Color', 'k');
    end

    % Legenda
    legend('Miscare cu tractiune', 'Miscare fara tractiune', 'Location', 'northwest');
    hold off;

    text(b/2/1e3, h/3/1e3, ['Vx = ', num2str(round(round(vx(index)))), ' m/s']);
    text((b/2-b/5)/1e3, h/3/1e3, ['t =', num2str(round(t(index))), ' s']);
    text((b/2+b/5)/1e3, h/3/1e3, ['Vy = ', num2str(round(round(vy(index)))), ' m/s']);
    pause(1e-3);
    simt = toc*simrate;
end
