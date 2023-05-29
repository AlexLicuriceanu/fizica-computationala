clc;
clear;
close all;

% Parametrii fizici ai sistemului.
l = 1.5;            % Lungimea tijei (m)
m = 2;              % Masa corpului (kg)
theta0 = 60;   % Unghiul initial in grade [-180, 180].

% Transforma din grade in radiani.
theta0 = theta0 * pi / 180;

% Selector timp real (1) / slow motion (0).
RT = 1;

% Acceleratia gravitationala terestra.
g = 9.80665;    % m/s^2

% Pulsatia proprie.
omega0 = sqrt(g/l);
% Perioada proprie.
T0 = 2 * pi / omega0;

% Timpul sistemului.
ti = 0;         % Timp initial.
tf = 5 * T0;    % Timp final.
N = 2000;       % Numarul de momente de timp.
t = linspace(ti, tf, N);
dt = t(2)-t(1); % Pasul.

theta = zeros(1, N);
theta(1) = theta0;
theta(2) = theta0;

aux = g / l * dt^2;

tic
% Simulare.
for i = 2:N-1
    theta(i+1) = 2*theta(i) - theta(i-1) - aux * sin(theta(i));
end
toc

% Alocare memorie.
omega = zeros(1, N);
T = zeros(1, N);
U = zeros(1, N);

% Calculeaza vitezele unghiulare.
omega = (theta(2:end) - theta(1:end-1)) / dt;
omega(1) = 0;
omega(N) = (theta(N) - theta(N-1)) / dt;

% Calculeaza energia cinetica.
T = 1/2 * m * (l * omega).^2;
T(1) = 0;
% Calculeaza energia potentiala.
U = m * g * l * (1 - cos(theta));
% Calculeaza energia totala.
H = T + U;

% Coordonatele carteziene ale corpurilor.
x = l * sin(theta);
y = -l * cos(theta);

% Semilatura cadrului grafic.
lmax = l + 0.2;
% Controleaza dimensiunile grafice ale corpului.
coef = 30;
% Raza grafica a corpului.
rg = coef * m^(1/3);


figure('units', 'normalized', 'outerposition', [0.25 0.25 0.5 0.5]);

% Grafic pentru legea unghiulara de miscare.
subplot(1, 2, 1);
plot(t, theta * 180 / pi, '-r');
grid;
xlim([ti tf]);
xlabel('t (s)');
ylabel('theta (grade)');
title('Legea unghiulara de miscare.');
axis square;

% Simularea dinamica.
tic;
simt = 0;
while simt <= tf

  % Cauta cel mai apropiat t din discretizare.
  j = abs(t-simt) == min(abs(t-simt));
    
  % Fereastra grafica.
  subplot(1, 2, 2);


  % Tija.
  plot([0 x(j)], [0 y(j)], '-g', 'LineWidth', 2);
  hold on;

  % Articulatia de suspensie.
  plot(0, 0, '.k', 'MarkerSize', 10);
  % Corpul.
  plot(x(j), y(j), '.r', 'MarkerSize', rg);

  % Seteaza axele.
  axis([-lmax lmax -lmax lmax]);
  axis square;
  title('Simularea dinamica.');

  % Simulare in timp real sau in slow motion.
  if RT==1

    % Actualizeaza timpul simularii cu ceasul sistemului.
    simt = toc;
    text(3/5 * lmax, 4/5 * lmax, ['t = ', num2str(round(t(j))), ' s']);
  else

    % Incrementeaza cu o centisecunda.
    simt=simt+1e-3;
    text(3/5 * lmax, 4/5 * lmax, ['t =', num2str(round(t(j)*100)), ' cs']);
  end


  % Afisare energie cinetica, potentiala, totala.
  text(3/5 * lmax, 3.5/5 * lmax, ['E = ', num2str(round(H(j))), ' J']);
  text(3/5 * lmax, 2.5/5 * lmax, ['T = ', num2str(round(T(j))), ' J']);
  text(3/5 * lmax, 2/5 * lmax, ['U = ', num2str(round(U(j))), ' J']);

  pause(1e-6);
  hold off
end
