clc;
clear;
close all;

% Parametrii fizici.
kB = 1.38e-23;  % Constanta lui Boltzmann.
u = 1.66e-27;   % Unitatea atomica de masa. 
T = 300;        % Temperatura.
vmax = 6e3;     % Viteza maxima in reprezentarile grafice.
m0 = 28*u;       % Masa moleculei.

% Parametrii sistemului.
numParticles = 10;  % Numarul de particule.
roomWidth = 10;     % Latimea camerei.
roomHeight = 10;    % Inaltimea camerei.

% Generez perechi aleatoare de coordonate initiale.
x = rand(1, numParticles) * roomWidth;
y = rand(1, numParticles) * roomHeight;

% Generez unghiuri aleatoare pentru fiecare particula.
angles = 2*pi*rand(1, numParticles);

% Generez viteze initiale aleatoare cu o distributie conform
% functiei de distributiei Maxwell. De asemenea returnez si
% functia de densitate a probabilitatii (utilizata la reprezentarile grafice).
[v, fNormalized] = generateMaxwellVelocities(kB, m0, T, numParticles, vmax);

% Pasul.
dt = 0.00004;
% Numarul momentelor de timp.
N = 1000;

% Initialize matrices to store particle positions
xPositions = zeros(N, numParticles);
yPositions = zeros(N, numParticles);

% Alocare de memorie pentru matricile care vor memora
% coordonatele fiecarei particule pe timpul miscarii.
xPositions(1, :) = x;
yPositions(1, :) = y;

% Simularea MRU.
for i = 2:N

    % Calculez urmatoarea pozitie.
    x = x + v .* cos(angles) * dt;
    y = y + v .* sin(angles) * dt;
    
    % Verific daca vreo particula s-a atins de perete.
    hitWallX = (x <= 0) | (x >= roomWidth);
    hitWallY = (y <= 0) | (y >= roomHeight);
    
    % Daca se loveste un perete vertical, se reflecta unghiul.
    angles(hitWallX) = pi - angles(hitWallX);
    % Daca se loveste un perete orizontal, unghiul este negat.
    angles(hitWallY) = -angles(hitWallY);
    
    % Pastreaza particula in camera.
    x(hitWallX) = max(0, min(roomWidth, x(hitWallX)));
    y(hitWallY) = max(0, min(roomHeight, y(hitWallY)));
    
    % Actualizeaza pozitiile.
    xPositions(i, :) = x;
    yPositions(i, :) = y;
end

% Reprezentarile grafice.
figure(1);
set(gcf, 'Position', [200, 200, 1400 700]);

% Figura cu pozitiile initiale.
subplot(2, 3, 1);
scatter(xPositions(1, :), yPositions(1, :), 'filled', 'r');

grid on;
axis([0 roomWidth 0 roomHeight]);
title('Pozitiile initiale.');
xlabel('X');
ylabel('Y');


% Figura cu distributia Maxwell si vitezele generate.
subplot(2, 3, 4);

% Maxwell.
temp = linspace(1, vmax, 1000);
plot(temp, fNormalized, 'b', 'LineWidth', 2);
hold on;

% Viteze.
temp = zeros(1, length(v));
scatter(v, temp, 'filled', 'r');

legend('', "Vitezele generate")
xlabel('Viteza');
ylabel('Densitatea de probabilitate');
title('Distributia Maxwell.');
grid on;

hold off;

% Figura cu simularea dinamica.
for i = 1:N
    hold off;
    subplot(2, 3, [2 3 5 6]);
    
    for j = 1:numParticles
        grid on;
        plot(xPositions(i, j), yPositions(i, j), '.', 'MarkerSize', 35);
        axis([0 roomWidth 0 roomHeight]);

        xlabel('X');
        ylabel('Y');
        title('Simularea dinamica.');
        hold on;
    end

    pause(1e-3);
end

