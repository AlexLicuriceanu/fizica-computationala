function [velocities, fNormalized] = generateMaxwellVelocities(kB, m0, T, numParticles, vmax)    
    
    % Genereaza un interval de viteze intre 0 si viteza maxima data.
    v = linspace(0, vmax, 1000);
    
    f = Maxwell(v, kB, m0, T);
    % Normalizeaza functia de repartitie.
    fNormalized = f / sum(f);
    
    % Genereaza numere aleatoare intre 0 si 1.
    rng('shuffle');
    randomNumbers = rand(numParticles, 1);
    
    cumulativeDistance = cumsum(fNormalized);
    index = sum(randomNumbers >= cumulativeDistance, 2) + 1;
    velocities = v(index);
end


function f = Maxwell(v, kB, m0, T)
    f = 4*pi*v.^2*(m0/(2*pi*kB*T))^(3/2).*exp(-m0*v.^2/(2*kB*T));
end
