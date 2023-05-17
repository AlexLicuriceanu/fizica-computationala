function MaxwellDistribLaw

%Studiul functiei de distributie Maxwell dupa modulul vitezei moleculelor unui gaz ideal

%1. Se studiaza dependenta de temperatura

%2. Se studiaza dependenta de natura gazului

%3. Se calculeaza probabilitati

%3a. Verificarea normarii functiei de distributie

%3b. Se calculeaza probabilitatea ca modulul vitezei sa fie cuprins intr-un interval


clc;clear all;close all

kB=1.38e-23; %constanta lui Boltzmann, J/K

u=1.66e-27; %unitatea atomica de masa, kg


m0=[2 28 32]*u;%masa unei molecule de H_2, N_2, O_2; kg

T=[300 400 500]; %temperatura termodinamica a gazului, K


%f=@(v,kB,m,T) 4*pi*v.^2*(m0/(2*pi*kB*T))^(3/2).*exp(-m0*v.^2/(2*kB*T));


vmax=4e3; %viteza maxima in reprezentarile grafice, m/s

v=linspace(0,vmax,201);


color=['b' 'g' 'r'];


%1. Se studiaza dependenta de temperatura

%Graficul functiei de distributie pentru m0(1) la temperaturile din vectorul T

subplot(211)

hold on

for i=1:3

    %y=f(v,kB,m0(1),T(i));

    f=Maxwell(v,kB,m0(1),T(i));

    plot(v,f,color(i))

end

legend('300 K','400 K','500 K')

title('Dependenta de temperatura a functiei de distributie Maxwell a hidrogenului')


%2. Se studiaza dependenta de natura gazului

%Graficul functiei de distributie pentru T(1) si gazele din vectorul m0

subplot(212)

hold on

for i=1:3

    %y=f(v,kB,m0(i),T(1));

    f=Maxwell(v,kB,m0(i),T(1));

    plot(v,f,color(i))

end

legend('H_2','N_2','O_2')

title('Dependenta de natura gazului a functiei de distributie Maxwell la {\it{T}} = 300 K')

xlabel('{\it{v}} / (m/s)')

ylabel('{\it{f}}({\it{v}}) / (m/s)^{-1}')


%3. Se calculeaza probabilitati

%3a. Verificarea normarii functiei de distributie

%Consideram gazul 1 aflat la temperatura T(1)

disp('Integrala de la 0 la infinit a functiei de distributie pentru m0(1) si T(1) este:')

integral(@(v) Maxwell(v,kB,m0(1),T(1)),0,inf)


%3b. Se calculeaza probabilitatea ca modulul vitezei sa fie cuprins intr-un interval

%Consideram gazul 1 aflat la temperatura T(1) si intervalul de viteze

%v1=500 m/s, v2=3000 m/s

v1=500;%m/s

v2=3000;%m/s

disp('Probabilitatea ca modulul vitezei unei molecule a gazului 1 aflat la temp T(1) sa fie cuprins intre v1 si v2 este:')

integral(@(v) Maxwell(v,kB,m0(1),T(1)),v1,v2)

end


function f=Maxwell(v,kB,m0,T)

f=4*pi*v.^2*(m0/(2*pi*kB*T))^(3/2).*exp(-m0*v.^2/(2*kB*T));

end