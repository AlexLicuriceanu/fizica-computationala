% Oscilator elastic dublu

clear; close all; clc;

% Selectori:

LM=0; % afisarea legilor de miscare eta(t) - 1

TR=0; % dinamica in timp real - 1

% Parametrii fizici ai sistemului:

m1=0.35; m2=0.15; % kg; masele corpurilor

ka=10; kb=5; kc=7; % N/m; constantele elastice ale resorturilor

La=0.5; Lb=0.5; Lc=0.5; % m; lungimile resorturilor nedeformate

% Deplasari si viteze initiale:

eta10=0.2; eta20=-0.2; % m

v10=-0.13; v20=0.13; % m/s

% Pulsatii caracteristice sistemului:

omega11=sqrt((ka+kb)/m1); omega12=sqrt(kb/m1);

omega21=sqrt(kb/m2); omega22=sqrt((kb+kc)/m2);

% Perioade caracteristice:

T11=2*pi/omega11; T12=2*pi/omega12; T21=2*pi/omega21; T22=2*pi/omega22;

Tmax=max([T11,T12,T21,T22]); % "timp" caracteristic sistemului

ti=0; tf=5*Tmax; N=100000; t=linspace(ti,tf,N); dt=t(2)-t(1); % timp discret

eta1=zeros(1,N); eta2=eta1; % prealocare deplasari

v1=zeros(1,N); v2=v1; % prealocare viteze

% Valori de start:

eta1(1)=eta10; eta2(1)=eta20; % deplasari initiale - valori de start pas 1

eta1(2)=eta10+v10*dt; eta2(2)=eta20+v20*dt; % deplasari la t=dt - valori de start pas 2

v1(1)=v10; v2(1)=v20; % viteze de start

for i=2:N-1

    % Recurente de ordinul II (vezi Curs 3):

    eta1(i+1)=2*eta1(i)-eta1(i-1)-dt^2*(omega11^2*eta1(i)-omega12^2*eta2(i));

    eta2(i+1)=2*eta2(i)-eta2(i-1)-dt^2*(omega22^2*eta2(i)-omega21^2*eta1(i));

    v1(i)=(eta1(i+1)-eta1(i))/dt;

    v2(i)=(eta2(i+1)-eta2(i))/dt;

end

v1(N)=v1(N-1); v2(N)=v2(N-1);

T1=1/2*m1*v1.^2; T2=1/2*m2*v2.^2;

T=T1+T2; % energia cinetica

Ua=1/2*ka*eta1.^2; Ub=1/2*kb*(eta2-eta1).^2; Uc=1/2*kc*eta2.^2;

U=Ua+Ub+Uc; % energia elastica

H=T+U; % hamiltoniana sistemului (energia totala)

if LM==1 % afisarea legilor de miscare

    figure(1);

    eta1min=min(eta1); eta2min=min(eta2); etamin=min(eta1min,eta2min);

    eta1max=max(eta1); eta2max=max(eta2); etamax=max(eta1max,eta2max);

    plot(t,eta1,'-r',t,eta2,'-b');

    xlabel('t / s'); ylabel('deplasari / m');

    grid; legend('oscilator 1','oscilator 2');

    title('Legile de miscare ale componentelor');

    axis([ti tf 1.1*etamin 1.5*etamax]); % [xmin xmax ymin ymax]

end

xs=0; x1=La+eta1; x2=La+Lb+eta2; xd=La+Lb+Lc; % coordonate x

% Grosimile resorturilor nedeformate (acelasi material):

ga=sqrt(ka*La); gb=sqrt(kb*Lb); gc=sqrt(kc*Lc);

% Grosimile resorturilor deformate (vectori linie):

ga=ga*La./(La+eta1); gb=gb*La./(La+eta2-eta1); gc=gc*La./(La-eta2);

coef=80; % controleaza dimensiunea grafica a corpurilor

rg1=coef*m1^(1/3); rg2=coef*m2^(1/3); % raze grafice

figure(2); % Simularea dinamica a oscilatiilor

tic; simt=0; % porneste cronometrul si initializeaza timpul simularii

while simt<=tf % ciclul grafic

    hold off; % sterge grafica anterioara

    j=abs(t-simt)==min(abs(t-simt)); % cauta cel mai apropiat t de simt

    plot([xs x1(j)],[0 0],'-m','LineWidth',ga(j)); hold on; % resort stanga (a)

    plot([x1(j) x2(j)],[0 0],'-g','LineWidth',gb(j)); % resort median (b)

    plot([x2(j) xd],[0 0],'-c','LineWidth',gc(j)); % resort dreapta (c)

    plot(x1(j),0,'.r','MarkerSize',rg1); % oscilator 1

    plot(x2(j),0,'.b','MarkerSize',rg2); % oscilator 2

    plot(La,0,'+k',La+Lb,0,'+k'); % pozitiile de echilibru

    xlabel('x / m');

    axis([xs xd -5 5]);

    if TR==1 % in timp real

        simt=toc; % actualizeaza timpul simularii cu ceasul sistemului

        text(0.8*xd,4,['t = ',num2str(round(t(j))),' s']);

    else

        simt=simt+1e-2; % incrementeaza timpul simularii

        text(0.1*xd,4.0,['t = ',num2str(round(t(j)*10)),' ds']);

        text(0.8*xd,4.0,['T = ',num2str(round(T(j)*1e3)),' mJ']);

        text(0.8*xd,3.5,['U = ',num2str(round(U(j)*1e3)),' mJ']);

        text(0.8*xd,3.0,['H = ',num2str(round(H(j)*1e3)),' mJ']);

        text(x1(j),1.0,['T1 = ',num2str(round(T1(j)*1e3)),' mJ'],'Color','r');

        text(x2(j),-1.0,['T2 = ',num2str(round(T2(j)*1e3)),' mJ'],'Color','b');

        text(0.1*xd,-2,['Ua = ',num2str(round(Ua(j)*1e3)),' mJ'],'Color','m');

        text(0.4*xd,-2,['Ub = ',num2str(round(Ub(j)*1e3)),' mJ'],'Color','g');

        text(0.75*xd,-2,['Uc = ',num2str(round(Uc(j)*1e3)),' mJ'],'Color','c');

    end

    pause(1e-6)

end