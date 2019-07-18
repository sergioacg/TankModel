% Controlador

%PID por Ziegler e Nichols
[z,P,k] = zpkdata(H);
G = zpk(z,P,k);

K=dcgain(H); %Ganho do sistema
L1=0.01;
L2=0.01;
tau1=-1/P{1,1};
tau2=-1/P{2,2};
K1=K(1,1);
K2=K(2,2);

% PID 1
Kc1=(0.9*tau1)/(K1*L1); 
ti1=3.33*L1; 
td1=0;
C1=tf(Kc1*[ti1*td1 ti1 1],[ti1 0]);

% PID 2
Kc2=(0.9*tau2)/(K2*L2); 
ti2=3.33*L2; 
td2=0;
C2=tf(Kc2*[ti2*td2 ti2 1],[ti2 0]);

%Malha Fechada (FT SISO)

GH1=feedback(C1*G(1,1),1);
GH2=feedback(C2*G(2,2),1);

figure
step(GH1)
figure
step(GH2)