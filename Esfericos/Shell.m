%Controle Coluna de Destilação da Shell MIMO 3x3
clc
clear all
close all
%% Sistema
Ps=[tf(4.05,[50 1]) tf(1.77,[60 1]) tf(5.88,[50 1]);...
    tf(5.39,[50 1]) tf(5.72,[60 1]) tf(6.9,[40 1]);...
    tf(4.38,[33 1]) tf(4.42,[44 1]) tf(7.2,[19 1]);];
Ps.iodelay=[27 28 27;18 14 15;20 22 0];
% Ps.iodelay=[7 8 7;8 4 5;2 3 1];

[~,polo,~] = zpkdata(Ps);

%% Matriz RGA
K=dcgain(Ps);
R=inv(K)';
%R=(K\eye(2))'
RGA=K.*R


L1=Ps.iodelay(1,1);
L2=Ps.iodelay(2,2);
L3=Ps.iodelay(3,3);

tau1=-1/polo{1,1};
tau2=-1/polo{2,2};
tau3=-1/polo{3,3};

K1=K(1,1);
K2=K(2,2);
K3=K(2,2);

% PID 1
Kc1=(0.9*tau1)/(K1*L1); 
Kc1=Kc1/2;
ti1=3.33*L1; 
td1=0;
C1=tf(Kc1*[ti1*td1 ti1 1],[ti1 0]);

% PID 2
Kc2=(0.9*tau2)/(K2*L2); 
Kc2=Kc2/2;
ti2=3.33*L2; 
td2=0;
C2=tf(Kc2*[ti2*td2 ti2 1],[ti2 0]);

% PID 3
Kc3=(0.9*tau3)/(K3*L3); 
Kc3=Kc3/2;
ti3=3.33*L3; 
td3=0;
C3=tf(Kc3*[ti3*td3 ti3 1],[ti3 0]);

C(1,1)=C1;
C(2,2)=C2;
C(3,3)=C3;

%Malha Fechada (FT SISO)

GH1=feedback(C1*Ps(1,1),1);
GH2=feedback(C2*Ps(2,2),1);
GH3=feedback(C3*Ps(3,3),1);

figure
step(GH1)
figure
step(GH2)
figure
step(GH3)

