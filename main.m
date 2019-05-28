% Modelagem dos Tanques em Cascata
clc
clear all
close all

%Parametros
k1=0.04;
k2=0.03;
k3=0.055;
a1=0:0.1:1;   %Abertura da Valvula entrada
a2=0:0.1:1;   %Abertura da Valvula saida
A1=1;
A2=1.5;
%Ponto de Equilibrio das Valvulas
a1s=0.5;
a2s=0.45;

par.k1=k1; 
par.k2=k2; 
par.k3=k3; 
par.A1=A1;
par.A2=A2;

u=[a1s,a2s];
x0=[0,0];
t=0:0.01:1000;

%Perfil do Estado Estacionario
for i=1:length(a1)
    u(1)=a1(i);
    X(i,:) = fsolve(@(x)tankmodel(t,x,u,par),x0); %X=[h1 h2]
end

%Grafica do Estado Estacionario
plot(a1,X(:,2),'Linewidth',2)
ylabel('Altura tanque 2 (m)');
xlabel('Abertura Válvula a1');
title('Curva Estática com a2=0.45');

%% Obtenção do Modelo linear

%Estado estacionario seleccionado
u=[a1s,a2s];
X = fsolve(@(x)tankmodel(t,x,u,par),x0);

%Estados Estacionario
h1s = X(1);
h2s = X(2);

a11=-k2/(2*A1*sqrt(h1s));
a12=0;
a21= k2/(2*A1*sqrt(h1s));
a22=-k3*a2s/(2*A2*sqrt(h2s));

b11=k1/A1;
b21=0;

A=[a11 a12;a21 a22]
B=[b11;b21]
C=[0 1];

%Função de Transferencia
H = ss(A,B,C,0);
[z,P,k] = zpkdata(H); %Obtem os Zeros, Polos e Ganhos de H
G = zpk(z,P,k)

du=0.001;   %Pequena Variação na Entrada do Sistema
u(1)=a1s+du;

%Solução Númerica da ODE
Eq=[h1s h2s];
[ts,X1] = ode45(@(t,x)tankmodel(t,x,u,par), t , Eq);


u2=du*ones(length(ts),1); %Vetor de Entrada do Sistema Linear
%Solução do Sistema Linear
[ylin]=lsim(G,u2,ts);

%Grafica
figure
plot(ts,X1(:,2),'-b',ts,ylin+Eq(2),'--r','LineWidth',2);
title('Comparação modelo Linear vs Não Linear');
xlabel('tempo (s)');
ylabel('h_2 (m)');
legend('Não linear','Linear');