%% Sergio Castaño
% Rio de Janeiro - 2019
% Dos tanques Esféricos proceso MIMO
%%
clc
close all
clear all

FS=14; %Tamanho da Fonte das Figuras

nit=300;
Ts=1;
tsim = 0:Ts:(nit-1)*Ts; %Tiempo de Simulación


%Parametros
R1 = 4;  %Radio 1 [m]
R2 = 6;  %Radio 2 [m] 
k1 = 3.2; %Coeficiente de la Válvula Tanque 1 [m^2/s]
k2 = 6.8; %Coeficiente de la Válvula Tanque 1 [m^2/s]
a1 = 0.4;
a2 = 0.5;
par.R1 = R1;  %Radio 1 [m]
par.R2 = R2;  %Radio 2 [m] 
par.k1 = k1; %Coeficiente de la Válvula Tanque 1 [m^2/s]
par.k2 = k2; %Coeficiente de la Válvula Tanque 1 [m^2/s]
par.a1 = a1;
par.a2 = a2;

%Entrada
F1=1.75;
F2=1;
u=[F1 F2];

%Condição inicial
x0=[0.1  0.1];

%Simulação do Sistema
[t,x] = ode45(@(t,x)TwoSphericalT(t,x,u,par), tsim , x0, u);

%Grafica do Sistema
figure
plot(t,x(:,1),[t(1) t(end)],[F1 F1],'linewidth',2)
title('Altura h_1');
ylabel('h_1');
xlabel('tempo (s)');
set(gca,'fontsize',FS)

figure
plot(t,x(:,2),[t(1) t(end)],[F2 F2],'linewidth',2)
title('Altura h_2');
ylabel('h_2');
xlabel('tempo (s)');
set(gca,'fontsize',FS)


%Grafico estacionario
F1x=0:0.1:3;
F2x=1;
H1x=(F1x/(k1*a1)).^2;
H2x=(F2x/(k2*a2)).^2+(F1x/(k2*a2)).^2+((2*F1x*F2x)/(a2*k2)^2);

figure
subplot(2,1,1);
plot(F1x,H1x),grid
title('Variación de nivel h1');
ylabel('Altura (h)');
xlabel('Variação F1 com F2=1');
set(gca,'fontsize',FS)
subplot(2,1,2);
plot(F1x,H2x),grid
title('Variación de nivel h2');
ylabel('Altura (h)');
xlabel('Variação F1 com F2=1');
set(gca,'fontsize',FS)

%Encuentra el estado Estacionario
X = fsolve(@(x)TwoSphericalT(t,x,u,par),x0); %X=[h1 h2]

%Establesco los Estados en el Estado Estacionario
h1s = X(1);
h2s = X(2);
Eq=[h1s h2s F1 F2];

%% Linealizacion del sistema por Jacobiana (Libreria SYMBOLICA)
% syms x1 x2 u1 u2
% 
% A1 = pi*(2*R1*x1-x1^2);
% A2 = pi*(2*R2*x2-x2^2);
% 
% % Funciones para realizar la jacobiana
% fx1 = 1/A1 * (u1 -a1*k1*sqrt(x1));
% fx2  = 1/A2 * (u2 + a1*k1*sqrt(x1) - a2*k2*sqrt(x2));
% 
% %Conforma vectores
% fx = [fx1;fx2];
% x = [x1;x2];
% uj = [u1;u2]; 
% %Linealiza por medio de la Jacobiana
% A = jacobian(fx,x);
% B = jacobian(fx,uj);
% C = eye(2); %Matriz de salida (Cb)
% 
% %Reemplazo los puntos de equilibrio en la Jacobiana
% display('Representación en Espacio de Estado Caso 1')
% Ad = double(subs(A,{x1,x2,u1,u2},{Eq}))
% Bd = double(subs(B,{x1,x2,u1,u2},{Eq}))



% Linealización usando Jacobiana
%Cálculo da matriz A
a11 = (-F1*(2*pi*R1-2*pi*h1s))/(2*pi*R1*h1s-pi*h1s^2)^2 ...
       + ( (-k1*a1)/(2*sqrt(h1s))*(2*pi*R1*h1s-pi*h1s^2)+k1*a1*sqrt(h1s)*(2*pi*R1-2*pi*h1s) ) / (2*pi*R1*h1s-pi*h1s^2)^2;
a12 = 0;

a21 = ( (k1*a1)/(2*sqrt(h1s)) ) / (2*pi*R2*h2s-pi*h2s^2);

a22 = (-F2*(2*pi*R2-2*pi*h2s))/(2*pi*R2*h2s-pi*h2s^2)^2 ...
       + ( (-k2*a2)/(2*sqrt(h2s))*(2*pi*R2*h2s-pi*h2s^2)+k2*a2*sqrt(h2s)*(2*pi*R2-2*pi*h2s) ) / (2*pi*R2*h2s-pi*h2s^2)^2 ...
       +( -k1*a1*sqrt(h1s)*(2*pi*R2-2*pi*h2s) ) / (2*pi*R2*h2s-pi*h2s^2)^2;


A = [a11 a12      
     a21 a22];      

%Cálculo da matriz B
B = [1/(pi*(2*R1*h1s-h1s^2))  0
     0                      1/(pi*(2*R2*h2s-h2s^2))];

%Cálculo da matriz C
C = [1 0 
     0 1];

%Cálculo da matriz D
D = [0 0
     0 0];

% %Determino la Funcion de Transferencia para la entrada 1 Cb
% display('Representación en Función de Transferencia para Cb')
% [num1,den1]=ss2tf(A,B,C,D,1);
% ft1=tf(num1,den1)

H = ss(A,B,C,D);
[z,P,k] = zpkdata(H);
K = zpk(z,P,k)

%Comprobar Modelo Lineal VS Modelo NO Lineal
du=0.1;
u(1)=F1+du;
[t,x] = ode45(@(t,x)TwoSphericalT(t,x,u,par), tsim , X, u);

in=zeros(length(tsim),2);
in(:,1)=du*ones(length(tsim),1);
in(:,2)=0*ones(length(tsim),1);
[ylin]=lsim(H,in,tsim);

%Grafica
figure
plot(t,x(:,1),'-b',tsim,ylin(:,1)+X(1),'--r','LineWidth',2)
ylabel('h_1','fontsize',FS);
titulo='Comparação modelo Linear vs Não Linear (h_1) para u= ','fontsize',FS;
titulo=strcat(titulo,num2str(F1),'h^{-1}');
title(titulo);
xlabel('tempo (s)','fontsize',FS);
legend('Não Linear','Linear','Location','SouthEast')
set(gca,'fontsize',FS);


figure
plot(t,x(:,2),'-b',tsim,ylin(:,2)+X(2),'--r','LineWidth',2)
ylabel('h_2','fontsize',FS);
titulo="Comparação modelo Linear vs Não Linear (h_2) para u= ";
titulo=strcat(titulo,num2str(F2),"h^{-1}");
title(titulo);
xlabel('tempo (s)','fontsize',FS);
legend('Não Linear','Linear','Location','SouthEast')
set(gca,'fontsize',FS);
