 function x = tankmodel(~,x,u,par)

 %Parâmetros do Modelo
k1  = par.k1; 
k2  = par.k2; 
k3  = par.k3; 
A1  = par.A1;
A2  = par.A2;

%Entradas do Modelo
a1  = u(1);
a2  = u(2);


%Variáveis de Estado
h1   = x(1);
h2   = x(2);

% Sistema de Equações
dh1dt = 1/A1*(k1*a1-k2*sqrt(h1));
dh2dt = 1/A2*(k2*sqrt(h1)-k3*a2*sqrt(h2));

%Vetor de Estados
x = [dh1dt;dh2dt];