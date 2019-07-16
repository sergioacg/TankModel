function dhdt = TwoSphericalT(t,x,u,par)

%% Datos   
% Entradas
Fi1  = u(1);         %h^-1  (Flujo)
Fi2  = u(2);          %h^-1  (Flujo)

%Parametros
R1 = par.R1;  %Radio 1 [m]
R2 = par.R2;  %Radio 2 [m]
k1 = par.k1; %Coeficiente de la Válvula Tanque 1 [m^2/s]
k2 = par.k2; %Coeficiente de la Válvula Tanque 1 [m^2/s]
a1 = par.a1;
a2 = par.a2;
% 
%  Notación para las Variables de Estado
%
   h1   = x(1);     %height 1 [m]
   h2   = x(2);     %height 2 [m]
   
%   Areas
A1 = pi*(2*R1*h1-h1^2);
A2 = pi*(2*R2*h2-h2^2);
%
%  Ecuaciones Diferenciales del Sistema
%
dhdt = zeros(2,1);
dhdt(1)  = 1/A1 * (Fi1 -a1*k1*sqrt(h1));
dhdt(2)  = 1/A2 * (Fi2 + a1*k1*sqrt(h1) - a2*k2*sqrt(h2));


