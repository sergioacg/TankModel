%Analise de Interações RGA

K=dcgain(H);
R=inv(K)';
%R=(K\eye(2))'
% Matriz RGA
RGA=K.*R

%% estabilidade tendo por base o índice de Niederlinski
for i=1:length(F1x)
    u(1)=F1x(i);
    %Encuentra el estado Estacionario
    X = fsolve(@(x)TwoSphericalT(t,x,u,par),x0); %X=[h1 h2]

    %Establesco los Estados en el Estado Estacionario
    h1s = X(1);
    h2s = X(2);

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

    % %Determino a Funcion de Transferencia
    H = ss(A,B,C,D);
    
    K1=dcgain(H);
    R=inv(K1)';
    % Matriz RGA
    RGA=K1.*R;
    lambda(i)=RGA(1,1);
    
    %estabilidade tendo por base o índice de Niederlinski
    
    Ni(i)=det(K1)/prod(diag(K1));
end

% Grafica de los Puntos de Equilibrio
figure
subplot(2,1,1)
plot(F1x,lambda,'Linewidth',2)
%axis([0 200 -1.4 2])
title('Relação entre \lambda vs F/V','FontSize',FS);
ylabel('\lambda','FontSize',FS);
set(gca,'fontsize',FS);
xlabel('F_1 (h^{-1})','FontSize',FS);
subplot(2,1,2)
plot(F1x,Ni,'Linewidth',2)
%axis([0 160 -1 2])
title('Relação entre Indice Niederlinski vs F/V','FontSize',FS);
ylabel('Indice Niederlinski','FontSize',FS);
xlabel('F_1 (h^{-1})','FontSize',FS);
set(gca,'fontsize',FS);