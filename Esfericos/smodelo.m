function [sys,x0]=smodelo(t,x,u,flag,par,h1s,h2s)
%O Flag 0 é o passo onde indicamos para a S-Function o que ela vai
%achar quando ler o modelo pela primeira vez. Isso é feito com un vetor de 6 elementos
%chamado nesse caso como sys
if flag==0
   %Elemento 1: quantidade de estados continuos (Equações diferenciaveis)= 2
   %Elemento 2: Numero de estados discretos: 0
   %Elemento 3: Numero de saídas do modelo: 2 (h1 e h2)
   %Elemento 4: Numero de Entradas do modelo: 2 entradas (a1 e a2)
   %Elemento 5: Parametro de control, colcar 0
   %Elemento 6: Tempo de Amostragem
   [sys]=[2,0,2,2,0,0];
   %Incluir condições iniciai
   x0=[h1s h2s];
elseif flag==1
    %Flag 1 chama modelo
   [sys]=TwoSphericalT(t,x,u,par);
elseif flag==3
    %Flag 3 indica a resposta que se debe obter, neste caso, um vetor con as 2
    %variables de estado
   [sys]=x;
else
    %Fecha o loop
   [sys]=[];
end
end