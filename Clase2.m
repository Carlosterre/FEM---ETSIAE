E=1;
A=1;
Le=1;
Je=Le/2;
L=1;
a0=L/20;

%Punto x de integracion 1
xi=-1/sqrt(3);  w1=1;
B1=1/Je*[xi-0.5 -2*xi xi+0.5];
Ke1=B1'*E*A*B1*Je*w1;

%Punto x de integracion 2
xi=1/sqrt(3);  w2=1;
B2=1/Je*[xi-0.5 -2*xi xi+0.5];
Ke2=B2'*E*A*B2*Je*w2;

%Matriz de rigidez del elemento
Ke=Ke1+Ke2;

%% Bucle para integrar el elemento
%Cuadratura de orden 2
w2=[1 1];
ipoint2=[-1 1]*1/sqrt(3) ;
Ke=zeros(3);

for ip=1:2
   
 xi=ipoint2(ip);
 w=w2(ip);
 Je=Le/2; 
 B=1/Je*[xi-0.5 -2*xi xi+0.5];
      
 Ke=Ke+B'*E*A*B*Je*w;

%% Para un elemento barra de 3 nodos
%Cuadratura de orden 3
w3=[5/9 8/9 5/9];
ipoint3=[-sqrt(3/5) 0 sqrt(3/5)];
Ke=zeros(3);

for ip=1:3
   
 xi=ipoint3(ip);
 w=w3(ip);
 Je=Le/2;
    
 x1=0;
 x2=L/2;
 x3=L;
    
 A1=(a0^2)*exp(-(4*x1)/(3*L));
 A2=(a0^2)*exp(-(4*x2)/(3*L));
 A3=(a0^2)*exp(-(4*x3)/(3*L));
    
 N1=0.5*xi*(xi-1);
 N2=1-xi^2;
 N3=0.5*xi*(xi+1);
    
 B=1/Je*[xi-0.5,-2*xi, xi+0.5];
 N=[N1, N2, N3]; 
 An=[A1 A2 A3]';  %An=[A1; A2; A3] tambien valdria para introducir vector columna
 A_=N*An;
      
 Ke=Ke+B'*E*A_*B*Je*w;
    
end


%% Vector de fuerzas nodales equivalentes
f0=1;  %Con una fuerza constante.
fe=zeros(3,1);

for ip=1:3
   
 xi = ipoint3(ip);
 w = w3(ip);
 Je=Le/2;
   
 N1=0.5*xi*(xi-1);
 N2=1-xi^2;
 N3=0.5*xi*(xi+1);
 N=[N1 N2 N3];  
 fn= [1 1 1]'*f0;
 f_= N*fn; %f virgulilla
    
 fe=fe+N'*f_*Je*w
 
end 
