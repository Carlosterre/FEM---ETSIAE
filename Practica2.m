%% DATOS DE ENTRADA
E=1000; L=1; f0=20; a0=L/40; P=1;

%% SOL ANALITICA FUERTE (1500 elementos)
u_strong=24.6473
sigma_max_strong=41557

%% DATOS DE MALLA
nele=9;  %9 elementos
nnod=2*nele+1;
coordx=linspace(0,L,nnod);
ngdl=nnod;  %Numero de grados de libertad

%Conectividad
conectividad=zeros(nele,3);

for e=1:nele  
 conectividad(e,:)=[(e-1)*2+1 (e-1)*2+2 (e-1)*2+3];
end

K=zeros(ngdl);  %Matriz de rigidez global
f=zeros(ngdl,1);  %Vector de fzas nodales
u=zeros(ngdl,1);  %Vector de desplazam

%% BUCLES PARA DETERMINAR LA RIGIDEZ GLOBAL Y LA FZA NODAL GLOBAL
%Ptos de integracion
w1=2;  %1 pto de integracion (cuadratura de orden 1)
ipoint1=0;

w2=[1 1];  %2 ptos de integracion (cuadratura de orden 2)
ipoint2=[-1/sqrt(3) 1/sqrt(3)];

w3=[5/9 8/9 5/9]; %3 ptos de integracion (cuadratura de orden 3)
ipoint3=[-sqrt(3/5) 0 sqrt(3/5)]; 

%Integracion numerica
K_almacen=zeros(3,3,nele);

for e=1:nele  
 nodo1=conectividad(e,1);
 nodo2=conectividad(e,2);
 nodo3=conectividad(e,3);
 
 x1=coordx(nodo1);
 x2=coordx(nodo2);
 x3=coordx(nodo3);
 
 Le=x3-x1;
 Je=Le/2;
 
 gdl1=nodo1;  %Grados de libertad para el ensamblaje
 gdl2=nodo2;
 gdl3=nodo3;
 
 A1=a0^2*(4*L-3*x1)^2/16*L^2;  %Areas nodales
 A2=a0^2*(4*L-3*x2)^2/16*L^2;
 A3=a0^2*(4*L-3*x3)^2/16*L^2;
 An=[A1; A2; A3];

 f1=f0*(x1/L)^3;  %Fzas nodales
 f2=f0*(x2/L)^3;
 f3=f0*(x3/L)^3;
 fn=[f1; f2; f3];

 Ke=zeros(3);
 
 for ip=1:3
  xi=ipoint3(ip);  %Integ point
  w=w3(ip);
 
  N1=0.5*xi*(xi-1);  %Fciones  de forma
  N2=1-xi^2;
  N3=0.5*xi*(xi+1);
  N=[N1 N2 N3];

  A_=N*An;

  B1=(xi-0.5);  %Matriz de compatibilidad
  B2=-2*xi;
  B3=(xi+0.5);
  B=1/Je*[B1 B2 B3];

  Ke=Ke+B'*E*A_*B*Je*w;
 end

 K_almacen(:,:,e)=K_almacen(:,:,e)+Ke(:,:);
 index_gdl=conectividad(e,1:3);
 
 %Ensamblaje
 K(index_gdl,index_gdl)=K(index_gdl,index_gdl)+Ke;
  
 fe=zeros(3,1);
  
 for ip=1:3
  xi=ipoint3(ip);  %Integ point
  w=w3(ip);
  
  N1=0.5*xi*(xi-1);
  N2=1-xi^2;
  N3=0.5*xi*(xi+1);
  N=[N1 N2 N3];

  f_=N*fn;
  fe=N'*f_*Je*w+fe;
 end
 
 %Ensamblaje
 f(index_gdl)=f(index_gdl)+fe;
end 
     
gradosL=(2:ngdl);
gradosR=(1);
gLoad=ngdl;

f(gLoad)=f(gLoad)+P;

%Reduccion (condensacion estatica por aplicacion de las condiciones de contorno)
KLL=K(gradosL,gradosL);
KLR=K(gradosL,gradosR);
fL=f(gradosL);

uL=KLL\fL;  %uL=KLL^-1*fL;

K_quest=K  %Matriz de rigidez global
f_quest=f  %Vector global de fzas nodales

u=zeros(ngdl,1);
u(gradosL)=uL;

fR=KLR'*uL;  %f=K*u

u_quest=u
error_u=(abs((u_strong-u(19,1))/u_strong))*100;
error_u=round(error_u,3)  %Error (%) en desplazam entre la solucion  fuerte y debil

%%  VECTOR DE DEFORM Y TENSION EN LOS PTOS DE INTEGRACION
epsilon=zeros(nele,3);
sigma=zeros(nele,3);

for e=1:nele
 gdl1=conectividad(e,1);
 gdl2=conectividad(e,2);
 gdl3=conectividad(e,3);
 ue=u([gdl1 gdl2 gdl3]);
 
 x1=coordx(gdl1);
 x2=coordx(gdl2);
 x3=coordx(gdl3);
 
 Le=x3-x1;
 Je=Le/2;
 
 for ip=1:3
  xi=ipoint3(ip);  %Integ point

  B1=(xi-0.5);
  B2=-2*xi;
  B3=xi+0.5;
  B=1/Je*[B1 B2 B3];
  
  epsilon(e,ip)=B*ue;
  sigma(e,ip)=epsilon(e,ip)*E; 
 end
end

%sigma_max, error_smax
sigma_max=max(sigma(:));
sigma_max=round(sigma_max,0)  %Tension max calculada en los ptos de integracion

error_smax=(abs((sigma_max_strong-sigma_max)/sigma_max_strong))*100;
error_smax=round(error_smax,3)  %Error (%) en tensiones max entre la solucion fuerte y debil

%sigma_end
epsilon_end=zeros(nele,3);
sigma_end=zeros(nele,3);
ipoint_end=[-1 0 1];

for e=1:nele
 gdl1=conectividad(e,1);
 gdl2=conectividad(e,2);
 gdl3=conectividad(e,3);
 ue=u([gdl1 gdl2 gdl3]);
 
 x1=coordx(gdl1);
 x2=coordx(gdl2);
 x3=coordx(gdl3);
 
 Le=x3-x1;
 Je=Le/2;
 
 for ip=1:3
  xi=ipoint_end(ip);  %Integ point, nodos

  B1=(xi-0.5);
  B2=-2*xi;
  B3=xi+0.5;
  B=1/Je*[B1 B2 B3];
  
  epsilon_end(e,ip)=B*ue;
  sigma_end(e,ip)=epsilon_end(e,ip)*E; 
 end
end

sigma_end=sigma_end(nele,3);
sigma_end=round(sigma_end,0)  %Tension en el extremo libre

% %sigma_end obtenido con interpolacion lineal
% xpto1=x_i(1,2);  
% xpto2=x_i(1,3);
% sigmapto1=sigma(9,2);
% sigmapto2=sigma(9,3);
% m=(sigmapto2-sigmapto1)/(xpto2-xpto1);
% n=sigmapto1-m*xpto1;
% sigmapto3=m*L+n;
% sigma_end=sigmapto3;
% sigma_end=round(sigma_end,0)  %Tension en el extremo libre

%% REPRESENTACION GRAFICA DE LA DISTRIB DE TENSIONES EN LOS PTOS DE INTEGRACION
for e=1:nele    
 nodo1=conectividad(e,1);
 nodo2=conectividad(e,2);
 nodo3=conectividad(e,3);
 
 x1=coordx(nodo1);
 x2=coordx(nodo2);
 x3=coordx(nodo3);
 x_=[x1; x2; x3];
 
 for ip=1:3    
  xi=ipoint3(ip);  %Integ point
  N1=0.5*xi*(xi-1);
  N2=1-xi^2;
  N3=0.5*xi*(xi+1);
  N=[N1,N2,N3];       
  x_i(ip)=N*x_;  %Coord global del pto de integracion
 end
 
 grafstress=plot(x_i,sigma(e,:),'bs-','LineWidth',1);
 hold on
end

grid on  %Grafica de tensiones
title('Tensiones calculadas en puntos de integración')
xlabel('Posición en longitud L')
ylabel('Tensión en elementos')