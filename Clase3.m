E=1; A=1; L=1; a0=L/20; nele=3; nnod=2*nele+1;  %Numero de nodos para elementos de 3 nodos
coordx=linspace(0,L,nnod);
f0=1;
P=1;

ngdl=nnod;

%Conectividad
conectividad=zeros(nele,3);

for e=1:nele    
 conectividad(e,:)=[(e-1)*2+1 (e-1)*2+2 (e-1)*2+3]; %Grado de libertad 1 del elem e (e-1)*2+1
end

%% Matriz de rigidez

K=zeros(ngdl);  %Matriz de rigidez
f=zeros(ngdl,1);  %Se declara el vector de fuerzas
u=zeros(ngdl,1);

%% Bucle para integrar el elemento
w1=2;
ipoint1=0;

%Cuadratura de orden 2
w2=[1 1];
ipoint2=[-1/sqrt(3) 1/sqrt(3)];

%Cuadratura de orden 3
w3=[5/9 8/9 5/9];
ipoint3=[-sqrt(3/5) 0 sqrt(3/5)];

%% Integracion de la matriz de rigidez
K_almacen=zeros(3,3,nele);

for e=1:nele    
 nodo1=conectividad(e,1);
 nodo2=conectividad(e,2);
 nodo3=conectividad(e,3);
 x1=coordx(nodo1);
 x2=coordx(nodo2);
 x3=coordx(nodo3);
 Le=x3-x1;
 J=Le/2;
 gdl1=nodo1;  gdl2=nodo2;  gdl3=nodo3;  %Grados de libertad para el ensamblaje

 A1=a0^2*exp((-4*x1)/(3*L));  %L del solido
 A2=a0^2*exp((-4*x2)/(3*L));
 A3=a0^2*exp((-4*x3)/(3*L));
 An=[A1; A2; A3];

 f1=f0*(x1/L)^2;  % Valor de fza en el nodol
 f2=f0*(x2/L)^2;  % no depende puntos integracion, solo de los ptos
 f3=f0*(x3/L)^2;
 fn=[f1; f2; f3];

 Ke=zeros(3);

 for ip = 1:3

  xi=ipoint3(ip); %Integ point

  w=w3(ip);

  N1=0.5*xi*(xi-1);
  N2=1-xi^2;
  N3=0.5*xi*(xi+1);
  N=[N1 N2 N3];

  A_ = N*An;

  B1=(xi-0.5);
  B2=-2*xi;
  B3=(xi+0.5);

  B=J^-1*[B1 B2 B3];

  Ke=Ke+B'*E*A_*B*J*w;
 end

  K_almacen(:,:,e)=K_almacen(:,:,e)+Ke(:,:);
  index_gdl= conectividad(e,1:3);
  K(index_gdl,index_gdl)=K(index_gdl,index_gdl)+Ke; %Ensamblaje

  %Calculo del vector completo de fzas nodales
  fe=zeros(3,1);
  
 for ip=1:3
    
  xi=ipoint3(ip);  %Integ point
  w=w3(ip);
  
  N1=0.5*xi*(xi-1);
  N2=1-xi^2;
  N3=0.5*xi*(xi+1);
  N=[N1 N2 N3];

  f_=N*fn;
  fe=N'*f_*J*w+fe;
 end

 f(index_gdl)=f(index_gdl)+fe; %Ensamblaje
end

%% Calculo de desplazamientos
gradosL=(2:ngdl);  gradosR=(1);  gLoad = ngdl;
f(gLoad)=f(gLoad)+P;

%Reduccion por aplicacion de las condiciones de contorno
KLL=K(gradosL,gradosL);
KLR=K(gradosL,gradosR);
fL=f(gradosL);

%uL=KLL^-1*fL;
uL=KLL\fL;

%u=zeros(ngdl,1);
u(gradosL)=uL;

%f=K*u
fR=KLR'*uL;

%Tensiones se calculan con deformacion
ue1=u([1 2 3]);

for e=1:nele

 gdl1=conectividad(e,1);  gdl2=conectividad(e,2);  gdl3=conectividad(e,3);
 ue=u([gdl1 gdl2 gdl3]);
 x1=coordx(gdl1);
 x2=coordx(gdl2);
 x3=coordx(gdl3);
 Le=x3-x1;
 J=Le/2;

 for ip=1:3
 
  xi = ipoint3(ip);

  B1=(xi-0.5);
  B2=-2*xi;
  B3=(xi+0.5);
  B=J^-1*[B1 B2 B3];
  epsilon(e,ip)=B*ue;
  sigma(e,ip)=epsilon(e,ip)*E;
    
 end

end
