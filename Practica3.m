%% DEFINICION DEL MATERIAL, E[=]MPa, L[=]mm, b[=]mm
E=70e3; L=1000; a0=L/10; b=1;

%Cargas
Px=1000; Py=-50; %Mz=0

%% DEFINICION GEOMETRICA
gdln=3;  %Grados de libertad (gdl) por nodo
nele=10;  %Numero de elementos
nnod=nele+1;  %Numero de nodos
coord=zeros(nnod,2); coordx=coord(:,1); coordy=coord(:,2);
for n=0:nele
    coordx(n+1)=n*(L/nele);
    coord(n+1,1)=coordx(n+1);
    coord(n+1,2)=coordy(n+1);
end

GDL=nnod*gdln;
conectividad_e=zeros(nele,2);
gdl_e=zeros(nele,2*gdln);  %Almacen gdl por elemento
for e=1:nele
    conectividad_e(e,:)=[e e+1];
    index=conectividad_e(e,:);
    gdl_e(e,:)=[index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...  %[Axil Flecha Giro
                index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln];
end

%% CUADRATURAS DE INTEGRACION
ipoint1=0; w1=2;  %1 pto de integracion
ipoint2=1/sqrt(3)*[-1 1]; w2=[1 1];  %2 ptos de integracion

%% MATRIZ DE RIGIDEZ Y FZAS NODALES EQUIVALENTES
K=zeros(GDL);  %Matriz de rigidez completa vacia
for e=1:nele
    index=conectividad_e(e,:);  %Nodos del elemento
    x1=coordx(index(1));  %Coord de los nodos
    x2=coordx(index(2));
    Le=x2-x1;
    Je=Le/2; iJe=Je^-1;
    
    A1=b*a0*exp(-2*x1/L);  %Areas nodales
    A2=b*a0*exp(-2*x2/L);
    An=[A1;A2];
    
    I1=b*a0^3*exp(3*(-2*x1/L))/12;  %Momentos de inercia
    I2=b*a0^3*exp(3*(-2*x2/L))/12;
    In=[I1;I2];
    
    %Integracion de la rigidez axil
    Kea=zeros(2);
    for ip=1:2
        xi=ipoint2(ip);  %Ptos de integracion
        
        n1=(1-xi)/2;  %Fciones de forma, 2 nodos
        n2=(1+xi)/2;
        N=[n1 n2];
        
        A_=N*An;
        
        b1=-1/2;  %B axil, derivada de las fciones de forma (Ba)
        b2=1/2;
        Ba=iJe*[b1 b2];
        
        Kea=Kea+Ba'*E*A_*Ba*Je*w2(ip);
    end
    
    %Integracion de la rigidez flectora
    Kef=zeros(4);
    for ip=1:2
        xi=ipoint2(ip);  %Ptos de integracion
        
        k1=xi*3/2; k2=xi*3/2-1/2; k3=-xi*3/2; k4=xi*3/2+1/2;  %Curvatura, derivada segunda de las fciones de forma (Bf)
        K0f=iJe^2*[k1 Je*k2 k3 Je*k4];
        
        n1=(1-xi)/2;  %Fciones de forma, 2 nodos
        n2=(1+xi)/2;
        N=[n1 n2];
        
        I_=N*In;
        
        Kef=Kef+K0f'*E*I_*K0f*Je*w2(ip);
    end
    
    gdle=gdl_e(e,:);  %Gdl totales de cada elemento
    gdle1=gdle([1 4]);  %Gdl axiles de cada elemento
    gdle2=gdle([2 3 5 6]);  %Gdl en flexion de cada elemento
    
    K(gdle1,gdle1)=K(gdle1,gdle1)+Kea; %Se sustituye la rigidez axil
    K(gdle2,gdle2)=K(gdle2,gdle2)+Kef; %Se sustituye la rigidez a flexion  
end

K_quest=round(K,0)

%% COND DE CONTORNO (CC)
U=zeros(GDL,1);  %Se declara vector de desplazam vacio
uR=[1 2 3];
uL=setdiff(1:length(U),uR);
F=zeros(GDL,1);  %Vector de fzas nodales nulo
F(GDL-2)=F(GDL-2)+Px; F(GDL-1)=F(GDL-1)+Py; %F(GDL)=F(GDL)+Mz

%Reduccion de las CC
FL=F(uL);
KLL=K(uL,uL);
KLR=K(uL,uR);
UL=KLL\FL;
FR=KLR'*UL;
U(uL)=UL;

%Desplazam nodo de carga, U[=]mm
ux_end=U(GDL-2,1);
uy_end=U(GDL-1,1);
uxy_end=U(GDL,1);
U_end=[ux_end; uy_end; uxy_end];

U_quest=round(U_end,1)

%% REPRESENTACION GRAFICA DEFORMADA ELASTICA
U_round=round(U,1);
Ux=U_round(1:3:GDL); Uy=U_round(2:3:GDL); Uxy=U_round(3:3:GDL);

escala=0.1*L/(max(abs(Uy)));
hold on
fn=1; figure(fn);
grafdefelas=plot(coordx,coordy,'b*-',coordx+escala*Ux,coordy+escala*Uy,'r*--','LineWidth',2)
grid on
axis equal
title('Deformada elástica','Fontsize',10);
xlabel('Posiciones nodales X','Fontsize',10)  %Etiqueta eje horizontal
ylabel('Posiciones nodales Y','Fontsize',10)  %Etiqueta eje vertical

%% TENSION EN LOS NODOS
epsilon=zeros(nele,2);
sigmaa=zeros(nele,2);
sigmaf=zeros(nele,2);
sigma=zeros(nele,2);
ipointn=[-1 1];
for e=1:nele
    index=conectividad_e(e,:);  %Nodos del elemento
    x1=coordx(index(1));  %Coord de los nodos
    x2=coordx(index(2));
    Le=x2-x1;
    Je=Le/2; iJe=Je^-1;
          
    gdle=gdl_e(e,:);  %Gdl totales de cada elemento
    gdle1=gdle([1 4]);  %Gdl axiles de cada elemento
    gdle2=gdle([2 3 5 6]);  %Gdl en flexion de cada elemento
    Uea=U(gdle1);
    Uef=U(gdle2);
    
    z1=-b*a0*exp(-2*x1/L)/2;  %Zmax
    z2=-b*a0*exp(-2*x2/L)/2;
    Z=[z1;z2];
    
    %Tension axil en los nodos
    for ip=1:2
    xi=ipointn(ip);  %Ptos de integracion, nodos
     
    b1=-1/2;
    b2=1/2;
    Ba=iJe*[b1 b2];
    
    epsilon(e,ip)=Ba*Uea;
    sigmaa(e,ip)=epsilon(e,ip)*E;
    end
     
    %Momento flector en los nodos
    for ip=1:2
        xi=ipointn(ip);  %Ptos de integracion, nodos
        
        n1=(1-xi)/2;  %Fciones de forma, 2 nodos
        n2=(1+xi)/2;
        N=[n1 n2];
        
        Z_=N*Z;
              
        k1=xi*3/2; k2=xi*3/2-1/2; k3=-xi*3/2; k4=xi*3/2+1/2;  %Curvatura (Bf)
        K0f=iJe^2*[k1 Je*k2 k3 Je*k4];
        
        K_(e,ip)=K0f*Uef;
        
        sigmaf(e,ip)=E*K_(e,ip)*Z_; 
    end
    
    sigma=sigmaa+sigmaf;
end

sigma;

%% REPRESENTACION GRAFICA DE LA DISTRIB DE TENSIONES NODALES
for e=1:nele    
    nodo1=conectividad_e(e,1);
    nodo2=conectividad_e(e,2);
 
    x1=coordx(nodo1);
    x2=coordx(nodo2);
    x_=[x1;x2];
 
    for ip=1:2    
        xi=ipointn(ip);  %Integ point, nodos
  
        N1=(1-xi)/2;  %Fciones de forma
        N2=(1+xi)/2;
        N=[N1,N2];
  
        x_n(ip)=N*x_;  %Coord global de los nodos
    end
    
    hold on
    fn=2; figure(fn);
    grafstressnod=plot(x_n,sigma(e,:),'bs-','LineWidth',2);
end

grid on
title('Tensiones calculadas en nodos','Fontsize',10)
xlabel('Posición en longitud L','Fontsize',10)
ylabel('Tensión en elementos','Fontsize',10)

%% TENSION EN EL EMPOTRAMIENTO, sigma_first
sigma_first=sigma(1,1);
sigma1_quest=round(sigma_first,1)
%Error
sigma_first_strong=40;  %Obtenida para 1500 elem
err_first=abs(100*(sigma_first_strong-sigma_first)/sigma_first_strong);
err1=round(err_first,0)

%% TENSION MAX EN EL VOLADIZO, sigma_end
sigma_end=sigma(nele,2);
sigmaend_quest=round(sigma_end,1)
%Error
sigma_end_strong=73.8;  %Obtenida para 1500 elem
err_end=abs(100*(sigma_end_strong-sigma_end)/sigma_end_strong);
err2=round(err_end,0)

%% TENSION EN LOS PTOS DE INTEGRACION
epsilonpi=zeros(nele,1);
sigmaapi=zeros(nele,1);
sigmafpi=zeros(nele,1);
sigmapi=zeros(nele,1);
for e=1:nele
    index=conectividad_e(e,:);  %Nodos del elemento
    x1=coordx(index(1));  %Coord de los nodos
    x2=coordx(index(2));
    Le=x2-x1;
    Je=Le/2; iJe=Je^-1;
          
    gdle=gdl_e(e,:);  %Gdl totales de cada elemento
    gdle1=gdle([1 4]);  %Gdl axiles de cada elemento
    gdle2=gdle([2 3 5 6]);  %Gdl en flexion de cada elemento
    Uea=U(gdle1);
    Uef=U(gdle2);
    
    z1=-a0*exp(-2*x1/L)/2;  %Zmax
    z2=-a0*exp(-2*x2/L)/2;
    Z=[z1;z2];
    
    %Tension axil en los ptos de integracion
    for ip=1:2
    xi=ipoint2(ip);  %Ptos de integracion
     
    b1=-1/2;
    b2=1/2;
    Ba=iJe*[b1 b2];
    
    epsilonpi(e,ip)=Ba*Uea;
    sigmaapi(e,ip)=epsilonpi(e,ip)*E;
    end
     
    %Momento flector en los ptos de integracion
    for ip=1:2
        xi=ipoint2(ip);  %Ptos de integracion
        
        n1=(1-xi)/2;  %Fciones de forma, 2 nodos
        n2=(1+xi)/2;
        N=[n1 n2];
        
        Z_=N*Z;
              
        k1=(xi*3/2); k2=(xi*3/2-1/2); k3=(-xi*3/2); k4=(xi*3/2+1/2);  %Curvatura (Bf)
        K0f=iJe^2*[k1 Je*k2 k3 Je*k4];
        
        K_(e,ip)=K0f*Uef;
        
        sigmafpi(e,ip)=E*K_(e,ip)*Z_; 
    end
    
    sigmapi=sigmaapi+sigmafpi;
end

sigmapi;

%% TENSION MAX EN PTOS DE INTEGRACION, sigmaint_max
sigmaint_max=round(max(max(sigmapi)),0);
sigmamax_quest=round(sigmaint_max,1)
%Error
sigmaint_max_strong=197.1;  %Obtenida para 1500 elementos
err_sigmaint_max=abs(100*(sigmaint_max_strong-sigmaint_max)/sigmaint_max_strong);
err3=round(err_sigmaint_max,0)