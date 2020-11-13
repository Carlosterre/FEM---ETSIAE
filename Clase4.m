%% Codigo MEF para elementos viga hermitica
clear all

%% Definicion del material
E=1000; L=1; a0=L/20; A0=a0^2; I0=a0^4/12; zmax=-a0/2;
%Cargas
fa0=10; ff0=-10; Px=1; Py=1; Mz=1;

%% Definicion geometrica
gdln=3;  %Grados de libertad por nodo
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
    gdl_e(e,:)=[index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln ...  %[Axil Flecha Giro
                index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln];
end

%% Cuadraturas de integracion
ipoint1=0; w1=2;  %1 punto de integracion
ipoint2=1/sqrt(3)*[-1 1]; w2=[1 1];  %2 puntos de integracion (elemento viga hermitica de seccion cte)
ipoint3=sqrt(3/5)*[-1 0 1]; w3=[5/9 8/9 5/9];  %3 puntos de integracion

%% Matriz de rigidez
K=zeros(GDL);  %Se declara la matriz de rigidez completa vacia
Fn=zeros(GDL,1);  %Fuerza nodal equivalente de la viga
for e=1:nele
    index=conectividad_e(e,:);  %Nodos del elemento
    x1=coordx(index(1));
    x2=coordx(index(2));  %Coordenadas de los nodos
    Le=x2-x1;  %deltax=x2-x1; deltay=y2-y1;  Le=sqrt(deltax^2+deltay^2) para obtener cosenos directores
    Je=Le/2; iJe=Je^-1;
    
    %Integracion de la rigidez axil
    Kea=zeros(2);
    for ip=1:1 %Parte axil para seccion variable cambiamos el 1:1 a 1:2
        xi=ipoint1(ip);  %Para seccion variable xi=ipoint2
        n1=(1-xi)/2;  %Fciones de forma que se interpolan entre 2 nodos
        n2=(1+xi)/2;
        N=[n1 n2];
        b1=-1/2;
        b2=1/2;
        Ba=iJe*[b1 b2];  %B axil
        Kea=Kea+Ba'*E*A0*Ba*Je*w1(ip);
    end
    
    %Integracion de la rigidez flectora
    Kef=zeros(4);
    for ip=1:2
        xi=ipoint2(ip);  %Para seccion variable xi=ipoint3
        %h1=;...  %Polinomios hermiticos
        k1=(xi*3/2); k2=(xi*3/2-1/2); k3=(-xi*3/2); k4=(xi*3/2+1/2);  %Curvatura, derivada segunda de los polinomios hermiticos
        K0f=iJe^2*[k1 Je*k2 k3 Je*k4];  %Correccion en los terminos de giro con el jacobiano
        Kef=Kef+K0f'*E*I0*K0f*Je*w2(ip);
    end
    
    gdle=gdl_e(e,:);  %Gdl de cada elemento
    gdle1=gdle([1 4]);  %Gdl axiles de cada elemento
    gdle2=gdle([2 3 5 6]);  %Gdl de cada elemento en flexion
   
    K(gdle1,gdle1)=K(gdle1,gdle1) + Kea; %Se sustituye la rigidez axil
    K(gdle2,gdle2)=K(gdle2,gdle2) + Kef; %Se sustituye la rigidez a flexion
    
    %Fuerzas nodales equivalentes
    fe=zeros(6,1);
    fa1=fa0*(x1/L);  %f(x)=fa0*x/L triangular
    fa2=fa0*(x2/L);
    ff1=ff0*(x1/L);  %f(x)=ff0*x/L triangular
    ff2=ff0*(x2/L);
    for ip=1:3
        xi=ipoint3(ip);
        n1=(1-xi)/2;
        n2=(1+xi)/2;
        N=[n1 n2];
        h1=(2-3*xi+xi^3)/4;
        h2=(1-xi-xi^2+xi^3)/4;
        h3=(2+3*xi-xi^3)/4;
        h4=(-1-xi+xi^2+xi^3)/4;
        H=[h1 h2 h3 h4];
       
        fa_=N*[fa1;fa2];  %Fuerza axil en los nodos. El comando ; genera vector columna
        fa=N'*fa_*w3(ip)*Je;
       
        ff_=N*[ff1;ff2];  %Fuerza cortante en los nodos
        ff=H'*ff_*w3(ip)*Je;
       
        fe([1 4])=fe([1 4])+fa;  %Ensamblaje de la fuerza axil
        fe([2 3 5 6])=fe([2 3 5 6])+ff;  %Ensamblaje de la fuerza flectora
    end
    
    Fn(gdle)=Fn(gdle)+fe;       
end

K
Fn

%% Condiciones de contorno
U=zeros(GDL,1);
uR=[1 2 3 (GDL-1)];  %Restringidos, [1 2 3] si el extremo esta en voladizo
uL=setdiff((1:length(U)),uR);  %Libres
F=Fn;
F(GDL-2)=F(GDL-2)+Px; F(GDL-1)=F(GDL-1)+Py; F(GDL)=F(GDL)+Mz;

%Reduccion de las CC
FL=F(uL);
KLL=K(uL,uL);
KLR=K(uL,uR);
UL=KLL\FL;
FR=KLR'*UL;
U(uL)=UL