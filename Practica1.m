%INICIO
%% DEFINICIONES
E=100; A=100; I=5; EA=E*A; EI=E*I; L=1; 
q0=1000;

%% DEFINICION GEOMETRICA DEL MODELO
%Grados de libertad por nodo
GDLn=3;

%Coord de los nodos de la estructura por filas, coord x primera columna, coord y segunda columna
coord_nodos=[0 0; L/2 L/2; L L; 1.5*L L; 2*L L; 2*L L/2; 2*L 0];

%Coord de los nodos
xx=coord_nodos(:,1);  yy=coord_nodos(:,2);

%Numero de nodos
nnodos=size(xx,1);

%Numero de elementos
nelementos=size(coord_nodos,1)-1;

%Grados de libertad de toda la estructura
GDL=nnodos*GDLn; 

%% VECTOR DE FZAS NODALES EQUIVALENTES
%Vector inicial vacio
q=zeros(GDL,1);

%Solo hay fzas nodales equivalentes en nodos 1, 2, 3.
L13=L*sqrt(2);
L1=L13/2;  L2=L13/2;
alfa=atan(1);

%Distribución 12 triangular
F11=(3/20)*(q0/2*L1)*sin(alfa);  F21t=(7/20)*(q0/2*L1)*sin(alfa);
F12=(3/20)*(-q0/2*L1)*cos(alfa);  F22t=(7/20)*(-q0/2*L1)*cos(alfa);
F13=(-q0/2*L1^2)/30;  F23t=(-q0/2*L1^2)/20;

%Distribución 23 trapezoidal
F21tr=((3/20)*(q0*L2)+(7/20)*(q0/2*L2))*sin(alfa);  F31=((7/20)*(q0*L2)+(3/20)*(q0/2*L2))*sin(alfa);
F22tr=((3/20)*(-q0*L2)+(7/20)*(-q0/2*L2))*cos(alfa);  F32=((7/20)*(-q0*L2)+(3/20)*(-q0/2*L2))*cos(alfa);
F23tr=(-1/30)*(q0*L2^2)+(1/20)*(q0/2*L2^2);  F33=(1/20)*(q0*L2^2)+(1/30)*(q0/2*L2^2);
F21=F21t+F21tr;  F22=F22t+F22tr;  F23=F23t+F23tr;

q(1:9)=[F11;F12;F13;F21;F22;F23;F31;F32;F33]

%% MATRIZ DE RIGIDEZ
%Declarar y crear matriz de rigidez global vacia
K=zeros(GDL);
L3=L/2;L4=L/2;L5=L/2;L6=L/2;

%Elemento 1
K1=[EA/L1 0 0 -EA/L1 0 0; 0 12*EI/L1^3 6*EI/L1^2 0 -12*EI/L1^3 6*EI/L1^2;...
    0 6*EI/L1^2 4*EI/L1 0 -6*EI/L1^2 2*EI/L1; -EA/L1 0 0 EA/L1 0 0;...
    0 -12*EI/L1^3 -6*EI/L1^2 0 12*EI/L1^3 -6*EI/L1^2; 0 6*EI/L1^2 2*EI/L1 0 -6*EI/L1^2 4*EI/L1];

alfa1=alfa;
R1=[cos(alfa1) sin(alfa1) 0; -sin(alfa1) cos(alfa1) 0; 0 0 1];
ceros=zeros(3);
R1=[R1 ceros; ceros R1];
K1g=R1'*K1*R1;
Lg1=[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
 
K1g=Lg1*K1g*Lg1';

%Elemento 2
K2=[EA/L2 0 0 -EA/L2 0 0; 0 12*EI/L2^3 6*EI/L2^2 0 -12*EI/L2^3 6*EI/L2^2;...
    0 6*EI/L2^2 4*EI/L2 0 -6*EI/L2^2 2*EI/L2; -EA/L2 0 0 EA/L2 0 0;...
    0 -12*EI/L2^3 -6*EI/L2^2 0 12*EI/L2^3 -6*EI/L2^2;0 6*EI/L2^2 2*EI/L2 0 -6*EI/L2^2 4*EI/L2];

alfa2=alfa;
R2=[cos(alfa2) sin(alfa2) 0; -sin(alfa2) cos(alfa2) 0; 0 0 1];
ceros=zeros(3);
R2=[R2 ceros; ceros R2];
K2g=R2'*K2*R2;
Lg2=[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0;...
     0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
 
K2g=Lg2*K2g*Lg2';

%Elemento 3 
K3=[EA/L3 0 0 -EA/L3 0 0; 0 12*EI/L3^3 6*EI/L3^2 0 -12*EI/L3^3 6*EI/L3^2;...
    0 6*EI/L3^2 4*EI/L3 0 -6*EI/L3^2 2*EI/L3; -EA/L3 0 0 EA/L3 0 0;...
    0 -12*EI/L3^3 -6*EI/L3^2 0 12*EI/L3^3 -6*EI/L3^2; 0 6*EI/L3^2 2*EI/L3 0 -6*EI/L3^2 4*EI/L3];

alfa3=0;
R3=[cos(alfa3) sin(alfa3) 0; -sin(alfa3) cos(alfa3) 0; 0 0 1];
ceros=zeros(3);
R3=[R3 ceros; ceros R3];
K3g=R3'*K3*R3;
Lg3=[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
 
K3g=Lg3*K3g*Lg3';

%Elemento 4
K4=[EA/L4 0 0 -EA/L4 0 0; 0 12*EI/L4^3 6*EI/L4^2 0 -12*EI/L4^3 6*EI/L4^2;...
    0 6*EI/L4^2 4*EI/L4 0 -6*EI/L4^2 2*EI/L4; -EA/L4 0 0 EA/L4 0 0;...
    0 -12*EI/L4^3 -6*EI/L4^2 0 12*EI/L4^3 -6*EI/L4^2; 0 6*EI/L4^2 2*EI/L4 0 -6*EI/L4^2 4*EI/L4];

alfa4=0;
R4=[cos(alfa4) sin(alfa4) 0; -sin(alfa4) cos(alfa4) 0; 0 0 1];
ceros=zeros(3);
R4=[R4 ceros; ceros R4];
K4g=R4'*K4*R4;
Lg4=[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0;...
     0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
 
K4g=Lg4*K4g*Lg4';

%Elemento 5
K5=[EA/L5 0 0 -EA/L5 0 0; 0 12*EI/L5^3 6*EI/L5^2 0 -12*EI/L5^3 6*EI/L5^2;...
    0 6*EI/L5^2 4*EI/L5 0 -6*EI/L5^2 2*EI/L5; -EA/L5 0 0 EA/L5 0 0;...
    0 -12*EI/L5^3 -6*EI/L5^2 0 12*EI/L5^3 -6*EI/L5^2; 0 6*EI/L5^2 2*EI/L5 0 -6*EI/L5^2 4*EI/L5];
alfa5=-pi/2;
R5=[cos(alfa5) sin(alfa5) 0; -sin(alfa5) cos(alfa5) 0; 0 0 1];
ceros=zeros(3);
R5=[R5 ceros; ceros R5];
K5g=R5'*K5*R5;
Lg5=[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
 
K5g=Lg5*K5g*Lg5';

%Elemento 6
K6=[EA/L6 0 0 -EA/L6 0 0; 0 12*EI/L6^3 6*EI/L6^2 0 -12*EI/L6^3 6*EI/L6^2;...
    0 6*EI/L6^2 4*EI/L6 0 -6*EI/L6^2 2*EI/L6; -EA/L6 0 0 EA/L6 0 0;...
    0 -12*EI/L6^3 -6*EI/L6^2 0 12*EI/L6^3 -6*EI/L6^2; 0 6*EI/L6^2 2*EI/L6 0 -6*EI/L6^2 4*EI/L6];
alfa6=-pi/2;
R6=[cos(alfa6) sin(alfa6) 0; -sin(alfa6) cos(alfa6) 0; 0 0 1];
ceros=zeros(3);
R6=[R6 ceros; ceros R6];
K6g=R6'*K6*R6;
Lg6=[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;...
     0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0;...
     0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
 
K6g=Lg6*K6g*Lg6'; 

K=K1g+K2g+K3g+K4g+K5g+K6g

%% MATRIZ DE RIGIDEZ, METODO DE ENSAMBLAJE
% K=zeros(GDL);

% %Ensamblaje del elemento 1
% K(1:3,1:3)=K1_(1:3,1:3);
% K(4:6,4:6)=K1_(4:6,4:6);
% K(1:3,4:6)=K1_(1:3,4:6);
% K(4:6,1:3)=K1_(4:6,1:3);

% %Ensamblaje del elemento 2
% K(4:6,4:6)=K(4:6,4:6)+K2_(1:3,1:3);
% K(7:9,7:9)=K2_(4:6,4:6);
% K(4:6,7:9)=K2_(1:3,4:6);
% K(7:9,4:6)=K2_(4:6,1:3);

% %Ensamblaje del elemento 3
% K(7:9,7:9)=K(7:9,7:9)+K3_(1:3,1:3);
% K(10:12,10:12)=K3_(4:6,4:6);
% K(7:9,10:12)=K3_(1:3,4:6);
% K(10:12,7:9)=K3_(4:6,1:3);

% %Ensamblaje del elemento 4
% K(10:12,10:12)=K(10:12,10:12)+K4_(1:3,1:3);
% K(13:15,13:15)=K4_(4:6,4:6);
% K(10:12,13:15)=K4_(1:3,4:6);
% K(13:15,10:12)=K4_(4:6,1:3);

% %Ensamblaje del elemento 5
% K(13:15,13:15)=K(13:15,13:15)+K5_(1:3,1:3);
% K(16:18,16:18)=K5_(4:6,4:6);
% K(13:15,16:18)=K5_(1:3,4:6);
% K(16:18,13:15)=K5_(4:6,1:3);

% %Ensamblaje del elemento 6
% K(16:18,16:18)= K(16:18,16:18) +K6_(1:3,1:3);
% K(19:21,19:21)=K6_(4:6,4:6);
% K(16:18,19:21)=K6_(1:3,4:6);
% K(19:21,16:18)=K6_(4:6,1:3) 

% %Otra forma de ensamblar la matriz K, se puede hacer asi porque los elementos van encadenados en linea
% %K(1:6,1:6)=K1_;
% %K(4:9,4:9)=K(4:9,4:9)+K2_;
% %K(7:12,7:12)=K(7:12,7:12)+K3_;
% %K(10:15,10:15)=K(10:15,10:15)+K4_;
% %K(13:18,13:18)=K(13:18,13:18)+K5_;
% %K(16:21,16:21)=K(16:21,16:21)+K6_

%% COND DE CONTORNO (CC)
%Crear vector de desplazam vacio
U=zeros(GDL,1);
%Nodos de apoyo de la estructura
u_contorno=[1:GDLn,(nnodos*GDLn-GDLn+1):nnodos*GDLn];
%Nodos libres, eliminando nodos de contorno
u_libres=setdiff((1:length(U))',u_contorno);
%Crear vector de fzas nodales vacio
f=zeros(GDL,1);
%Definir vector de fzas nodales exteriores
L=f(u_libres);

%% DESPLAZAM DE NODOS CENTRALES 2, 4, 6
%Matriz reducida por las CC
KLL=K(u_libres,u_libres);

%Vector de fzas nodales generalizadas: P=f+q
FL=f(u_libres)+q(u_libres);

%Calculo del vector de desplazam en nodos libres
u_lib=U(4:18);
u_lib=KLL\FL;

%Vector de desplazam completo
U=[0; 0; 0; u_lib; 0; 0; 0];

%Vector de desplazam del nodo 2
u_2=U(4:6);

%Vector de desplazam del nodo 4
u_4=U(10:12);

%Vector de desplazam del nodo 6
u_6=U(16:18);

%Vector inicial vacio
u_quest=zeros(GDLn*3,1);

u_quest(1:3,1)=u_2;
u_quest(4:6,1)=u_4;
u_quest(7:9,1)=u_6

%% VALOR DE LA REACCION EN EL APOYO IZQUIERDO
f=q;
f=K*U;

%Vector inicial vacío
r_quest=zeros(GDLn*3,1);

r_quest=f(1:3,1)-q(1:3,1)

%FIN