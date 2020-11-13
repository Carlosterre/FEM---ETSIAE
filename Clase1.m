AE=1; L=1;

%Matriz de rigidez local del elemento 1
L1=L*sqrt(2);  %Por trigonometria se saca la longitud del elemento 1
K1=AE/L1*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];  %Declaracion de matriz de rigidez de la barra en ejes locales 2D
%Matriz de cambio de ejes
alfa1=-pi/4;
R1=[cos(alfa1) sin(alfa1); -sin(alfa1) cos(alfa1)];  %Matriz de giro
Zeros=[0 0; 0 0];  %zeros(2): Esta es otra forma de declarar la matriz de 2X2 de ceros
R1_=[R1 Zeros; Zeros R1];  %Se construye esta matriz 4X4 poniendo en la diagonal la matriz de giro y el resto ceros

K1_= R1_*K1*R1_';  %El apostrofe '  da la matriz traspuesta

%Se repite el proceso para los elementos 2 y 3 cambiando las longitudes y los angulos
%Matriz de rigidez local del elemento 2
L2=L;
K2=AE/L2*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];

%Matriz de cambio de ejes
alfa2=0;
R2=[cos(alfa2) sin(alfa2); -sin(alfa2) cos(alfa2)];
Zeros=[0 0; 0 0];  %zeros(2)
R2_=[R2 Zeros; Zeros R2];

K2_=R2_*K2*R2_';

%Matriz de rigidez local del elemento 3
L3=L*sqrt(5)/2;
K3=AE/L3*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];

%Matriz de cambio de ejes
alfa3=atan(1/2);
R3=[cos(alfa3) sin(alfa3); -sin(alfa3) cos(alfa3)];
Zeros=[0 0; 0 0];  %zeros(2)
R3_=[R3 Zeros; Zeros R3];

K3_=R3_*K3*R3_';

%% Ensamblaje de la matriz de rigidez. Este metodo no se va a a utilizar, solo se muestra para ver como seria
% L_1=[1 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0;...  3 puntos + Enter para bajar de linea
%      0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1]
% K_1_=L_1'*K1_*L_1
 
K=zeros(8);  % Se define una matriz K 8X8 la cual se llena inicialmente de ceros y mediante el ensamblaje
             % de cada elemento se completa con los numeros correspondientes
               
%Ensamblaje del elemento 1
K(1:2,1:2)=K1_(1:2,1:2);  %Fila 1 a 2, columna 1 a 2 de K se rellenan con fila 1 a 2, columna 1 a 2 de K1
K(7:8,7:8)=K1_(3:4,3:4);
K(1:2,7:8)=K1_(1:2,3:4);
K(7:8,1:2)=K1_(3:4,1:2);

%Ensamblaje del elemento 2
K(3:4,3:4)=K(3:4,3:4)+K2_(1:2,1:2);  %Hay que hacer esta suma de elementos para no perder los que ya se hayan colocado
K(7:8,7:8)=K(7:8,7:8)+K2_(3:4,3:4);
K(3:4,7:8)=K(3:4,7:8)+K2_(1:2,3:4);
K(7:8,3:4)=K(7:8,3:4)+K2_(3:4,1:2);

%Ensamblaje del elemento 3
K(5:6,5:6)=K(5:6,5:6)+K3_(1:2,1:2);
K(7:8,7:8)=K(7:8,7:8)+K3_(3:4,3:4);
K(5:6,7:8)=K(5:6,7:8)+K3_(1:2,3:4);
K(7:8,5:6)=K(7:8,5:6)+K3_(3:4,1:2);

%% Calculo de desplazamientos
fL= [1; -1]*100/(sqrt(2));  %Vector de fuerzas aplicadas
KLL=K(7:8,7:8);  %Matriz de K reducida eliminando filas y columnas cuyos desplazamientos son nulos
uL=KLL\fL;  %El comando \ da la matriz inversa de la matriz que este delante
u=zeros(8,1);  % Se construye el vector de desplazamientos de 8X1 y se rellena con los desplazamientos
               % calculados en las posicioens 7,8 que son las unicas que se desplazan
u(7:8)=uL

%% Calculo de reacciones
f=K*u;  %Vector de fuerzas totales: Reacciones + Fuerzas aplicadas
KRL=K(7:8,1:6);
fR=KRL'*uL  %Vector unicamente de reacciones
