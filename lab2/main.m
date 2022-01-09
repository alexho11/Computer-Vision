clear all;
close all;
clc;
% Kameramatrix und aeussere Orientierungen
K=[-3933.36363636	     -0.00000000	   2143.50000000;
	     -0.00000000	   3933.36363636	   1423.50000000;
	      0.00000000	      0.00000000	      1.00000000];
      
     
R1=[ -0.867233884963	 -0.496431897214	 -0.038219892721;
	  0.497756374147	 -0.862573995849	 -0.090579764190;
	  0.011999198600	 -0.097578036019	  0.995155538657];
  X1=[  513010.498970	;  5427654.064230	  ;    514.392000];
  
  R2=[-0.839539757811	 -0.528716554759	 -0.125027196158;
	  0.538291662815	 -0.840647830782	 -0.059609649767;
	 -0.073587232584	 -0.117345768246	  0.990360989678];
 X2=[512996.502990;	  5427678.079250	  ;    513.677870]; 

%  figure;
%  img=imread('R0020849.jpg');
%  imshow(img);
%  hold on;
%  [xi,yi] = getpts;
% xf=[xi';yi';ones(1,3)];

load('xf.mat');
m=3; % Anzahl der Punkte

b=-R1*(X2-X1); % Basis
R12=R1*R2';
B=[0,-b(3),b(2);b(3),0,-b(1);-b(2),b(1),0];
E=B*R12; % Essentialmatrix
E=E';
F=inv(K)'*E*inv(K); % Fundamentalmatrix

l=F*xf; % Epipolarlinie

x=0:2848*2; % X-Koordinaten

y_ep=(-l(3,:)'-l(1,:)'*x)./l(2,:)'; % Y-Koordinaten

% Linien ploten
figure;
img=imread('R0020850.jpg');
imshow(img)
hold on;
plot(x,y_ep,'r','LineWidth',5)

%% Alternative Berechnung aus Projektionsmatrizen
P1=K*[R1,-R1*X1];
P2=K*[R2,-R2*X2];

Ps=P1'*inv(P1*P1');% pseudo-inverse
P2X1=P2*[X1;1];
P2X1_x=[0,-P2X1(3),P2X1(2);
    P2X1(3),0,-P2X1(1);
    -P2X1(2),P2X1(1),0];

F2=P2X1_x*P2*Ps;
l2=F2*xf;
y_ep2=(-l2(3,:)'-l2(1,:)'*x)./l2(2,:)';
figure;
img=imread('R0020850.jpg');
imshow(img)
hold on;
plot(x,y_ep2,'r','LineWidth',5)
