%Eigenvalues analisys closed-loop system

close all
% clear all
clc


%% Infinite-Dimensional Analisys

P1 = [0,1;1,0];
I2 = eye(2);
R = 1/sqrt(2)*[P1,-P1;I2,I2];
W = 1/sqrt(2)*[-1,0,0,1;0,1,1,0];
Wt = 1/sqrt(2)*[0,1,-1,0;1,0,0,1];











%% Finite-dimensonal Analisys
    
N = 100;      % Total of tate variables
long = 1;   % length of the string
rho = 1;    % mass densitu
T = 1;      % Young's modulus
Dis = 0;    % Dissipation along the string
%Create a model (ABCD) in which the state x = [qd;pd] where qd is the
%discretized strain and pd is the discretized momentum
[A,B,C,D,Q,h,np,nq] = VibratingStringModelWithoutVi(N,long,rho,T,Dis);

% With this matrix we can obtain the displacement of the string
Cw = [h*tril(ones(N/2,N/2)),zeros(N/2,N/2)];

Ce = h*[ones(1,N/2),zeros(1,N/2)];
Cp = h*[zeros(1,N/2),ones(1,N/2)];


dc1 = 1000*0;
dc2 = 1000;
qc1 = 10*0;
qc2 = 10;

Qc = diag([qc1,qc2]);
Dc = diag([dc1,dc2]);

K = inv(eye(2)+Qc*h/2*diag([rho,T]))*(Dc*C+Qc*[Cp;Ce]);
ALC = A-B*K

eig(ALC)