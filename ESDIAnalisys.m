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
[A,B,C,D,Q,h,np,nq] = VibratingStringModel(N,long,rho,T,Dis);

A = A(1:end-1,1:end-1);
B = B(1:end-1,:);
C = C(:,1:end-1);



% With this matrix we can obtain the displacement of the string
Cw = [h*tril(ones(N/2,N/2)),zeros(N/2,N/2)];

Ce = h*[ones(1,N/2),zeros(1,N/2)];
Cp = h*[zeros(1,N/2),ones(1,N/2)];
Cpe = [Cp;Ce];


gamma1 = 10;
gamma2 = 10;
gamma = diag([gamma2,gamma1]);


Dc = h/2*diag([rho,1/T]);

Cea = [1,zeros(1,N-1)];
Cpb = [zeros(1,N-1),1];

theta1 = 1;
theta2 = 10;
Theta = diag([theta1*(T+gamma1*long),-theta2*(1/rho+gamma2*long)]);

K = inv(I2 + gamma*Dc)*(gamma*Cpe-Theta*[Cea;Cpb])

ACL = A-B*K;
eig(ACL)



