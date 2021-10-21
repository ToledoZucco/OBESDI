%This function create the model of an attached string with a force actuator
%at the oposite side.
%N: Number of state variables to describe the system (even number)
%long: Length of the string
%rho_c: mass density (constant)
%T_c: Youngs modulus (constant)
%Dis: Dissipation
%% Model
%dx = Ax+Bu
% y = Cx+Du
% u = [u1;u2]. u1 is the velocity at a and u2 is the force at b
% y = [y1;y2]. y1 is the force at a and y2 is the velocity at b
% x = [qd;pd;vi] with qd the discretized strain, pd the discretized
% momentum and vi = int_0^t [v(t,0)] dt
function [A,B,C,D,Q,h,zp,zq] = VibratingStringModel(N,long,rho_c,T_c,Dis)

np = N/2;
a = 0;
b = long;

Zn = zeros(np);
nq = np;   %Number of Strain Variables
n = np+nq;
d = long/(n+1);
h = 2*d;

zq = (a+d):h:(b-h); zq = zq';
zp = (a+h):h:(b-d); zp = zp';
% zw = zp;

T = T_c*zq./zq;
rho = rho_c*zp./zp;



Qq = h*diag(T);
Qp = h*diag(1./rho);

D1 =  -eye(nq);
D2 = [zeros(nq,1),[eye(nq-1);zeros(1,nq-1)]];
D = 1/h^2*(D1+D2);
D = -D';

Disn = 1/h*Dis*eye(np);

b1 = 1/h*[-1;zeros(nq-1,1)];
b2 = 1/h*[zeros(np-1,1);1];
J = [zeros(np),D;-D',zeros(nq)];
R = [Zn,Zn;Zn,Disn];
Q = [Qq,Zn;Zn,Qp];
B = [[b1,zeros(nq,1)];[zeros(nq,1),b2]];

A = (J-R)*Q;
C = B'*Q;
D = 0;


%Include one state vi = int_0^t [v(t,0)] dt
%This new state is for rebuilding w(t,zeta)
bi = [1,0];
A = [A,zeros(n,1);zeros(1,N),0];
B = [B;bi];
C = [C,zeros(2,1)];
D = D*1;

end