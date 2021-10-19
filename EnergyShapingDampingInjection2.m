%This code shows the Energy shaping and daming injection Simulation

close all
% clear all
clc
    
N = 100;      % Total of tate variables
long = 1;   % length of the string
rho = 1;    % mass densitu
T = 1;      % Young's modulus
Dis = 0;    % Dissipation along the string
%Create a model (ABCD) in which the state x = [qd;pd] where qd is the
%discretized strain and pd is the discretized momentum
[A,B,C,D,Q,h,np,nq] = VibratingStringModel(N,long,rho,T,Dis);

% With this matrix we can obtain the displacement of the string
Cw = [h*tril(ones(N/2,N/2)),zeros(N/2,N/2),ones(N/2,1)];

C_int_e = h*[ones(1,N/2),zeros(1,N/2),0];
C_int_p = h*[zeros(1,N/2),ones(1,N/2),0];


dc1 = 100;
dc2 = 100;
qc1 = 100;
qc2 = 100;

Qc = diag([qc1,qc2]);
Dc = diag([dc1,dc2]);

K = Dc*C+Qc*[C_int_p;C_int_e];
ALC = A-B*K

eig(ALC)

%%

t0 = 0;
dt = 1e-5;
t = t0:dt:400000*dt;
Nt = length(t);
[Ad,Bd,Cd,Dd] = MidPointTimeDiscretization(A,B,C,D,dt);

%% Damping Injection

k1 = 5;
k2 = 5;
K = diag([k1,k2]);
% K = inv(eye(2)+K*Dd)*K*Cd;
K = K*Cd;

%% Energy shaping

Sig1 = 100;
Sig2 = 100;
Sig = diag([Sig1,Sig2]);

%% Controller
Jc = [0,1;-1,0];
Qc = diag([Sig1,Sig2]);
Dc = diag([k1,k2]);
Ac = Jc*Qc;
Bc = eye(2);
Cc = Bc'*Qc;


[Acd,Bcd,Ccd,Dcd] = MidPointTimeDiscretization(Ac,Bc,Cc,Dc,dt);


%Initial condition
w00 = 0;
q0 = 1*ones(N/2,1);
p0 = zeros(N/2,1);
vi0 = 0;
z0 = [q0;p0;vi0];

z = zeros(N+1,Nt);
z(:,1) = z0;
u = zeros(2,Nt);
y = zeros(2,Nt);

xc0 = [0;0];
xc = zeros(2,Nt);
xc(:,1) = xc0;
uc = zeros(2,Nt);
yc = zeros(2,Nt);
for k = 1:Nt
    
    y(:,k) = Cd*z(:,k) + Dd*u(:,k);
    uc(:,k) = y(:,k);
    
    yc(:,k) = Ccd*xc(:,k) + Dcd*uc(:,k);
    
    u(:,k) = -yc(:,k);
    
    
    z(:,k+1) = Ad*z(:,k) + Bd*u(:,k);
    xc(:,k+1) = Acd*xc(:,k) + Bcd*uc(:,k);
end
z = z(:,1:end-1);
xc = xc(:,1:end-1);

w = Cw*z+w00;

%% Figures
x0screen=100;y0screen=50;width=1000;height=600;font=35;lw=4;ms = 15;

%Output and input
figure
subplot(2,1,1)
hold on
plot(t,u,'LineWidth',lw)
legend({'$u(t)$'},'Interpreter','latex','FontSize',font)
grid on
set(gca,'FontSize',font);

subplot(2,1,2)
hold on
plot(t,y,'LineWidth',lw)
legend({'$y(t)$'},'Interpreter','latex','FontSize',font)
grid on
set(gca,'FontSize',font);



% End-tip position
figure
hold on
plot(t,w(end,:),'LineWidth',lw)
plot(t,zeros(1,Nt),'--','LineWidth',lw)
legend({'$w(b,t)$','$w_{desired}(b,t)$'},'Interpreter','latex','FontSize',font)
grid on
set(gca,'FontSize',font);



figureDirectory = 'Damp';
filename = 'Sim';
LegendSim1 = '$w(t,\zeta)$';
saveFigures = false;
N_vid = 100;
dk_vid = (Nt-1)/N_vid;
k_vid = 1:dk_vid:Nt-1;
zeta = linspace(0,long,N/2);

Save_Figures_1plot(t,np,w,k_vid,figureDirectory,filename,saveFigures,LegendSim1)