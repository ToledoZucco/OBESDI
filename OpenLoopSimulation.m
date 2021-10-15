%This code shows the open loop simulation

close all
% clear all
clc
    
N = 200;      % Total of tate variables
long = 1;   % length of the string
rho = 1;    % mass densitu
T = 1;      % Young's modulus
Dis = 0.01;    % Dissipation along the string
%Create a model (ABCD) in which the state x = [qd;pd] where qd is the
%discretized strain and pd is the discretized momentum
[A,B,C,D,Q,h,np,nq] = VibratingStringModel(N,long,rho,T,Dis);

% With this matrix we can obtain the displacement of the string
Cw = [h*tril(ones(N/2,N/2)),zeros(N/2,N/2),ones(N/2,1)];


t0 = 0;
dt = 1e-3;
t = t0:dt:10000*dt;
Nt = length(t);
[Ad,Bd,Cd,Dd] = MidPointTimeDiscretization(A,B,C,D,dt);


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
for k = 1:Nt
    
    y(:,k) = Cd*z(:,k) + Dd*u(:,k);
    z(:,k+1) = Ad*z(:,k) + Bd*u(:,k);
end
z = z(:,1:end-1);

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



figureDirectory = 'OpenLoopSimulation';
filename = 'Sim';
LegendSim1 = '$w(t,\zeta)$';
saveFigures = false;
N_vid = 100;
dk_vid = (Nt-1)/N_vid;
k_vid = 1:dk_vid:Nt-1;
zeta = linspace(0,long,N/2);

Save_Figures_1plot(t,np,w,k_vid,figureDirectory,filename,saveFigures,LegendSim1)


