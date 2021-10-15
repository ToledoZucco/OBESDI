%This function create a time-discretized model of an LTI (A,B,C,D)
%using a mid-point method

function [Ad,Bd,Cd,Dd] = MidPointTimeDiscretization(A,B,C,D,dt)


I = eye(length(A));
Ad = inv(I-dt/2*A)*(I+dt/2*A);
Bd = inv(I-dt/2*A)*dt*B;
Cd = C;
Dd = D;

end