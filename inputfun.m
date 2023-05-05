function [dx_final,dy_final,M,N,P,tn]=inputfun(D1A,D1B,D3A,D3B,DxA,DxB,C1,C2,tfinal)
%Calculate parameters based on inputs
% L, horizontal length (m)

global L dt T

% %formula: x=Ft=> dx=-Dx*grad(C)*dt=-Dx*((C2-C1)/L)*dt, 
% where F is veclocity flux
Dxavg=(DxA+DxB)/2;
dx=-Dxavg*((C2-C1)/L)*dt;

%formula: dx/dy=sqrt(k1/k3)
D1avg=(D1A+D1B)/2;
D3avg=(D3A+D3B)/2;
dz=dx/(sqrt(D1avg/D3avg));

%formula: dx=L/(N-1), dy=T/(M-1);
n=round(L/dx)+1;%number of nodes in x-direction
m=round(T/dz)+1;%%number of nodes in z-direction

% to get symmetry about layers no. of nodes need to be even
if mod(n,2)==0
    N=n;
else
    N=n+1;
end

if mod(m,2)==0
    M=m;
else
    M=m+1;
end
P=M*N;%total number of nodes
dx_final=L/(N-1);%recalculate dx to consistent with N
dy_final=T/(M-1);%recalculate dy to conistent with  M
tn=round(tfinal/dt);%total number of time steps
end