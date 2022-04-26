clc
clear 

c0=3e8;                %light speed
mu0=4*pi*10^(-7);      %the permeability of vacuum
eps0=1/(c0^2*mu0);     %the permitimeivity of cacuum

MAX_Z=200;
step=800;

Ex(1:MAX_Z)=0;
Hy(1:MAX_Z)=0;
CFNL=1;
dz=0.03;
dt=CFNL/(c0*sqrt(1/dz^2+1/dz^2+1/dz^2));
CA=1;
CB=(dt/2/eps0);
CP=1;
CQ=dt/2/mu0;

a=-dt*dt/(4*mu0*eps0*dz*dz);
b=1-2*a;
c=-dt*dt/(4*mu0*eps0*dz*dz);
% A(1,1)=0;
for m=1:MAX_Z-3
    A(m+1,m) = a;
    A(m,m) = b;
    A(m,m+1) = c;
end

A(m+1,m+1) = b;
% d(1)=0;
% d(MAX_Z)=0;
JA=sparse(A);
A_inv=inv(A);
clear A

Tw=1.0618e-9;
t0=4*Tw;
for n=1:step
    n
    pulse=exp(-(n*dt-t0)^2/(Tw)^2);
    Ex(1,50)=Ex(1,50)+pulse;  
    
    
    Ex1=Ex;Hy1=Hy;
    k=2:MAX_Z-1;
    d(k-1)=Ex1(1,k)-CB*(Hy1(1,k)-Hy1(1,k-1))/dz;
    Ex(1,2:MAX_Z-1)=JA\d';
%     
%     pulse=exp(-(n*dt-t0)^2/(Tw)^2);
%     Ex(1,50)=pulse;  
    
    Hy(1,1:MAX_Z-1)=Hy1(1,1:MAX_Z-1)-CQ*(Ex(1,(1:MAX_Z-1)+1)-Ex(1,1:MAX_Z-1))/dz;
    
    
    Ex1=Ex;Hy1=Hy;
    Ex(1,2:MAX_Z-1)=Ex1(1,2:MAX_Z-1)-CB*(Hy1(1,2:MAX_Z-1)-Hy1(1,(2:MAX_Z-1)-1))/dz;
    
    Hy(1,1:MAX_Z-1)=Hy1(1,1:MAX_Z-1)-CQ*(Ex1(1,(1:MAX_Z-1)+1)-Ex1(1,1:MAX_Z-1))/dz;
    
    figure(1)
    plot(Ex)
%     axis([1 MAX_Z,-0.05 0.05]);
    AAA(n)=Ex(1,100);
end
