clc;
clear;

mu0=pi*4e-7; %permability of free space
eps0=8.854E-12; %permittivity of free space
c0=1/sqrt(eps0*mu0); %velocity of light
% 
% Lx=0.5;
% Ly=0.5;
% Lz=0.5;
%设置网格大小
nx1=1;
ny1=nx1;
nz1=nx1;
nx2=30;
ny2=35;
nz2=40;

ic=nx2/2;
jc=ny2/2;
kc=nz2/2;
% f_0=8e9;
% tw=2/(pi*f_0);
% t0=4*tw;


% dx=Lx/(nx2-1);
% dy=Ly/(ny2-1);
% dz=Lz/(nz2-1);
dx=0.05;
dy=0.05;
dz=0.05;
dt=3/(c0*sqrt(1/dx^2+1/dy^2+1/dz^2));  %time step

Tw=1.0618e-9;
t0=4*Tw;

Hx(nx1:nx2+1,ny1:ny2,nz1:nz2)=0;
Hy(nx1:nx2,ny1:ny2+1,nz1:nz2)=0;
Hz(nx1:nx2,ny1:ny2,nz1:nz2+1)=0;
Ex(nx1:nx2,ny1:ny2+1,nz1:nz2+1)=0;
Ey(nx1:nx2+1,ny1:ny2,nz1:nz2+1)=0;
Ez(nx1:nx2+1,ny1:ny2+1,nz1:nz2)=0;

CA=1;
CB=(dt/2/eps0);
CP=1;
CQ=dt/2/mu0;
mur_x=(dx-c0*dt)/(dx+c0*dt);
mur_y=(dy-c0*dt)/(dy+c0*dt);
mur_z=(dz-c0*dt)/(dz+c0*dt);
% x方向
a1=-dt*dt/(4*mu0*eps0*dx*dx);
b1=1-2*a1;
c1=a1;

for m=nx1:nx2-2
A(m+1,m) = a1;
A(m,m) = b1;
A(m,m+1) = c1;
end
A(1,1)=b1-a1*mur_x;
A(m+1,m+1)=b1-c1*mur_x;
JAx=sparse(A);
clear A;
a_x=a1;   c_x=c1;

%y方向
a1=-dt*dt/(4*mu0*eps0*dy*dy);
b1=1-2*a1;
c1=a1;

for m=ny1:ny2-2
A(m+1,m) = a1;
A(m,m) = b1;
A(m,m+1) = c1;
end
A(1,1)=b1-a1*mur_y;
A(m+1,m+1)=b1-a1*mur_y;
JAy=sparse(A);
clear A;
a_y=a1;   c_y=c1;

%z方向
a1=-dt*dt/(4*mu0*eps0*dy*dy);
b1=1-2*a1;
c1=a1;
for m=nz1:nz2-2
A(m+1,m) = a1;
A(m,m) = b1;
A(m,m+1) = c1;
end
A(1,1)=b1-a1*mur_z;
A(m+1,m+1)=b1-a1*mur_z;
JAz=sparse(A);
clear A;
a_z=a1;   c_z=c1;

for n=1:70
    n
%激励
pulse=exp(-(n*dt-t0)^2/(Tw)^2);
Ez(10,10,10)=pulse;    
%第一步
Ex1=Ex; Ey1=Ey; Ez1=Ez; 
Hx1=Hx; Hy1=Hy; Hz1=Hz;
%----------EZ---------------
     
for j=ny1+1:ny2
    for k=nz1:nz2     
        i=nx1+1:nx2;
        d_z1(i-1)=CA*Ez1(i,j,k)-CQ*CB/dx/dz*(Ex1(i,j,k+1)-Ex1(i,j,k)-Ex1(i-1,j,k+1)+Ex1(i-1,j,k))...
            -CB/dy*(Hx1(i,j,k)-Hx1(i,j-1,k))+CB/dx*(Hy1(i,j,k)-Hy1(i-1,j,k));
        d_z1(1)=d_z1(1)-a_x*dx/(dx+c0*dt)*(Ez1(1,j,k)+Ez1(2,j,k));
        d_z1(nx2-1)=d_z1(nx2-1)-c_x*dx/(dx+c0*dt)*(Ez1(nx2+1,j,k)+Ez1(nx2,j,k));
        Ez(i,j,k)=JAx\d_z1';
    end
end
% k=nx1:nx2-1;
% i=ny1+1:ny2; 
Ez(1,:,:)=-mur_x*Ez(2,:,:)+dx/(dx+c0*dt)*(Ez1(1,:,:)+Ez1(2,:,:));
Ez(nx2+1,:,:)=-mur_x*Ez(nx2,:,:)+dx/(dx+c0*dt)*(Ez1(nx2+1,:,:)+Ez1(nx2,:,:));

Ez(:,1,:)=-mur_y*Ez(:,2,:)+dy/(dy+c0*dt)*(Ez1(:,1,:)+Ez1(:,2,:));
Ez(:,ny2+1,:)=-mur_y*Ez(:,ny2,:)+dy/(dy+c0*dt)*(Ez1(:,ny2+1,:)+Ez1(:,ny2,:));

%--------HY----------------

    i=nx1:nx2;
    j=ny1:ny2+1;
    k=nz1:nz2;
    Hy(i,j,k)=Hy1(i,j,k)+CQ/dx*(Ez(i+1,j,k)-Ez(i,j,k))-CQ/dz*(Ex1(i,j,k+1)-Ex1(i,j,k));
    
%--------EX------------------

for i=nx1:nx2
    for k=nz1+1:nz2
        j= ny1+1:ny2; 
        d_x1(j-1)=CA*Ex1(i,j,k)-CQ*CB/dx/dy*(Ey1(i+1,j,k)-Ey1(i,j,k)-Ey1(i+1,j-1,k)+Ey1(i,j-1,k))...
            -CB/dz*(Hy1(i,j,k)-Hy1(i,j,k-1))+CB/dy*(Hz1(i,j,k)-Hz1(i,j-1,k));  
        d_x1(1)=d_x1(1)- a_y *  dy/(dy+c0*dt) *( Ex1(i,1,k) + Ex1(i,2,k) );
        d_x1(ny2-1)=d_x1(ny2-1)- c_y *  dy/(dy+c0*dt) *(Ex1(i,ny2+1,k) + Ex1(i,ny2,k) );
        
        Ex(i,j,k)=JAy\d_x1';     
    end 
end
% i=nx1:nx2-1;
% j=ny1+1:ny2; 
Ex(:,1,:)=-mur_y*Ex(:,2,:)+dy/(dy+c0*dt)*(Ex1(:,1,:)+Ex1(:,2,:));
Ex(:,ny2+1,:)=-mur_y*Ex(:,ny2,:)+dy/(dy+c0*dt)*(Ex1(:,ny2+1,:)+Ex1(:,ny2,:));

Ex(:,:,1)=-mur_z*Ex(:,:,2)+dz/(dz+c0*dt)*(Ex1(:,:,1)+Ex1(:,:,2));
Ex(:,:,nz2+1)=-mur_z*Ex(:,:,nz2)+dz/(dz+c0*dt)*(Ex1(:,:,nz2+1)+Ex1(:,:,nz2));

%--------HZ-----------------

    i=nx1:nx2;
    j=ny1:ny2;
    k=nz1:nz2+1;
    Hz(i,j,k)=Hz1(i,j,k)+CQ/dy*(Ex(i,j+1,k)-Ex(i,j,k))-CQ/dx*(Ey1(i+1,j,k)-Ey1(i,j,k));
    
%---------EY----------------  

for i=nx1+1:nx2
    for j=ny1:ny2
        k=nz1+1:nz2;
        d_y1(k-1)=CA*Ey1(i,j,k)-CQ*CB/dy/dz*(Ez1(i,j+1,k)-Ez1(i,j,k)-Ez1(i,j+1,k-1)+Ez1(i,j,k-1))...
            -CB/dx*(Hz1(i,j,k)-Hz1(i-1,j,k))+CB/dz*(Hx1(i,j,k)-Hx1(i,j,k-1));
        d_y1(1)=d_y1(1)-a_z*+dz/(dz+c0*dt)*(Ey1(i,j,1)+ Ey1(i,j,2));
        d_y1(nz2-1)=d_y1(nz2-1)-c_z*+dz/(dz+c0*dt)*(Ey1(i,j,nz2+1)+ Ey1(i,j,nz2));
        Ey(i,j,k)=JAz\d_y1';        
    end
end
% k=nz1+1:nz2;
% j=ny1:ny2-1;
Ey(:,:,1)=-mur_z*Ey(:,:,2)+dz/(dz+c0*dt)*(Ey1(:,:,1)+Ey1(:,:,2));
Ey(:,:,nz2+1)=-mur_z*Ey(:,:,nz2)+dz/(dz+c0*dt)*(Ey1(:,:,nz2+1)+Ey1(:,:,nz2));

Ey(1,:,:)=-mur_x*Ey(2,:,:)+dx/(dx+c0*dt)*(Ey1(1,:,:)+Ey1(2,:,:));
Ey(nx2+1,:,:)=-mur_x*Ey(nx2,:,:)+dx/(dx+c0*dt)*(Ey1(nx2,:,:)+Ey1(nx2+1,:,:));
%----------HX----------------  

    i=nx1:nx2+1;
    j=ny1:ny2;
    k=nz1:nz2;
    Hx(i,j,k)=Hx1(i,j,k)+CQ/dz*(Ey(i,j,k+1)-Ey(i,j,k))-CQ/dy*(Ez1(i,j+1,k)-Ez1(i,j,k));
    
    
%第二步
Ex1=Ex; Ey1=Ey; Ez1=Ez; 
Hx1=Hx; Hy1=Hy; Hz1=Hz;
%--------Ey----------------   
for j=ny1:ny2
    for k=nz1+1:nz2
        i=nx1+1:nx2;
        d_y2(i-1)=Ey1(i,j,k)-CQ*CB/dx/dy*(Ex1(i,j+1,k)-Ex1(i,j,k)-Ex1(i-1,j+1,k)+Ex1(i-1,j,k))...
            -CB/dx*(Hz1(i,j,k)-Hz1(i-1,j,k))+CB/dz*(Hx1(i,j,k)-Hx1(i,j,k-1));
        d_y2(1)=d_y2(1)-a_x*dx/(dx+c0*dt)*(Ey1(1,j,k)+Ey1(2,j,k));
        d_y2(nx2-1)=d_y2(nx2-1)-c_x*dx/(dx+c0*dt)*(Ey1(nx2+1,j,k)+Ey1(nx2,j,k));
        Ey(i,j,k)=JAx\d_y2';
    end
end

Ey(:,:,1)=-mur_z*Ey(:,:,2)+dz/(dz+c0*dt)*(Ey1(:,:,1)+Ey1(:,:,2));
Ey(:,:,nz2+1)=-mur_z*Ey(:,:,nz2)+dz/(dz+c0*dt)*(Ey1(:,:,nz2)+Ey1(:,:,nz2+1));

Ey(1,:,:)=-mur_x*Ey1(2,:,:)+dx/(dx+c0*dt)*(Ey1(1,:,:)+Ey1(2,:,:));
Ey(nx2+1,:,:)=-mur_x*Ey1(nx2,:,:)+dx/(dx+c0*dt)*(Ey1(nx2,:,:)+Ey1(nx2+1,:,:));
%---------hz---------------  
i=nx1:nx2;
j=ny1:ny2;
k=nz1:nz2+1;
Hz(i,j,k)=Hz1(i,j,k)-CQ/dx*(Ey(i+1,j,k)-Ey(i,j,k))+CQ/dy*(Ex1(i,j+1,k)-Ex1(i,j,k));
%----------ez---------------  
for i=nx1+1:nx2
    for k=nz1:nz2
        j=ny1+1:ny2;
        d_z2(j-1)=Ez1(i,j,k)-CQ*CB/dy/dz*(Ey1(i,j,k+1)-Ey1(i,j,k)-Ey1(i,j-1,k+1)+Ey1(i,j-1,k))...
            -CB/dy*(Hx1(i,j,k)-Hx1(i,j-1,k))+CB/dx*(Hy(i,j,k)-Hy(i-1,j,k));
        d_z2(1)=d_z2(1)-a_y*dy/(dy+c0*dt) *(Ez1(i,1,k)+Ez1(i,2,k));
        d_z2(ny2-1)=d_z2(ny2-1)-c_y*dy/(dy+c0*dt) *(Ez1(i,ny2+1,k)+Ez1(i,ny2,k));
        Ez(i,j,k)=JAy\d_z2';
    end
end
% 
Ez(1,:,:)=-mur_x*Ez(2,:,:)+dx/(dx+c0*dt)*(Ez1(1,:,:)+Ez1(2,:,:));
Ez(nx2+1,:,:)=-mur_x*Ez(nx2,:,:)+dx/(dx+c0*dt)*(Ez1(nx2,:,:)+Ez1(nx2+1,:,:));

Ez(:,1,:)=-mur_y*Ez(:,2,:)+dy/(dy+c0*dt)*(Ez1(:,1,:)+Ez1(:,2,:));
Ez(:,ny2+1,:)=-mur_y*Ez(:,ny2,:)+dy/(dy+c0*dt)*(Ez1(:,ny2,:)+Ez1(:,ny2+1,:));
%----------hx--------------- 
i=nx1:nx2+1;
j=ny1:ny2;
k=nz1:nz2;
Hx(i,j,k)=Hx1(i,j,k)-CQ/dy*(Ez(i,j+1,k)-Ez(i,j,k))+CQ/dz*(Ey1(i,j,k+1)-Ey1(i,j,k));
%----------ex--------------- 
for i=nx1:nx2
    for j=ny1+1:ny2
        k=nz1+1:nz2;
        d_x2(k-1)=Ex1(i,j,k)-CQ*CB/dx/dz*(Ez1(i+1,j,k)-Ez1(i,j,k)-Ez1(i+1,j,k-1)+Ez1(i,j,k-1))...
             -CB/dz*(Hy1(i,j,k)-Hy1(i,j,k-1))+CB/dy*(Hz1(i,j,k)-Hz1(i,j-1,k));
        d_x2(1)= d_x2(1)-a_z*dz/(dz+c0*dt)*(Ex1(i,j,1)+Ex1(i,j,2));
        d_x2(nz2-1)= d_x2(nz2-1)-c_z*dz/(dz+c0*dt)*(Ex1(i,j,nz2+1)+Ex1(i,j,nz2));
        Ex(i,j,k)=JAz\ d_x2';      
    end
end

Ex(:,1,:)=-mur_y*Ex(:,2,:)+dy/(dy+c0*dt)*(Ex1(:,1,:)+Ex1(:,2,:));
Ex(:,ny2+1,:)=-mur_y*Ex(:,ny2,:)+dy/(dy+c0*dt)*(Ex1(:,ny2,:)+Ex1(:,ny2+1,:));

Ex(:,:,1)=-mur_z*Ex(:,:,2)+dz/(dz+c0*dt)*(Ex1(:,:,1)+Ex1(:,:,2));
Ex(:,:,nz2+1)=-mur_z*Ex(:,:,nz2)+dz/(dz+c0*dt)*(Ex1(:,:,nz2)+Ex1(:,:,nz2+1));
%----------hy---------------- 

i=nx1:nx2;
j=ny1:ny2+1;
k=nz1:nz2;
Hy(i,j,k)=Hy1(i,j,k)-CQ/dz*(Ex(i,j,k+1)-Ex(i,j,k))+CQ/dx*(Ez1(i+1,j,k)-Ez1(i,j,k));

%---------------------------- 

Ez_time1(n)=Ez(2,7,2);
Ez_time2(n)=Ez(28,20,20);
Ez_time3(n)=Ez(15,15,15);
figure(1)
mesh((Ez(:,:,7)))
% plot((Ez(:,7,7)))
axis([0 35 0 30 -0.005 0.005])
end


