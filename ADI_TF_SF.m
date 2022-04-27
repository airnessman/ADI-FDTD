clc;
clear;

mu0=pi*4e-7; %permability of free space
eps0=8.854E-12; %permittivity of free space
c0=1/sqrt(eps0*mu0); %velocity of light
nx1=1;
ny1=nx1;
nz1=nx1;

ja=20;

nx2=10;
ny2=200;
nz2=10;


je=200;
 CFLN=3;
dx=0.05;
dy=0.05;
dz=0.05;
dt=CFLN/(c0*sqrt(1/dx^2+1/dy^2+1/dz^2));
dt1=1/(c0*sqrt(1/dx^2+1/dy^2+1/dz^2));

Tw=0.8*1.0618e-9;
t0=2*Tw;

Hx(nx1:nx2+1,ny1:ny2,nz1:nz2)=0;
Hy(nx1:nx2,ny1:ny2+1,nz1:nz2)=0;
Hz(nx1:nx2,ny1:ny2,nz1:nz2+1)=0;
Ex(nx1:nx2,ny1:ny2+1,nz1:nz2+1)=0;
Ey(nx1:nx2+1,ny1:ny2,nz1:nz2+1)=0;
Ez(nx1:nx2+1,ny1:ny2+1,nz1:nz2)=0;

ez_inc=zeros(1,je+1);
hx_inc=zeros(1,je);

CA=1;
CB=(dt/2/eps0);
CP=1;
CQ=dt/2/mu0;

CB1=(dt1/2/eps0);
CQ1=dt1/2/mu0;

mur_y=(dy-c0*dt)/(dy+c0*dt);
mur_y1=(dy-c0*dt1)/(dy+c0*dt1);
% x方向
a1=-dt*dt/(4*mu0*eps0*dx*dx);
b1=1-2*a1;
c1=a1;

for m=nx1:nx2-2
A(m+1,m) = a1;
A(m,m) = b1;
A(m,m+1) = c1;
end
A(m+1,m+1)=b1;
A(1,m+1)=a1;
A(m+1,1)=c1;
JAx=sparse(A);
clear A;

%z方向
a1=-dt*dt/(4*mu0*eps0*dz*dz);
b1=1-2*a1;
c1=a1;

for m=nz1:nz2-2
A(m+1,m) = a1;
A(m,m) = b1;
A(m,m+1) = c1;
end
A(m+1,m+1)=b1;
A(1,m+1)=a1;
A(m+1,1)=c1;
JAz=sparse(A);
clear A;

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

%一维运算矩阵
a = -dt1^2/(4*mu0*eps0*dy^2);
b = 1-2*a;
c = a;
for i = 2:je
    A(i,i-1) = a;
    A(i,i) = b;
    A(i,i+1) = c;
end

A(1,1) = 1;
A(1,2) = mur_y1;
A(je+1,je) = mur_y1;
A(je+1,je+1) = 1;
A1 = sparse(A);



for n= 1:400
    n
    pulse=10*exp(-(n*dt1-t0)^2/(Tw)^2);

    ez_inc(10)=ez_inc(10)+(1)^0.5* pulse;
%     if n==10
%         Ez(:,10,:)=1;
%     end
    %第一步
    ez_pre=ez_inc;     hx_pre=hx_inc;
    Ex1=Ex; Ey1=Ey; Ez1=Ez; 
    Hx1=Hx; Hy1=Hy; Hz1=Hz;
    
    %-----1D
    k=2:je;
    d(1) = (dy/(c0*dt1+1*dy))*ez_pre(1)+(dy/(c0*dt1+1*dy))*ez_pre(2);
    d(k)=ez_pre(k)-CB1*(hx_pre(k)-hx_pre(k-1))/dy;
    d(je+1) = (dy/(c0*dt1+1*dy))*ez_pre(je+1)+(dy/(c0*dt1+1*dy))*ez_pre(je);
    ez_inc(1:je+1)=A1(1:je+1,1:je+1)\d(1:je+1)';
    
    hx_inc(1:je)=hx_pre(1:je)-CQ1*(ez_inc((1:je)+1)-ez_inc(1:je))/dy;
    
    %----------EZ---------------
     temp(ny1:ny2+1)=0; temp(ja)=CB*hx_inc(ja-1)/dz;

for k=nz1:nz2  
    for j=ny1+1:ny2   
        i=nx1+1:nx2;
        d_z1(i-1)=CA*Ez1(i,j,k)-CQ*CB/dx/dz*(Ex1(i,j,k+1)-Ex1(i,j,k)-Ex1(i-1,j,k+1)+Ex1(i-1,j,k))...
            -CB/dy*(Hx1(i,j,k)-Hx1(i,j-1,k))+CB/dx*(Hy1(i,j,k)-Hy1(i-1,j,k))+temp(j)';
        Ez(i,j,k)=JAx\d_z1';
    end
end

% Ez(nx1+1:nx2,ja,nz1:nz2)=Ez(nx1+1:nx2,ja,nz1:nz2)-CB*hx_inc(ja-1)/dz;

Ez(:,1,:)=-mur_y*Ez(:,2,:)+dy/(dy+c0*dt)*(Ez1(:,1,:)+Ez1(:,2,:));
Ez(:,ny2+1,:)=-mur_y*Ez(:,ny2,:)+dy/(dy+c0*dt)*(Ez1(:,ny2+1,:)+Ez1(:,ny2,:));

Ez(nx1,:,:)=Ez(nx1+1,:,:);
Ez(nx2+1,:,:)=Ez(nx2,:,:);
Ez(:,:,nz2)=Ez(:,:,nz1);

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
Ex(:,1,:)=-mur_y*Ex(:,2,:)+dy/(dy+c0*dt)*(Ex1(:,1,:)+Ex1(:,2,:));
Ex(:,ny2+1,:)=-mur_y*Ex(:,ny2,:)+dy/(dy+c0*dt)*(Ex1(:,ny2+1,:)+Ex1(:,ny2,:));

Ex(nx1,:,:)=Ex(nx2,:,:);
Ex(:,:,nz2+1)=Ex(:,:,nz2);
Ex(:,:,nz1)=Ex(:,:,nz1+1);
     
%---------EY----------------     
for i=nx1+1:nx2
    for j=ny1:ny2
        k=nz1+1:nz2;
        d_y1(k-1)=CA*Ey1(i,j,k)-CQ*CB/dy/dz*(Ez1(i,j+1,k)-Ez1(i,j,k)-Ez1(i,j+1,k-1)+Ez1(i,j,k-1))...
            -CB/dx*(Hz1(i,j,k)-Hz1(i-1,j,k))+CB/dz*(Hx1(i,j,k)-Hx1(i,j,k-1));        
        Ey(i,j,k)=JAz\d_y1';        
    end
end 
Ey(nx1,:,:)=Ey(nx1+1,:,:);
Ey(nx2+1,:,:)=Ey(nx2,:,:);
Ey(:,:,nz1)=Ey(:,:,nz1+1);
Ey(:,:,nz2+1)=Ey(:,:,nz2);

%----------HX----------------  
    i=nx1:nx2+1;
    j=ny1:ny2;
    k=nz1:nz2;
    Hx(i,j,k)=Hx1(i,j,k)+CQ/dz*(Ey(i,j,k+1)-Ey(i,j,k))-CQ/dy*(Ez1(i,j+1,k)-Ez1(i,j,k));
    
    Hx(nx1:nx2+1,ja-1,nz1:nz2)=Hx(nx1:nx2+1,ja-1,nz1:nz2)+CQ*ez_inc(ja)/dx;
    Hx(:,:,nz2)=Hx(:,:,nz1);
    
%--------HY----------------
    i=nx1:nx2;
    j=ny1:ny2+1;
    k=nz1:nz2;
    Hy(i,j,k)=Hy1(i,j,k)+CQ/dx*(Ez(i+1,j,k)-Ez(i,j,k))-CQ/dz*(Ex1(i,j,k+1)-Ex1(i,j,k));
    
    Hy(nx2,:,:)=Hy(nx1,:,:);
    Hy(:,:,nz2)=Hy(:,:,nz1);
    
%--------HZ-----------------
    i=nx1:nx2;
    j=ny1:ny2;
    k=nz1:nz2+1;
    Hz(i,j,k)=Hz1(i,j,k)+CQ/dy*(Ex(i,j+1,k)-Ex(i,j,k))-CQ/dx*(Ey1(i+1,j,k)-Ey1(i,j,k));
    
    Hy(nz2,:,:)=Hy(nz1,:,:);
    
    %第二步
    %1D
    
    ez_pre=ez_inc;     hx_pre=hx_inc;
    Ex1=Ex; Ey1=Ey; Ez1=Ez; 
    Hx1=Hx; Hy1=Hy; Hz1=Hz;
    
    ez_inc(2:je)=ez_pre(2:je)-CB1*(hx_pre(2:je)-hx_pre(1:je-1))/dy;
    ez_inc(1) = -mur_y1*ez_inc(2)+dy/(dy+c0*dt1)*(ez_pre(1)+ez_pre(2));
    ez_inc(je+1) = -mur_y1*ez_inc(je)+dy/(dy+c0*dt1)*(ez_pre(je)+ez_pre(je+1));
    
    hx_inc(1:je)=hx_pre(1:je)-CQ1*(ez_pre(2:je+1)-ez_pre(1:je))/dy;
    
    
    
    %--------Ey----------------   
for j=ny1:ny2
    for k=nz1+1:nz2
        i=nx1+1:nx2;
        d_y2(i-1)=Ey1(i,j,k)-CQ*CB/dx/dy*(Ex1(i,j+1,k)-Ex1(i,j,k)-Ex1(i-1,j+1,k)+Ex1(i-1,j,k))...
            -CB/dx*(Hz1(i,j,k)-Hz1(i-1,j,k))+CB/dz*(Hx1(i,j,k)-Hx1(i,j,k-1));
        Ey(i,j,k)=JAx\d_y2';
    end
end
Ey(nx1,:,:)=Ey(nx1+1,:,:);
Ey(nx2+1,:,:)=Ey(nx2,:,:);
Ey(:,:,nz1)=Ey(:,:,nz1+1);
Ey(:,:,nz2+1)=Ey(:,:,nz2);

%---EZ-------

temp(ny1:ny2+1)=0; 
temp(ja-1)=CB*CQ*ez_inc(ja-1)/dz^2; 
temp(ja)=-CB*CQ*ez_inc(ja)/dz^2+CB*hx_inc(ja)/dz;
% temp(ja)=CB*hx_inc(ja)/dz;
for i=nx1+1:nx2
    for k=nz1:nz2
        j=ny1+1:ny2;
        d_z2(j-1)=Ez1(i,j,k)-CQ*CB/dy/dz*(Ey1(i,j,k+1)-Ey1(i,j,k)-Ey1(i,j-1,k+1)+Ey1(i,j-1,k))...
            -CB/dy*(Hx1(i,j,k)-Hx1(i,j-1,k))+CB/dx*(Hy(i,j,k)-Hy(i-1,j,k))+temp(ny1+1:ny2);
%         d_z2(ja-1)=d_z2(ja-1)-c_y*ez_inc(ja-1);
%         d_z2(ja)=d_z2(ja)-c_y*ez_inc(je)+CB*hx_inc(ja)/dz;
        d_z2(1)=d_z2(1)-a_y*dy/(dy+c0*dt) *(Ez1(i,1,k)+Ez1(i,2,k));
        d_z2(ny2-1)=d_z2(ny2-1)-c_y*dy/(dy+c0*dt) *(Ez1(i,ny2+1,k)+Ez1(i,ny2,k));
        Ez(i,j,k)=JAy\d_z2';
    end
end

% Ez(nx1+1:nx2,ja,nz1:nz2)=Ez(nx1+1:nx2,ja,nz1:nz2)-CB*hx_inc(ja-1)/dz;

Ez(:,1,:)=-mur_y*Ez(:,2,:)+dy/(dy+c0*dt)*(Ez1(:,1,:)+Ez1(:,2,:));
Ez(:,ny2+1,:)=-mur_y*Ez(:,ny2,:)+dy/(dy+c0*dt)*(Ez1(:,ny2,:)+Ez1(:,ny2+1,:));

Ez(nx1,:,:)=Ez(nx1+1,:,:);
Ez(nx2+1,:,:)=Ez(nx2,:,:);
Ez(:,:,nz2)=Ez(:,:,nz1);

%----------ex--------------- 
for i=nx1:nx2
    for j=ny1+1:ny2
        k=nz1+1:nz2;
        d_x2(k-1)=Ex1(i,j,k)-CQ*CB/dx/dz*(Ez1(i+1,j,k)-Ez1(i,j,k)-Ez1(i+1,j,k-1)+Ez1(i,j,k-1))...
             -CB/dz*(Hy1(i,j,k)-Hy1(i,j,k-1))+CB/dy*(Hz1(i,j,k)-Hz1(i,j-1,k));
        Ex(i,j,k)=JAz\ d_x2';      
    end
end
Ex(:,1,:)=-mur_y*Ex(:,2,:)+dy/(dy+c0*dt)*(Ex1(:,1,:)+Ex1(:,2,:));
Ex(:,ny2+1,:)=-mur_y*Ex(:,ny2,:)+dy/(dy+c0*dt)*(Ex1(:,ny2,:)+Ex1(:,ny2+1,:));

Ex(nx1,:,:)=Ex(nx2,:,:);
Ex(:,:,nz2+1)=Ex(:,:,nz2);
Ex(:,:,nz1)=Ex(:,:,nz1+1);

%----------hx--------------- 
i=nx1:nx2+1;
j=ny1:ny2;
k=nz1:nz2;
Hx(i,j,k)=Hx1(i,j,k)-CQ/dy*(Ez(i,j+1,k)-Ez(i,j,k))+CQ/dz*(Ey1(i,j,k+1)-Ey1(i,j,k));

Hx(nx1:nx2+1,ja-1,nz1:nz2)=Hx(nx1:nx2+1,ja-1,nz1:nz2)+CQ*ez_inc(ja)/dx;

Hx(:,:,nz2)=Hx(:,:,nz1);

%----------hy---------------- 

i=nx1:nx2;
j=ny1:ny2+1;
k=nz1:nz2;
Hy(i,j,k)=Hy1(i,j,k)-CQ/dz*(Ex(i,j,k+1)-Ex(i,j,k))+CQ/dx*(Ez1(i+1,j,k)-Ez1(i,j,k));
    Hy(nx2,:,:)=Hy(nx1,:,:);
    Hy(:,:,nz2)=Hy(:,:,nz1);
%---------hz---------------  
i=nx1:nx2;
j=ny1:ny2;
k=nz1:nz2+1;
Hz(i,j,k)=Hz1(i,j,k)-CQ/dx*(Ey(i+1,j,k)-Ey(i,j,k))+CQ/dy*(Ex1(i,j+1,k)-Ex1(i,j,k));
Hy(nz2,:,:)=Hy(nz1,:,:);


figure(1)
mesh(squeeze(Ez(5,:,:)))
axis([1,nz2+1,1,ny2+1,-1,20]);
% plot(ez_inc)
% axis([1,je+1,-1,1]);
    
end










