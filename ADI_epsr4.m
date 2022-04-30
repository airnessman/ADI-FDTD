clc;
clear;

%% 初始设置
mu0=pi*4e-7; %permability of free space
eps0=8.854E-12; %permittivity of free space
epsr=4;
c0=1/sqrt(eps0*mu0); %velocity of light
%网格大小
nx1=1;
ny1=nx1;
nz1=nx1;
nx2=8;
ny2=500;
nz2=8;

je=200;

ja=20; %总场散射场位置

CFLN=3;
dx=0.05;
dy=0.05;
dz=0.05;
dt=CFLN/(c0*sqrt(1/dx^2+1/dy^2+1/dz^2));
dt1=1/(c0*sqrt(1/dx^2+1/dy^2+1/dz^2));

Tw=0.5*1.0618e-9;  %源参数
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

%媒质位置
my1=80;  my2=my1+29;

CBx(1:nx2,1:ny2+1,1:nz2+1)=CB;
CBx(:,my1:my2+1,:)=CB/epsr;

CBy(1:nx2+1,1:ny2,1:nz2+1)=CB;
CBy(:,my1:my2,:)=CB/epsr;

CBz(1:nx2+1,1:ny2+1,1:nz2)=CB;
CBz(:,my1:my2+1,:)=CB/epsr;
%% 矩阵设置

%一维矩阵
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
clear A

%————————————————
%x方向 自由空间

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

%x方向 媒质
a1=-dt*dt/(4*mu0*eps0*epsr*dx*dx);
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
JAx1=sparse(A);
clear A;

%z 方向
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

%z 方向媒质
a1=-dt*dt/(4*mu0*eps0*epsr*dz*dz);
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
JAz1=sparse(A);
clear A;

%y 方向
a1=-dt*dt/(4*mu0*eps0*dy*dy);
b1=1-2*a1;
c1=a1;

a_e=-dt*dt/(4*mu0*eps0*epsr*dy*dy);
b_e=1-2*a1;
c_e=a1;

for m=ny1:ny2-2
A(m+1,m) = a1;
A(m,m) = b1;
A(m,m+1) = c1;
end
A(1,1)=b1-a1*mur_y;
A(m+1,m+1)=b1-a1*mur_y;

for m=my1-1:my2-1
A(m,m-1) = a1/epsr;
A(m,m) = 1-2*a1/epsr;
A(m,m+1) = c1/epsr;
end
JAy=sparse(A);
clear A;
a_y=a1;   c_y=c1;
%%  main

for n= 1:400
    n
    pulse=10*exp(-(n*dt1-t0)^2/(Tw)^2);
    ez_inc(10)=(1)^0.5* pulse;   
    %% 第一步       
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
 %----———Ez———————————   
    temp(ny1:ny2+1)=0; temp(ja)=CB*hx_inc(ja-1)/dz;

for j=ny1+1:ny2
    for k=nz1:nz2         
        i=nx1+1:nx2;
        
        d_z1(i-1)=CA*Ez1(i,j,k)-CQ*CBz(i,j,k)/dx/dz.*(Ex1(i,j,k+1)-Ex1(i,j,k)-Ex1(i-1,j,k+1)+Ex1(i-1,j,k))...
            -CBz(i,j,k)/dy.*(Hx1(i,j,k)-Hx1(i,j-1,k))+CBz(i,j,k)/dx.*(Hy1(i,j,k)-Hy1(i-1,j,k))+temp(j)';
        
        if(j<my1||j>my2)
            Ez(i,j,k)=JAx\d_z1';
        else
            Ez(i,j,k)=JAx1\d_z1';
        end
    end
end
    Ez(:,1,:)=-mur_y*Ez(:,2,:)+dy/(dy+c0*dt)*(Ez1(:,1,:)+Ez1(:,2,:));
    Ez(:,ny2+1,:)=-mur_y*Ez(:,ny2,:)+dy/(dy+c0*dt)*(Ez1(:,ny2+1,:)+Ez1(:,ny2,:));

    Ez(nx1,:,:)=Ez(nx2,:,:);
    Ez(nx2+1,:,:)=Ez(nx1+1,:,:);
%     Ez(:,:,nz2)=Ez(:,:,nz1);
 %----———Ex———————————   
 for i=nx1:nx2
    for k=nz1+1:nz2
        j= ny1+1:ny2; 
        d_x1(j-1)=CA*Ex1(i,j,k)-CQ*CBx(i,j,k)/dx/dy.*(Ey1(i+1,j,k)-Ey1(i,j,k)-Ey1(i+1,j-1,k)+Ey1(i,j-1,k))...
            -CBx(i,j,k)/dz.*(Hy1(i,j,k)-Hy1(i,j,k-1))+CBx(i,j,k)/dy.*(Hz1(i,j,k)-Hz1(i,j-1,k)); 
        
        d_x1(1)=d_x1(1)- a_y *  dy/(dy+c0*dt) *( Ex1(i,1,k) + Ex1(i,2,k) );
        d_x1(ny2-1)=d_x1(ny2-1)- c_y *  dy/(dy+c0*dt) *(Ex1(i,ny2+1,k) + Ex1(i,ny2,k) );
        
        Ex(i,j,k)=JAy\d_x1';     
    end 
 end
    Ex(:,1,:)=-mur_y*Ex(:,2,:)+dy/(dy+c0*dt)*(Ex1(:,1,:)+Ex1(:,2,:));
    Ex(:,ny2+1,:)=-mur_y*Ex(:,ny2,:)+dy/(dy+c0*dt)*(Ex1(:,ny2+1,:)+Ex1(:,ny2,:));


    Ex(:,:,nz2+1)=Ex(:,:,nz1+1);
    Ex(:,:,nz1)=Ex(:,:,nz2);
 
 %----———Ey———————————  
 
for j=ny1:ny2   
    for i=nx1+1:nx2
        k=nz1+1:nz2;
        
        d_y1(k-1)=CA*Ey1(i,j,k)-CQ*CBy(i,j,k)/dy/dz.*(Ez1(i,j+1,k)-Ez1(i,j,k)-Ez1(i,j+1,k-1)+Ez1(i,j,k-1))...
            -CBy(i,j,k)/dx.*(Hz1(i,j,k)-Hz1(i-1,j,k))+CBy(i,j,k)/dz.*(Hx1(i,j,k)-Hx1(i,j,k-1));
        
        if(j<my1||j>my2)
            Ey(i,j,k)=JAz\d_y1'; 
        else
            Ey(i,j,k)=JAz1\d_y1'; 
        end
    end   
end
    Ey(nx1,:,:)=Ey(nx2,:,:);
    Ey(nx2+1,:,:)=Ey(nx1+1,:,:);
    Ey(:,:,nz1)=Ey(:,:,nz2);
    Ey(:,:,nz2+1)=Ey(:,:,nz1+1);
 
 %----------HX----------------  
    i=nx1:nx2+1;
    j=ny1:ny2;
    k=nz1:nz2;
    Hx(i,j,k)=Hx1(i,j,k)+CQ/dz*(Ey(i,j,k+1)-Ey(i,j,k))-CQ/dy*(Ez1(i,j+1,k)-Ez1(i,j,k));
    
    Hx(nx1:nx2+1,ja-1,nz1:nz2)=Hx(nx1:nx2+1,ja-1,nz1:nz2)+CQ*ez_inc(ja)/dx;

    
%--------HY----------------
    i=nx1:nx2;
    j=ny1:ny2+1;
    k=nz1:nz2;
    Hy(i,j,k)=Hy1(i,j,k)+CQ/dx*(Ez(i+1,j,k)-Ez(i,j,k))-CQ/dz*(Ex1(i,j,k+1)-Ex1(i,j,k));
    

    
%--------HZ-----------------
    i=nx1:nx2;
    j=ny1:ny2;
    k=nz1:nz2+1;
    Hz(i,j,k)=Hz1(i,j,k)+CQ/dy*(Ex(i,j+1,k)-Ex(i,j,k))-CQ/dx*(Ey1(i+1,j,k)-Ey1(i,j,k));
    
    
    
    %% 第二步
        ez_pre=ez_inc;     hx_pre=hx_inc;
    Ex1=Ex; Ey1=Ey; Ez1=Ez; 
    Hx1=Hx; Hy1=Hy; Hz1=Hz;
   %———1D————— 
    ez_inc(2:je)=ez_pre(2:je)-CB1*(hx_pre(2:je)-hx_pre(1:je-1))/dy;
    ez_inc(1) = -mur_y1*ez_inc(2)+dy/(dy+c0*dt1)*(ez_pre(1)+ez_pre(2));
    ez_inc(je+1) = -mur_y1*ez_inc(je)+dy/(dy+c0*dt1)*(ez_pre(je)+ez_pre(je+1));
    
    hx_inc(1:je)=hx_pre(1:je)-CQ1*(ez_pre(2:je+1)-ez_pre(1:je))/dy;
    
    %--------Ey----------------   
for j=ny1:ny2
    for k=nz1+1:nz2
        i=nx1+1:nx2;
        d_y2(i-1)=Ey1(i,j,k)-CQ*CBy(i,j,k)/dx/dy.*(Ex1(i,j+1,k)-Ex1(i,j,k)-Ex1(i-1,j+1,k)+Ex1(i-1,j,k))...
            -CBy(i,j,k)/dx.*(Hz1(i,j,k)-Hz1(i-1,j,k))+CBy(i,j,k)/dz.*(Hx1(i,j,k)-Hx1(i,j,k-1));
        if(j<my1||j>my2)
            Ey(i,j,k)=JAx\d_y2';
        else
            Ey(i,j,k)=JAx1\d_y2';
        end
    end
end
Ey(nx1,:,:)=Ey(nx2,:,:);
Ey(nx2+1,:,:)=Ey(nx1+1,:,:);
Ey(:,:,nz1)=Ey(:,:,nz2);
Ey(:,:,nz2+1)=Ey(:,:,nz1+1);

%---EZ-------

temp(ny1:ny2+1)=0; 
temp(ja-1)=-CB*CQ*ez_inc(ja)/dz^2; 
temp(ja)=CB*CQ*ez_inc(ja-1)/dz^2+CB*hx_inc(ja-1)/dz;

for i=nx1+1:nx2
    for k=nz1:nz2
        j=ny1+1:ny2;
        d_z2(j-1)=Ez1(i,j,k)-CQ*CBy(i,j,k)/dy/dz.*(Ey1(i,j,k+1)-Ey1(i,j,k)-Ey1(i,j-1,k+1)+Ey1(i,j-1,k))...
            -CBy(i,j,k)/dy.*(Hx1(i,j,k)-Hx1(i,j-1,k))+CBy(i,j,k)/dx.*(Hy(i,j,k)-Hy(i-1,j,k))+temp(ny1+1:ny2);

        d_z2(1)=d_z2(1)-a_y*dy/(dy+c0*dt) *(Ez1(i,1,k)+Ez1(i,2,k));
        d_z2(ny2-1)=d_z2(ny2-1)-c_y*dy/(dy+c0*dt) *(Ez1(i,ny2+1,k)+Ez1(i,ny2,k));
        Ez(i,j,k)=JAy\d_z2';
    end
end

Ez(:,1,:)=-mur_y*Ez(:,2,:)+dy/(dy+c0*dt)*(Ez1(:,1,:)+Ez1(:,2,:));
Ez(:,ny2+1,:)=-mur_y*Ez(:,ny2,:)+dy/(dy+c0*dt)*(Ez1(:,ny2,:)+Ez1(:,ny2+1,:));

Ez(nx1,:,:)=Ez(nx2,:,:);
Ez(nx2+1,:,:)=Ez(nx1+1,:,:);


 %----------ex--------------- 

for j=ny1+1:ny2       
    for i=nx1:nx2
        k=nz1+1:nz2;
        d_x2(k-1)=Ex1(i,j,k)-CQ*CBx(i,j,k)/dx/dz.*(Ez1(i+1,j,k)-Ez1(i,j,k)-Ez1(i+1,j,k-1)+Ez1(i,j,k-1))...
             -CBx(i,j,k)/dz.*(Hy1(i,j,k)-Hy1(i,j,k-1))+CBx(i,j,k)/dy.*(Hz1(i,j,k)-Hz1(i,j-1,k));
         if(j<my1||j>my2)
           Ex(i,j,k)=JAz\ d_x2';    
         else
            Ex(i,j,k)=JAz1\ d_x2';  
         end
    end
end
Ex(:,1,:)=-mur_y*Ex(:,2,:)+dy/(dy+c0*dt)*(Ex1(:,1,:)+Ex1(:,2,:));
Ex(:,ny2+1,:)=-mur_y*Ex(:,ny2,:)+dy/(dy+c0*dt)*(Ex1(:,ny2,:)+Ex1(:,ny2+1,:));

Ex(:,:,nz2+1)=Ex(:,:,nz1+1);
Ex(:,:,nz1)=Ex(:,:,nz2);   
    
    %----------hx--------------- 
i=nx1:nx2+1;
j=ny1:ny2;
k=nz1:nz2;
Hx(i,j,k)=Hx1(i,j,k)-CQ/dy*(Ez(i,j+1,k)-Ez(i,j,k))+CQ/dz*(Ey1(i,j,k+1)-Ey1(i,j,k));

Hx(nx1:nx2+1,ja-1,nz1:nz2)=Hx(nx1:nx2+1,ja-1,nz1:nz2)+CQ*ez_inc(ja)/dx;


%----------hy---------------- 

i=nx1:nx2;
j=ny1:ny2+1;
k=nz1:nz2;
Hy(i,j,k)=Hy1(i,j,k)-CQ/dz*(Ex(i,j,k+1)-Ex(i,j,k))+CQ/dx*(Ez1(i+1,j,k)-Ez1(i,j,k));

%---------hz---------------  
i=nx1:nx2;
j=ny1:ny2;
k=nz1:nz2+1;
Hz(i,j,k)=Hz1(i,j,k)-CQ/dx*(Ey(i+1,j,k)-Ey(i,j,k))+CQ/dy*(Ex1(i,j+1,k)-Ex1(i,j,k));



figure(1)
mesh(squeeze(Ez(:,:,5)))
axis([1,ny2,1,nx2,-10,10]);
 
 
end



