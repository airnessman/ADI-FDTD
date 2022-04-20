clc
clear
%-----定义常数-------------
epsz=1/pi/36e9;uz=pi*(4e-7);ita=sqrt(uz/epsz);
cc=3e8;
delt=1e-3;
lamt=20*delt;
omiga=cc/lamt;
tao=2/omiga;
ta=0.8*tao;

dx=delt;dy=delt;dz=0.2*delt;
dt=dx/(0.3*cc);
%-----定义参数-------------
nx1=1;nx2=nx1+200;nljx=nx1+100;nmx=nljx+20;
ny1=1;ny2=ny1+14;
nz1=1;nz2=nz1+30; 

nfk=10;
nfg=5;
nyy1=ny1+2;nyy2=ny1+12;%10mm
nzz1=nz1+10;nzz2=nz1+15;%1mm
%---------定义电参数变量--------------------
a1=cc*dt/2;
a2=cc*dt/2;
a3=a1*a2;
aa1=(cc*dt/2-dx)/(cc*dt/2+dx);

%-----定义变量-------------
Hx(nx1:nx2,ny1:ny2,nz1:nz2)=0;
Hy(nx1:nx2,ny1:ny2,nz1:nz2)=0;
Hz(nx1:nx2,ny1:ny2,nz1:nz2)=0;
Ex(nx1:nx2,ny1:ny2,nz1:nz2)=0;
Ey(nx1:nx2,ny1:ny2,nz1:nz2)=0;
Ez(nx1:nx2,ny1:ny2,nz1:nz2)=0;
%-----定义稀疏矩阵-------------
%计算用矩阵
a=-a3/dz/dz;
b=1+2*a3/dz/dz;
c=-a3/dz/dz;
%Z方向

for i=1:nz2-1  
    if i<nz2-1
        j=i+1;
    else
        j=1;
    end
    B(j,i)=a;    
    B(i,i)=b;
    B(i,j)=c;
end
JA3=sparse(B);
clear B;

%%%--Z方向 for Ey
for i=1:nfg-2  

    B(i+1,i)=a;    
    B(i,i)=b;
    B(i,i+1)=c;
end
B(i+1,i+1)=b;
JA3y=sparse(B);
clear B;

%Y方向
a=-a3/dy/dy;
b=1+2*a3/dy/dy;
c=-a3/dy/dy;

for i=ny1:ny2-1  
    if i<ny2-1
        j=i+1;
    else
        j=1;
    end
    B(j,i)=a;    
    B(i,i)=b;
    B(i,j)=c;
end
JA2=sparse(B);
clear B;


%%%--Y方向 for Ez
for i=1:nfk-2  

    B(i+1,i)=a;    
    B(i,i)=b;
    B(i,i+1)=c;
end
B(i+1,i+1)=b;
JA2z=sparse(B);
clear B;

%X方向
a=-a3/dx/dx;
b=1+2*a3/dx/dx;
b1=1+2*a3/dx/dx-aa1*a3/dx/dx;
c=-a3/dx/dx;
for i=1:nx2-3  

    B(i+1,i)=a;    
    B(i,i)=b;
    B(i,i+1)=c;
end
B(1,1)=b1;
B(i+1,i+1)=b1;
JA1=sparse(B);
clear B;

%X方向 for Ey and Ez
for i=1:nx2-3 
    if i==nmx-1
        B(i,i)=1;
    elseif i==nmx-2
        B(i,i)=b;  
    else   
    B(i+1,i)=a;
    B(i,i)=b;
    B(i,i+1)=c;
    end
  
end
B(1,1)=b1;
B(i+1,i+1)=b1;
JA1x=sparse(B);
clear  a b c b1;

%-------------计算主体--------------------
tic
for n=1:1024*6
    n
    l=5*dx;

    t1=n*dt-(l-0.5*dx)/cc;
    t11=(n+0.5)*dt-(l-0.5*dx)/cc;
    t2=(n+0.5)*dt-l/cc;
    t3=(n+0.5)*dt-(l-1*dx)/cc;
    if t1<0
        Hy0=0;
    else
        Hy0=-100*exp(-4*pi*(t1-ta)^2/tao/tao);
    end
    if t11<0
        Hy1=0;
    else
        Hy1=-100*exp(-4*pi*(t11-ta)^2/tao/tao);
    end
    if t2<0
        Ez0=0;
    else
        Ez0=100*exp(-4*pi*(t2-ta)^2/tao/tao);
    end
    if t3<0
        Ez1=0;
    else
        Ez1=100*exp(-4*pi*(t3-ta)^2/tao/tao);
    end
%------------------------- 
%计算第一步
%------------------------- 

Ezl=Ez;
temp(nx1:nx2)=0; temp(nljx-1)=-a3/dx/dx*Ez0; temp(nljx)=a3/dx/dx*Ez1-a1/dx*Hy0;
for j=ny1+1:ny2
    for k=nz1:nz2-1
        JBB1(nx1:nx2-2)=(Ezl(nx1+1:nx2-1,j,k)-a3/dx/dz*(Ex(nx1+1:nx2-1,j,k+1)-Ex(nx1+1:nx2-1,j,k)-Ex(nx1:nx2-2,j,k+1)+Ex(nx1:nx2-2,j,k))...
            -a1/dy*(Hx(nx1+1:nx2-1,j,k)-Hx(nx1+1:nx2-1,j-1,k))+a1/dx*(Hy(nx1+1:nx2-1,j,k)-Hy(nx1:nx2-2,j,k))+temp(nx1+1:nx2-1)');
        JBB1(nmx-1)=0; 
        JBB1(nx1)=JBB1(nx1)+a3/dx/dx*(Ezl(nx1+1,j,k)-aa1*Ezl(nx1,j,k));   
        JBB1(nx2-2)=JBB1(nx2-2)+a3/dx/dx*(Ezl(nx2-1,j,k)-aa1*Ezl(nx2,j,k));
        Ez(nx1+1:nx2-1,j,k)=JA1x\JBB1(nx1:nx2-2)';%---x
    end
end

%-----模型-----
for j=nyy1+1:nyy2-1
    for k=nzz1:nzz2-1
          JBBB1(nx1:nx2-2)=(Ezl(nx1+1:nx2-1,j,k)-a3/dx/dz*(Ex(nx1+1:nx2-1,j,k+1)-Ex(nx1+1:nx2-1,j,k)-Ex(nx1:nx2-2,j,k+1)+Ex(nx1:nx2-2,j,k))...
            -a1/dy*(Hx(nx1+1:nx2-1,j,k)-Hx(nx1+1:nx2-1,j-1,k))+a1/dx*(Hy(nx1+1:nx2-1,j,k)-Hy(nx1:nx2-2,j,k))+temp(nx1+1:nx2-1)');
          JBBB1(nx1)=JBBB1(nx1)+a3/dx/dx*(Ezl(nx1+1,j,k)-aa1*Ezl(nx1,j,k));
          JBBB1(nx2-2)=JBBB1(nx2-2)+a3/dx/dx*(Ezl(nx2-1,j,k)-aa1*Ezl(nx2,j,k));
        Ez(nx1+1:nx2-1,j,k)=JA1\JBBB1(nx1:nx2-2)';
    end
end

Ez(nx1,:,:)=Ezl(nx1+1,:,:)-aa1*Ezl(nx1,:,:)+aa1*Ez(nx1+1,:,:);
Ez(nx2,:,:)=Ezl(nx2-1,:,:)-aa1*Ezl(nx2,:,:)+aa1*Ez(nx2-1,:,:);

Ez(:,ny1,:)=Ez(:,ny2,:);
Ez(:,:,nz2)=Ez(:,:,nz1);
%------------------------- 
Hyl=Hy;
for i=nx1:nx2-1
    for j=ny1:ny2
        for k=nz1:nz2-1
            Hy(i,j,k)=Hyl(i,j,k)+a2/dx*(Ez(i+1,j,k)-Ez(i,j,k))-a2/dz*(Ex(i,j,k+1)-Ex(i,j,k));
        end
    end
end
Hy(nljx-1,:,:)=Hy(nljx-1,:,:)-a2/dx*Ez0;
Hy(:,:,nz2)=Hy(:,:,nz1);
%------------------------- 
Exl=Ex;
for i=nx1:nx2-1
    for k=nz1+1:nz2
        Ex(i,ny1+1:ny2,k)=JA2\(Exl(i,ny1+1:ny2,k)-a3/dx/dy*(Ey(i+1,ny1+1:ny2,k)-Ey(i,ny1+1:ny2,k)-Ey(i+1,ny1:ny2-1,k)+Ey(i,ny1:ny2-1,k))...
            -a1/dz*(Hyl(i,ny1+1:ny2,k)-Hyl(i,ny1+1:ny2,k-1))+a1/dy*(Hz(i,ny1+1:ny2,k)-Hz(i,ny1:ny2-1,k)))';
    end
end
Ex(:,ny1,:)=Ex(:,ny2,:);
Ex(:,:,nz1)=Ex(:,:,nz2);
%------------------------- 
Hzl=Hz;
for i=nx1:nx2-1
    for j=ny1:ny2-1
        for k=nz1:nz2
            Hz(i,j,k)=Hzl(i,j,k)+a2/dy*(Ex(i,j+1,k)-Ex(i,j,k))-a2/dx*(Ey(i+1,j,k)-Ey(i,j,k));
        end
    end
end
Hz(:,ny2,:)=Hz(:,ny1,:);
%------------------------- 
Eyl=Ey;
for i=nx1+1:nx2-1
    for j=ny1:ny2-1
        JB1(nz1:nz2-1)=Eyl(i,j,nz1+1:nz2)-a3/dy/dz*(Ezl(i,j+1,nz1+1:nz2)-Ezl(i,j,nz1+1:nz2)-Ezl(i,j+1,nz1:nz2-1)+Ezl(i,j,nz1:nz2-1))...
            -a1/dx*(Hzl(i,j,nz1+1:nz2)-Hzl(i-1,j,nz1+1:nz2))+a1/dz*(Hx(i,j,nz1+1:nz2)-Hx(i,j,nz1:nz2-1));
        Ey(i,j,nz1+1:nz2)=JA3\JB1';
    end
end

%------模型------
Ey(nmx,:,:)=0;
for j=nyy1:nyy2-1
    JB11(1:nfg-1)=Eyl(nmx,j,nzz1+1:nzz2-1)-a3/dy/dz*(Ezl(nmx,j+1,nzz1+1:nzz2-1)-Ezl(nmx,j,nzz1+1:nzz2-1)-Ezl(nmx,j+1,nzz1:nzz2-2)+Ezl(nmx,j,nzz1:nzz2-2))...
        -a1/dx*(Hzl(nmx,j,nzz1+1:nzz2-1)-Hzl(nmx-1,j,nzz1+1:nzz2-1))+a1/dz*(Hx(nmx,j,nzz1+1:nzz2-1)-Hx(nmx,j,nzz1:nzz2-2));
    Ey(nmx,j,nzz1+1:nzz2-1)=JA3y\JB11';
end
Ey(nx1,:,:)=Eyl(nx1+1,:,:)-aa1*Eyl(nx1,:,:)+aa1*Ey(nx1+1,:,:);
Ey(nx2,:,:)=Eyl(nx2-1,:,:)-aa1*Eyl(nx2,:,:)+aa1*Ey(nx2-1,:,:);

Ey(:,ny2,:)=Ey(:,ny1,:);
Ey(:,:,nz1)=Ey(:,:,nz2);
%------------------------- 
Hxl=Hx;
for i=nx1:nx2
    for j=ny1:ny2-1
        for k=nz1:nz2-1
            Hx(i,j,k)=Hxl(i,j,k)+a2/dz*(Ey(i,j,k+1)-Ey(i,j,k))-a2/dy*(Ezl(i,j+1,k)-Ezl(i,j,k));
        end
    end
end
Hx(:,ny2,:)=Hx(:,ny1,:);
Hx(:,:,nz2)=Hx(:,:,nz1);
%------------------------- 
%计算第二步
%------------------------- 
Eyl1=Ey;
for j=ny1:ny2-1
    for k=nz1+1:nz2
        JBB2(nx1:nx2-2)=(Eyl1(nx1+1:nx2-1,j,k)-a3/dx/dy*(Ex(nx1+1:nx2-1,j+1,k)-Ex(nx1+1:nx2-1,j,k)-Ex(nx1:nx2-2,j+1,k)+Ex(nx1:nx2-2,j,k))...
            -a1/dx*(Hz(nx1+1:nx2-1,j,k)-Hz(nx1:nx2-2,j,k))+a1/dz*(Hx(nx1+1:nx2-1,j,k)-Hx(nx1+1:nx2-1,j,k-1)));
        JBB2(nmx-1)=0;
          JBB2(nx1)=JBB2(nx1)+a3/dx/dx*(Eyl1(nx1+1,j,k)-aa1*Eyl1(nx1,j,k));
          JBB2(nx2-2)=JBB2(nx2-2)+a3/dx/dx*(Eyl1(nx2-1,j,k)-aa1*Eyl1(nx2,j,k));
        Ey(nx1+1:nx2-1,j,k)=JA1x\JBB2(nx1:nx2-2)';%---x
    end
end

%------模型------
for j=nyy1:nyy2-1
    for k=nzz1+1:nzz2-1
          JBBB2(nx1:nx2-2)=(Eyl1(nx1+1:nx2-1,j,k)-a3/dx/dy*(Ex(nx1+1:nx2-1,j+1,k)-Ex(nx1+1:nx2-1,j,k)-Ex(nx1:nx2-2,j+1,k)+Ex(nx1:nx2-2,j,k))...
            -a1/dx*(Hz(nx1+1:nx2-1,j,k)-Hz(nx1:nx2-2,j,k))+a1/dz*(Hx(nx1+1:nx2-1,j,k)-Hx(nx1+1:nx2-1,j,k-1)));
          JBBB2(nx1)=JBBB2(nx1)+a3/dx/dx*(Eyl1(nx1+1,j,k)-aa1*Eyl1(nx1,j,k));
          JBBB2(nx2-2)=JBBB2(nx2-2)+a3/dx/dx*(Eyl1(nx2-1,j,k)-aa1*Eyl1(nx2,j,k));
        Ey(nx1+1:nx2-1,j,k)=JA1\JBBB2(nx1:nx2-2)';
    end
end

Ey(nx1,:,:)=Eyl1(nx1+1,:,:)-aa1*Eyl1(nx1,:,:)+aa1*Ey(nx1+1,:,:);
Ey(nx2,:,:)=Eyl1(nx2-1,:,:)-aa1*Eyl1(nx2,:,:)+aa1*Ey(nx2-1,:,:);

Ey(:,ny2,:)=Ey(:,ny1,:);
Ey(:,:,nz1)=Ey(:,:,nz2);
%------------------------- 
Hzl1=Hz;
for i=nx1:nx2-1
    for j=ny1:ny2-1
        for k=nz1:nz2
            Hz(i,j,k)=Hzl1(i,j,k)-a2/dx*(Ey(i+1,j,k)-Ey(i,j,k))+a2/dy*(Ex(i,j+1,k)-Ex(i,j,k));
        end
    end
end
Hz(:,ny2,:)=Hz(:,ny1,:);
%------------------------- 
temp(nx1:nx2)=0;temp(nljx)=-a1/dx*Hy1;
Ezl1=Ez;
for i=nx1+1:nx2-1
    for k=nz1:nz2-1
        Ez(i,ny1+1:ny2,k)=JA2\(Ezl1(i,ny1+1:ny2,k)-a3/dy/dz*(Eyl1(i,ny1+1:ny2,k+1)-Eyl1(i,ny1:ny2-1,k+1)-Eyl1(i,ny1+1:ny2,k)+Eyl1(i,ny1:ny2-1,k))...
            -a1/dy*(Hx(i,ny1+1:ny2,k)-Hx(i,ny1:ny2-1,k))+a1/dx*(Hy(i,ny1+1:ny2,k)-Hy(i-1,ny1+1:ny2,k))+temp(i))';
    end
end

%------模型------
Ez(nmx,:,:)=0;
for k=nzz1:nzz2-1
    Ez(nmx,nyy1+1:nyy2-1,k)=JA2z\(Ezl1(nmx,nyy1+1:nyy2-1,k)-a3/dy/dz*(Eyl1(nmx,nyy1+1:nyy2-1,k+1)-Eyl1(nmx,nyy1:nyy2-2,k+1)-Eyl1(nmx,nyy1+1:nyy2-1,k)+Eyl1(nmx,nyy1:nyy2-2,k))...
        -a1/dy*(Hx(nmx,nyy1+1:nyy2-1,k)-Hx(nmx,nyy1:nyy2-2,k))+a1/dx*(Hy(nmx,nyy1+1:nyy2-1,k)-Hy(nmx-1,nyy1+1:nyy2-1,k)))';
end

Ez(nx1,:,:)=Ezl1(nx1+1,:,:)-aa1*Ezl1(nx1,:,:)+aa1*Ez(nx1+1,:,:);
Ez(nx2,:,:)=Ezl1(nx2-1,:,:)-aa1*Ezl1(nx2,:,:)+aa1*Ez(nx2-1,:,:);

Ez(:,ny1,:)=Ez(:,ny2,:);
Ez(:,:,nz2)=Ez(:,:,nz1);
%-------------------------
Hxl1=Hx;
for i=nx1:nx2
    for j=ny1:ny2-1
        for k=nz1:nz2-1
            Hx(i,j,k)=Hxl1(i,j,k)-a2/dy*(Ez(i,j+1,k)-Ez(i,j,k))+a2/dz*(Eyl1(i,j,k+1)-Eyl1(i,j,k));
        end
    end
end
Hx(:,ny2,:)=Hx(:,ny1,:);
Hx(:,:,nz2)=Hx(:,:,nz1);
%-------------------------
Exl1=Ex;
for i=nx1:nx2-1
    for j=ny1+1:ny2
        JB2(nz1:nz2-1)=Exl1(i,j,nz1+1:nz2)-a3/dx/dz*(Ezl1(i+1,j,nz1+1:nz2)-Ezl1(i+1,j,nz1:nz2-1)-Ezl1(i,j,nz1+1:nz2)+Ezl1(i,j,nz1:nz2-1))...
            -a1/dz*(Hy(i,j,nz1+1:nz2)-Hy(i,j,nz1:nz2-1))+a1/dy*(Hzl1(i,j,nz1+1:nz2)-Hzl1(i,j-1,nz1+1:nz2));
        Ex(i,j,nz1+1:nz2)=JA3\JB2';
    end
end
Ex(:,:,nz1)=Ex(:,:,nz2);
Ex(:,ny1,:)=Ex(:,ny2,:);
%-------------------------
Hyl1=Hy;
for i=nx1:nx2-1
    for j=ny1:ny2
        for k=nz1:nz2-1
            Hy(i,j,k)=Hyl1(i,j,k)-a2/dz*(Ex(i,j,k+1)-Ex(i,j,k))+a2/dx*(Ezl1(i+1,j,k)-Ezl1(i,j,k));
        end
    end
end
Hy(nljx-1,:,:)=Hy(nljx-1,:,:)-a2/dx*Ez0;
Hy(:,:,nz2)=Hy(:,:,nz1);
%-------------------------
e00(n,:,:)=Ez(nmx+20,:,:);
h00(n,:,:)=Hy(nmx+20,:,:);
e50(n,:,:)=Ez(nljx-5,:,:);
h50(n,:,:)=Hy(nljx-5,:,:);
ezz(n)=Ez0;
hyy(n)=Hy0;

% figure(1)
% for z=nz1:nz2
%     for x=nx1:nx2
%     en(x,z)=Ez(x,7,z);
%     end
% end
% mesh(en)
% colorbar

% if mod(n,5)==0
%     figure(1)
%     mesh(squeeze(Ez(:,3,:)))
%     axis([0 51 0 201 -100 100]);
%     pause(0.001)

% end
EZ_time(n)=Ez(3,7,26);
end

toc




