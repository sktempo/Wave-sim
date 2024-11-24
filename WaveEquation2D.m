clc
clear

Lx=100;
Ly=100;
dx=0.1;
dy=dx;
nx=fix(Lx/dx);
ny=fix(Ly/dy);
x=linspace(0,Lx,nx);
y=linspace(0,Ly,ny);

T=10; %time of total run

wn=zeros(nx,ny); %wave elavation at time n 
wnminus1=wn; %wave elavation at time n-1
wnplus1=wn; %wave elavation at time n+1 
%improtant to set the first two time steps to zero as intial condition

CFL=0.5;
c=1;
dt=CFL*dx/c;
f=0.1; %Hz
wavelength=c/f;

visualize=zeros(fix(T/dt),nx,ny);
visualizeindex=1;
t=0;
while(t<T)
    
    %Reflecting B.C
    %wn(:,[1 end])=0; % put 0 over all x and in y=1 & y=end (last y)
    %wn([1 end],:)=0;
    wn((Lx/2):(Lx/2+(0.1*Lx)),:)=0;
%     wn((Lx/2),(Ly/5):(Ly:2+(0.1*Ly)))=0;
%     wn((Lx/2+(0.1*Lx)),(Ly/2):(Ly:2+0.1*Ly))=0;
    
    %Absorbing B.C by mur's absorption approxmation
    wnplus1(1,:)=wn(2,:)+((CFL-1)/(CFL+1))*(wnplus1(2,:)-wn(1,:));
    wnplus1(end,:)=wn(end-1,:)+((CFL-1)/(CFL+1))*(wnplus1(end-1,:)-wn(end,:));
    wnplus1(:,1)=wn(:,2)+((CFL-1)/(CFL+1))*(wnplus1(:,2)-wn(:,1));
    wnplus1(:,end)=wn(:,end-1)+((CFL-1)/(CFL+1))*(wnplus1(:,end-1)-wn(:,end));
    
    % single slit  
%     slit=wavelength*20;
%     wnplus1(nx/2,(1:(ny/2)-slit))=0;
%     wnplus1(nx/2,((ny/2)+slit: end))=0;
    
    % double slit  
%     slit=wavelength;
%     wnplus1(nx/2,(1:(ny/2)-1.5*slit))=0;
%     wnplus1(nx/2,((ny/2)-0.5*slit:(ny/2)+0.5*slit))=0;
%     wnplus1(nx/2,((ny/2)+1.5*slit:end))=0;
    
    
    %solution
    t=t+dt;
    wnminus1=wn;
    wn=wnplus1;
    
    %source 
    wn(nx/10,ny/2)=(dt^2)*20*sin(2*pi*f*t); %source wave sin with amp of 20 at middle
    

    
    for i=2:nx-1
        for j=2:ny-1
            wnplus1(i,j)=2*wn(i,j)-wnminus1(i,j)+(CFL^2)*(wn(i+1,j)+wn(i,j+1)-4*wn(i,j)+wn(i-1,j)+wn(i,j-1));
            
        end
    end
    visualize(visualizeindex,:,:)=wn;
    visualizeindex=visualizeindex+1;


end

for k=1:visualizeindex-1
    visualizemat(:,:)=visualize(k,:,:);
    surf(x,y,visualizemat)
    shading interp
    %axis off
    zlim ([-0.1 0.1])
    pause(dt/100000)
end
