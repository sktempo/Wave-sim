% the 2-d wave equation - simulation of a square drum hit in the center
clear
nx=350;
nt=1500;
dt = .005;

dx = 2.5/(nx-1);
dy = 2.5/(nx-1);

c=1;
c1=dt^2*c^2/dx^2;

un=zeros(nx,nx);
unm1=zeros(nx,nx);
u=zeros(nx,nx);
CFL=0.5;

for n=1:nt
    u(100:140,150)=0;
    u(100,150:220)=0;
    u(140,150:220)=0;
    if(n<200)
        u(100,2)=1*sin(n*2*pi/100);
    end
%     for k=1:nx
%         u(1,k)=u(2,k);
%         u(k,1)=u(k,2);
%         u(k,nx)=u(k,nx-1);
%         u(nx,k)=u(nx-1,k);
%     end
    u(1,:)=un(2,:)+((CFL-1)/(CFL+1))*(u(2,:)-un(1,:));
    u(end,:)=un(end-1,:)+((CFL-1)/(CFL+1))*(u(end-1,:)-un(end,:));
    u(:,1)=un(:,2)+((CFL-1)/(CFL+1))*(u(:,2)-un(:,1));
    u(:,end)=un(:,end-1)+((CFL-1)/(CFL+1))*(u(:,end-1)-un(:,end));
    
    unm1 = un;   % un becomes unm1, u becomes un
    un = u;      % and we're ready to calculate un+1
    
    if(1==1)     % use Laplace operator
        D=[0 1 0; 1 -4 1; 0 1 0]; % 2d Laplace operator
        u=2*un - unm1 + c1*conv2(un,D,'same'); % 
    else
      

      for i = 2:nx-1     % do explicit calculation
          for j = 2:nx-1
             u(i,j) = 2*un(i,j) - unm1(i,j) ...
             + c1*(un(i+1,j) -2*un(i,j) +  un(i-1,j) ...
              + un(i,j+1) -2*un(i,j) +  un(i,j-1) );
          end
       end
    end
        
    if mod(n,5)==0     % animate graph
        mesh(un)   
        axis([0 nx 0 nx -.12 .12])
        shg
        pause(0.01)
    end
end