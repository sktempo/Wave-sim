
% the 2-d wave equation - simulation of a square drum hit in the center
clear
nx=350;
nt=350;
dt = .005;

dx = 2.5/(nx-1);
dy = 2.5/(nx-1);

c=1;
c1=dt^2*c^2/dx^2;

un=zeros(nx,nx);
unm1=zeros(nx,nx);
u=zeros(nx,nx);

R = 20;
for i=1:nx      % set up initial pulse in center of grid
    for j=1:nx
        d = sqrt((i-125)^2 + (j-125)^2);
        if  d<R
            u(i,j) = .2*sin(pi/2+d/R*pi/2);
            un(i,j) = u(i,j);
        end
    end
end
 

for n=1:nt
    u(100:140,150)=0;
    u(100,150:220)=0;
    u(140,150:220)=0;
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