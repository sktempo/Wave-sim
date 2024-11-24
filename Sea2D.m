clc
clear

Hs=2.8;
Tp=5;
dw=0.05;
w=0.2:dw:2;
s=487*(Hs^2)./((Tp^4).*(w.^5)).*exp(-1948./((Tp^4).*(w.^4)));
x=0:0.1:20;
dt=0.05;
t=0:dt:3;
K=(w.^2)./9.81;
phix=(2*pi)*rand(1,length(w));
phiy=(2*pi)*rand(1,length(w));
amp=sqrt(2*dw*s);
sinx=zeros(length(x),length(w));
siny=zeros(length(x),length(w));
sinmult=zeros(length(t),length(x),length(x),length(w));
wave=zeros(length(t),length(x),length(x));

% [T,X,Y,W]=ndgrid(t,x,x,w);
% K=(W.^2)./9.81;
% phix=(2*pi)*rand(length(x),length(x),length(w));
% phiy=(2*pi)*rand(length(x),length(x),length(w));
% PHIX=zeros(length(t),length(x),length(x),length(w));
% PHIY=zeros(length(t),length(x),length(x),length(w));
% AMP=zeros(length(t),length(x),length(x),length(w));
% for o=1:length(t)
%     PHIX(o,:,:,:)=phix;
%     PHIY(o,:,:,:)=phiy;
%     for i=1:length(x)
%         for j=1:length(x)
%             AMP(o,i,j,:)=amp;
%         end
%     end
% end

% sinmult=(AMP.*sin(W.*T+K.*X+PHIX)).*(AMP.*(sin(W.*T+K.*Y+PHIY)));
%  for l=1:length(w)
%      wave=wave+sinmult(:,:,:,l);
%  end

for i=1:length(t)
    for j=1:length(x)
        for k=1:length(x)
     
            sinx(j,:)=amp.*sin(w*t(i)+K*x(j)+phix);
            siny(k,:)=amp.*sin(w*t(i)+K*x(k)+phiy);
            sinmult(i,j,k,:)=sinx(j,:).*siny(k,:);
           
        end   
    end
    for l=1:length(w)
        wave(i,:,:)=wave(i,:,:)+sinmult(i,:,:,l);
    end
    
end

for v=1:length(t)
    wavemat(:,:)=wave(v,:,:);
    surf(x,x,wavemat)
    shading interp
    %axis off
    zlim ([-inf inf])
    pause(dt)
end
