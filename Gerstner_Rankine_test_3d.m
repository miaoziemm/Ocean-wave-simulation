%Gerstner-Rankine Model
%the analyze of 3D x-y-z
clc;
clear all;
g=9.8;
v=8;
m=30;
df=pi/3;
xmin=1;
ymin=1;
xmax=100;
ymax=100;
dw=pi/100;
x1=[xmin:xmax];
y1=[ymin:ymax];
w=[dw:dw:2*pi];
wm=8.565/v;
dtheta=pi/m;
z0=0;
t=1;
for i=1:m
    theta(1,i)=df-pi/2+dtheta*(i-0.5);
end
Sw=(((8.1*10^(-3)*g.^2))./(w.^5)).*exp(-0.74.*(g./(v.*w)).^4);
p=(0.5+0.82*exp(-0.5*(w./wm).^4));
q=0.32*exp(-0.5*(w./wm).^4);
for i=1:length(w)
    for j=1:length(theta)
        G(i,j)=(1./pi).*(1+p(1,i).*cos(2*theta(1,j))+q(1,i).*cos(4*theta(1,j)));
    end
end
for i=1:length(w)
    for j=1:length(theta)
        Swu(i,j)=Sw(1,i)*G(i,j);
    end
end
a=(2*Swu*dw*dtheta).^0.5;

for t=1:1
    xx1=zeros(length(x1),length(y1));
    yy1=zeros(length(x1),length(y1));
    zz1=zeros(length(x1),length(y1));
    for i=1:length(x1)
        for j=1:length(y1)
            for k=1:length(w)
                for l=1:length(theta)
                    xx1(i,j)=(xx1(i,j)+cos(theta(1,l))*a(k,l)*sin(w(1,k)^2/g*(x1(1,i)*...
                        cos(theta(1,l))+y1(1,j)*sin(theta(1,l)))-w(1,k)*t));
                    yy1(i,j)=(yy1(i,j)+sin(theta(l))*a(k,l)*sin(w(1,k)^2/g*(x1(1,i)*...
                        cos(theta(1,l))+y1(1,j)*sin(theta(1,l)))-w(1,k)*t));
                    zz1(i,j)=(zz1(i,j)+a(k,l)*cos(w(1,k)^2/g*(x1(1,i)*...
                        cos(theta(1,l))+y1(1,j)*sin(theta(1,l)))-w(1,k)*t));
                end
            end
            xx1(i,j)=x1(1,i)-xx1(i,j);
            yy1(i,j)=y1(1,j)-yy1(i,j);
            zz1(i,j)=z0+zz1(i,j);
        end
    end
    figure(1)
    surf(xx1,yy1,zz1);axis([-10 150 -10 150 -20 20]);
    saveas(gcf,['/home/hanfx/Desktop/vm_marine/data/out/gerstner/',int2str(t*1000),'.jpg']);
end

% k=1;
% for i=1:100
%     for j=1:100
%         zs(1,k)=xx1(i,j);
%         zs(2,k)=yy1(i,j);
%         zs(3,k)=zz1(i,j);
%         k=k+1;
%     end
% end

% theta=[2*pi/200:dtheta:2*pi];
% e=rand(length(w),length(theta))*2*pi;
% [T,M]=zm_MSS3d(9.8,2*pi/200,2*pi/200,1,100,1,100,12.5,0.05,e);
e=zeros(length(w),length(theta));
for t=1:2
    [T,M]=zm_MSS3d(9.8,dw,dtheta,1,100,1,100,v,t,e,30);
    figure(2)
    surf(T);axis([-10 150 -10 150 -20 20]);
    saveas(gcf,['/home/hanfx/Desktop/vm_marine/data/out/lh/',int2str(t*1000),'.jpg']);
    
end

% for i=1:20
%     L(:,:)=X(i,:,:);
%     surf(L);
%     pause(0.05);
% end





