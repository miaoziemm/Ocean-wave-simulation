function [T,M] = zm_MSS3d(g,dw,dtheta,xmin,xmax,ymin,ymax,v,t,e,m)

x1=[xmin:xmax];
y1=[ymin:ymax];
w=[dw:dw:2*pi];
Sw=(((8.1*10^(-3)*g.^2))./(w.^5)).*exp(-0.74.*(g./(v.*w)).^4);
% theta=[dtheta:dtheta:2*pi];
df=pi/3;
for i=1:m
    theta(1,i)=df-pi/2+dtheta*(i-0.5);
end
wm=8.565/v;
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
T=zeros(length(x1),length(y1));
for x=1:length(x1)
    for y=1:length(y1)
        for i=1:length(w)
            for j=1:length(theta)
            T(x,y)=T(x,y)+(2*dw*dtheta*2*Swu(i,j)).^0.5*cos(w(1,i)*t-...
                w(1,i)^2/g*(x*cos(theta(1,j))+y*sin(theta(1,j)))+e(i,j));
            end
        end
    end
end
k=1;
for x=xmin:xmax
    for y=ymin:ymax
        M(1,k)=x;
        M(2,k)=y;
        M(3,k)=T(x,y);
        k=k+1;
    end
end

end

