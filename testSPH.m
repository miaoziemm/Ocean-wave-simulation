
%% 参数初始化
clc;
clear;
g=9.8;
w=[2*pi/100:4*pi/100:2*pi]; 
xx=[1:50];
v=15;
e=rand(1,length(w))*2*pi;
base=[0:-1:-10];

k=w.^2/g;
Sw=(((8.1*10^(-3)*g.^2))./(w.^5)).*exp(-0.74.*(g./(v.*w)).^4);
r=(Sw.*2*pi/100).^0.5;

z=zeros(length(base),length(xx));
x=zeros(length(base),length(xx));

%% 利用Gerstner波生成粒子位置
for t=1:1
    for b=1:length(base)
        for i=1:length(x)
            for j=1:length(w)
                z(b,i)=z(b,i)+(-r(1,j)*(1/b)*cos(k(1,j)*xx(1,i)-w(1,j)*t+e(1,j)));
                x(b,i)=(r(1,j)*(1/b)*sin(k(1,j)*xx(1,i)-w(1,j)*t)+e(1,j));
            end
        end
         z(b,:)=base(1,b)+z(b,:);
         x(b,:)=xx-x(b,:);
    end
end

%% 计算密度参数
% x_rou=[1:100];
% y_rou=[1:100];
m_water=3*10^(-16)*20;
h_rou=8;

rou=zeros(length(xx),length(base));
for xxx=1:length(xx)
    for yyy=1:length(base)    
        for i=1:length(w)
            for j=1:length(base)
                cha=(xxx-x(j,i))^2+(-yyy/2-z(j,i))^2;
                if cha<=h_rou^2
                    rou(xxx,yyy)=rou(xxx,yyy)+m_water*(4/pi*h_rou^8)*(h_rou^2-cha)^3;
                end
            end
        end
    end
end

rou_t=rou(15:35,3:8);
T=zeros(length(rou_t(:,1))*length(rou_t(1,:)),5);
k=1;
for i=1:length(rou_t(:,1))
    for j=1:length(rou_t(1,:))
        T(k,1)=i;
        T(k,2)=-j;
        T(k,3)=rou_t(i,j);
        k=k+1;
    end
end

%% 计算声速
for i=1:length(T(:,3))
    T(i,4)=(T(i,3)-1)/T(i,3);
end
for i=1:length(T(:,3))
    T(i,5)=shengsu(25,100*T(i,4),abs(T(i,2)));
end

%% 保存数据
colums={'x','y','rou','s','v'};
data=table(T(:,1),T(:,2),T(:,3),T(:,4),T(:,5),'VariableNames',colums);
writetable(data,'result1.csv');
 

         
