%关于SPH方法的实现（2D）
%参数初始化
clc;
clear;
g=9.8; %重力值
dw=2*pi/100;
w=[2*pi/100:dw:2*pi]; %频率细分
x=[1:100]; %二维模拟范围
v=8; %风速
e=rand(1,length(w))*2*pi; %设计随机矩阵
z0=0; %设计海面基准
z11=-0.5
z2=-1
z3=-1.5
z4=-2
z5=-2.5
k1=3.5;
k=w.^2/g;
Sw=(((8.1*10^(-3)*g.^2))./(w.^5)).*exp(-0.74.*(g./(v.*w)).^4);
r=(Sw.*dw).^0.5;
t=1;
z1=zeros(5,100);
x1=zeros(5,100);

for i=1:length(x)
    for j=1:length(w)
        z1(1,i)=z0+(z1(1,i)-r(1,j)*cos(k1*x(1,i)-w(1,j)*t+e(1,j)));
        x1(1,i)=x1(1,i)+(x(1,i)+r(1,j)*sin(k1*x(1,i)-w(1,j)*t+e(1,j)));
    end
end
subplot(5,1,1)
scatter(x1(1,:),z1(1,:));

t=1
for i=1:length(x)
    for j=1:length(w)
        z1(2,i)=z1(2,i)+z11-r(1,j)*(9/10)*cos(k1*x(1,i)-w(1,j)*t+e(1,j));
        x1(2,i)=x1(2,i)+x(1,i)+r(1,j)*(9/10)*sin(k1*x(1,i)-w(1,j)*t+e(1,j));
    end
end
subplot(5,1,2)
scatter(x1(2,:),z1(2,:));

for i=1:length(x)
    for j=1:length(w)
        z1(3,i)=z1(3,i)+z2-r(1,j)*(8/10)*cos(k1*x(1,i)-w(1,j)*t+e(1,j));
        x1(3,i)=x1(3,i)+x(1,i)+r(1,j)*(8/10)*sin(k1*x(1,i)-w(1,j)*t+e(1,j));
    end
end
subplot(5,1,3)
scatter(x1(3,:),z1(3,:));

for i=1:length(x)
    for j=1:length(w)
        z1(4,i)=z1(4,i)+z3-r(1,j)*(7/10)*cos(k1*x(1,i)-w(1,j)*t+e(1,j));
        x1(4,i)=x1(4,i)+x(1,i)+r(1,j)*(7/10)*sin(k1*x(1,i)-w(1,j)*t+e(1,j));
    end
end
subplot(5,1,4)
scatter(x1(4,:),z1(4,:));

for i=1:length(x)
    for j=1:length(w)
        z1(5,i)=z1(5,i)+z4-r(1,j)*(6/10)*cos(k1*x(1,i)-w(1,j)*t+e(1,j));
        x1(5,i)=x1(5,i)+x(1,i)+r(1,j)*(6/10)*sin(k1*x(1,i)-w(1,j)*t+e(1,j));
    end
end
subplot(5,1,5)
scatter(x1(5,:),z1(5,:));
for m=1:5
    for i=1:length(x)
        figure(2)
        scatter(x1(m,i),z1(m,i));hold on
    end
end



    