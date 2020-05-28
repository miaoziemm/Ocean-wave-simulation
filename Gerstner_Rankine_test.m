%Gerstner-Rankine Model
%the analyze of 2D x-z
clear;
clc;
k1=4.5;
k2=3.5;
k3=2.5;
k4=1.5;
r=0.2;
w=pi/3;
t=1;
z0=0;
x0=[1:0.01:10];
for i=1:length(x0)  
    z1(1,i)=z0-r*cos(k1*x0(1,i)-w*t);
    z2(1,i)=z0-r*cos(k2*x0(1,i)-w*t);
    z3(1,i)=z0-r*cos(k3*x0(1,i)-w*t);
    z4(1,i)=z0-r*cos(k4*x0(1,i)-w*t);
    
    x1(1,i)=x0(1,i)+r*sin(k1*x0(1,i)-w*t);
    x2(1,i)=x0(1,i)+r*sin(k2*x0(1,i)-w*t);
    x3(1,i)=x0(1,i)+r*sin(k3*x0(1,i)-w*t);
    x4(1,i)=x0(1,i)+r*sin(k4*x0(1,i)-w*t);
end
xx=x1+x2+x3+x4;
zz=z1+z2+z3+z4;
figure(1)
plot(x1,z1);hold on
plot(x2,z2);hold on
plot(x3,z3);hold on
plot(x4,z4);legend('k1=4.5','k2=3.5','k3=2.5','k4=1.5');
figure(2)
plot(xx,zz);