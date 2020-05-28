%波动方程法模拟海面
clc;
clear;
tic;
load t1.mat
load t2.mat
g=9.8;
p=0.05;
dt=1;
dx=1;
x=[1:1:100];
y=[1:1:100];
h=zeros(100,100,101);
h(:,:,1)=T1;
h(:,:,2)=T2;
for k=2:400
    for i=2:99
        for j=2:99
            h(i,j,k+1)=2*h(i,j,k)-h(i,j,k-1)+(g*p*dt^2/dx^2)*(h(i+1,j,k)-2*h(i,j,k)...
                +h(i-1,j,k))+(g*p*dt^2/dx^2)*(h(i,j+1,k)-2*h(i,j,k)+h(i,j-1,k));
        end
    end
end
toc;
for i=1:400
    surf(h(:,:,i));axis([20 100 20 100 -10 10]);
    pause(0.1)
    saveas(gcf,['C:\Users\MiaoZHY\Desktop\毕设相关\vm_marine\data\wave/',int2str(i*1000),'.jpg']);
end