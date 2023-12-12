clc;
clear;
% % %% 测试Z
% E = 35.5e9;
% G = 14.8e9;
% A = 0.29;
% rho = 2500;
% m = rho*A;
% % 坐标参数(形心与剪心相对坐标)
% ay =0.046;
% az =0.032;
% % 几何参数
% L = 10;
% Iy = 1.68E-2;
% Iz = 2.89E-2;
% Iw = 1.47E-3;
% J =  1.02E-3;
% rz = sqrt(Iz/A);
% ry = sqrt(Iy/A);
% rw = sqrt(sqrt(Iw/A));  
% r = sqrt(J/A);
% % 荷载
% qy = 2e3;
% qz = 2e3;
% y0 = 0.3; % 与荷载qy形成扭矩
% z0= 0.2;
% % 荷载频率与速度、时间
% 
% % omega = 77.53;% omega_n1
% % omega = 109.9;%omega_n2
% % omega = 336.1;%omega_n3
% omega= 336.1;%任意荷载频率
% 
% c = 200;
%% 测试H
E = 35.5e9;
G = 14.8e9;
A = 0.14;
rho = 2500;
m = rho*A;

% 坐标参数(形心与剪心相对坐标)
ay =0.006;
az =0.037;

% 几何参数
L = 10;% 原来速度10
Iy = 4.95E-3;
Iz = 7.32E-3;
Iw = 5.37E-5;
J =  4.96E-4;
rz = sqrt(Iz/A);
ry = sqrt(Iy/A);
rw = sqrt(sqrt(Iw/A));  
r = sqrt(J/A);
% 荷载
qy = 2e3;
qz = 2e3;
y0 = 0.3; % 与荷载qy形成扭矩
z0=0.2;


% % 荷载频率与速度、时间
% omega1 = 69.75; % omega_n1
% omega1 = 81.34; % omega_n2
% omega1 = 250.5; % omega_n3
omega = 80; %任意荷载频率

c = 1200;% 速度
%% 中间参数
dispsum = 0;
for n=1:1:1
    omega_vn = sqrt( E*Iz*(n*pi/L)^4/(m+m*(n*pi*rz/L)^2) );
    omega_wn = sqrt( E*Iy*(n*pi/L)^4/(m+m*(n*pi*ry/L)^2) );
    omega_tn = sqrt( (E*Iw+G*J)*(n*pi/L)^4/(m*r^2+m*(n*pi*rw^2)^2) );
    alpha_n = az/(1+(n*pi*rz/L)^2);
    beta_n = ay/(1+(n*pi*ry/L)^2);
    gamma_vn = az/(r^2+(n*pi*rw^2/L)^2);
    gamma_wn = ay/(r^2+(n*pi*rw^2/L)^2);
    alpha_fy = qy/(L*m)/(1+(n*pi*rz/L)^2);
    alpha_fz = qz/(L*m)/(1+(n*pi*ry/L)^2);
    alpha_ft = (qy*y0+qz*z0)/(m*L)/(r^2+(n*pi*rw^2/L)^2);
    Omega_cn = c*n*pi/L; % 速度参数
    %% 刚度矩阵，质量矩阵，阻尼矩阵，荷载矩阵
    M = [1,0,alpha_n;0,1,-beta_n;gamma_vn,-gamma_wn,1];
    C= zeros(3,3);
    K = [omega_vn^2,0,0;0,omega_wn^2,0;0,0,omega_tn^2] ;
    dt = 0.001;% 时间步长
    t = 0:dt:L/c+0;
%     disp(max(t));
    R1 = 2*alpha_fy*sin(omega*t).*sin(Omega_cn*t).*(heaviside(t)-heaviside(t-L/c));
    R2 = 2*alpha_fz*sin(omega*t).*sin(Omega_cn*t).*(heaviside(t)-heaviside(t-L/c));
    R3 = 2*alpha_ft*sin(omega*t).*sin(Omega_cn*t).*(heaviside(t)-heaviside(t-L/c));
    R= [R1;R2;R3];
    %% 求解
    initDisp = [0,0,0];
    initVelocity  = [0,0,0];
    gamma = 1/2;
    beta = 1/4;
    % [displ,velocity,acceleration,T] = Newmark(M,C,K,R,initDisp,initVelocity,dt,gamma,beta);
    [displ,accl,velo,T] = precise(M,C,K,R,initDisp,initVelocity,dt);
    dispsum = dispsum+displ;
    [m,~] = size(t);
end
max(dispsum,[],2);
% writematrix([T;dispsum],'E:\Project\TriplyCoupledBeam\resonance.xlsx')
%% 绘图
[tttt,ssss] = title('Response Curve--Triply coupled system','Velocity='+string(c),'Color','blue'); % 字体设置
tttt.FontSize = 16;% 字体大小
ssss.FontAngle = 'italic'; %字体
set(gcf, 'Position', [100, 100, 800,700]);
movegui(gcf, 'center');
% 竖向振动v
subplot(2,2,1)
plot(T,dispsum(1,:),'r');
title('竖向振动V')
xlabel('时间/s');
ylabel('位移/m');
grid on;
% 横向振动w
subplot(2,2,2)
plot(T,dispsum(2,:),'k');
title('横向振动W')
xlabel('时间/s');
ylabel('位移/m');
grid on;
% 纵向扭转ψ
subplot(2,2,3)
plot(T,dispsum(3,:),'b');
title('扭转振动ψ')
xlabel('时间/s');
ylabel('扭转角/rad');
grid on;
% 总图
subplot(2,2,4)
plot(T,dispsum(1,:) ,'r--',T,dispsum(2,:),'k',T,dispsum(3,:),'b-.');
legend('位移v','位移w','位移ψ');
title('ω='+string(omega)+'|v='+string(c));
grid on;






