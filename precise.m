function [displ,accl,velo,T] = precise(M,C,K,F,init_disp,init_velo,dt)
% M、C、K，F分别为动力学方程 对应的矩阵
% init_disp 初始位移
% init_velo 初始速度
% dt 时间间隔 F=[Ft1;Ft1....;Ftn]
H1 = -1/2*M\C;
Ey = eye(size(M));
H2 = M\Ey;
H3 = -( K-1/4*C'*H2*C );
H4 = 1/2*C/M;
H = [H1 H2;H3 H4];
I = eye(size(H));
Ht = H*dt/2^20;
Ht2 = Ht*Ht;
Ha = Ht+Ht2*(I+Ht/3+Ht2/12)/2;
for i = 1:1:20
    Ha = 2*Ha + Ha*Ha;
end
E =Ha+I;
%% 位移 速度 加速度
displ = zeros(size(F));
velo = zeros(size(F));
accl = zeros(size(F));
displ(:,1) = init_disp';
velo(:,1) = init_velo';
Z0 = M*velo(:,1) + 1/2*C*displ(:,1);
Y0 = init_disp;
[m,n] = size(F);
X = zeros(m+m,n); % 用于存储Z(i)和Y(i) 竖着存储[Y(0) Y(1)...;Z(0) Z(1)...]
X(:,1) = [Y0'; Z0]; % 计算
for i = 2:1:n
    Ra = [zeros(size(F(:,i)));F(:,i-1)] ;
    Rb = [zeros(size(F(:,i)));( F(:,i)-F(:,i-1) )/dt] ;
    X(:,i) = E*( X(:,i-1) + H\(Ra+H\Rb) ) - H\(Ra+H\Rb+dt*Rb);
end
%% 提取位移和对偶向量
displ = X(1:m,:);
zCouple = X(m+1:m+m,:);% 提取对偶向量
%% 提取速度和加速度
for j =1:1:n
    velo(:,j) = M\zCouple(:,j) - 1/2*M\C*displ(:,j);
    accl(:,j) = M\( F(:,j) - K*displ(:,j) - C*velo(:,j) );
end
% 返回时间轴
T = 0:dt:dt*(n-1);
end