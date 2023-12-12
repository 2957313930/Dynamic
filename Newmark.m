function [displament,velocity,acceleration,T] = Newmark(M,C,K,R,initDisp,initVelocity,dt,gamma,beta)
% 参数全部采用行向量 initDisp=[0,0,0],initVelocity=[0,0,0] 
% 输出结果为ti时刻的各项位移，速度，加速度
% displament = [第一项位移;第二项位移;第三项位移] ,velocity与accelaration同样的数据排列形式（行排列）
% gamma = 1/2   beta = 1/4 经验取值
    % 基本参数
    c1 = 1/(beta*dt^2);
    c2 = gamma/(beta*dt);
    c3 = 1/(beta*dt);
    c4 = 1/(2*beta)-1 ;
    c5 = gamma/beta-1;
    c6 = (gamma/2/beta-1)*dt;
    c7 = dt*(1-gamma);
    c8 = gamma*dt;
    % 初始化参数
    [row,col] = size(R); % 迭代次数col
    displament = zeros(row,col); % 位移存储
    velocity = zeros(row,col); % 速度存储
    acceleration = zeros(row,col); % 加速度存储
    % 计算初始加速度
    initAccleartion =  M\(R(:,1)  -K*initDisp' - C*initVelocity');
    acceleration(:,1) = initAccleartion ;% 存储初始加速度
    displament(:,1) = initDisp';
    velocity(:,1) = initVelocity';
    acceleration(:,1) = initAccleartion';
    Kbar = c1*M+c2*C+K;
    for i = 1:1:col
        Rbar = R(:,i) + M*(c1*displament(:,i)+c3*velocity(:,i)+c4*acceleration(:,i)) +...
               C*( c2*displament(:,i) + c5*velocity(:,i) + c6*acceleration(:,i) );
        displament(:,i+1) = Kbar\Rbar;
        acceleration(:,i+1) = c1*( displament(:,i+1) - displament(:,i)) - c3*velocity(:,i) - c4*acceleration(:,i);
        velocity(:,i+1) = velocity(:,i) + c7*acceleration(:,i) + c8*acceleration(:,i+1);
    end
 T = 0:dt:dt*col;
end

