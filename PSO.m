%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 粒子群优化算法设计
% 压缩粒子群优化算法设计
% 《说明》
% c1,c2低的值使粒子在目标区域外徘徊，而高的值导致粒子越过目标区域。 推荐取值范围：[0,4];
% 典型取值 : c1=c2=2           ;c1=1.6、c2=1.8         ;c1=1.6、c2 =2
% 针对不同的问题有不同的取值，一般通过在一个区间内试凑来调整这两个值
% 速度一般取边界范围的10%-20%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [groupBest,groupBestX,convergence] = PSO(func,xbound,velo_ratio,iternum,N)
% fit：目标函数
% bound: 变量约束[a,b;c,d]表示a<x1<b、c<x2<d
% iternum : 迭代次数
% vbound :速度约束[最小速度，最大速度]
% N 种群数量
convergence = zeros(1,iternum); % 收敛曲线
figure; % 创建动图窗口
scatterPlot = scatter(nan,nan,'filled',MarkerEdgeColor='k',MarkerFaceColor='#ffff66',LineWidth=1); % 散点动图
%% 速度变化幅度
vbound = zeros(size(xbound));
vbound(:,1) = -( xbound(:,2) - xbound(:,1) )*velo_ratio; % 速度下限
vbound(:,2) =  ( xbound(:,2) - xbound(:,1) )*velo_ratio; % 速度上限
%% 学习率 参数
c1u =2.5;c1l =0.5;
c2u =3.5;c2l =0.8;
%% 初始化位置和速度
[m,~]=size(xbound);
individualsPos = xbound(:,1) + ( xbound(:,2) - xbound(:,1) ).*rand(m,N) ; % 初始化位置 X_i,k
individualsVel = zeros(size(individualsPos)); % 初始化速度V_i,k
individualsP = func(individualsPos);% 计算个体适应度P_i,k
individualsPbest = individualsP; % 个体最优适应度p_best_i,k
[minv,ind] = min(individualsP); % 求解种群最优适应度g_best_i
groupBest = minv; % 种群最优适应度g_best_i
groupBestX = individualsPos(:,ind); % 种群最优位置g_ix
individualsBestX = individualsPos; % 个体最优位置p_xik
for i = 1:1:iternum
    %% 自适应学习因子《一种加权变异的粒子群算法》
    w = 0.8*exp(-2.5*i/iternum);
    c1 = c1u - rand()*(1-exp(-iternum/(iternum-i)))*(c1u-c1l);
    c2 = c2l + rand()*(1-exp(-iternum/(iternum-i)))*(c2u-c2l);
    % 速度更新公式
    individualsVel = w*individualsVel + ...
        c1*rand()*(individualsBestX - individualsPos) ...
        + c2*rand()*(groupBestX - individualsPos);% 速度更新
    % 速度约束ref function:constraint
    [individualsVel] = constraint(individualsVel,vbound,m,N);% 速度约束

    individualsPos = individualsPos + individualsVel; % 更新位置

    % 位置约束ref function:constraint
    [individualsPos] = constraint(individualsPos,xbound,m,N);% 位置约束

    % 继续
    individualsP = func(individualsPos); % 个体适应度
    for k=1:1:N
        if individualsP(k) < individualsPbest(k)
            individualsPbest(k) = individualsP(k); % 更新个体最优适应度
            individualsBestX(:,k) = individualsPos(:,k); % 更新个体最优位置
        end
    end
    [minv,ind] = min(func(individualsBestX)); % 计算种群最优适应度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 法一：小概率变异
%     if groupBest > minv && rand()<0.85 % 百分之90的变异概率，跳出局部最优
%             groupBest = minv; % 更新种群最优适应度
%             groupBestX = individualsBestX(:,ind);% 更新种群最优位置
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    groupBest = minv; % 更新种群最优适应度
    groupBestX = individualsBestX(:,ind);% 更新种群最优位置
    convergence(i) = groupBest; % 收敛曲线图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 法二：高斯震荡
%     if i>=20 && all(convergence(i-10)==convergence(i))
%         individualsPos = individualsPos.*randn(m,N);
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% 动图收敛查看
    set(scatterPlot,'XData',individualsBestX(1,:),'YData',individualsBestX(2,:));
    drawnow;
end
disp('全局最优解:');
disp(groupBestX);
disp('全局最优适应度:');
disp(groupBest);
%% 收敛曲线
plot(convergence,'b',LineWidth=1);
title('Convergence Curve');
%% 参考文献
%《一种加权变异的粒子群算法》
%
%% 测试
% f = @(x,y) 1+ 1/4000*(x.^2+y.^2) - cos(x).* cos(y/2);
% x = -10:0.5:10;
% y = -10:0.5:10;
% [xx,yy] = meshgrid(x,y);
% zz = f(xx,yy);
% surf(xx,yy,zz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% func = @(x)   -x(1,:).*sin(sqrt(abs(x(1,:)))) - x(2,:).*sin(sqrt(abs(x(2,:)))) ;
% func = @(x) 1+ 1/4000*(x(1,:).^2+x(2,:).^2) - cos(x(1,:)).* cos(x(2,:)/2);
% xbound = [-10,10;-10,10];
% iternum = 50;
% N = 50;
% velo_ratio = 0.10;
% [groupBest,groupBestX,convergence] = PSO(func,xbound,velo_ratio,iternum,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%