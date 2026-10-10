%% 注水算法 Water-filling Algorithm
clear; clc; close all;

% ===== 参数设置 =====
N = 8;               % 子信道数目
P_total = 10;        % 总发射功率
N0 = 1;              % 噪声功率谱密度
% 信道增益 |H_k|^2
H_sq = [0.2, 0.5, 1.2, 2.0, 0.8, 3.0, 0.1, 1.5]; 

% 噪声基底: noise_floor(k) = N0 / |H_k|^2
noise_floor = N0 ./ H_sq;

% ===== 注水算法核心：搜索拉格朗日乘子 λ（水面高度）=====
% 目标：sum( max(1/λ - noise_floor, 0) ) = P_total
% 采用二分法求解λ
lambda_low = 1/(max(noise_floor)+0.01);
lambda_high = 1/min(noise_floor);
tol = 1e-6;    % 迭代精度
max_iter = 100;

for iter = 1:max_iter
    lambda_mid = (lambda_low + lambda_high)/2;
    Pk = max(1/lambda_mid - noise_floor, 0);
    P_sum = sum(Pk);

    if abs(P_sum - P_total) < tol
        break;
    elseif P_sum > P_total
        lambda_low = lambda_mid;
    else
        lambda_high = lambda_mid;
    end
end
lambda_opt = lambda_mid;
Pk = max(1/lambda_opt - noise_floor, 0);
C = sum( log2(1 + Pk.*H_sq / N0) ); % 总信道容量 bit/s/Hz

%% 输出结果
disp('子信道序号：'); disp(1:N);
disp('信道增益 |H|^2：'); disp(H_sq);
disp('噪声基底：'); disp(noise_floor);
disp('各子信道分配功率 Pk：'); disp(Pk);
disp(['总分配功率：',num2str(sum(Pk))]);
disp(['系统总容量：',num2str(C),' bit/s/Hz']);
disp(['最优注水线 1/λ：',num2str(1/lambda_opt)]);

%% 绘图：注水示意图
figure('Color','w');
bar(1:N, noise_floor,'FaceColor',[0.7 0.7 0.7]); hold on;
bar(1:N, Pk,'FaceColor',[0.2 0.6 0.8]);
yline(1/lambda_opt,'r--','LineWidth',1.5);
text(N-1,1/lambda_opt+0.05,'注水线 1/\lambda','Color','r');
xlabel('子信道编号'); ylabel('功率/噪声基底');
legend('噪声基底','分配功率','注水线');
title('注水算法功率分配');
grid on;
