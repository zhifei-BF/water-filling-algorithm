%% 注水算法 + CVX 对比仿真，Capacity vs Transmit Power
clear; clc; close all;

%====参数设置====
N = 8;                  %子信道数量，沿用之前例子
H_sq = [0.2, 0.5, 1.2, 2.0, 0.8, 3.0, 0.1, 1.5]; %信道增益 |Hk|^2
N0 = 1;                 %噪声功率谱密度
noise_floor = N0 ./ H_sq;

% 发射功率dB范围：50 ~ 300 dB，和图中横坐标保持一致
Pt_dB_list = 50:50:300;
Pt_linear_list = 10.^(Pt_dB_list/10);  % dB转线性功率

C_waterfilling = zeros(size(Pt_dB_list)); %注水容量
C_cvx = zeros(size(Pt_dB_list));          %CVX容量

%% 循环遍历每个总发射功率
for idx = 1:length(Pt_dB_list)
    P_total = Pt_linear_list(idx);

    % ========== 方法1：二分法注水算法（Water Filling）==========
    lambda_low = 1/(max(noise_floor)+0.01);
    lambda_high = 1/min(noise_floor);
    tol = 1e-6;
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
    Pk_wf = max(1/lambda_opt - noise_floor, 0);
    C_waterfilling(idx) = sum( log2(1 + Pk_wf .* H_sq / N0) );

    % ========== 方法2：CVX凸优化求解同一问题 ==========
    cvx_begin quiet
    variable Pk_cvx(1,N) nonnegative;
    maximize sum( log2( 1 + Pk_cvx .* H_sq / N0 ) )
    subject to
    sum(Pk_cvx) <= P_total;
    cvx_end
    C_cvx(idx) = cvx_optval;
end

%% 绘图，和你截图样式完全对齐
figure('Color','w');
plot(Pt_dB_list, C_waterfilling, '-k','LineWidth',1.5); hold on; %Water Filling实线
plot(Pt_dB_list, C_cvx, 'k*','MarkerSize',8); %CVX用*标记

xlabel('Transmit Power in [dB]');
ylabel('Sum of Capacity');
title('Water Filling and CVX Capacity Curve');
legend('Water Filling','CVX','Location','best');
grid on;
