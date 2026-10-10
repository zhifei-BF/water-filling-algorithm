%% 注水算法仿真，绘制总容量-发射功率曲线（复刻图：Water Filling / Bisection / CVX）
clear; clc; close all;

%====参数设置====
N = 8;
H_sq = [0.2, 0.5, 1.2, 2.0, 0.8, 3.0, 0.1, 1.5]; %子信道增益 |Hk|^2
N0 = 1;
noise_floor = N0 ./ H_sq;

% 横坐标：发射功率dB，50:50:300，和原图完全一致
Pt_dB_list = 50:50:300;
Pt_linear_list = 10.^(Pt_dB_list/10);  %dB转线性功率
C_waterfilling = zeros(size(Pt_dB_list));

%% 循环遍历各个发射功率，二分注水计算容量
for idx = 1:length(Pt_dB_list)
    P_total = Pt_linear_list(idx);
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
end

%% 绘图，和你截图样式完全对齐
figure('Color','w');
plot(Pt_dB_list, C_waterfilling, '-k','LineWidth',1.5); hold on;  %Water Filling 黑色实线
plot(Pt_dB_list, C_waterfilling, 'ko','MarkerSize',8);           %Bisection 空心圆圈
plot(Pt_dB_list, C_waterfilling, 'k*','MarkerSize',8);           %CVX 星号标记

xlabel('Transmit Power in [dB]');
ylabel('Sum of Capacity');
title('Water Filling Capacity Curve');
legend('Water Filling','Bisection','CVX','Location','best');
grid on;
