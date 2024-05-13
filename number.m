%% 计算单个个体的路线规划、人数
%输入chrom：           单个个体
%输入R：               订单
%输入B：               公交车
%输出Number：          单个个体的人数
function [number1, number2] = number(chrom, R, B)
n_o = size(R, 1);          % 订单的个数       
n_B = size(B, 1);          % 公交车站点的个数
% taxi = 80;
taxi = n_o;
n_p = taxi * 9;         % 出租车的个数
number1 = zeros(1, n_p);
number2 = zeros(1, n_p);
for i = 1 : n_o
    if chrom(i, 1) ~= 0 && chrom(i, 2) ~= 0    % 第一程
        number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
    end
    if chrom(i, 2) == 0     % 第二程
        number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B)+R(i, size(R, 2));
        number2(chrom(i, 1) - n_o * 2 - n_B) = number2(chrom(i, 1) - n_o * 2 - n_B)+R(i, size(R, 2));
    end
    if chrom(i, 4) ~= 0   % 第三程
        number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));
    end
end
end

