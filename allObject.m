%% 计算种群中每个个体的总价钱，总路程
%输入SelCh：               种群
%输入R：                   订单
%输入B：                   公交车
%输入D：                   距离矩阵
%输出functionvalue：       种群中每个个体的总价钱，总路程（，共享出租车的单位公里收益）
function functionvalue = allObject(SelCh, R, B, D)
NIND = size(SelCh, 3);  % 种群大小
functionvalue = zeros(NIND, 3);
for i = 1 : NIND
    [price, distance, profit] = chromObj(SelCh(:, :, i), R, B, D);
    functionvalue(i, :) = [price, distance, profit];
end
end