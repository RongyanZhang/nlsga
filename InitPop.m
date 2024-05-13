%% 初始化种群
%输入NIND：                 种群大小
%输入R：                    订单
%输入B，rB：                公交车
%输入D：                    距离矩阵
%输入station：              订单-站点信息表
%输入Region                 订单终点分区信息
%输出Chrom                  初始种群
function Chrom = InitPop(NIND, R, B, rB, D, station, Region)
n_o = size(R, 1);   % 订单的个数
Chrom = zeros(n_o, 4, NIND);   % 用于存储种群
%% 调用encode函数进行初始化
for i = 1 : NIND
    Chrom(:, :, i) = encode(R, B, rB, D, station, Region);     % 随机生成初始种群
end
end
