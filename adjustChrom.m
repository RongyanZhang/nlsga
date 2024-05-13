%% 调整种群染色体，将不满足约束条件的染色体进行调整
% 输入Chrom：      种群
% 输入R：          订单
% 输入B，rB：      公交车
% 输入D：          距离矩阵
% 输入Region       订单终点分区信息
% 输出newChrom：   调整后的染色体，全部满足约束条件
function newChrom = adjustChrom(Chrom, R, B, rB, D, Region)
NIND = size(Chrom, 3);   %NIND种群大小
newChrom = Chrom;
for i = 1 : NIND
    newChrom(:, :, i) = repair(Chrom(:, :, i), R, B, rB, D, Region);
end
end

