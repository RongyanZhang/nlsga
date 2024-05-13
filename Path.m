%% 计算单个个体的路线规划、标记
%输入chrom：       单个个体
%输入R,B,D：       订单,公交车,距离矩阵
%输出Path：        单个个体的路线
%输出Mark：        标记
function [Path, Mark] = Path(chrom, R, B, D)
n_o = size(R,1);          % 订单的个数
n_B = size(B,1);          % 公交车站点的个数
% taxi = 80;
taxi = n_o;
n_p = taxi * 9;           % 出租车的个数         

Path = cell(n_p, 1);
Mark = [];
for j = n_o * 2 + n_B + 1 : n_o * 2 + n_B + n_p    % 每辆车的行驶轨迹
    L = [j];
    marki = [];
    markj = [];
    for i = 1 : n_o
        if chrom(i, 1) == j
            if chrom(i,2) ~= 0
                L = [L, i, chrom(i, 2)];   % 起点和公交车上车点
            else
                L = [L, i, i + n_o];        % 起点和终点
            end
        end
        if chrom(i, 4) == j
            L = [L, chrom(i,3), i+n_o];   % 公交车下车点和终点
            markj = [markj; j, i, chrom(i,3)];
        end
    end

    if length(L) > 1
        L_best = Greedy(L, D, n_o, n_B);     % 最优路径
    else
        L_best = L;
    end
    
    Path{j - n_o * 2 - n_B} = L_best;
    Mark = [Mark; marki];
end
end

