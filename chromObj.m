%% 计算单个个体的总价钱，总路程
%输入chrom：       单个个体
%输入R：           订单
%输入B：           公交车
%输入D：           距离矩阵
%输出price：       单个个体的总价钱
%输出distance：    单个个体的总路程
%输出profit：      共享出租车的单位公里收益
function [price, distance, profit] = chromObj(chrom, R, B, D)
n_o = size(R,1);
n_B = size(B,1);
% taxi = 80;
taxi = n_o;
n_p = taxi * 9;

pj = 2;    % 出租车每一公里的价钱 2元
a = 0.1;   % 打折费
b = 2;    % 座位费
c = 2;    % 公交车费
l = 0.5; % 走路折扣费

%% 所有订单花费的金钱
[path, ~]=Path(chrom, R, B, D);
price1 = 0;
for i = 1 : n_o
    if chrom(i, 1) ~= 0
        path1 = path{chrom(i,1) - n_o * 2 - n_B};
        if chrom(i,2) == 0           % 出租车（起点-终点）
            Dis=0;
            for k = find(path1 == i) : find(path1 == i + n_o) - 1
                Dis = Dis + D(path1(k), path1(k+1));
            end
            price1=price1+D(i,i+n_o)*pj-a*(Dis-D(i,i+n_o))+(R(i,size(R,2))-1) * b;
            
        elseif  chrom(i,2) ~= 0    % 第一程为出租车（起点-公交车上车点）
            Dis = 0;
            for k=find(path1==i):find(path1==chrom(i,2))-1
                Dis=Dis+D(path1(k),path1(k+1));
            end
            price1=price1+D(i,chrom(i,2))*pj-a*(Dis-D(i,chrom(i,2)))+(R(i,size(R,2))-1) * b;
        end
    end
    if chrom(i,4)~=0 && chrom(i,3)~=0     % 第三程为出租车（公交车下车点-终点）
        path1=path{chrom(i,4)-n_o*2-n_B};
        Dis=0;
        for k=find(path1==chrom(i,3)):find(path1==i+n_o)-1
            Dis=Dis+D(path1(k),path1(k+1));
        end
        price1=price1+D(chrom(i,3),i+n_o)*pj-a*(Dis-D(chrom(i,3),i+n_o))+(R(i,size(R,2))-1) * b;
    end
end
price=price1;

for i = 1 : n_o
    if chrom(i, 2) ~= 0
        price = price + R(i,size(R,2)) * c;     % 总花费+公交车的花费
        if chrom(i, 1) == 0
            price = price - R(i,size(R,2)) * l * D(i, chrom(i,2));
        end
        if chrom(i, 4) == 0
            price = price - R(i,size(R,2)) * l * D(chrom(i,3), i + n_o);
        end
    end
    
end

%% 出租车平均利润和出租车的行车距离
distance = 0;
distance2 = 0;
for j = 1 : n_p 
    path1 = path{j};
    if size(path1, 2) > 1
        n = size(path1, 2);
        for i = 1 : length(path1) - 1
            distance = distance + D(path1(1, i), path1(1, i + 1));
        end
        for i = 2 : length(path1) - 1
            distance2 = distance2 + D(path1(1, i), path1(1, i + 1));
        end
        D2 = D(path1(1, n) , n_o * 2 + n_B + 1 : n_o * 2 + n_B + 90);
        distance = distance + min(D2);
    end   
end

profit = price1 / distance2;   % 没有计算出租车没载客的那段距离
end