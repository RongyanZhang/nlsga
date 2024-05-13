clear
clc
tic
%% 创建数据
[R, B, rB, D, station, Region] = Get_data();

%% 参数设置
NIND = 100;           % 种群大小
Pc = 0.6;             % 交叉概率
Pc2 = 0.5;
Pm = 0.05;            % 变异概率
MAXGEN = 500;         % 迭代次数
Price = zeros(1, MAXGEN);
Distance =zeros(1, MAXGEN);
Profit = zeros(1, MAXGEN);
up = 1;
ud = 1;

%% 初始种群
Population.Particle = InitPop(NIND, R, B, rB, D, station, Region);
Population.PopObj = allObject(Population.Particle, R, B, D);
                          
%% 迭代优化
gen=0;
while gen<MAXGEN
    %% 交叉操作
    NewPopulation.Particle = Crossover(Population.Particle, Pc, R, B);
    NewPopulation.Particle = SelfCross(NewPopulation.Particle, Pc2, R, B);
    %% 变异操作
    NewPopulation.Particle = Mutate(NewPopulation.Particle, Pm, R, B);
    %% 越界处理
    NewPopulation.Particle = adjustChrom(NewPopulation.Particle, R, B, rB, D, Region);
    NewPopulation.PopObj = allObject(NewPopulation.Particle, R, B, D);
    %% 种群合并
    NewPopulation.Particle = cat(3, Population.Particle, NewPopulation.Particle);
    NewPopulation.PopObj = [Population.PopObj;NewPopulation.PopObj];
    NewPopulation.frontvalue = Non_DS(NewPopulation.PopObj);
    %% 父代选择（基于线性排名的）
    % 参数迭代取值
    if mod(gen, 50) == 0 && gen ~= 0
        up = up + 0.2;
        ud = ud - 0.2;
        if up > 2 || ud < 0
            up = 2;
            ud = 0;
        end
    end
    Population = Select(NewPopulation, up, ud, NIND, gen);
    Obj = Population.PopObj;

    plot(Obj(:, 1), Obj(:, 2), 'r*');
    xlabel('总价钱');
    ylabel('总路程');
    pause(0.01)
    %% 更新迭代次数
    gen = gen + 1
    min_price = min(Obj(:, 1));
    Price(gen) = min_price;
    min_distance = min(Obj(:, 2));
    Distance(gen) = min_distance;
    max_profit = max(Obj(:, 3));
    Profit(gen) = max_profit;
end

Population.frontvalue = Non_DS(Population.PopObj);
mean_Obj = Population.PopObj(Population.frontvalue == 1, :);
num=size(mean_Obj, 1);  % num个帕累托为1的解
mean_Obj2 = Population.PopObj;

mean_price = sum(mean_Obj(:, 1)) / num;    % 帕累托为1的解
mean_price2 = sum(mean_Obj2(:, 1)) / NIND;
min_price;
mean_distance = sum(mean_Obj(:, 2)) / num;
mean_distance2 = sum(mean_Obj2(:, 2)) / NIND;
min_distance;
mean_Profit = sum(mean_Obj(:, 3)) / num;
mean_Profit2 = sum(mean_Obj2(:, 3)) / NIND;
max_profit;
Mean = zeros(9, 1);
Mean(1, 1) = mean_price;
Mean(2, 1) = mean_price2;
Mean(3, 1) = min_price;
Mean(4, 1) = mean_distance;
Mean(5, 1) = mean_distance2;
Mean(6, 1) = min_distance;
Mean(7, 1) = mean_Profit;
Mean(8, 1) = mean_Profit2;
Mean(9, 1) = max_profit;

figure
plot(Price)
xlabel('迭代次数')
ylabel('目标函数值')
title('总价钱进化曲线')
figure
plot(Distance)
xlabel('迭代次数')
ylabel('目标函数值')
title('总路程进化曲线')
figure
plot(Profit)
xlabel('迭代次数')
ylabel('目标函数值')
title('共享出租车单位公里收益进化曲线')

toc