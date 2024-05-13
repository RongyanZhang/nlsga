%% 计算拥挤距离
% 输入functionvalue                    目标函数值，总价钱、总距离
% 输入frontvalue                       每个个体的前沿面编号
% 输入fnum                             当前前沿面编号
% 输出distancevalue                    第fnum个面上的个体的拥挤距离
function distancevalue = CD_Assign(functionvalue, frontvalue, fnum)
popu = find(frontvalue == fnum);                                            %popu记录第fnum个面上的个体编号
distancevalue = zeros(size(popu));                                        %popu各个体的拥挤距离
fmax = max(functionvalue(popu, :), [], 1);                                   %popu每维上的最大值
fmin = min(functionvalue(popu, :), [], 1);                                   %popu每维上的最小值
for i = 1 : size(functionvalue, 2)                                           %分目标计算每个目标上popu各个体的拥挤距离
    [~, newsite] = sortrows(functionvalue(popu, i));
    distancevalue(newsite(1)) = inf;
    distancevalue(newsite(end)) = inf;
    for j = 2 : length(popu) - 1
        distancevalue(newsite(j)) = distancevalue(newsite(j)) + (functionvalue(popu(newsite(j + 1)), i) - functionvalue(popu(newsite(j - 1)), i)) / (fmax(i) - fmin(i));
    end
end
end

