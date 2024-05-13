%% 贪婪算法求最短距离/路径
%输入L                 原路径所经过的点
%输入D                 距离矩阵
%输入n_o               订单个数
%输入n_B               公交车站点个数
%输出distance          最短距离
%输出path1             最短路径
function path1 = Greedy(L, D, n_o, n_B)
L1 = [];
L2 = [];
L3 = [];
L4 = [];
for k = 2 : length(L)
    if L(k) <= n_o
        L1 = [L1; L(k)];    % 起点
    elseif L(k) > n_o * 2 && mod(k, 2) ~= 0
        L2 = [L2; L(k)];    % 公交车上车点
    elseif L(k) > n_o * 2 && mod(k, 2) == 0
        L3 = [L3; L(k)];     % 公交车下车点
    elseif L(k) > n_o && L(k) <= n_o * 2
        L4 = [L4; L(k)];     % 终点
    
    end
end
% taxi = 80;
taxi = n_o;
park = 10;
best_o=L(1);
if mod((best_o - n_o * 2 - n_B), taxi) == 0
    kp = fix((best_o - n_o * 2 - n_B - 1) / taxi);
else
    kp = fix((best_o - n_o * 2 - n_B) / taxi);
end
best_o = [n_o * 2 + n_B + kp * park + 1 : n_o * 2 + n_B + kp * park + park]';

n_A = size(L1, 1);
n_B = size(L2, 1);
n_C = size(L3, 1);
n_D = size(L4, 1);

allLocations = [best_o; L1; L2; L3; L4]; % 将所有地点的坐标合并成一个矩阵
numLocations = size(allLocations, 1);    % 初始化
distance = zeros(numLocations, numLocations);
for i = 1 : numLocations
    for j = 1 : numLocations
        distance(i, j) = D(allLocations(i, 1), allLocations(j, 1));
    end
end
if n_A > 0   % 车-起点
    distance0 = distance(1 : 10, 11 : 10 + n_A);
    [~, linearIndex] = min(distance0(:));
    [rowIndex, ~] = ind2sub(size(distance0), linearIndex);   % 第一辆车
    visitedA = false(1, n_A);   % 将要遍历的点取0
    visitedB = false(1, n_B);   % 将要遍历的点取0
    visitedC = false(1, n_C);   % 将要遍历的点取0
    visitedD = false(1, n_D);   % 将要遍历的点取0
    path = zeros(1, 1 + n_A + n_B + n_C + n_D);    
    currentLocation = rowIndex; 
    path(1) = currentLocation;
else   % 车-下车点
    distance0 = distance(1 : 10, 11 + n_A + n_B + n_C : 10 + n_A + n_B + n_C + n_D);
    [~, linearIndex] = min(distance0(:));
    [rowIndex, ~] = ind2sub(size(distance0), linearIndex);   % 第一辆车
    visitedA = false(1, n_A);   % 将要遍历的点取0
    visitedB = false(1, n_B);   % 将要遍历的点取0
    visitedC = false(1, n_C);   % 将要遍历的点取0
    visitedD = false(1, n_D);   % 将要遍历的点取0
    
    path = zeros(1, 1 + n_A + n_B + n_C + n_D); 
    currentLocation = rowIndex; 
    path(1) = currentLocation;
end


% 贪婪算法
if n_A > 0
    for j = 1 : n_A
        % 选择下一个要访问的点
        unvisitedIndices = find(~visitedA);
        % 如果所有点都被访问过，则直接跳出循环
        if isempty(unvisitedIndices)
            break;
        end
        distance1 = distance(currentLocation, unvisitedIndices + 10);
        [~, minIndex] = min(distance1);
        nextIndex = unvisitedIndices(minIndex);
        
        % 更新路径和访问状态
        path(j + 1) = nextIndex + 10;
        visitedA(nextIndex) = true;
        currentLocation = nextIndex + 10;
    end
end
if n_B > 0
    for j = n_A + 1 : n_A + n_B

        % 选择下一个要访问的点
        unvisitedIndices = find(~visitedB);
        % 如果所有点都被访问过，则直接跳出循环
        if isempty(unvisitedIndices)
            break;
        end
        distance2 = distance(currentLocation, unvisitedIndices + 10 + n_A);
        [~, minIndex] = min(distance2);
        nextIndex = unvisitedIndices(minIndex);
        
        % 更新路径和访问状态
        path(j+1) = nextIndex+10 + n_A;
        visitedB(nextIndex) = true;
        currentLocation = nextIndex + 10 + n_A;
    end
end
if n_C > 0
    for j = n_A + n_B + 1 : n_A + n_B + n_C

        % 选择下一个要访问的点
        unvisitedIndices = find(~visitedC);
        % 如果所有点都被访问过，则直接跳出循环
        if isempty(unvisitedIndices)
            break;
        end
        distance3 = distance(currentLocation, unvisitedIndices + 10 + n_A + n_B);
        [~, minIndex] = min(distance3);
        nextIndex = unvisitedIndices(minIndex);
        
        % 更新路径和访问状态
        path(j + 1) = nextIndex + 10 + n_A + n_B;
        visitedC(nextIndex) = true;
        currentLocation = nextIndex + 10 + n_A + n_B;
    end
end
if n_D > 0
    for j=n_A+n_B+n_C+1:n_A+n_B+n_C+n_D

        % 选择下一个要访问的点
        unvisitedIndices = find(~visitedD);
        % 如果所有点都被访问过，则直接跳出循环
        if isempty(unvisitedIndices)
            break;
        end
        distance4 = distance(currentLocation,unvisitedIndices + 10 + n_A + n_B + n_C);
        [~, minIndex] = min(distance4);
        nextIndex = unvisitedIndices(minIndex);
        
        % 更新路径和访问状态
        path(j + 1) = nextIndex + 10 + n_A + n_B + n_C;
        visitedD(nextIndex) = true;
        currentLocation = nextIndex + 10 + n_A + n_B + n_C;
    end
end
allLocations = allLocations';
path1 = allLocations(path);
end

