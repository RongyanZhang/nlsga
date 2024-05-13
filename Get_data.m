%% 从excel中读取数据
%输出R：                 订单
%输出B：                 公交车
%输出rB：                返程公交车
%输出D2：                距离矩阵
%输出station：           订单-站点信息表
%输出Region              订单终点分区信息
function [R, B, rB, D2, station, Region] = Get_data()
% R=xlsread('data\Request1_5.xlsx');       % 订单 
% R=xlsread('data\Request6_10.xlsx');      % 订单
R = xlsread('data\Request1.xlsx');     % 订单
P = xlsread('data\Taxi.xlsx');             % 出租车
% D=importdata('data\Distance1.xlsx');    % 距离矩阵
B = xlsread('data\Bus1.xlsx');             % 公交车
rB = xlsread('data\Bus2.xlsx');            % 公交车


%% 得到订单，公交车矩阵
[mr, ~] = size(R);    % 订单信息
for i = 1 : mr
    num = datevec(R(i, 6));
    num1= num(4) * 60 + num(5) + num(6) / 60;
    R(i, 6) = num1;
end

[mb1, nb1] = size(B);   % 公交车信息
for i = 1 : mb1
    for j = 4 : nb1
        num = datevec(B(i, j));
        num1 = num(4) * 60 + num(5) + num(6) / 60;
        B(i, j) = num1;
    end
end
[mb2, nb2] = size(rB);   % 返程公交车信息
for i = 1 : mb2     
    for j = 4 : nb2
        num = datevec(rB(i, j));
        num1= num(4) * 60 + num(5) + num(6) / 60;
        rB(i, j) = num1;
    end
end



%% 距离信息
RIDUS = 6372.8;    % km
All = [R(:, 2:3); R(:, 4:5); B(:, 2:3); P(:, 2:3)];
[m, ~] = size(All);
D1 = zeros(m, m);
D2 = zeros(m, m);

for i = 1 : m
    for j = 1 : m
        rand_number = 1.2 + 0.6 * rand;    % 1.2 - 1.8
        lon1 = All(i, 1);
        lat1 = All(i, 2);
        lon2 = All(j, 1);
        lat2 = All(j, 2);
    
        dlat = (lat2 - lat1) .* pi /180;
        dlon = (lon2 - lon1) .* pi /180;
        lat1 = lat1 .* pi /180;
        lat2 = lat2 .* pi /180;
        
        a = (sin(dlat./2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
        distance1 = 2 .* RIDUS .* asin(sqrt(a));
        distance2 = distance1 * rand_number;
        D1(i, j) = distance1;
        D2(i, j) = distance2;
    end
end

%% 得到订单-站点信息表
n_o = size(R, 1);   % 订单的个数
n_B = size(B, 1);  % 公交车的站点个数
% 五辆公交车
n_B1 = sum(B(:) == 101);
n_B2 = sum(B(:) == 102);
n_B3 = sum(B(:) == 103);
n_B4 = sum(B(:) == 104);
n_B5 = sum(B(:) == 105);

% 获得合适的公交车信息表
r = 6371e3;  % 地球半径
lat_ref = 39.92;   % 定义参考经纬度
lon_ref = 116.27;
lat_ref_rad = deg2rad(lat_ref);
lon_ref_rad = deg2rad(lon_ref);
radius = zeros(n_o, 1);
center = zeros(n_o, 2);

for i = 1 : n_o
    lat = [R(i, 2)
         R(i, 4)];
    lon = [R(i, 3)
         R(i, 5)];
    lat_rad = deg2rad(lat);   % 将经纬度转换为弧度
    lon_rad = deg2rad(lon);
    d = r * acos(sin(lat_rad) .* sin(lat_ref_rad) + cos(lat_rad) .* cos(lat_ref_rad) .* cos(lon_rad-lon_ref_rad));   % 计算两点之间的距离
    x =d .* cos(lat_rad) .* sin(lon_rad-lon_ref_rad);
    y =d .* (cos(lat_ref_rad) .* sin(lat_rad) - sin(lat_ref_rad) .* cos(lat_rad) .* cos(lon_rad - lon_ref_rad));
    radius(i, 1) = sqrt((x(1) - x(2))^2 + (y(1)-y(2))^2) / 2;
    center(i, :) = mean([x y]);
end


lat = [B(:, 2)];
lon = [B(:, 3)];
lat_rad = deg2rad(lat);   % 将经纬度转换为弧度
lon_rad = deg2rad(lon);
d = r * acos(sin(lat_rad) .* sin(lat_ref_rad) + cos(lat_rad) .* cos(lat_ref_rad) .* cos(lon_rad - lon_ref_rad));   % 计算两点之间的距离
x = d .* cos(lat_rad) .* sin(lon_rad - lon_ref_rad);
y = d .* (cos(lat_ref_rad) .* sin(lat_rad) - sin(lat_ref_rad) .* cos(lat_rad) .* cos(lon_rad - lon_ref_rad));

B1 = [x(1 : n_B1) y(1 : n_B1)];
B2 = [x(n_B1 + 1 : n_B1 + n_B2) y(n_B1 + 1 : n_B1 + n_B2)];
B3 = [x(n_B1 + n_B2 + 1 : n_B1 + n_B2 + n_B3) y(n_B1 + n_B2 + 1 : n_B1 + n_B2 + n_B3)];
B4 = [x(n_B1 + n_B2 + n_B3 + 1 : n_B1 + n_B2 + n_B3 + n_B4) y(n_B1 + n_B2 + n_B3 + 1 : n_B1 + n_B2 + n_B3 + n_B4)];
B5 = [x(n_B1 + n_B2 + n_B3 + n_B4 + 1 : n_B1 + n_B2 + n_B3 + n_B4 + n_B5) y(n_B1 + n_B2 + n_B3 + n_B4 + 1 : n_B1 + n_B2 + n_B3 + n_B4 + n_B5)];

count = zeros(n_o, 5);
for i = 1 : n_o
    for j = 1 : n_B1
        distance = norm(B1(j, :) - center(i, :));
        if distance <= radius(i, 1)
            count(i, 1) = count(i, 1) + 1;
        end
    end
    for j = 1 : n_B2
        distance = norm(B2(j, :) - center(i, :));
        if distance <= radius(i, 1)
            count(i, 2) = count(i, 2) + 1;
        end
    end
    for j = 1 : n_B3
        distance = norm(B3(j, :) - center(i, :));
        if distance <= radius(i, 1)
            count(i, 3) = count(i, 3) + 1;
        end
    end
    for j = 1 : n_B4
        distance = norm(B4(j, :) - center(i, :));
        if distance <= radius(i, 1)
            count(i, 4) = count(i, 4) + 1;
        end
    end
    for j = 1 : n_B5
        distance = norm(B5(j, :) - center(i, :));
        if distance <= radius(i, 1)
            count(i, 5) = count(i, 5) + 1;
        end
    end

end

% 构造站点结构体
station.Don1 = [];
station.Doff1 = [];
station.Don2 = [];
station.Doff2 = [];
station.Don3 = [];
station.Doff3 = [];
station.Don4 = [];
station.Doff4 = [];
station.Don5 = [];
station.Doff5 = [];
station.D_on1 = [];
station.D_off1 = [];
station.D_on2 = [];
station.D_off2 = [];
station.D_on3 = [];
station.D_off3 = [];
station.D_on4 = [];
station.D_off4 = [];
station.D_on5 = [];
station.D_off5 = [];
station.index = [];
station = repmat(station, n_o, 1);    % 复制和平铺矩阵

for i = 1 : n_o

    D_index = [];
    D_on0 = zeros(1, n_B);  
    D_off0 = zeros(1, n_B);  
    for k = 1 : n_B
        D_on0(k) = D2(i, n_o * 2 + k);                % 订单起点到公交车站点的距离
        D_off0(k) = D2(n_o * 2 + k, i + n_o);            % 订单公交车站点到终点的距离
    end
    if count(i,1) > 10
        station(i).D_on1 = D_on0(1 : n_B1);  
        station(i).D_off1 = D_off0(1 : n_B1);    % 公交车B1
        [~, index_on1] = sort(station(i).D_on1);
        [~, index_off1] = sort(station(i).D_off1);
        station(i).Don1 = index_on1(1 : 3);
        station(i).Doff1 = index_off1(1 : 3);   % 选取距离起点/终点最近的三个站点
        D_index = [D_index, 1];
    end
    
    if count(i,2) > 10
        station(i).D_on2 = D_on0(n_B1 + 1 : n_B1 + n_B2);  
        station(i).D_off2 = D_off0(n_B1 + 1 : n_B1 + n_B2);    % 公交车B2
        [~, index_on2] = sort(station(i).D_on2);
        [~, index_off2] = sort(station(i).D_off2);
        station(i).Don2 = index_on2(1 : 3);
        station(i).Doff2 = index_off2(1 : 3);   % 选取距离起点/终点最近的三个站点
        D_index = [D_index, 2];
    end
    if count(i,3) > 10
        station(i).D_on3 = D_on0(n_B1 + n_B2 + 1 : n_B1 + n_B2 + n_B3);  
        station(i).D_off3 = D_off0(n_B1 + n_B2 + 1 : n_B1 + n_B2 + n_B3);    % 公交车B3
        [~, index_on3] = sort(station(i).D_on3);
        [~, index_off3] = sort(station(i).D_off3);
        station(i).Don3 = index_on3(1 : 3);
        station(i).Doff3 = index_off3(1 : 3);   % 选取距离起点/终点最近的三个站点
        D_index = [D_index, 3];
    end
    if count(i, 4) > 10
        station(i).D_on4 = D_on0(n_B1 + n_B2 + n_B3 + 1 : n_B1 + n_B2 + n_B3 + n_B4);  
        station(i).D_off4 = D_off0(n_B1 + n_B2 + n_B3 + 1 : n_B1 + n_B2 + n_B3 + n_B4);    % 公交车B4
        [~, index_on4] = sort(station(i).D_on4);
        [~, index_off4] = sort(station(i).D_off4);
        station(i).Don4 = index_on4(1 : 3);
        station(i).Doff4 = index_off4(1 : 3);   % 选取距离起点/终点最近的三个站点
        D_index = [D_index, 4];
    end
    if count(i,5)>10
        station(i).D_on5 = D_on0(n_B1 + n_B2 + n_B3 + n_B4 + 1 : n_B1 + n_B2 + n_B3 + n_B4 + n_B5);  
        station(i).D_off5 = D_off0(n_B1 + n_B2 + n_B3 + n_B4 + 1 : n_B1 + n_B2 + n_B3 + n_B4 + n_B5);    % 公交车B5
        [~, index_on5] = sort(station(i).D_on5);
        [~, index_off5] = sort(station(i).D_off5);
        station(i).Don5 = index_on5(1 : 3);
        station(i).Doff5 = index_off5(1 : 3);   % 选取距离起点/终点最近的三个站点
        D_index = [D_index,5];
    end

    station(i).index = D_index;
end

% 构建Region,判断订单起点和终点位于哪些区域
Re = R(:, 4 : 5);
Region = cell(9, 1);
region1 = [];
region2 = [];
region3 = [];
region4 = [];
region5 = [];
region6 = [];
region7 = [];
region8 = [];
region9 = [];

for i = 1 : n_o
    if (Re(i, 1) > 116.27 && Re(i, 1) < 116.343)  &&  (Re(i, 2) > 40.07 && Re(i, 2) < 40.13)
        region1 = [region1; i];
    elseif (Re(i, 1) > 116.343 && Re(i, 1) < 116.426)  &&  (Re(i, 2) > 40.07 && Re(i, 2) < 40.13)
        region2 = [region2; i];
    elseif (Re(i, 1) > 116.416 && Re(i, 1) < 116.49)  &&  (Re(i, 2) > 40.07 && Re(i, 2) < 40.13)
        region3 = [region3; i];

    elseif (Re(i, 1) > 116.27 && Re(i, 1) < 116.343)  &&  (Re(i, 2) > 40.01 && Re(i, 2) < 40.07)
        region4 = [region4; i];
    elseif (Re(i, 1) > 116.343 && Re(i, 1) < 116.416)  &&  (Re(i, 2) > 40.01 && Re(i, 2) < 40.07)
        region5 = [region5; i];
    elseif (Re(i, 1) > 116.416 && Re(i, 1) < 116.49)  &&  (Re(i, 2) > 40.01 && Re(i, 2) < 40.07)
        region6 = [region6; i];

    elseif (Re(i, 1) > 116.27 && Re(i, 1) < 116.343)  &&  (Re(i, 2) > 39.95 && Re(i, 2) < 40.01)
        region7 = [region7; i];
    elseif (Re(i, 1) > 116.343 && Re(i, 1) < 116.416)  &&  (Re(i, 2) > 39.95 && Re(i, 2) < 40.01)
        region8 = [region8; i];
    elseif (Re(i, 1) > 116.416 && Re(i, 1) < 116.49)  &&  (Re(i, 2) > 39.95 && Re(i, 2) < 40.01)
        region9 = [region9; i];
    end
end

Region{1} = region1;
Region{2} = region2;
Region{3} = region3;
Region{4} = region4;
Region{5} = region5;
Region{6} = region6;
Region{7} = region7;
Region{8} = region8;
Region{9} = region9;
end

