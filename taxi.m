clear
clc

R = xlsread('data\request56_00.xlsx');     % 订单

%% 距离信息
RIDUS = 6372.8;    % km
All = [R(:, 2:3), R(:, 4:5)];
[m, ~] = size(All);
D2 = zeros(m, 1);

for i = 1 : m
    rand_number = 1.2 + 0.6 * rand;    % 1.2 - 1.8
    lon1 = All(i, 1);
    lat1 = All(i, 2);
    lon2 = All(i, 3);
    lat2 = All(i, 4);

    dlat = (lat2 - lat1) .* pi /180;
    dlon = (lon2 - lon1) .* pi /180;
    lat1 = lat1 .* pi /180;
    lat2 = lat2 .* pi /180;
    
    a = (sin(dlat./2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
    distance1 = 2 .* RIDUS .* asin(sqrt(a));
    distance2 = distance1 * rand_number;
    D2(i, 1) = distance2;
end
A = sum(D2)
B = 2.3 * A

