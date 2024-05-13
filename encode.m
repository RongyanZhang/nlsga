%% 编码，满足约束
%输入R：                    订单
%输入P：                    共享出租车
%输入B：                    公交车
%输入D：                    距离矩阵
%输出chrom                  初始个体
function chrom=encode(R,B,rB,D,station,Region)
n_o = size(R,1);   % 订单的个数
% n_p = size(P,1);   % 出租车的个数
n_B = size(B,1);   % 公交车站点的个数

% 五辆公交车
n_B1 = sum(B(:) == 101);
n_B2 = sum(B(:) == 102);
n_B3 = sum(B(:) == 103);
n_B4 = sum(B(:) == 104);

z = 1.0;        % 超参数（步行距离 1km）
% taxi = 80;
taxi = n_o;
n_p = taxi * 9;

%% 随机初始化个体
chrom = zeros(n_o, 4);   %个体初始化为n_o*4的0矩阵
number1 = zeros(1, n_p);
number2 = zeros(1, n_p);

for i = 1 : n_o
    Don1 = station(i).Don1;
    Doff1 = station(i).Doff1;
    Don2 = station(i).Don2;
    Doff2 = station(i).Doff2;
    Don3 = station(i).Don3;
    Doff3 = station(i).Doff3;
    Don4 = station(i).Don4;
    Doff4 = station(i).Doff4;
    Don5 = station(i).Don5;
    Doff5 = station(i).Doff5;
    D_on1 = station(i).D_on1;
    D_off1 = station(i).D_off1;
    D_on2 = station(i).D_on2;
    D_off2 = station(i).D_off2;
    D_on3 = station(i).D_on3;
    D_off3 = station(i).D_off3;
    D_on4 = station(i).D_on4;
    D_off4 = station(i).D_off4;
    D_on5 = station(i).D_on5;
    D_off5 = station(i).D_off5;
    D_index = station(i).index;
    
    if ~isempty(D_index)
        B_id = D_index(randi(length(D_index)));
    
        if B_id == 1      % 距离101公交车最近
            B_on1 = find(D_on1 <= z);
            B_off1 = find(D_off1 <= z);
    
            if isempty(B_on1) == 0   % 不为空
                if isempty(B_off1) == 0   % 不为空
                    get_on = B_on1(randi(length(B_on1))) + n_o * 2;  % 上车点
                    get_off = B_off1(randi(length(B_off1))) + n_o * 2;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                               % 公交车
                else
                    get_on = B_on1(randi(length(B_on1))) + n_o * 2;  % 上车点
                    get_off = Doff1(randi(length(Doff1))) + n_o * 2;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    
                    for j = 1 : 9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;                  
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o*2-n_B)+R(i,size(R,2));               % 公交车+出租车

                end
            else
                if isempty(B_off1) == 0   % 不为空
                    get_on = Don1(randi(length(Don1))) + n_o * 2;  % 上车点
                    get_off = B_off1(randi(length(B_off1))) + n_o * 2;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
    
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                                                       % 出租车+公交车
                else
                    get_on = Don1(randi(length(Don1))) + n_o * 2;  % 上车点
                    get_off = Doff1(randi(length(Doff1))) + n_o * 2;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
    
                    for j = 1 : 9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;                 
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));    % 出租车+公交车+出租车
                end
            end
    
        elseif B_id == 2      % 距离102公交车最近
            B_on2 = find(D_on2 <= z);
            B_off2 = find(D_off2 <= z);
    
            if isempty(B_on2) == 0   % 不为空
                if isempty(B_off2) == 0   % 不为空
                    get_on = B_on2(randi(length(B_on2))) + n_o * 2 + n_B1;  % 上车点
                    get_off = B_off2(randi(length(B_off2))) + n_o * 2 + n_B1;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                                                        % 公交车
                else
                    get_on = B_on2(randi(length(B_on2))) + n_o * 2 + n_B1;  % 上车点
                    get_off = Doff2(randi(length(Doff2))) + n_o * 2 + n_B1;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    
                    for j=1:9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;              
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));     % 公交车+出租车
                end
            else
                if isempty(B_off2) == 0   % 不为空
                    get_on = Don2(randi(length(Don2))) + n_o * 2 + n_B1;  % 上车点
                    get_off = B_off2(randi(length(B_off2))) + n_o * 2 + n_B1;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                                                      % 出租车+公交车
                else
                    get_on = Don2(randi(length(Don2))) + n_o * 2 + n_B1;  % 上车点
                    get_off = Doff2(randi(length(Doff2))) + n_o * 2 + n_B1;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    for j = 1 : 9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (j - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;                  
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (j - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));   % 出租车+公交车+出租车
                end
            end
    
        elseif B_id == 3      % 距离103公交车最近
            B_on3 = find(D_on3 <= z);
            B_off3 = find(D_off3 <= z);
    
            if isempty(B_on3) == 0   % 不为空
                if isempty(B_off3) == 0   % 不为空
                    get_on = B_on3(randi(length(B_on3))) + n_o * 2 + n_B1 + n_B2;  % 上车点
                    get_off = B_off3(randi(length(B_off3))) + n_o * 2 + n_B1 + n_B2;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                               % 公交车
                else
                    get_on = B_on3(randi(length(B_on3))) + n_o * 2 + n_B1 + n_B2;  % 上车点
                    get_off = Doff3(randi(length(Doff3))) + n_o * 2 + n_B1 + n_B2;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    for j = 1 : 9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;                 
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));     % 公交车+出租车
    
                end
            else
                if isempty(B_off3) == 0   % 不为空
                    get_on = Don3(randi(length(Don3))) + n_o * 2 + n_B1 + n_B2;  % 上车点
                    get_off = B_off3(randi(length(B_off3))) + n_o * 2 + n_B1 + n_B2;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                                                      % 出租车+公交车
                else
                    get_on = Don3(randi(length(Don3))) + n_o * 2 + n_B1 + n_B2;  % 上车点
                    get_off = Doff3(randi(length(Doff3))) + n_o * 2 + n_B1 + n_B2;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
                    
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    for j = 1 : 9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;                  
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));     % 出租车+公交车+出租车
                end
            end
    
        elseif B_id == 4      % 距离104公交车最近
            B_on4 = find(D_on4 <= z);
            B_off4 = find(D_off4 <= z);
    
            if isempty(B_on4) == 0   % 不为空
                if isempty(B_off4) == 0   % 不为空
                    get_on = B_on4(randi(length(B_on4))) + n_o * 2 + n_B1 + n_B2 + n_B3;  % 上车点
                    get_off = B_off4(randi(length(B_off4))) + n_o * 2 + n_B1 + n_B2 + n_B3;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                                                    % 公交车
                else
                    get_on = B_on4(randi(length(B_on4))) + n_o * 2 + n_B1 + n_B2 + n_B3;  % 上车点
                    get_off = Doff4(randi(length(Doff4))) + n_o * 2 + n_B1 + n_B2 + n_B3;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    for j = 1 : 9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;                  
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));      % 公交车+出租车
     
                end
            else
                if isempty(B_off4) == 0   % 不为空
                    get_on = Don4(randi(length(Don4))) + n_o * 2 + n_B1 + n_B2 + n_B3;  % 上车点
                    get_off = B_off4(randi(length(B_off4))) + n_o * 2 + n_B1 + n_B2 + n_B3;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                                               % 出租车+公交车
                else
                    get_on = Don4(randi(length(Don4))) + n_o * 2 + n_B1 + n_B2 + n_B3;  % 上车点
                    get_off = Doff4(randi(length(Doff4))) + n_o * 2 + n_B1 + n_B2 + n_B3;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B)+R(i,size(R,2));
                    
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    for j=1:9
                        if ismember(i,Region{j})~=0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi+randi([1,taxi])+n_o*2+n_B;            
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));    % 出租车+公交车+出租车
                end     
            end
           
        elseif B_id == 5      % 距离105公交车最近
            B_on5 = find(D_on5 <= z);
            B_off5 = find(D_off5 <= z);
    
            if isempty(B_on5) == 0   % 不为空
                if isempty(B_off5) == 0   % 不为空
                    get_on = B_on5(randi(length(B_on5))) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;  % 上车点
                    get_off = B_off5(randi(length(B_off5))) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                                                    % 公交车
                else
                    get_on = B_on5(randi(length(B_on5))) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;  % 上车点
                    get_off = Doff5(randi(length(Doff5))) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;  % 下车点
                    chrom(i, 1) = 0;
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    for j = 1 : 9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;                 
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i,4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));   % 公交车+出租车
    
                end
            else
                if isempty(B_off5) == 0   % 不为空
                    get_on = Don5(randi(length(Don5))) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;  % 上车点
                    get_off = B_off5(randi(length(B_off5))) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点
                    chrom(i, 4) = 0;                                                                   % 出租车+公交车
                else 
                    get_on = Don5(randi(length(Don5))) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;  % 上车点
                    get_off = Doff5(randi(length(Doff5))) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;  % 下车点
                    chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i,1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
                    end
                    number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i,1)-n_o*2-n_B)+R(i,size(R,2));
    
                    chrom(i, 2) = get_on;  % 上车点
                    chrom(i, 3) = get_off;  % 下车点

                    for j = 1 : 9
                        if ismember(i, Region{j}) ~= 0
                            taxip = j;
                            chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;                 
                        end
                    end
                    while number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                        chrom(i, 4) = (taxip - 1) * taxi + randi([1, taxi]) + n_o * 2 + n_B;
                    end
                    number2(chrom(i, 4) - n_o * 2 - n_B) = number2(chrom(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2));    % 出租车+公交车+出租车

                end
            end
        end

        % 判断是否符合条件 
        if abs(get_off - get_on) < 5      % 乘坐出租车大于5站
            
            chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
            while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5 || number2(chrom(i, 1)- n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
            end
            number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
            number2(chrom(i, 1) - n_o * 2 - n_B) = number2(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
            chrom(i, 2) = 0;  
            chrom(i, 3) = 0;
            chrom(i, 4) = 0;                                                      % 出租车
        elseif D(i, get_on) + D(get_off, i + n_o) > D(i, i + n_o)      % 出租车乘坐距离小于起始点距离
            chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
            while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5 || number2(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
            end
            number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
            number2(chrom(i, 1) - n_o * 2 - n_B) = number2(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
            chrom(i, 2) = 0;  
            chrom(i, 3) = 0;
            chrom(i, 4) = 0;                                                      % 出租车
        end
    else
        chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
        while number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5 || number2(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
            chrom(i, 1) = randi([1, taxi]) + n_o * 2 + n_B + taxi;
        end
        number1(chrom(i, 1) - n_o * 2 - n_B) = number1(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
        number2(chrom(i, 1) - n_o * 2 - n_B) = number2(chrom(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2));
        chrom(i, 2) = 0;  
        chrom(i, 3) = 0;
        chrom(i, 4) = 0;
    end
end
chrom = repair(chrom, R, B, rB, D, Region);
end
