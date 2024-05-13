%% 修复不满足约束条件的个体
% 输入chrom：          个体
% 输入R：              订单
% 输入B：              公交车
% 输入rB：             返程公交车
% 输入D：              距离矩阵
% 输入Region           订单终点分区信息
% 输出chromR：         修复的个体
function chromR = repair(chrom, R, B, rB, D, Region)
n_o = size(R, 1);   % 订单的个数 
n_B = size(B, 1);          % 公交车站点的个数
% taxi = 80;
taxi = n_o;
n_p = taxi * 9;   % 出租车的个数
park = 10;

T2 = 610;      % 右时间窗口 （605 610 615 620 625)
vp = 0.5;      % 出租车的速度    30km/h=0.5km/minute
vr = 0.09;     % 人的速度     1.5m/s(5.4km/h)=0.09km/minute
tp = 10;       % 10分钟内上出租车
tp2 = 10;      % 第三程公交车下车点最多等10分钟上出租车

tt = 0.5;      % 上下车时间

%% 判断第一程是否在10分钟内上车
[path, ~] = Path(chrom, R, B, D);
for i = 1 : n_o
    if chrom(i, 1) ~= 0
        path1 = path{chrom(i, 1) - n_o * 2 - n_B};
        t_table = [0];
        for k = 2 : size(path1, 2)
            t = D(path1(k - 1), path1(k)) / vp + tt;
            % t = D(path1(k-1),path1(k)) / vp;
            t_table = [t_table, t];   % 出租车到达每一个点的时刻
        end
        time_table = [];
        for b = 1 : length(t_table)
            time_table = [time_table, sum(t_table(1 : b))];
        end
        time_table = time_table + T2;     % 每辆出租车到达某个点的时刻
        ts = time_table(path1 == i);  % 找到上出租车的时刻
        if ts - T2 > tp
            for k = 1 : taxi
                if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                    chrom(i, 1) = n_o * 2 + n_B + taxi;
                    break
                end
            end

        end
    end

end

%% 判断是否超过规定时间
[path, mark] = Path(chrom, R, B, D);
J = [];
for j = n_o * 2 + n_B + 1 : n_o * 2 + n_B + n_p    % 每辆车的行驶轨迹
    if ismember(j, chrom) ~= 0
        J = [J, j];
    end
end
for e = 1 : size(J, 2)     % 每辆车的行驶轨迹
    j = J(e);
    path1 = path{j - n_o * 2 - n_B};
    t_table = [0];
    for k = 2 : size(path1, 2)
        
        % t = D(path1(k-1),path1(k))/vp  ;
        t = D(path1(k - 1), path1(k)) / vp + tt;
        t_table = [t_table, t];   % 出租车到达每一个点的时刻
    end
    time_table = [];
    for b = 1 : length(t_table)
        time_table = [time_table, sum(t_table(1 : b))];
    end
    time_table = time_table + T2;     % 每辆出租车到达某个点的时刻  

    if isempty(mark) == 0    % 不为空，存在下车点
        mark1 = mark([find(mark(:, 1) == j)], 2 : 3);
        if isempty(mark1) == 0    % 不为空，存在下车点

            for i = 1 : size(mark1, 1)
                k = find(path1 == mark1(i, 2));
                if length(k) > 1
                    k = k(2);
                end
                if mark1(i, 2) == chrom(mark1(i, 1), 3) && chrom(mark1(i, 1), 1) == 0  && chrom(mark1(i, 1), 4) == j    % 公交车+出租车     
                    % tr= R(mark1(i,1),6) + D(mark1(i,1),chrom(mark1(i,1),2)) /vr;   % 走到公交车上车点的时刻
                    tr = R(mark1(i, 1), 6) + D(mark1(i, 1), chrom(mark1(i, 1),2)) / vr + tt;   % 走到公交车上车点的时刻
                    u = chrom(mark1(i, 1), 2) - n_o * 2;
                    f = chrom(mark1(i, 1), 3) - n_o * 2;
                    if u - f < 0
                        B_t = B(u, 4 : size(B, 2)) - tr;    % 等公交车的时间
                        Bb = B_t(B_t >= 0);
                        if isempty(Bb) == 0   % 不为空
                            t1 = tr + Bb(1) + abs(B(f, 4) - B(u, 4)) + tt;   % 到达下车点的时刻
                        else
                            t1 = Inf;
                        end
                    else
                        B_t = rB(u, 4 : size(rB, 2)) - tr;    % 等公交车的时间
                        Bb = B_t(B_t >= 0);
                        if isempty(Bb) == 0   % 不为空
                            t1 = tr + Bb(1) + abs(rB(f, 4) - rB(u, 4)) + tt;   % 到达下车点的时刻
                        else
                            t1 = Inf;
                        end
                    end
                    if t1 < time_table(k)
                        if time_table(k) - t1 > tp2
                            for x = 1 : 9
                                if ismember(mark1(i, 1), Region{x}) ~= 0
                                    
                                    for y = 1 : taxi
                                        if ismember(n_o * 2 + n_B + (x-1) * taxi + y, chrom) == 0 % 该车未被安排行程
                                            if chrom(i, 4) ~= 0
                                                chrom(i, 4) = n_o * 2 + n_B + (x - 1) * taxi + y;
                                            end
                                            break
                                        end
                                    end
                                end
                            end
                            
                        end
                    else 
                        ty = t1 - time_table(k);    % 出租车在下车点等待延迟的时间
                        for c = k : size(time_table, 2)
                            time_table(c) = time_table(c) + ty;   % 加上延迟的时间
                        end
                    end

                    
                elseif mark1(i, 2) == chrom(mark1(i, 1), 3) && chrom(mark1(i, 1), 1) ~= 0  && chrom(mark1(i, 1),4) == j      % 出租车+公交车+出租车     
                    if chrom(mark1(i, 1), 1) == j                     % 是同一辆出租车
                        tk = find(path1 == chrom(mark1(i, 1), 2));
                        if length(tk) > 1
                            tk = tk(1);
                        end
                        tr = time_table(tk); 
                    else      % 不是同一辆出租车
                        L1 = [chrom(mark1(i, 1), 1)];  
                        for a = 1 : n_o
                            if chrom(a, 1) == chrom(mark1(i, 1), 1)
                                if chrom(a, 2) ~= 0
                                    L1 = [L1, a, chrom(a, 2)];  % 起点和公交车上车点
                                end
                            end
                        end
                        if length(L1) > 1
                            path11 = Greedy(L1, D, n_o, n_B);     % 最优路径
                        else
                            path11 = L1;
                        end              

                        t_table1 = [0];
                        for s = 2 : size(path11, 2)
                            % t = D(path11(s),path11(s-1))/vp ;
                            t = D(path11(s), path11(s - 1)) / vp + tt;
                            t_table1 = [t_table1, t];   % 出租车到达每一个点的时刻
                        end
                        time_table1 = [];
                        for b = 1 : length(t_table1)
                            time_table1 = [time_table1, sum(t_table1(1 : b))];
                        end
                        time_table1 = time_table1 + T2; 
                        tk = find(path11 == chrom(mark1(i, 1), 2) );
                        if length(tk) > 1
                            tk = tk(1);
                        end
                        tr = time_table1(tk);             % 找到到达上车点的时刻 
                    end

                    u = chrom(mark1(i, 1), 2) - n_o * 2;
                    f = chrom(mark1(i, 1), 3) - n_o * 2;
                    if u - f < 0
                        
                        B_t = B(u, 4 : size(B, 2)) - tr(1);    % 等公交车的时间
                        Bb = B_t(B_t >= 0);
                        if isempty(Bb)==0   % 不为空
                            % t1=tr+ Bb(1) + abs(B(f,4) - B(u,4));   % 到达下车点的时刻
                            t1 = tr(1) + Bb(1) + abs(B(f, 4) - B(u, 4)) + tt;   % 到达下车点的时刻
                        else
                            t1 = Inf;
                        end
                    else
                        B_t = rB(u, 4 : size(rB, 2)) - tr(1);    % 等公交车的时间
                        Bb = B_t(B_t >= 0);
                        if isempty(Bb)==0   % 不为空
                            t1 = tr(1) + Bb(1) + abs(rB(f, 4) - rB(u, 4)) + tt;   % 到达下车点的时刻
                        else
                            t1 = Inf;
                        end
                    end
                    if t1 < time_table(k)
                        if time_table(k) - t1 > tp2
                            for x = 1 : 9
                                if ismember(mark1(i, 1), Region{x}) ~= 0
                                    
                                    for y = 1 : taxi
                                        if ismember(n_o * 2 + n_B + (x-1) * taxi + y, chrom) == 0 % 该车未被安排行程
                                            if chrom(i, 4) ~= 0
                                                chrom(i, 4) = n_o * 2 + n_B + (x - 1) * taxi + y;
                                            end
                                            break
                                        end
                                    end
                                end
                            end
                            
                        end
                    else 
                        ty = t1 - time_table(k);    % 出租车在下车点等待延迟的时间
                        for c = k : size(time_table, 2)
                            time_table(c) = time_table(c) + ty;   % 加上延迟的时间
                        end
                    end    
                end  
            end
        end
    end


    for i = 1 : n_o
        time = 0;
        if chrom(i, 1) == j && chrom(i, 2) ~= 0 && chrom(i, 4) == 0   % 出租车+公交车
            tk = find(path1 == chrom(i, 2));
            if ~isempty(tk)
                tk = tk(1);
            end
            td = time_table(tk);  % 找到到达上车点的时刻
            time = time + (td - R(i, 6)) + D(chrom(i, 3), i + n_o) / vr + tt;  % 出租车时间（等出租车时间+出租车行驶时间）+ 公交车下车步行时间
            u = chrom(i, 2) - n_o * 2;
            f = chrom(i, 3) - n_o * 2;
            if u - f < 0
                B_t = B(u, 4 : size(B, 2)) - td;    % 等公交车的时间
                Bb = B_t(B_t >= 0);
                if isempty(Bb) == 0   % 不为空
                    time = time + Bb(1) + abs(B(f, 4) - B(u, 4)) + tt;    % 出租车运行时间+等公交车时间+公交车运行时间+步行时间
                    % time= time + Bb(1) + abs(B(f,4)-B(u,4));    % 出租车运行时间+等公交车时间+公交车运行时间+步行时间
                else
                    time = Inf;
                end
            else
                B_t = rB(u, 4 : size(rB, 2)) - td;    % 等公交车的时间
                Bb = B_t(B_t >= 0);
                if isempty(Bb) == 0   % 不为空
                    % time= time + Bb(1) + abs(rB(f,4)-rB(u,4));    % 出租车运行时间+等公交车时间+公交车运行时间+步行时间
                    time = time + Bb(1) + abs(rB(f, 4) - rB(u, 4)) + tt;    % 出租车运行时间+等公交车时间+公交车运行时间+步行时间
                else
                    time = Inf;
                end
            end
            
        elseif chrom(i, 1) == j  && chrom(i, 2) == 0 && chrom(i, 4) == 0   % 出租车
            time = time_table( path1 == i + n_o ) - R(i, 6);  % 出租车直接运送到终点的时间

        elseif chrom(i, 1) == 0  && chrom(i, 2) ~= 0 && chrom(i, 4) == j   % 公交车+出租车
            time = time_table( path1 == i + n_o ) - R(i, 6);   % 出租车直接运送到终点的时间
           
        elseif chrom(i, 1) ~= 0  && chrom(i,4) == j     % 出租车+公交车+出租车
            time = time_table( path1 == i + n_o ) - R(i, 6);  % 出租车终点时刻-订单发出时刻
        end

        % 判断是否超过规定时间
        if time > R(i, 7)            
            if chrom(i, 1) ~= 0 && chrom(i, 2) == 0    % 出租车
                for k = 1 : taxi
                    if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                        chrom(i, 1) = n_o * 2 + n_B + taxi + k;   % 出租车-------->出租车 
                        break
                    end
                end
                    
            else
                if chrom(i, 1) ~= 0 && chrom(i, 2) ~= 0 && chrom(i, 4) == 0      % 出租车+公交车
                    tr = time_table( path1 == chrom(i, 2) ) + tt;  % 到公交车上车点的时刻
                    u = chrom(i, 2) - n_o * 2;
                    f = chrom(i, 3) - n_o * 2;
                    if u - f < 0
                        B_t = B(u, 4 : size(B, 2)) - tr(1);    % 等公交车的时间
                        Bb = B_t(B_t >= 0);
                        if isempty(Bb) == 0   % 不为空
                            t1 = tr(1) + Bb(1) + abs(B(f, 4) - B(u, 4)) + tt;   % 到达下车点的时刻
                        else
                            t1 = Inf;
                        end
                    else
                        B_t = rB(u, 4 : size(rB, 2)) - tr(1);    % 等公交车的时间
                        Bb = B_t(B_t >= 0);
                        if isempty(Bb) == 0   % 不为空
                            t1 = tr(1) + Bb(1) + abs(rB(f, 4) - rB(u, 4)) + tt;   % 到达下车点的时刻
                        else
                            t1 = Inf;
                        end
                    end
                    for x = 1 : 9
                        if ismember(i, Region{x}) ~= 0
                            for y = 1 : taxi
                                if ismember(n_o * 2 + n_B + (x - 1) * taxi + y, chrom) == 0 % 该车未被安排行程
                                    chrom(i, 4) = n_o * 2 + n_B + (x - 1) * taxi + y;  % 出租车+公交车-------->出租车+公交车+出租车
                                    break
                                end
                            end
                            
                        end
                    end 
                    time2 = t1 - R(i, 6) + D(chrom(i, 3), i + n_o) / vp + tt;
                    if time2 > R(i, 7)     % 仍超过规定时间
                        for k = 1 : taxi
                            if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                                chrom(i, 1) = n_o * 2 + n_B + taxi + k; 
                                break
                            end
                        end 
                        
                        if chrom(i, 1) ~= 0  && chrom(i, 2) ~= 0   % 找到满足条件的车辆
%                             if chrom(i,1)>n_o*2+n_B  &&  chrom(i,1)<=n_o*2+n_B+taxi1
%                                 best_o=[n_o*2+n_B+1:n_o*2+n_B+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+taxi1  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2
%                                 best_o=[n_o*2+n_B+10*1+1:n_o*2+n_B+10*1+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+taxi1+taxi2  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2*2
%                                 best_o=[n_o*2+n_B+10*2+1:n_o*2+n_B+10*2+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+taxi1+taxi2*2  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2*3
%                                 best_o=[n_o*2+n_B+10*3+1:n_o*2+n_B+10*3+10]';
%                             end
%                             if chrom(i,1) > n_o * 2 + n_B  &&  chrom(i,1) <= n_o * 2 + n_B + taxi
%                                 best_o=[n_o * 2 + n_B + 1 : n_o * 2 + n_B + 10]';
%                             elseif chrom(i,1) > n_o * 2 + n_B + taxi  && chrom(i,1) <= n_o * 2 + n_B + 2 * taxi
%                                 best_o=[n_o*2+n_B+10*1+1:n_o*2+n_B+10*1+10]';
%                             elseif chrom(i,1) > n_o*2+n_B+ 2*taxi  && chrom(i,1)<=n_o*2+n_B+ 3*taxi
%                                 best_o=[n_o*2+n_B+10*2+1:n_o*2+n_B+10*2+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+3* taxi  && chrom(i,1)<=n_o*2+n_B+4* taxi
%                                 best_o=[n_o*2+n_B+10*3+1:n_o*2+n_B+10*3+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+4*taxi  && chrom(i,1)<=n_o*2+n_B+5*taxi
%                                 best_o=[n_o*2+n_B+10*4+1:n_o*2+n_B+10*4+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+5*taxi  && chrom(i,1)<=n_o*2+n_B+6*taxi
%                                 best_o=[n_o*2+n_B+10*5+1:n_o*2+n_B+10*5+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+6*taxi  && chrom(i,1)<=n_o*2+n_B+7*taxi
%                                 best_o=[n_o*2+n_B+10*6+1:n_o*2+n_B+10*6+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+7*taxi  && chrom(i,1)<=n_o*2+n_B+8*taxi
%                                 best_o=[n_o*2+n_B+10*7+1:n_o*2+n_B+10*7+10]';
%                             elseif chrom(i,1)>n_o*2+n_B+8*taxi  && chrom(i,1)<=n_o*2+n_B+9*taxi
%                                 best_o=[n_o*2+n_B+10*8+1:n_o*2+n_B+10*8+10]';
%                             end
                            if mod((chrom(i, 1) - n_o * 2 - n_B), taxi) == 0
                                kp = fix((chrom(i,1) - n_o * 2 - n_B - 1) / taxi);
                            else
                                kp = fix((chrom(i,1) - n_o * 2 - n_B) / taxi);
                            end
                            % kp = fix((chrom(i,1) - n_o * 2 - n_B) / taxi);
                            best_o = [n_o * 2 + n_B + kp * park + 1 : n_o * 2 + n_B + kp * park + park]';
                            distance0 = D(best_o, i);
                            [~, linearIndex] = min(distance0(:));
                            [rowIndex, ~] = ind2sub(size(distance0), linearIndex);   % 第一辆车
                            taxip = best_o(rowIndex);
                            tr = T2 + ( D(taxip, i) + D(i, chrom(i, 2)) ) / vp + tt;   % 到公交车上车点的时刻
                            % tr= T2 + ( D(chrom(i,1),i)+D(i,chrom(i,2)) ) /vp + tt;   % 到公交车上车点的时刻
                            u = chrom(i, 2) - n_o * 2;
                            f = chrom(i, 3) - n_o * 2;
                            if u - f < 0
                                B_t = B(u, 4 : size(B, 2)) - tr;    % 等公交车的时间
                                Bb = B_t(B_t >= 0);
                                if isempty(Bb) == 0   % 不为空
                                    t1 = tr + Bb(1) + abs(B(f, 4) - B(u, 4)) + tt;   % 到达下车点的时刻
                                else
                                    t1 = Inf;
                                end
                            else
                                B_t = rB(u, 4 : size(rB, 2)) - tr;    % 等公交车的时间
                                Bb = B_t(B_t >= 0);
                                if isempty(Bb) == 0   % 不为空
                                    t1 = tr + Bb(1) + abs(rB(f, 4) - rB(u, 4)) + tt;   % 到达下车点的时刻
                                else
                                    t1 = Inf;
                                end
                            end
                            time2 = t1 - R(i, 6) + D(chrom(i, 3), i + n_o) / vp + tt;   % 出租车+公交车+出租车-------->出租车+公交车+出租车
                        end
                        
                    end
                elseif chrom(i, 1) == 0 && chrom(i, 2) ~= 0 && chrom(i, 4) ~= 0    % 公交车+出租车
                    for k = 1 : taxi
                        if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                            chrom(i, 1) = n_o * 2 + n_B + taxi + k;   
                            break
                        end
                    end
                    
%                     if chrom(i,1)>n_o*2+n_B  &&  chrom(i,1)<=n_o*2+n_B+taxi1
%                         best_o=[n_o*2+n_B+1:n_o*2+n_B+10]';
%                     elseif chrom(i,1)>n_o*2+n_B+taxi1  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2
%                         best_o=[n_o*2+n_B+10*1+1:n_o*2+n_B+10*1+10]';
%                     elseif chrom(i,1)>n_o*2+n_B+taxi1+taxi2  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2*2
%                         best_o=[n_o*2+n_B+10*2+1:n_o*2+n_B+10*2+10]';
%                     elseif chrom(i,1)>n_o*2+n_B+taxi1+taxi2*2  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2*3
%                         best_o=[n_o*2+n_B+10*3+1:n_o*2+n_B+10*3+10]';
%                     end
%                         if chrom(i,1) > n_o * 2 + n_B  &&  chrom(i,1) <= n_o * 2 + n_B + taxi
%                             best_o=[n_o * 2 + n_B + 1 : n_o * 2 + n_B + 10]';
%                         elseif chrom(i,1) > n_o * 2 + n_B + taxi  && chrom(i,1) <= n_o * 2 + n_B + 2 * taxi
%                             best_o=[n_o*2+n_B+10*1+1:n_o*2+n_B+10*1+10]';
%                         elseif chrom(i,1) > n_o*2+n_B+ 2*taxi  && chrom(i,1)<=n_o*2+n_B+ 3*taxi
%                             best_o=[n_o*2+n_B+10*2+1:n_o*2+n_B+10*2+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+3* taxi  && chrom(i,1)<=n_o*2+n_B+4* taxi
%                             best_o=[n_o*2+n_B+10*3+1:n_o*2+n_B+10*3+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+4*taxi  && chrom(i,1)<=n_o*2+n_B+5*taxi
%                             best_o=[n_o*2+n_B+10*4+1:n_o*2+n_B+10*4+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+5*taxi  && chrom(i,1)<=n_o*2+n_B+6*taxi
%                             best_o=[n_o*2+n_B+10*5+1:n_o*2+n_B+10*5+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+6*taxi  && chrom(i,1)<=n_o*2+n_B+7*taxi
%                             best_o=[n_o*2+n_B+10*6+1:n_o*2+n_B+10*6+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+7*taxi  && chrom(i,1)<=n_o*2+n_B+8*taxi
%                             best_o=[n_o*2+n_B+10*7+1:n_o*2+n_B+10*7+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+8*taxi  && chrom(i,1)<=n_o*2+n_B+9*taxi
%                             best_o=[n_o*2+n_B+10*8+1:n_o*2+n_B+10*8+10]';
%                         end
                    if mod((chrom(i, 1) - n_o * 2 - n_B), taxi) == 0
                        kp = fix((chrom(i,1) - n_o * 2 - n_B - 1) / taxi);
                    else
                        kp = fix((chrom(i,1) - n_o * 2 - n_B) / taxi);
                    end
                    % kp = fix((chrom(i,1) - n_o * 2 - n_B) / taxi);
                    best_o = [n_o * 2 + n_B + kp * park + 1 : n_o * 2 + n_B + kp * park + park]';
                    distance0=D(best_o,i);
                    [~, linearIndex] = min(distance0(:));
                    [rowIndex, ~] = ind2sub(size(distance0), linearIndex);   % 第一辆车
                    taxip=best_o(rowIndex);
                    % tr= R(i,6) + ( D(chrom(i,1),i)+D(i,chrom(i,2)) ) /vp+tt;   % 到公交车上车点的时刻
                    tr = T2 + ( D(taxip, i) + D(i, chrom(i, 2)) ) / vp + tt;   % 到公交车上车点的时刻
                    u = chrom(i, 2) - n_o * 2;
                    f = chrom(i, 3) - n_o * 2;
                    if u - f < 0
                        B_t = B(u, 4 : size(B, 2)) - tr;    % 等公交车的时间
                        Bb = B_t(B_t >= 0);
                        if isempty(Bb) == 0   % 不为空
                            t1 = tr + Bb(1) + abs(B(f, 4) - B(u, 4)) + tt;   % 到达下车点的时刻
                        else
                            t1 = Inf;
                        end
                    else
                        B_t = rB(u, 4 : size(rB, 2)) - tr;    % 等公交车的时间
                        Bb = B_t(B_t >= 0);
                        if isempty(Bb) == 0   % 不为空
                            t1 = tr + Bb(1) + abs(rB(f, 4) - rB(u, 4)) + tt;   % 到达下车点的时刻
                        else
                            t1 = Inf;
                        end 
                    end
                    for x = 1 : 9
                        if ismember(i, Region{x}) ~= 0                                
                            for y = 1 : taxi
                                if ismember(n_o * 2 + n_B + (x - 1) * taxi + y, chrom) == 0 % 该车未被安排行程
                                    chrom(i, 4) = n_o * 2 + n_B + (x - 1) * taxi + y;  % 公交车+出租车-------->出租车+公交车+出租车
                                    break
                                end
                            end
                            
                        end
                    end
                    time2 = t1 - R(i, 6) + D(chrom(i, 3), i + n_o) / vp + tt;
                elseif chrom(i, 1) ~= 0 && chrom(i, 2) ~= 0  && chrom(i, 4) ~= 0   % 出租车+公交车+出租车                   
                    for k = 1 : taxi
                        if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                            chrom(i, 1) = n_o * 2 + n_B + taxi + k; 
                            break
                        end
                    end 
                    
                    if chrom(i, 1) ~= 0  && chrom(i, 2) ~= 0   % 找到满足条件的车辆
%                         if chrom(i,1)>n_o*2+n_B  &&  chrom(i,1)<=n_o*2+n_B+taxi1
%                             best_o=[n_o*2+n_B+1:n_o*2+n_B+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+taxi1  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2
%                             best_o=[n_o*2+n_B+10*1+1:n_o*2+n_B+10*1+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+taxi1+taxi2  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2*2
%                             best_o=[n_o*2+n_B+10*2+1:n_o*2+n_B+10*2+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+taxi1+taxi2*2  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2*3
%                             best_o=[n_o*2+n_B+10*3+1:n_o*2+n_B+10*3+10]';
%                         end
%                         if chrom(i,1) > n_o * 2 + n_B  &&  chrom(i,1) <= n_o * 2 + n_B + taxi
%                             best_o=[n_o * 2 + n_B + 1 : n_o * 2 + n_B + 10]';
%                         elseif chrom(i,1) > n_o * 2 + n_B + taxi  && chrom(i,1) <= n_o * 2 + n_B + 2 * taxi
%                             best_o=[n_o*2+n_B+10*1+1:n_o*2+n_B+10*1+10]';
%                         elseif chrom(i,1) > n_o*2+n_B+ 2*taxi  && chrom(i,1)<=n_o*2+n_B+ 3*taxi
%                             best_o=[n_o*2+n_B+10*2+1:n_o*2+n_B+10*2+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+3* taxi  && chrom(i,1)<=n_o*2+n_B+4* taxi
%                             best_o=[n_o*2+n_B+10*3+1:n_o*2+n_B+10*3+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+4*taxi  && chrom(i,1)<=n_o*2+n_B+5*taxi
%                             best_o=[n_o*2+n_B+10*4+1:n_o*2+n_B+10*4+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+5*taxi  && chrom(i,1)<=n_o*2+n_B+6*taxi
%                             best_o=[n_o*2+n_B+10*5+1:n_o*2+n_B+10*5+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+6*taxi  && chrom(i,1)<=n_o*2+n_B+7*taxi
%                             best_o=[n_o*2+n_B+10*6+1:n_o*2+n_B+10*6+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+7*taxi  && chrom(i,1)<=n_o*2+n_B+8*taxi
%                             best_o=[n_o*2+n_B+10*7+1:n_o*2+n_B+10*7+10]';
%                         elseif chrom(i,1)>n_o*2+n_B+8*taxi  && chrom(i,1)<=n_o*2+n_B+9*taxi
%                             best_o=[n_o*2+n_B+10*8+1:n_o*2+n_B+10*8+10]';
%                         end
                        if mod((chrom(i, 1) - n_o * 2 - n_B), taxi) == 0
                            kp = fix((chrom(i,1) - n_o * 2 - n_B - 1) / taxi);
                        else
                            kp = fix((chrom(i,1) - n_o * 2 - n_B) / taxi);
                        end
                        % kp = fix((chrom(i,1) - n_o * 2 - n_B) / taxi);
                        best_o = [n_o * 2 + n_B + kp * park + 1 : n_o * 2 + n_B + kp * park + park]';
                        distance0 = D(best_o, i);
                        [~, linearIndex] = min(distance0(:));
                        [rowIndex, ~] = ind2sub(size(distance0), linearIndex);   % 第一辆车
                        taxip = best_o(rowIndex);
                        tr = T2 + ( D(taxip, i) + D(i, chrom(i, 2)) ) / vp + tt;   % 到公交车上车点的时刻
                        % tr= T2 + ( D(chrom(i,1),i)+D(i,chrom(i,2)) ) /vp + tt;   % 到公交车上车点的时刻
                        u = chrom(i, 2) - n_o * 2;
                        f = chrom(i, 3) - n_o * 2;
                        if u - f < 0
                            B_t = B(u, 4 : size(B, 2)) - tr;    % 等公交车的时间
                            Bb = B_t(B_t >= 0);
                            if isempty(Bb) == 0   % 不为空
                                t1 = tr + Bb(1) + abs(B(f, 4) - B(u, 4)) + tt;   % 到达下车点的时刻
                            else
                                t1 = Inf;
                            end
                        else
                            B_t = rB(u, 4 : size(rB, 2)) - tr;    % 等公交车的时间
                            Bb = B_t(B_t >= 0);
                            if isempty(Bb) == 0   % 不为空
                                t1 = tr + Bb(1) + abs(rB(f, 4) - rB(u, 4)) + tt;   % 到达下车点的时刻
                            else
                                t1 = Inf;
                            end
                        end
                        for x = 1 : 9
                            if ismember(i, Region{x}) ~= 0                                
                                for y = 1 : taxi
                                    if ismember(n_o * 2 + n_B + (x - 1) * taxi + y, chrom) == 0 % 该车未被安排行程
                                        chrom(i, 4) = n_o * 2 + n_B + (x - 1) * taxi + y;  % 出租车+公交车/公交车+出租车/出租车+公交车+出租车-------->出租车+公交车+出租车
                                        break
                                    end
                                end
                                
                            end
                        end
                        time2 = t1 - R(i, 6) + D(chrom(i, 3), i + n_o) / vp + tt;
                    end
                end
                if time2 > R(i, 7)                            
                    for k = 1 : taxi
                        if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                            chrom(i, 1) = n_o * 2 + n_B + taxi + k;  
                            break
                        end
                    end 
                    chrom(i, 2) = 0;
                    chrom(i, 3) = 0;
                    chrom(i, 4) = 0;     % 出租车+公交车+出租车-------->出租车
                end
            end  
        end
    end
end

% 公交车
for i = 1 : n_o
    time = 0;
    if chrom(i, 1) == 0  && chrom(i, 4) == 0 && chrom(i, 2) ~= 0     % 公交车
        tr = R(i, 6) + D(i, chrom(i, 2)) / vr + tt;   % 走到公交车上车点的时刻

        time = D(i, chrom(i, 2)) / vr + D(chrom(i, 3), i + n_o) / vr;  %  步行时间

        u = chrom(i, 2) - n_o * 2;
        f = chrom(i, 3) - n_o * 2;
        if u - f < 0
            B_t = B(u, 4 : size(B, 2)) - tr;    % 等公交车的时间
            Bb = B_t(B_t >= 0);
            if isempty(Bb) == 0   % 不为空
                time = time + Bb(1) + abs( B(f, 4) - B(u, 4)) + tt;    % 等公交车时间+公交车行驶时间
            else
                time = Inf;
            end            
        else
            B_t = rB(u, 4 : size(rB, 2)) - tr;    % 等公交车的时间
            Bb = B_t(B_t >= 0);
            if isempty(Bb) == 0   % 不为空
                time = time + Bb(1) + abs( rB(f, 4) - rB(u, 4)) + tt;    % 等公交车时间+公交车行驶时间
            else
                time = Inf;
            end
        end
    end

    if time > R(i, 7)        
        for k = 1 : taxi
            if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                chrom(i, 1) = n_o * 2 + n_B + taxi + k;  % 公交车-------->出租车+公交车
                break
            end
        end 
        

        if chrom(i, 1) ~= 0   % 找到满足条件的出租车
%             if chrom(i,1)>n_o*2+n_B  &&  chrom(i,1)<=n_o*2+n_B+taxi1
%                 best_o=[n_o*2+n_B+1:n_o*2+n_B+10]';
%             elseif chrom(i,1)>n_o*2+n_B+taxi1  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2
%                 best_o=[n_o*2+n_B+10*1+1:n_o*2+n_B+10*1+10]';
%             elseif chrom(i,1)>n_o*2+n_B+taxi1+taxi2  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2*2
%                 best_o=[n_o*2+n_B+10*2+1:n_o*2+n_B+10*2+10]';
%             elseif chrom(i,1)>n_o*2+n_B+taxi1+taxi2*2  && chrom(i,1)<=n_o*2+n_B+taxi1+taxi2*3
%                 best_o=[n_o*2+n_B+10*3+1:n_o*2+n_B+10*3+10]';
%             end
%             if chrom(i,1) > n_o * 2 + n_B  &&  chrom(i,1) <= n_o * 2 + n_B + taxi
%                 best_o=[n_o * 2 + n_B + 1 : n_o * 2 + n_B + 10]';
%             elseif chrom(i,1) > n_o * 2 + n_B + taxi  && chrom(i,1) <= n_o * 2 + n_B + 2 * taxi
%                 best_o=[n_o*2+n_B+10*1+1:n_o*2+n_B+10*1+10]';
%             elseif chrom(i,1) > n_o*2+n_B+ 2*taxi  && chrom(i,1)<=n_o*2+n_B+ 3*taxi
%                 best_o=[n_o*2+n_B+10*2+1:n_o*2+n_B+10*2+10]';
%             elseif chrom(i,1)>n_o*2+n_B+3* taxi  && chrom(i,1)<=n_o*2+n_B+4* taxi
%                 best_o=[n_o*2+n_B+10*3+1:n_o*2+n_B+10*3+10]';
%             elseif chrom(i,1)>n_o*2+n_B+4*taxi  && chrom(i,1)<=n_o*2+n_B+5*taxi
%                 best_o=[n_o*2+n_B+10*4+1:n_o*2+n_B+10*4+10]';
%             elseif chrom(i,1)>n_o*2+n_B+5*taxi  && chrom(i,1)<=n_o*2+n_B+6*taxi
%                 best_o=[n_o*2+n_B+10*5+1:n_o*2+n_B+10*5+10]';
%             elseif chrom(i,1)>n_o*2+n_B+6*taxi  && chrom(i,1)<=n_o*2+n_B+7*taxi
%                 best_o=[n_o*2+n_B+10*6+1:n_o*2+n_B+10*6+10]';
%             elseif chrom(i,1)>n_o*2+n_B+7*taxi  && chrom(i,1)<=n_o*2+n_B+8*taxi
%                 best_o=[n_o*2+n_B+10*7+1:n_o*2+n_B+10*7+10]';
%             elseif chrom(i,1)>n_o*2+n_B+8*taxi  && chrom(i,1)<=n_o*2+n_B+9*taxi
%                 best_o=[n_o*2+n_B+10*8+1:n_o*2+n_B+10*8+10]';
%             end
            if mod((chrom(i, 1) - n_o * 2 - n_B), taxi) == 0
                kp = fix((chrom(i,1) - n_o * 2 - n_B - 1) / taxi);
            else
                kp = fix((chrom(i,1) - n_o * 2 - n_B) / taxi);
            end
            % kp = fix((chrom(i,1) - n_o * 2 - n_B) / taxi);
            best_o = [n_o * 2 + n_B + kp * park + 1 : n_o * 2 + n_B + kp * park + park]';

            distance0 = D(best_o, i);
            [~, linearIndex] = min(distance0(:));
            [rowIndex, ~] = ind2sub(size(distance0), linearIndex);   % 第一辆车
            taxip = best_o(rowIndex);
            tr = T2 + ( D(taxip, i) + D(i, chrom(i, 2)) ) / vp + tt;   % 到公交车上车点的时刻
            % tr= T2 + ( D(chrom(i,1),i)+D(i,chrom(i,2)) ) /vp+tt;   % 到公交车上车点的时刻
            u = chrom(i, 2) - n_o * 2;
            f = chrom(i, 3) - n_o * 2;
            if u - f < 0
                B_t = B(u, 4 : size(B, 2)) - tr;    % 等公交车的时间
                Bb = B_t(B_t >= 0);
                if isempty(Bb) == 0   % 不为空
                    t1 = tr + Bb(1) + abs(B(f, 4) - B(u, 4)) + tt;   % 到达下车点的时刻
                else
                    t1 = Inf;
                end 
            else
                B_t = rB(u, 4 : size(rB, 2)) - tr;    % 等公交车的时间
                Bb = B_t(B_t >= 0);
                if isempty(Bb) == 0   % 不为空
                    t1 = tr + Bb(1) + abs(rB(f, 4) - rB(u, 4)) + tt;   % 到达下车点的时刻
                else
                    t1 = Inf;
                end
            end
            time2 = t1 - R(i, 6) + D(chrom(i, 3), i + n_o) / vr;

            if time2 > R(i, 7)
                for x = 1 : 9
                    if ismember(i, Region{x}) ~= 0                        
                        for y = 1 : taxi
                            if ismember(n_o * 2 + n_B + (x - 1) * taxi + y, chrom) == 0       % 该车未被安排行程
                                chrom(i, 4) = n_o * 2 + n_B + (x - 1) * taxi + y;    % 出租车+公交车-------->出租车+公交车+出租车
                                break
                            end
                        end
                        
                    end
                end                
                time2 = t1 - R(i, 6) + D(chrom(i, 3), i + n_o) / vp + tt;
                if time2 > R(i, 7)                    
                    for k = 1 : taxi
                        if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                            chrom(i, 1) = n_o * 2 + n_B + taxi + k;
                            break
                        end
                    end 
                    
                    chrom(i, 2) = 0;       
                    chrom(i, 3) = 0;
                    chrom(i, 4) = 0;     % 出租车+公交车+出租车-------->出租车
                end
            end
        else            
            for k = 1 : taxi
                if ismember(n_o * 2 + n_B + taxi + k, chrom) == 0 % 该车未被安排行程
                    chrom(i, 1) = n_o * 2 + n_B + taxi + k;
                    break
                end
            end 
            
            chrom(i, 2) = 0;       
            chrom(i, 3) = 0;
            chrom(i, 4) = 0;     % 出租车+公交车+出租车-------->出租车
        end        
    end
end
chromR = chrom;
end