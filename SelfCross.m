%% 交叉操作2
%输入SelCh：                   所有种群
%输入Pc：                      交叉概率
%输入R,B：                     订单公交车                 
%输出SelCh：                   交叉后的种群
function SelCh = SelfCross(SelCh, Pc, R, B)
[NSel1, ~, NSel3] = size(SelCh);
n_o = size(R, 1);          % 订单的个数
n_B = size(B, 1);          % 公交车站点的个数

% 五辆公交车
n_B1 = sum(B(:) == 101);
n_B2 = sum(B(:) == 102);
n_B3 = sum(B(:) == 103);
n_B4 = sum(B(:) == 104);
n_B5 = sum(B(:) == 105);
B1_stop = (1 : n_B1) + n_o * 2;                               % 公交车B1的站点
B2_stop = (1 : n_B2) + n_o * 2 + n_B1;                          % 公交车B2的站点
B3_stop = (1 : n_B3) + n_o * 2 + n_B1 + n_B2;                     % 公交车B3的站点
B4_stop = (1 : n_B4) + n_o * 2 + n_B1 + n_B2 + n_B3;                % 公交车B4的站点
B5_stop = (1 : n_B5) + n_o * 2 + n_B1 + n_B2 + n_B3 + n_B4;           % 公交车B5的站点

cross = 0.5;

for j = 1 : NSel3
    cross_Selch = SelCh(:, :, j);                    % 第j个进行交叉操作的个体
    if Pc >= rand    % 交叉概率Pc2
        n1 = randi([1, NSel1]);
        n2 = randi([1, NSel1]);
        for i = min(n1, n2) : max(n1, n2) - 1
            for k = i + 1 : max(n1, n2)
                if cross_Selch(i, 2) == 0 && cross_Selch(k, 2) == 0    % 出租车
                    if cross_Selch(k, 1) ~= cross_Selch(i, 1)
                        [number1, ~] = number(cross_Selch, R, B);
                        if cross > rand
                            if number1(cross_Selch(i, 1) - n_o * 2 - n_B) + R(k, size(R, 2)) < 5
                                cross_Selch(k, 1) = cross_Selch(i, 1);
                            end
                        else
                            if number1(cross_Selch(k, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) < 5
                                cross_Selch(i, 1) = cross_Selch(k, 1);
                            end
                        end
                    end
                end

                if (ismember(cross_Selch(i, 2), B1_stop) == 1 && ismember(cross_Selch(k, 2), B1_stop) == 1)...
                || (ismember(cross_Selch(i, 2), B2_stop) == 1 && ismember(cross_Selch(k, 2), B2_stop) == 1)...
                || (ismember(cross_Selch(i, 2), B3_stop) == 1 && ismember(cross_Selch(k, 2), B3_stop) == 1)...
                || (ismember(cross_Selch(i, 2), B4_stop) == 1 && ismember(cross_Selch(k, 2), B4_stop) == 1)...
                || (ismember(cross_Selch(i, 2), B5_stop) == 1 && ismember(cross_Selch(k, 2), B5_stop) == 1)
                    if abs(cross_Selch(k, 2) - cross_Selch(i, 2)) < 3 && cross_Selch(i, 2) ~= 0 && cross_Selch(k, 2) ~= 0    
                        if cross >= rand    % 50%的概率交换
                            cross_Selch(k, 2) = cross_Selch(i, 2);
                        else
                            cross_Selch(i, 2) = cross_Selch(k, 2);
                        end
                    end
                    if abs(cross_Selch(k, 3) - cross_Selch(i, 3)) < 3 && cross_Selch(i, 3) ~= 0 && cross_Selch(k, 3) ~= 0 
                        if cross >= rand    % 50%的概率交换
                            cross_Selch(k, 3) = cross_Selch(i, 3);
                        else
                            cross_Selch(i, 3) = cross_Selch(k, 3);
                        end
                    end
                    if cross_Selch(i, 1) ~= 0 && cross_Selch(k, 1) ~= 0 && cross_Selch(i, 2) ~= 0 && cross_Selch(k, 2) ~= 0 &&  abs(cross_Selch(k, 2) - cross_Selch(i, 2)) < 3   % 上车点相同，出租车1改成一样的
                        [number1, ~] = number(cross_Selch, R, B); 
                        if cross >= rand    % 50%的概率交换
                            if number1(cross_Selch(i, 1) - n_o * 2 - n_B) + R(k, size(R, 2)) < 5
                                cross_Selch(k, 1) = cross_Selch(i, 1);
                                cross_Selch(k, 2) = cross_Selch(i, 2);
                            end
                        else
                            if number1(cross_Selch(k, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) < 5
                                cross_Selch(i, 1) = cross_Selch(k, 1);
                                cross_Selch(i, 2) = cross_Selch(k, 2);
                            end
                        end
                    end
                    if cross_Selch(i, 4) ~= 0 && cross_Selch(k, 4) ~= 0 && cross_Selch(i, 3) ~= 0 && cross_Selch(k, 3) ~= 0 && abs(cross_Selch(k, 3) - cross_Selch(i, 3)) < 3       % 下车点相同，出租车2改成一样的
                        [~, number2] = number(cross_Selch, R, B);
                        if cross >= rand    % 50%的概率交换
                            if number2(cross_Selch(i, 4) - n_o * 2 - n_B) + R(k, size(R, 2)) < 5
                                cross_Selch(k, 4) = cross_Selch(i, 4);
                                cross_Selch(k, 3) = cross_Selch(i, 3);
                            end
                        else
                            if number2(cross_Selch(k, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) < 5
                                cross_Selch(i, 4) = cross_Selch(k, 4);
                                cross_Selch(i, 3) = cross_Selch(k, 3);
                            end
                        end
                    end     
                end

            end
        end
        SelCh(:, :, j) = cross_Selch;                    % 更新第i个个体
    end
end
end