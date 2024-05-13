%% 变异操作
%输入SelCh：              所有种群
%输入Pm：                 变异概率
%输入R,B：                订单,公交车
%输出SelCh：              变异后的个体
function SelCh = Mutate(SelCh, Pm, R,B)
[NSel1, ~, NSel3] = size(SelCh);
n_o = size(R, 1);             % 订单的个数
n_B = size(B, 1);             % 公交车站点个数

cross = 0.5;
for j = 1 : NSel3
    mutate_Selch = SelCh(:, :, j);                    % 第i个进行变异操作的个体
    if Pm >= rand                                   % 变异概率Pm
        a = [1, 2];
        Prob = [0.3, 0.7];
        s = randsrc(1, 1, [a; Prob]);
        switch s
            case 1
                if rand > 0.5
                    sp = mutate_Selch(:, 1);
                    p1 = sp(randi(length(sp)));
                    p2 = sp(randi(length(sp)));
                    if p1 ~= p2 && p1 ~= 0 && p2 ~= 0
                        p_1 = mutate_Selch == p1;
                        mutate_Selch(p_1) = p2;
                        p_2 = mutate_Selch == p2;
                        mutate_Selch(p_2) = p1;
                    end
                else
                    sp = mutate_Selch(:, 4);
                    p1 = sp(randi(length(sp)));
                    p2 = sp(randi(length(sp)));
                    if p1 ~= p2 && p1 ~= 0 && p2 ~= 0
                        p_1 = mutate_Selch == p1;
                        mutate_Selch(p_1) = p2;
                        p_2 = mutate_Selch == p2;
                        mutate_Selch(p_2) = p1;
                    end
                end
            case 2
                n1 = randi([1, NSel1]);
                n2 = randi([1, NSel1]);
                for i = min(n1, n2) : max(n1, n2) - 1
                    for k = i + 1 : max(n1, n2)
                        if cross > rand
                           if mutate_Selch(i, 1) ~= 0 && mutate_Selch(k, 1) ~= 0
                                [number1, ~] = number(mutate_Selch, R, B);
                                if cross >= rand    % 50%的概率交换
                                    if number1(mutate_Selch(i, 1) - n_o * 2 - n_B) + R(k, size(R, 2)) < 5
                                        mutate_Selch(k, 1) = mutate_Selch(i, 1);
                                    end
                                else
                                    if number1(mutate_Selch(k, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) < 5
                                        mutate_Selch(i, 1) = mutate_Selch(k, 1);
                                    end
                                end
                           end
                        end
                        if cross > rand
                            if mutate_Selch(i, 4) ~= 0 && mutate_Selch(k, 4) ~= 0
                                [~, number2] = number(mutate_Selch, R, B);
                                if cross >= rand    % 50%的概率交换
                                    if number2(mutate_Selch(i, 4) - n_o * 2 - n_B) + R(k, size(R, 2)) < 5
                                        mutate_Selch(k, 4) = mutate_Selch(i, 4);
                                    end
                                else
                                    if number2(mutate_Selch(k, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) < 5
                                        mutate_Selch(i, 4) = mutate_Selch(k, 4);
                                    end
                                end
                            end 
                        end
                    end
                end
        end 
    end
    SelCh(:, :, j) = mutate_Selch;                    % 更新第i个个体
end
end