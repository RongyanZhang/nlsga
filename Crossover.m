%% 交叉操作
%输入SelCh：                 所有种群
%输入Pc：                    交叉概率
%输入R,B：                   订单，公交车
%输出SelCh：                 交叉后的种群
function SelCh = Crossover(SelCh, Pc, R, B)
[NSel1, ~, NSel3] = size(SelCh);
n_o = size(R, 1);          % 订单的个数
n_B = size(B, 1);          % 公交车站点的个数

for j = 1 : 2 : NSel3
    cross_Selch1 = SelCh(:, :, j);                    % 第j个进行交叉操作的个体
    cross_Selch2 = SelCh(:, :, j + 1);                 % 第j+1个进行交叉操作的个体
    if Pc >= rand    % 交叉概率Pc
        n1 = randi([1, NSel1]);
        n2 = randi([1, NSel1]);
        for i = min(n1, n2) : max(n1, n2)
            flag = 0;
            [number1_1, number1_2] = number(cross_Selch1, R, B);
            [number2_1, number2_2] = number(cross_Selch2, R, B);
            if cross_Selch2(i, 1) ~= 0 && cross_Selch2(i, 2) ~= 0 
                if number1_1(cross_Selch2(i, 1) - n_o * 2 - n_B) + R(i,size(R, 2)) > 5
                    flag = 1;
                end
            end
            if cross_Selch1(i, 1) ~= 0 && cross_Selch1(i, 2) ~= 0
                if number2_1(cross_Selch1(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                    flag = 1;
                end
            end

            if cross_Selch2(i, 4) ~= 0 
                if number1_2(cross_Selch2(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                    flag = 1;
                end
            end
            if cross_Selch1(i, 4) ~= 0
                if number2_2(cross_Selch1(i, 4) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                    flag = 1;
                end
            end
            
            if cross_Selch2(i, 2) == 0
                if number1_1(cross_Selch2(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5 || number1_2(cross_Selch2(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5
                    flag = 1;
                end
            end
            if cross_Selch1(i, 2) == 0
                if number2_1(cross_Selch1(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5 || number2_2(cross_Selch1(i, 1) - n_o * 2 - n_B) + R(i, size(R, 2)) > 5 
                    flag = 1;
                end
            end
            if flag ~= 1
                for k = 1 : 4
                    temp = cross_Selch1(i, k);
                    cross_Selch1(i, k) = cross_Selch2(i, k);
                    cross_Selch2(i, k) = temp;
                end
            end
        end
    end   
    SelCh(:, :, j) = cross_Selch1;                    % 更新第i个个体
    SelCh(:, :, j + 1) = cross_Selch2;                  % 更新第i+1个个体

end
end