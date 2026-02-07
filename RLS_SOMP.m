function [pos_array, aug_x] = RLS_SOMP(y, T_Mat, G, L, gamma)

    [m, ~] = size(y);
    s = L; %  测量值维数
    
    % hat_x=zeros(G,n); %  待重构的谱域(变换域)向量
    Aug_t = [];        %  增量矩阵(初始值为空矩阵)
    r_n = y;  %  残差值
    pos_array = zeros(s, 1);
    TT_Mat = T_Mat;
    product = zeros(G, 1);
    for n = 1 : s %  迭代次数(稀疏度是测量的1/4)
        pro = TT_Mat' * r_n;
        for i = 1 : G
            product(i) = sum(abs(pro(i, :)));
        end
        
        [~, pos] = max(product);   %最大投影系数对应的位置
        Aug_t = [Aug_t, TT_Mat(:, pos)];   %矩阵扩充
        TT_Mat(:,pos) = zeros(m, 1); %选中的列置零
        aug_x = (Aug_t' * Aug_t + gamma * eye(n)) \ Aug_t' * y;
        %   aug_x=pinv(Aug_t)*y;
        r_n = y - Aug_t * aug_x;   %残差

        pos_array(n) = pos;   %纪录最大投影系数的位置
    end

end