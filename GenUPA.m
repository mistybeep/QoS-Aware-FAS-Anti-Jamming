function t = GenUPA(N, lambda)
    factors = [];
    for i = 1:sqrt(N)
        if mod(N, i) == 0
            factors = [factors; i, N/i];
        end
    end

    % 选择 N_L2 最小的分解（即 N_L1 最大的情况）
    if ~isempty(factors)
        [N1, N2] = deal(factors(end, 1), factors(end, 2));
    else
        % 如果 N_L 是质数，退化成 ULA（线性阵列）
        N1 = 1;
        N2 = N;
    end
    % x_pos = linspace(-(N1-1)*lambda/4, (N1-1)*lambda/4, N1); % x 方向
    % y_pos = linspace(-(N2-1)*lambda/4, (N2-1)*lambda/4, N2); % y 方向
    x_pos = linspace(0, (N1-1)*lambda/2, N1); % x 方向
    y_pos = linspace(0, (N2-1)*lambda/2, N2); % y 方向
    % 生成网格坐标
    [X, Y] = meshgrid(x_pos, y_pos);
    % 组合成 2×N_L 的矩阵 [x; y]
    t = [X(:)'; Y(:)'];
end