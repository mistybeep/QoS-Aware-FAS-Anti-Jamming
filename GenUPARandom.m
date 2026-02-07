function t = GenUPARandom(N, lambda, rhoSize)
% 生成满足最小间距约束的随机稀疏 UPA 阵列
% 输入：
%   N        : 天线总数
%   lambda   : 波长
%   rhoSize  : 稀疏度参数（区域缩放因子）
% 输出：
%   t (2×N)  : [x; y] 坐标

    % === 参数设置 ===
    d_min = lambda / 2;  % 最小间距（半波长）
    L = lambda / 2 * rhoSize;  % 正方形半边长 → 区域 [-L, L] × [-L, L]
    maxAttempts = 1000;   % 单点最大尝试次数（防死循环）
    safetyFactor = 1.05;  % 实际用 d_min * safetyFactor 避免浮点误差

    % 理论最大点数估计（六边形密铺）→ 防止 N 过大导致无限循环
    area = (2*L)^2;
    hexDensity = 2 / (sqrt(3) * d_min^2);  % 最密排布密度
    N_max_est = floor(area * hexDensity * 0.9); % 90% 密度上限
    if N > N_max_est
        warning('Requested N=%d exceeds estimated maximum (%d) for rhoSize=%.2f. May fail or take long.', ...
                N, N_max_est, rhoSize);
    end

    % === 初始化 ===
    points = zeros(2, N);  % 预分配
    count = 0;             % 已生成点数

    % === 主循环：逐个生成满足约束的点 ===
    attempt = 0;
    while count < N && attempt < maxAttempts * N  % 总尝试上限
        % 1. 随机生成候选点
        % x_cand = (rand() * 2 - 1) * L;
        % y_cand = (rand() * 2 - 1) * L;
        cand = [(rand() * 2 - 1) * L; (rand() * 2 - 1) * L];

        % 2. 检查是否与已有所有点满足最小距离
        valid = true;
        if count > 0
            dists = vecnorm(points(:, 1:count) - cand, 2, 1);  % 欧氏距离
            if any(dists < d_min * safetyFactor)
                valid = false;
            end
        end

        % 3. 接受有效点
        if valid
            count = count + 1;
            points(:, count) = cand;
            attempt = 0;  % 重置尝试计数器（避免因连续失败提前退出）
        else
            attempt = attempt + 1;
        end
    end

    % === 输出结果 ===
    if count < N
        warning('Only %d/%d points generated (max attempts reached). Consider reducing N or increasing rhoSize.', count, N);
        t = points(:, 1:count);  % 返回实际生成的点
    else
        t = points;
    end
end