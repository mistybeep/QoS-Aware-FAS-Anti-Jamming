function coord = GenCoord(rhoSize, n, lambda)
    lambda = 10 * lambda;
    % 参数设置
    halfLambda = lambda / 2;  % 半波长
    maxBoundary = rhoSize * halfLambda;  % 正方形边界长度
    TOL = 1e-10;  % 提高数值稳定性

    % 初始猜测：均匀分布圆心，半径取边界长度的合理初始值
    x0 = zeros(2*n + 1, 1);
    gridSize = ceil(sqrt(n));  % 按网格初步分布
    spacing = 2*maxBoundary / (gridSize + 1);
    idx = 1;
    for i = 1:gridSize
        for j = 1:gridSize
            if idx <= n
                x0(2*idx-1) = -maxBoundary + j*spacing;  % x 坐标
                x0(2*idx) = -maxBoundary + i*spacing;    % y 坐标
                idx = idx + 1;
            end
        end
    end
    x0(end) = maxBoundary / (2*sqrt(n));  % 初始半径猜测

    % 目标函数：最大化半径（最小化负半径）
    obj_fun = @(z) -z(end);

    % 非线性约束
    nonlcon = @(z) nonlinear_constraints(z, n, maxBoundary, TOL);

    % 上下界
    lb = [-maxBoundary*ones(2*n, 1); 0];  % x, y 下界和 r >= 0
    ub = [maxBoundary*ones(2*n, 1); maxBoundary];  % r 有上界

    % 优化选项
    options = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 2000, ...
        'TolFun', 1e-10, ...
        'TolCon', 1e-10, ...
        'Display', 'off');  % 关闭显示
    
    % 求解优化问题
    [sol, ~, ~] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlcon, options);

    % 提取结果
    coord = zeros(2, n);
    for i = 1:n
        coord(1,i) = real(sol(2*i - 1));
        coord(2,i) = real(sol(2*i));
    end
    % movable_coordinates = sol(1:end-1);  % x 和 y 坐标
    r_opt = sol(end);  % 最优半径

    figure;
    hold on;
    axis equal;
    theta = linspace(0, 2*pi, 100);
    for i = 1 : n
        x_center = coord(2*i-1);
        y_center = coord(2*i);
        x_circle = x_center + r_opt*cos(theta);
        y_circle = y_center + r_opt*sin(theta);
        plot(x_circle, y_circle, 'b-', 'LineWidth', 1.5);
        plot(x_center, y_center, 'ro');  % 标记圆心
    end
    xlim([-maxBoundary maxBoundary]);
    ylim([-maxBoundary maxBoundary]);
    rectangle('Position', [-maxBoundary -maxBoundary 2*maxBoundary 2*maxBoundary], ...
              'LineStyle', '--', 'LineWidth', 1);
    title(sprintf('Optimal Antenna Placement (r = %.4f)', r_opt));
    grid on;
    hold off;

    coord = coord / 10;

end

function [c, ceq] = nonlinear_constraints(z, num_movable_antennas, max_boundary, TOL)
    % 非线性约束函数
    x = z(1:2:end-1);  % x 坐标
    y = z(2:2:end-1);  % y 坐标
    r = z(end);  % 半径

    % 1. 每个圆必须在正方形内
    c_boundary = [
        x + r - max_boundary;    % 右侧边界
        y + r - max_boundary;    % 上侧边界
        -x + r - max_boundary;   % 左侧边界
        -y + r - max_boundary    % 下侧边界
    ];

    % 2. 每两个圆之间不能重叠
    c_overlap = [];
    for i = 1:num_movable_antennas
        for j = i+1:num_movable_antennas
            dist = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2 + TOL);  % 避免零距离
            c_overlap = [c_overlap; 2*r - dist];
        end
    end

    % 所有不等式约束 (c <= 0)
    c = [c_boundary; c_overlap];
    ceq = [];  % 无等式约束
end