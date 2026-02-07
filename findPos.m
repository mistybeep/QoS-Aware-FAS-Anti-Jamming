function r_An = findPos(r_A_opt, r_Astar, lambda, normA, rad)

    D = 0.5 * lambda;
    x_min = -0.5*normA*lambda;
    x_max = 0.5*normA*lambda;
    y_min = -0.5*normA*lambda;
    y_max = 0.5*normA*lambda;

    r_Astar_2 = [
        max(x_min, min(x_max, r_Astar(1))), ...
        max(y_min, min(y_max, r_Astar(2)))
        ];

    % 已经优化的天线数量
    % [~, N] = size(r_A_opt);

    kdtree = createns(r_A_opt, 'NSMethod', 'kdtree');
    % 检查 P0 是否可行
    distances = sqrt(sum((r_A_opt - r_Astar_2).^2, 2));
    if all(distances >= D)
        % fprintf('最优点 P0 = (%.4f, %.4f) 已满足约束\n', P0(1), P0(2));
        r_An = r_Astar_2;
    else
        % 活跃约束筛选
        range = rad + D; % 搜索范围
        [active_idx, ~] = rangesearch(kdtree, r_Astar, range);
        active_idx = active_idx{1};
        M = length(active_idx);
        % active_P = r_A_opt(active_idx, :);

        % 候选点集合
        candidates = [];

        % 单圆投影
        for i = 1:M
            idx = active_idx(i);
            dir = r_Astar - r_A_opt(idx,:);
            dist = norm(dir);
            if dist < D
                cand = r_A_opt(idx,:) + D * (dir / dist);
                candidates = [candidates; cand];
            end
        end

        % 双圆交点（仅活跃点对）
        for i = 1:M
            for j = i+1:M
                idx_i = active_idx(i); idx_j = active_idx(j);
                c1 = r_A_opt(idx_i,:); c2 = r_A_opt(idx_j,:);
                d = norm(c1 - c2);
                if d < 2*D && d > 0
                    a = (d^2 + D^2 - D^2) / (2*d);
                    h = sqrt(D^2 - a^2);
                    mid = c1 + a * (c2 - c1) / d;
                    perp = [-(c2(2) - c1(2)), c2(1) - c1(1)] / d;
                    p1 = mid + h * perp;
                    in_bounds = p1(1) >= x_min && p1(1) <= x_max && p1(2) >= y_min && p1(2) <= y_max;
                    p2 = mid - h * perp;
                    candidates = [candidates; p1; p2];
                end
            end
        end

        % 正方形边界交点（仅活跃点）
        edges = {[x_min, y_min, x_min, y_max], [x_max, y_min, x_max, y_max], ...
            [x_min, y_min, x_max, y_min], [x_min, y_max, x_max, y_max]};
        for i = 1:M
            idx = active_idx(i);
            for e = 1:4
                x1 = edges{e}(1); y1 = edges{e}(2);
                x2 = edges{e}(3); y2 = edges{e}(4);
                [x, y] = circle_line_intersection(r_A_opt(idx,:), D, [x1, y1], [x2, y2]);
                if ~isempty(x) && ~isempty(y)
                    intersect_points = [x(:), y(:)];
                    candidates = [candidates; intersect_points];
                end
            end
        end

        % 筛选满足所有约束的点（用 KD 树加速）
        valid_candidates = [];
        for k = 1:size(candidates, 1)
            cand = candidates(k, :);
            in_bounds = cand(1) >= x_min && cand(1) <= x_max && cand(2) >= y_min && cand(2) <= y_max;
            if in_bounds
                [~, min_dist] = knnsearch(kdtree, cand, 'K', 1);
                if abs(min_dist - D) < 1e-4 % min_dist >= D
                    valid_candidates = [valid_candidates; cand];
                end
            end
        end

        % 选择离 P0 最近的点
        if isempty(valid_candidates)
            % fprintf('未找到满足所有约束的点，可能区域无解\n');
            r_An = [NaN NaN];
        else
            distances_to_P0 = sqrt(sum((valid_candidates - r_Astar).^2, 2));
            [~, idx] = min(distances_to_P0);
            P_new = valid_candidates(idx, :);
            % fprintf('找到满足约束的点 P_new = (%.4f, %.4f)\n', P_new(1), P_new(2));
            % fprintf('到 P0 的距离 = %.4f\n', min_dist);
            r_An = P_new;
        end
    end

end

function [x, y] = circle_line_intersection(center, r, p1, p2)
    x0 = center(1); y0 = center(2);
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    
    dx = x2 - x1; dy = y2 - y1;
    dr = sqrt(dx^2 + dy^2);
    D = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    delta = r^2 * dr^2 - D^2;
    
    if delta < 0
        x = []; y = [];
    else
        sgn = sign(dy); if sgn == 0, sgn = 1; end
        x = [(D * dy + sgn * dx * sqrt(delta)) / (dr^2) + x0, ...
             (D * dy - sgn * dx * sqrt(delta)) / (dr^2) + x0];
        y = [(-D * dx + abs(dy) * sqrt(delta)) / (dr^2) + y0, ...
             (-D * dx - abs(dy) * sqrt(delta)) / (dr^2) + y0];
        valid = (x >= min(x1, x2) & x <= max(x1, x2) & y >= min(y1, y2) & y <= max(y1, y2));
        x = x(valid); y = y(valid);
    end
end