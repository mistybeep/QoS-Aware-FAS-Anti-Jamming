function pos = GenPosCurrent(rhoSize, Lambda, N)
    % ==== rhoSize : 边界大小 Lambda ： 波长 N : 天线数量
    W = rhoSize * Lambda;          % 正方形区域边长
    borderMargin = - W / 2 + 0.25 * W;  % 边界间距
    n = sqrt(N);

    %% 计算点间距
    % effW = W - 2 * borderMargin;  % 有效区域大小
    % pointSpacing = effW / (n - 1); % 点间距
    
    %% 生成网格点
    x_coords = linspace(borderMargin, W / 2  - 0.05*W, n);
    y_coords = linspace(borderMargin, W / 2  - 0.05*W, n);
    [X, Y] = meshgrid(x_coords, y_coords);

    pos = [X(:), Y(:)]';

end