function [F, H, G] = GenChannel(t, r, lambda, thetaTransmit, phiTransmit, thetaReceive, phiReceive, Sigma)
    Nt = size(t, 2); % 发射端天线数
    Nr = size(r, 2); % 接收端天线数
    L = size(thetaReceive, 1);

    % 初始化发射端场响应矩阵
    if Nt == 1
        G = ones(L, 1);
    else
        G = zeros(L, Nt); % 发射端场响应矩阵 G (L x Nt)

        % 计算发射端波矢量和场响应向量
        for m = 1:L
            % 发射端归一化波矢量 n_t
            n_t = [cos(thetaTransmit(m)) * cos(phiTransmit(m)); sin(thetaTransmit(m))];
            for nt = 1:Nt
                % 计算信号传播距离差 rho_t
                rho_t = n_t' * t(:, nt);
                % 计算相位差并填入场响应向量
                G(m, nt) = exp(1j * 2 * pi * rho_t / lambda);
            end
        end
    end

    % 初始化接收端场响应矩阵
    F = zeros(L, Nr); % 接收端场响应矩阵 F (M x Nr)

    for m = 1:L
        % 接收端归一化波矢量 n_r
        n_r = [cos(thetaReceive(m)) * cos(phiReceive(m)); sin(thetaReceive(m))];
        for nr = 1:Nr
            % 计算信号传播距离差 rho_r
            rho_r = n_r' * r(:, nr);
            % 计算相位差并填入场响应向量
            F(m, nr) = exp(1j * 2 * pi * rho_r / lambda);
        end
    end

    % 计算信道矩阵 H = F^H * Sigma * G
    H = F' * Sigma * G;

end