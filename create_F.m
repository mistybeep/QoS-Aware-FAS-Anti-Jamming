function F = create_F(r, theta, phi, L, N, lambda)
    F = zeros(L, N);
    for l = 1 : L
        n_r = [cos(theta(l)) * cos(phi(l)); sin(theta(l))];
        for nr = 1 : N
            % 计算信号传播距离差 rho_r
            rho_r = n_r' * r(:, nr);
            % 计算相位差并填入场响应向量
            F(l, nr) = exp(1j * 2 * pi * rho_r / lambda);
        end
    end

end