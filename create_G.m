function G = create_G(t, theta, phi, L, N, lambda)
    G = zeros(L, N);
    for l = 1 : L
        n_t = [cos(theta(l)) * cos(phi(l)); sin(theta(l))];
        for nt = 1 : N
            % 计算信号传播距离差 rho_r
            rho_t = n_t' * t(:, nt);
            % 计算相位差并填入场响应向量
            G(l, nt) = exp(1j * 2 * pi * rho_t / lambda);
        end
    end

end