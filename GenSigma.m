function Sigma = GenSigma(L, rho, alpha, d, sigma)
    variance = db2pow(rho) * d^(-alpha);
    Sigma = zeros(L, L);
    diag_elements = sqrt(variance) * (randn(L, 1) + 1j * randn(L, 1)) / sqrt(2); % CSCG 分布
    Sigma(1, 1) = diag_elements(1); % 第一个元素为 LoS 分量
    for m = 2:L
        Sigma(m, m) = diag_elements(m) / (L - 1); % 非 LoS 分量除以 L-1
    end

    Sigma = Sigma / sigma;

end

