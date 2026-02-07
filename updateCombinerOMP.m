function [vCurrent, C] = updateCombinerOMP(vCurrent, HatH1, HatH2, B, C, wCurrent, Pt, P, Gr, Nr, R, Q1, Q2, sigma2)
    % F1 维度 L * Gr * K F2 维度 L * Gr * K * Q
    
    K = size(wCurrent, 2);

    Gamma = sigma2 * trace(wCurrent*wCurrent')/Pt;
    HatJamPow = sqrt(Gamma * P / sigma2);
    
    E = eye(K + R * Q1 * Q2);

    T_C = zeros(Gr, Nr, K);

    for k = 1 : K
        e = E(:, k);
        Phi = [HatH1(:, :, k) * B * wCurrent, HatJamPow * reshape(squeeze(HatH2(:, k, :, :)), Gr, R*Q1*Q2) / sqrt(Q1 * Q2)];

        Phi = Phi';
        [pos, VTemp] = RLS_SOMP(e, Phi, Gr, Nr, Gamma);
        % [pos, VTemp] = Row_IHT(Phi, e, vInit, Gamma, Nr, optsIHT);
        % [pos, VTemp] = Nonconvex_ProxGrad(Phi, e, Nr, Gamma);
        
        for m = 1 : Nr
            T_C(pos(m), m, k) = 1;
            % vCurrent(m, k) = VTemp(pos(m), :);
        end
        vCurrent(:,k) = VTemp / norm(VTemp, 2);
        % vCurrent(:,k) = vCurrent(:,k) / norm(vCurrent(:,k), 2);
        C(:, :, k) = T_C(:, :, k);
    end
    % J
end