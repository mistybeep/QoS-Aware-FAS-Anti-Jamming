function [wCurrent, Flag] = updatePrecoder(wCurrent, vCurrent, H1, UserRecJamPow, bCurrent, ACurrent, Nt, K, Pt, Gamma, sigma2)
    Flag = 1;
    % === 等效信道 ===
    G = zeros(Nt, K);
    for k = 1 : K
        G(:, k) = H1(:, :, k)' * vCurrent(:, k);
    end
    
    W = W_MM(vCurrent, H1, ACurrent, K, Nt, UserRecJamPow, Pt, sigma2);
    PowGamma = norm(vec(W), 2)^2 / Pt;
    FlagQoS = 1;
    for k = 1 : K
        lhsA = 2 * real(ACurrent(1, 2, k) * G(:, k)' * W(:, k)) ...
                + real(ACurrent(2, 2, k)) * (sum_square_abs(vec(G(:, k)' * W)) + PowGamma * (UserRecJamPow(k) + sigma2*norm(vCurrent(:, k), 2)^2));

        rhsA = bCurrent(k) - Gamma - real(ACurrent(1, 1, k));

        if lhsA > rhsA
            FlagQoS = 0;
        end
    end
    if FlagQoS==0
        cvx_begin quiet
            cvx_solver mosek
            variable W(Nt, K) complex
            MinMaxObj = 0;
            PowGamma = sum_square_abs(vec(W)) / Pt;
            for k = 1 : K
                MinMaxObj = MinMaxObj + 2 * real(ACurrent(1, 2, k) * G(:, k)' * W(:, k)) ...
                    + real(ACurrent(2, 2, k)) * (sum_square_abs(vec(G(:, k)' * W)) + PowGamma * (UserRecJamPow(k) + sigma2*norm(vCurrent(:, k), 2)^2));
            end

            minimize MinMaxObj;

            subject to
                for k = 1 : K

                    lhsA = 2 * real(ACurrent(1, 2, k) * G(:, k)' * W(:, k)) ...
                    + real(ACurrent(2, 2, k)) * (sum_square_abs(vec(G(:, k)' * W)) + PowGamma * (UserRecJamPow(k) + sigma2*norm(vCurrent(:, k), 2)^2));

                    rhsA = bCurrent(k) - Gamma - real(ACurrent(1, 1, k));

                    lhsA <= rhsA;

                end

        cvx_end
    end
    % 
    if ~any(isnan(W))
        wCurrent = W;
    else
        wCurrent = W_MM(vCurrent, H1, ACurrent, K, Nt, UserRecJamPow, Pt, sigma2);
        % wCurrent = sqrt(Pt) * W / norm(W, 'fro');
        Flag = 0;
    end
end

function W = W_MM(vCurrent, H1, ACurrent, K, Nt, UserRecJamPow, Pt, sigma2)
    A11 = zeros(K, K);
    A12 = zeros(K, K);

    GammaMinMax = 0;

    T_H = zeros(K, Nt);
    
    for i = 1:K
        A11(i, i) = real(ACurrent(1, 1, i));
        A12(i, i) = - ACurrent(1, 2, i);
        T_H(i, :) = vCurrent(:, i)' * H1(:, :, i);
        GammaMinMax = GammaMinMax + real(ACurrent(2, 2, i) * (UserRecJamPow(i) + sigma2*norm(vCurrent(:, i), 2)^2) / Pt);
    end
    
    Phi = (A11)^(-0.5) * A12 * T_H;
    W = (Phi' * Phi + GammaMinMax * eye(Nt))^(-1) * Phi' * (A11)^(0.5);

end