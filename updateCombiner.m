function vCurrent = updateCombiner(vCurrent, wCurrent, H1, K_r, thetaJammerDiscretization, phiJammerDiscretization, Sigma2, Pt, P, Nr, K, L, R, Q1, Q2, lambda, sigma2)
    Q = Q1 * Q2;
    normPow = trace(wCurrent*wCurrent') / Pt;
    G2 = ones(L, 1);
    for k = 1 : K
        % mu = 1e3;
        % nuOpt = -1;
        % nuCurrent = inf;
        % while abs(nuOpt - nuCurrent) > 1e-3
        %     nuCurrent = nuOpt;
        %     cvx_begin quiet
        %         cvx_solver mosek;
        %         variable V(Nr, 1) complex
        %         variable nuOpt
        % 
        %         maximize nuOpt;
        %         subject to
        %         lhsA = 2 * real(vCurrent(:, k)' * H1(:, :, k) * wCurrent(:, k) * (H1(:, :, k) * wCurrent(:, k))' * V);
        %         rhsA = mu*nuOpt + abs(vCurrent(:, k)' * H1(:, :, k) * wCurrent(:, k))^2;
        % 
        %         lhsA >= rhsA;
        % 
        %         for q = 1 : Q1 * Q2
        %             F = create_F(K_r{k}, thetaJammerDiscretization(:, k, q), phiJammerDiscretization(:, k, q), L, Nr, lambda);
        %             lhsB = mu;
        %             rhsB = 0.5 * (lhsB - 1);
        %             lhsB = 0.5 * (lhsB + 1);
        % 
        %             rhsB = [rhsB, V' * H1(:, :, k) * wCurrent(:, setdiff(1:K, k)), sqrt(P) * V' * F' * Sigma2(:, :, k) * G2(:, :, k), V'];
        % 
        %             lhsB >= norm(rhsB, 2);
        % 
        %         end
        %     cvx_end
        %     vCurrent(:, k) = V;
        % end
        ComTemp = zeros(Nr, Nr, Q);
        for r = 1 : R
            for q = 1 : Q
                F = create_F(K_r{k}, thetaJammerDiscretization(:, k, r, q), phiJammerDiscretization(:, k, r, q), L, Nr, lambda);
                ComTemp(:, :, q) = ComTemp(:, :, q) + F' * Sigma2(:, :, k, r) * G2 * (F' * Sigma2(:, :, k, r) * G2)' / Q;
            end
        end
        vCurrent(:, k) = (H1(:, :, k) * (wCurrent * wCurrent') * H1(:, :, k)' + normPow * ( P * sum(ComTemp, 3) + sigma2 * eye(Nr)))^(-1) * H1(:, :, k) * wCurrent(:, k);
        vCurrent(:, k) = vCurrent(:, k) / norm(vCurrent(:, k), 2);
    end

end