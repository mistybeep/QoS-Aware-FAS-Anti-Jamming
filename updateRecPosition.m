function [K_r, F1, H1] = updateRecPosition(data, K_r, wCurrent, vCurrent, Sigma1, G1, Sigma2, thetaJammerDiscretization, phiJammerDiscretization, Pt, P, Nt, Nr, K, L, R, Q1, Q2, lambda, rhoSize, sigma2)
    
    HalfLambda = lambda / 2;
    normPow = trace(wCurrent*wCurrent') / Pt;
    % === 合法用户角度信息 ===
    thetaRec = data.thetaReceive;
    phiRec = data.phiReceive;

    F1 = zeros(L, Nr, K);
    H1 = zeros(Nr, Nt, K);

    G2 = ones(L, 1);

    for k = 1 : K
        varphi = cos(thetaRec(:, k)) .* cos(phiRec(:, k));
        vartheta = sin(thetaRec(:, k));
        kappaIter = 0;
        kappaCurrent = 0.001;
        IterInner = 0;
        while (1)
            IterInner = IterInner + 1;
            if abs(kappaCurrent - kappaIter) < 1e-4 || IterInner > 100
                break;
            end
            kappaIter = kappaCurrent;
            F1Temp = create_F(K_r{k}, thetaRec(:, k), phiRec(:, k), L, Nr, lambda);
            JamTemp = zeros(Q1*Q2, 1);
            for r = 1 : R
                for q = 1 : Q1 * Q2
                    F = create_F(K_r{k}, thetaJammerDiscretization(:, k, r, q), phiJammerDiscretization(:, k, r, q), L, Nr, lambda);
                    JamTemp(q) = JamTemp(q) + P * abs(vCurrent(:, k)' * F' * Sigma2(:, :, k, r) * G2)^2 / (Q1 * Q2);
                end
            end
            MaxRecJamPow = sum(JamTemp);
            muCurrent = norm(vCurrent(:, k)' * F1Temp' * Sigma1(:, :, k) * G1(:, :, k) * wCurrent(:, setdiff(1:K, k)), 2)^2 + normPow * (MaxRecJamPow + sigma2 * norm(vCurrent(:, k), 2)^2);
            kappaCurrent = abs(vCurrent(:, k)' * F1Temp' * Sigma1(:, :, k) * G1(:, :, k) * wCurrent(:, k))^2 / muCurrent;
            % disp(kappaCurrent)
            % === 分子项 ===
            [UNumTotal, nuNumTotal, ~] = CalAbsRLower(K_r{k}, varphi, vartheta, vCurrent(:, k), Sigma1(:, :, k) * G1(:, :, k) * wCurrent(:, k), Nr, L, lambda);
            % F1 = create_F(K_r{k}, thetaRec(:, k), phiRec(:, k), L, Nr, lambda);
            UDomTotal = zeros(2*Nr, 2*Nr);
            nuDomTotal = zeros(2*Nr, 1);
            % dDomTotal = 0;
            for i = setdiff(1:K, k)
                [U, nu, ~] = CalAbsRUpper(K_r{k}, varphi, vartheta, vCurrent(:, k), Sigma1(:, :, k) * G1(:, :, k) * wCurrent(:, i), Nr, L, lambda);
                UDomTotal = UDomTotal + U;
                nuDomTotal = nuDomTotal + nu;
                % dDomTotal = dDomTotal + d;
            end
            for r = 1 : R
                for q = 1 : Q1 * Q2
                    varJamphi = cos(thetaJammerDiscretization(:, k, r, q)) .* cos(phiJammerDiscretization(:, k, r, q));
                    varJamtheta = sin(thetaJammerDiscretization(:, k, r, q));
                    [U, nu, ~] = CalAbsRUpper(K_r{k}, varJamphi, varJamtheta, vCurrent(:, k), Sigma2(:, :, k, r) * G2, Nr, L, lambda);
                    UDomTotal = UDomTotal + P * U / (Q1 * Q2);
                    nuDomTotal = nuDomTotal + P * nu / (Q1 * Q2);
                    % dDomTotal = dDomTotal + P * d / (Q1 * Q2);
                end
            end

            UCurrent = kappaCurrent * UDomTotal - UNumTotal;
            nuCurrent = kappaCurrent * nuDomTotal - nuNumTotal;
            % dCurrent = kappaCurrent * dDomTotal - dNumTotal;
    
            rCurrentVec = vec(K_r{k}');
            rOpt = (UCurrent + UCurrent') \ nuCurrent;
            % rOpt = max(-0.5 * rhoSize * lambda, min(0.5 * rhoSize * lambda, rOpt));
            Flag = CheckFeasiable(rOpt, rCurrentVec, HalfLambda, rhoSize, Nr);

            if Flag == 0

                UCurrent = UCurrent / 1e4;
                nuCurrent = nuCurrent / 1e2;
                rCurrentVec = 1e2 * rCurrentVec;

                cvx_begin quiet
                cvx_solver mosek
                variable rOpt(2*Nr, 1);
                variable t;
                minimize t;
                subject to

                    rOpt' * UCurrent * rOpt - rOpt' * nuCurrent <= t;

                    1e1 * rOpt >= - HalfLambda * rhoSize * 1e3;
                    1e1 * rOpt <= HalfLambda * rhoSize * 1e3;
    
                    for n = 1 : Nr - 1
                        for m = n + 1 : Nr
                            lhsC = (rCurrentVec([n, Nr+n]) - rCurrentVec([m, Nr+m]))'*(rOpt([n, Nr+n]) - rOpt([m, Nr+m]));
                            rhsC = HalfLambda * norm(rCurrentVec([n, Nr+n]) - rCurrentVec([m, Nr+m]), 2) * 1e2;
                            lhsC >= rhsC;
                        end
                    end

                cvx_end
                if any(isnan(rOpt))
                    F1(:, :, k) = create_F(K_r{k}, thetaRec(:, k), phiRec(:, k), L, Nr, lambda);
                    H1(:, :, k) = F1(:, :, k)' * Sigma1(:, :, k) * G1(:, :, k);
                    break;
                else
                    rOpt = rOpt / 1e2;
                end
            end
            
            K_r{k} = reshape(rOpt, [Nr, 2])';
            F1(:, :, k) = create_F(K_r{k}, thetaRec(:, k), phiRec(:, k), L, Nr, lambda);
            H1(:, :, k) = F1(:, :, k)' * Sigma1(:, :, k) * G1(:, :, k);
        end
        
    end
    
end

function Flag = CheckFeasiable(rOpt, rCurrentVec, HalfLambda, rhoSize, Nr)
    Flag = 1;
    for n = 1:Nr - 1
        for m = n+1 : Nr
            lhsC = (rCurrentVec([n, Nr+n]) - rCurrentVec([m, Nr+m]))'*(rOpt([n, Nr+n]) - rOpt([m, Nr+m]));
            rhsC = HalfLambda * norm(rCurrentVec([n, Nr+n]) - rCurrentVec([m, Nr+m]), 2);
            if lhsC < rhsC
                Flag = 0;
            end
        end
    end
    if any(rOpt < - HalfLambda * rhoSize) || any(rOpt > HalfLambda * rhoSize)
        Flag = 0;
    end

end

function [U_total, nu_total, d_total] = CalAbsRLower(R, varphi, vartheta, a, b, M, L, lambda)

    % 计算 f0 和 f
    f0 = zeros(L, L, M, M);
    tau = zeros(L, L, M, M);
    for l1 = 1:L
        for l2 = 1:L
            for m1 = 1:M
                for m2 = 1:M
                    % 计算相位项
                    phase_term = angle(conj(b(l1)) * b(l2) * a(m1) * conj(a(m2)));

                    % f0 在初始点
                    rho_l1_r0 = R(1, m1)*varphi(l1) + R(2, m1)*vartheta(l1);
                    rho_l2_r0 = R(1, m2)*varphi(l2) + R(2, m2)*vartheta(l2);
                    f0(l1, l2, m1, m2) = 2*pi/lambda*(rho_l1_r0 - rho_l2_r0) + phase_term;
                    tau(l1, l2, m1, m2) = 2*pi/lambda*(rho_l1_r0 - rho_l2_r0);
                end
            end
        end
    end

    a_mag = abs(a);
    kappa = sum(a_mag);

    U_total = zeros(2*M, 2*M);
    nu_total = zeros(2*M, 1);
    d_total = 0;

    for l1 = 1:L
        for l2 = 1:L
            bl1_mag = abs(b(l1)); bl2_mag = abs(b(l2));
            weight = bl1_mag * bl2_mag;

            U11 = (varphi(l1)^2 + varphi(l2)^2)*diag(kappa*a_mag) - 2*varphi(l1)*varphi(l2)*(a_mag*a_mag');
            U12 = (varphi(l1)*vartheta(l1) + varphi(l2)*vartheta(l2))*diag(kappa*a_mag) - (varphi(l1)*vartheta(l2) + varphi(l2)*vartheta(l1))*(a_mag*a_mag');
            U22 = (vartheta(l1)^2 + vartheta(l2)^2)*diag(kappa*a_mag) - 2*vartheta(l1)*vartheta(l2)*(a_mag*a_mag');

            U_l1l2 = -(2*pi^2/lambda^2) * [U11, U12; U12, U22];

            e = zeros(M, 1);
            p = zeros(M, 1);

            d_l1l2 = 0;
            for n = 1:M
                e_ele = 0;
                p_ele = 0;
                for m = 1:M
                    e_ele = e_ele + a_mag(m)*(sin(f0(l1, l2, n, m)) - tau(l1, l2, n, m));
                    p_ele = p_ele + a_mag(m)*(sin(f0(l1, l2, m, n)) - tau(l1, l2, m, n));

                    d_l1l2 = d_l1l2 + a_mag(n)*a_mag(m)*(cos(f0(l1, l2, n, m)) - 0.5*tau(l1, l2, n, m)^(2) + tau(l1, l2, n, m)*sin(f0(l1, l2, n, m)));
                end
                e(n) = a_mag(n)*e_ele;
                p(n) = a_mag(n)*p_ele;
            end


            nu_l1l2 = zeros(2*M, 1);

            nu_l1l2(1:M) = (2*pi/lambda) * (varphi(l1) * e - varphi(l2) * p);
            nu_l1l2(M+1:end) = (2*pi/lambda) * (vartheta(l1) * e - vartheta(l2) * p);

            U_total = U_total + weight * U_l1l2;
            nu_total = nu_total + weight * nu_l1l2;
            d_total = d_total + weight * d_l1l2;
        end
    end


end


function [U_total, nu_total, d_total] = CalAbsRUpper(R, varphi, vartheta, a, b, M, L, lambda)

    % 计算 f0 和 f
    f0 = zeros(L, L, M, M);
    tau = zeros(L, L, M, M);
    for l1 = 1:L
        for l2 = 1:L
            for m1 = 1:M
                for m2 = 1:M
                    % 计算相位项
                    phase_term = angle(conj(b(l1)) * b(l2) * a(m1) * conj(a(m2)));

                    % f0 在初始点
                    rho_l1_r0 = R(1, m1)*varphi(l1) + R(2, m1)*vartheta(l1);
                    rho_l2_r0 = R(1, m2)*varphi(l2) + R(2, m2)*vartheta(l2);
                    f0(l1, l2, m1, m2) = 2*pi/lambda*(rho_l1_r0 - rho_l2_r0) + phase_term;
                    tau(l1, l2, m1, m2) = 2*pi/lambda*(rho_l1_r0 - rho_l2_r0);
                end
            end
        end
    end

    a_mag = abs(a);
    kappa = sum(a_mag);

    U_total = zeros(2*M, 2*M);
    nu_total = zeros(2*M, 1);
    d_total = 0;

    for l1 = 1:L
        for l2 = 1:L
            bl1_mag = abs(b(l1)); bl2_mag = abs(b(l2));
            weight = bl1_mag * bl2_mag;

            U11 = (varphi(l1)^2 + varphi(l2)^2)*diag(kappa*a_mag) - 2*varphi(l1)*varphi(l2)*(a_mag*a_mag');
            U12 = (varphi(l1)*vartheta(l1) + varphi(l2)*vartheta(l2))*diag(kappa*a_mag) - (varphi(l1)*vartheta(l2) + varphi(l2)*vartheta(l1))*(a_mag*a_mag');
            U22 = (vartheta(l1)^2 + vartheta(l2)^2)*diag(kappa*a_mag) - 2*vartheta(l1)*vartheta(l2)*(a_mag*a_mag');

            U_l1l2 = (2*pi^2/lambda^2) * [U11, U12; U12, U22];

            e = zeros(M, 1);
            p = zeros(M, 1);

            d_l1l2 = 0;
            for n = 1:M
                e_ele = 0;
                p_ele = 0;
                for m = 1:M
                    e_ele = e_ele + a_mag(m)*(sin(f0(l1, l2, n, m)) + tau(l1, l2, n, m));
                    p_ele = p_ele + a_mag(m)*(sin(f0(l1, l2, m, n)) + tau(l1, l2, m, n));

                    d_l1l2 = d_l1l2 + a_mag(n)*a_mag(m)*(cos(f0(l1, l2, n, m)) + 0.5*tau(l1, l2, n, m)^(2) + tau(l1, l2, n, m)*sin(f0(l1, l2, n, m)));
                end
                e(n) = a_mag(n)*e_ele;
                p(n) = a_mag(n)*p_ele;
            end


            nu_l1l2 = zeros(2*M, 1);

            nu_l1l2(1:M) = (2*pi/lambda) * (varphi(l1) * e - varphi(l2) * p);
            nu_l1l2(M+1:end) = (2*pi/lambda) * (vartheta(l1) * e - vartheta(l2) * p);

            U_total = U_total + weight * U_l1l2;
            nu_total = nu_total + weight * nu_l1l2;
            d_total = d_total + weight * d_l1l2;
        end
    end


end