function [t, H1, G1] = updateTransPosition(data, t, wCurrent, vCurrent, Sigma1, F1, UserRecJamPow, ACurrent, bCurrent, Pt, Nt, Nr, K, L, lambda, rhoSize, Flag, Gamma)
    
    HalfLambda = lambda / 2;
    normPow = trace(wCurrent*wCurrent') / Pt;
    
    G1 = zeros(L, Nt, K);
    H1 = zeros(Nr, Nt, K);
    AlphaIter = 0;
    AlphaCurrent = - 0.01;
    IterInner = 0;
    while(1)
        IterInner = IterInner + 1;
        % disp(AlphaCurrent);
        if AlphaIter - AlphaCurrent < 1e-4 || IterInner > 100
            break;
        end
        AlphaIter = AlphaCurrent;
        UCurrent = zeros(2*Nt, 2*Nt, K);
        nuCurrent = zeros(2*Nt, K);
        dCurrent = zeros(K, 1);
    
        thetaTrans = data.thetaTransmit;
        phiTrans = data.phiTransmit;
        for k = 1 : K
            varphi = cos(thetaTrans(:, k)) .* cos(phiTrans(:, k));
            vartheta = sin(thetaTrans(:, k));
            % === real ===
            [U, nu, d] = CalRealT(t, varphi, vartheta, wCurrent(:, k), 2 * (ACurrent(1, 2, k) * vCurrent(:, k)' * F1(:, :, k)' * Sigma1(:, :, k))', L, Nt, lambda);
            UCurrent(:, :, k) = UCurrent(:, :, k) + U;
            nuCurrent(:, k) = nuCurrent(:, k) + nu;
            dCurrent(k) = dCurrent(k) + d;
    
            % G1 = create_G(t, thetaTrans(:, k), phiTrans(:, k), L, Nt, lambda);
            % 2*real(ACurrent(1, 2, k) * vCurrent(:, k)' * F1(:, :, k)' * Sigma1(:, :, k) * G1 * wCurrent(:, k))
            % vec(t')' * U * vec(t') - vec(t')' * nu + d
    
            % === quadratic ===
            for i = 1:K
    
                [U, nu, d] = CalAbsT(t, varphi, vartheta,  wCurrent(:, i), sqrt(real(ACurrent(2, 2, k)))*(vCurrent(:, k)' * F1(:, :, k)' * Sigma1(:, :, k))', Nt, L, lambda);
                UCurrent(:, :, k) = UCurrent(:, :, k) + U;
                nuCurrent(:, k) = nuCurrent(:, k) + nu;
                dCurrent(k) = dCurrent(k) + d;
    
                % G1 = create_G(t, thetaTrans(:, k), phiTrans(:, k), L, Nt, lambda);
                % real(ACurrent(2, 2, k)) * abs(vCurrent(:, k)' * F1(:, :, k)' * Sigma1(:, :, k) * G1 * wCurrent(:, i))^2
                % vec(t')' * U * vec(t') - vec(t')' * nu + d
    
            end
            % vec(t')' * UCurrent(:, :, k) * vec(t') - vec(t')' * nuCurrent(:, k) + dCurrent(k)
            % bCurrent(k) - Gamma - real(ACurrent(2, 2, k))*(UserRecJamPow(k) + norm(vCurrent(:, k), 2)^2) - real(ACurrent(1, 1, k))
        end
    
        tVecCurrent = 1e1 * vec(t');
        UCurrent = UCurrent / 1e2;
        nuCurrent = nuCurrent / 1e1;
        % dCurrent = dCurrent / 1e1;
        if Flag == 1
            tOpt = (sum(UCurrent, 3)+sum(UCurrent, 3)') \ sum(nuCurrent, 2);
            tOpt = tOpt / 1e1;
            [FlagTrans, alphaOpt] = CheckFeasiable(K, tOpt, tVecCurrent / 1e1, HalfLambda, rhoSize, Nt, bCurrent, ACurrent, UserRecJamPow, vCurrent, 1e2 * UCurrent, 1e1 * nuCurrent, dCurrent, Gamma, normPow);
            if FlagTrans == 0
                cvx_begin quiet
                    cvx_solver mosek
                    variable alphaOpt(K, 1)
                    variable tOpt(2*Nt, 1)

                    minimize sum(alphaOpt)

                    subject to
                        for k = 1 : K
                            lhsA = quad_form(tOpt, UCurrent(:, :, k)) - tOpt'*nuCurrent(:, k) + dCurrent(k);
                            rhsA = alphaOpt(k);

                            lhsA <= rhsA;

                            lhsB = bCurrent(k) - Gamma - normPow*real(ACurrent(2, 2, k))*(UserRecJamPow(k) + norm(vCurrent(:, k), 2)^2) - real(ACurrent(1, 1, k));
                            rhsB = alphaOpt(k);

                            lhsB >= rhsB;
                        end

                        1e1 * tOpt >= - HalfLambda * rhoSize * 1e2;
                        1e1 * tOpt <= HalfLambda * rhoSize * 1e2;

                        for n = 1 : Nt - 1
                            for m = n + 1 : Nt
                                lhsC = (tVecCurrent([n, Nt+n]) - tVecCurrent([m, Nt+m]))'*(tOpt([n, Nt+n]) - tOpt([m, Nt+m]));
                                rhsC = HalfLambda * norm(tVecCurrent([n, Nt+n]) - tVecCurrent([m, Nt+m]), 2) * 1e1;
                                lhsC >= rhsC;
                            end
                        end

                cvx_end
                if any(isnan(tOpt))
                    for k = 1 : K
                        G1(:, :, k) = create_G(t, thetaTrans(:, k), phiTrans(:, k), L, Nt, lambda);
                        H1(:, :, k) = F1(:, :, k)' * Sigma1(:, :, k) * G1(:, :, k);
                    end
                    break;
                else
                    tOpt = tOpt / 1e1;
                end
            end
            
            % disp(1)
        else
            [tOpt, alphaOpt] = CVXSolDis(UCurrent, nuCurrent, dCurrent, tVecCurrent, K, Nt, HalfLambda, rhoSize);
            if any(isnan(tOpt))
                for k = 1 : K
                    G1(:, :, k) = create_G(t, thetaTrans(:, k), phiTrans(:, k), L, Nt, lambda);
                    H1(:, :, k) = F1(:, :, k)' * Sigma1(:, :, k) * G1(:, :, k);
                end
                break;
            else
                tOpt = tOpt / 1e1;
            end
        end
    
        t = reshape(tOpt, [Nt, 2])';
        AlphaCurrent = sum(alphaOpt);
        for k = 1 : K
            G1(:, :, k) = create_G(t, thetaTrans(:, k), phiTrans(:, k), L, Nt, lambda);
            H1(:, :, k) = F1(:, :, k)' * Sigma1(:, :, k) * G1(:, :, k);
        end
    end

end

function [FlagTrans, SumAlpha] = CheckFeasiable(K, tOpt, tCurrentVec, HalfLambda, rhoSize, Nt, bCurrent, ACurrent, UserRecJamPow, vCurrent, UCurrent, nuCurrent, dCurrent, Gamma, normPow)
    FlagTrans = 1;
    SumAlpha = zeros(K, 1);
    for n = 1 : Nt - 1
        for m = n + 1 : Nt
            lhsC = (tCurrentVec([n, Nt+n]) - tCurrentVec([m, Nt+m]))'*(tOpt([n, Nt+n]) - tOpt([m, Nt+m]));
            rhsC = HalfLambda * norm(tCurrentVec([n, Nt+n]) - tCurrentVec([m, Nt+m]), 2);
            if lhsC < rhsC
                FlagTrans = 0;
            end
        end
    end
    for k = 1 : K
        lhsA = tOpt' * UCurrent(:, :, k) * tOpt - tOpt'*nuCurrent(:, k) + dCurrent(k);
        SumAlpha(k) = lhsA;
        rhsA = bCurrent(k) - Gamma - normPow*real(ACurrent(2, 2, k))*(UserRecJamPow(k) + norm(vCurrent(:, k), 2)^2) - real(ACurrent(1, 1, k));
        if lhsA > rhsA
            FlagTrans = 0;
        end
    end
    if any(tOpt < - HalfLambda * rhoSize) || any(tOpt > HalfLambda * rhoSize)
        FlagTrans = 0;
    end

end

function [tOpt, alphaOpt] = CVXSolDis(UCurrent, nuCurrent, dCurrent, tVecCurrent, K, Nt, HalfLambda, rhoSize)
    cvx_begin quiet
        cvx_solver mosek
        variable alphaOpt(K, 1)
        variable tOpt(2*Nt, 1)
    
        minimize sum(alphaOpt)
    
        subject to
            for k = 1:K
                lhsA = quad_form(tOpt, UCurrent(:, :, k)) - tOpt'*nuCurrent(:, k) + dCurrent(k);
                rhsA = alphaOpt(k);

                lhsA <= rhsA;
            end
        
            1e1 * tOpt >= - HalfLambda * rhoSize * 1e2;
            1e1 * tOpt <= HalfLambda * rhoSize * 1e2;
            for n = 1 : Nt - 1
                for m = n : Nt
                    lhsC = (tVecCurrent([n, Nt+n]) - tVecCurrent([m, Nt+m]))'*(tOpt([n, Nt+n]) - tOpt([m, Nt+m]));
                    rhsC = HalfLambda * norm(tVecCurrent([n, Nt+n]) - tVecCurrent([m, Nt+m]), 2) * 1e1;
                    lhsC >= rhsC;
                end
            end
        
    cvx_end

end


function [U_total, nu_total, d_total] = CalRealT(T, varphi, vartheta, a, b, L, N, lambda)
    
    % 计算f
    f0 = zeros(L, N);
    tau = zeros(L, N);
    for l = 1:L
        for n = 1:N
            phase_term = angle(conj(b(l)) * a(n));
            rho_l_t0 = T(1, n)*varphi(l) + T(2, n)*vartheta(l);
            f0(l, n) = 2*pi/lambda*rho_l_t0 + phase_term;
            tau(l, n) = 2*pi/lambda*rho_l_t0;
        end
    end

    a_mag = abs(a);
    
    % 初始化
    U_total = zeros(2*N, 2*N);
    nu_total = zeros(2*N, 1);
    d_total = 0;

    for l = 1 : L
        weight = abs(b(l));

        U11 = (varphi(l)^2 )*diag(a_mag);
        U12 = (varphi(l)*vartheta(l))*diag(a_mag);
        U22 = (vartheta(l)^2)*diag(a_mag);
        U = (2*pi^2/lambda^2) * [U11, U12; U12, U22];

        e = zeros(N, 1);
        d = 0;
        for n = 1:N
            d = d + a_mag(n)*(cos(f0(l, n)) + 0.5*tau(l, n)^(2) + tau(l, n)*sin(f0(l, n)));
            e(n) = a_mag(n)*(sin(f0(l, n)) + tau(l, n));
        end

        nu = zeros(2*N, 1);

        nu(1:N) = (2*pi/lambda) * (varphi(l) * e );
        nu(N+1:end) = (2*pi/lambda) * (vartheta(l) * e);

        U_total = U_total + weight * U;
        nu_total = nu_total + weight * nu;
        d_total = d_total + weight * d;
    end

end

function [U_total, nu_total, d_total] = CalAbsT(T, varphi, vartheta, a, b, N, L, lambda)

    % 计算 f0 和 f
    f0 = zeros(L, L, N, N);
    tau = zeros(L, L, N, N);
    for l1 = 1:L
        for l2 = 1:L
            for n1 = 1:N
                for n2 = 1:N
                    % 计算相位项
                    phase_term = angle(conj(b(l1)) * b(l2) * a(n1) * conj(a(n2)));

                    % f0 在初始点
                    rho_l1_t0 = T(1, n1)*varphi(l1) + T(2, n1)*vartheta(l1);
                    rho_l2_t0 = T(1, n2)*varphi(l2) + T(2, n2)*vartheta(l2);
                    f0(l1, l2, n1, n2) = 2*pi/lambda*(rho_l1_t0 - rho_l2_t0) + phase_term;
                    tau(l1, l2, n1, n2) = 2*pi/lambda*(rho_l1_t0 - rho_l2_t0);
                end
            end
        end
    end

    a_mag = abs(a);
    kappa = sum(a_mag);

    U_total = zeros(2*N, 2*N);
    nu_total = zeros(2*N, 1);
    d_total = 0;

    for l1 = 1:L
        for l2 = 1:L
            bl1_mag = abs(b(l1)); bl2_mag = abs(b(l2));
            weight = bl1_mag * bl2_mag;

            U11 = (varphi(l1)^2 + varphi(l2)^2)*diag(kappa*a_mag) - 2*varphi(l1)*varphi(l2)*(a_mag*a_mag');
            U12 = (varphi(l1)*vartheta(l1) + varphi(l2)*vartheta(l2))*diag(kappa*a_mag) - (varphi(l1)*vartheta(l2) + varphi(l2)*vartheta(l1))*(a_mag*a_mag');
            U22 = (vartheta(l1)^2 + vartheta(l2)^2)*diag(kappa*a_mag) - 2*vartheta(l1)*vartheta(l2)*(a_mag*a_mag');

            U_l1l2 = (2*pi^2/lambda^2) * [U11, U12; U12, U22];

            e = zeros(N, 1);
            p = zeros(N, 1);

            d_l1l2 = 0;
            for n = 1:N
                e_ele = 0;
                p_ele = 0;
                for m = 1:N
                    e_ele = e_ele + a_mag(m)*(sin(f0(l1, l2, n, m)) + tau(l1, l2, n, m));
                    p_ele = p_ele + a_mag(m)*(sin(f0(l1, l2, m, n)) + tau(l1, l2, m, n));

                    d_l1l2 = d_l1l2 + a_mag(n)*a_mag(m)*(cos(f0(l1, l2, n, m)) + 0.5*tau(l1, l2, n, m)^(2) + tau(l1, l2, n, m)*sin(f0(l1, l2, n, m)));
                end
                e(n) = a_mag(n)*e_ele;
                p(n) = a_mag(n)*p_ele;
            end


            nu_l1l2 = zeros(2*N, 1);

            nu_l1l2(1:N) = (2*pi/lambda) * (varphi(l1) * e - varphi(l2) * p);
            nu_l1l2(N+1:end) = (2*pi/lambda) * (vartheta(l1) * e - vartheta(l2) * p);

            U_total = U_total + weight * U_l1l2;
            nu_total = nu_total + weight * nu_l1l2;
            d_total = d_total + weight * d_l1l2;
        end
    end


end