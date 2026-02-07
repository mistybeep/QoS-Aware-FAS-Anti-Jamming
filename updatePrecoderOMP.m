function [wCurrent, B] = updatePrecoderOMP(wCurrent, vCurrent, BCurrent, HatH1, C, Pt, UserRecJamPow, bCurrent, ACurrent, Gt, Gamma, Nt, K, sigma2)
    
    nuTempJudge = zeros(K, 1);
    [WTemp, B] = wCurrentOMP(vCurrent, wCurrent, BCurrent, HatH1, C, nuTempJudge, ACurrent, K, Gt, Nt, UserRecJamPow, Pt, sigma2);
    rate = KUpdateSINR(vCurrent, WTemp, HatH1, B, C, UserRecJamPow, Pt, Nt, K, sigma2);

    IterInner = 0;

    nu = zeros(K, 1);
    RES = 0;
    wCurrent = WTemp;

    if rate <= K * (Gamma + 1)
        wCurrent = WTemp;
    else
        while(IterInner < 100)
            IterInner = IterInner + 1;
            RES2 = RES;
            nuTempJudge = nu;
            for k = 1:K
                % 构建公共矩阵
                nuTempJudge(k) = 0;
                [wCurrent, B] = wCurrentOMP(vCurrent, wCurrent, BCurrent, HatH1, C, nuTempJudge, ACurrent, K, Gt, Nt, UserRecJamPow, Pt, sigma2);
                res = DualFunc(HatH1(:, :, k), B, C(:, :, k), vCurrent, wCurrent, k, UserRecJamPow, Pt, ACurrent, sigma2);
                if res <= bCurrent(k) - Gamma
                    nu(k) = 0;
                    continue;
                end
                nuTempJudge(k) = 1e15;
                [wCurrent, B] = wCurrentOMP(vCurrent, wCurrent, BCurrent, HatH1, C, nuTempJudge, ACurrent, K, Gt, Nt, UserRecJamPow, Pt, sigma2);
                res = DualFunc(HatH1(:, :, k), B, C(:, :, k), vCurrent, wCurrent, k, UserRecJamPow, Pt, ACurrent, sigma2);
                if abs(res - (bCurrent(k) - Gamma)) < 1e-3
                    nuTempJudge(k) = 1e15;
                else
                    lower = 0;
                    upper = 1;

                    nuTempJudge(k) = upper;
                    [wCurrent, B] = wCurrentOMP(vCurrent, wCurrent, BCurrent, HatH1, C, nuTempJudge, ACurrent, K, Gt, Nt, UserRecJamPow, Pt, sigma2);
                    res = DualFunc(HatH1(:, :, k), B, C(:, :, k), vCurrent, wCurrent, k, UserRecJamPow, Pt, ACurrent, sigma2);
                    T_ITER = 0;
                    while res > bCurrent(k) - Gamma && T_ITER < 100
                        T_ITER = T_ITER + 1;
                        upper = upper * 2;
                        nuTempJudge(k) = upper;
                        [wCurrent, B] = wCurrentOMP(vCurrent, wCurrent, BCurrent, HatH1, C, nuTempJudge, ACurrent, K, Gt, Nt, UserRecJamPow, Pt, sigma2);
                        res = DualFunc(HatH1(:, :, k), B, C(:, :, k), vCurrent, wCurrent, k, UserRecJamPow, Pt, ACurrent, sigma2);
                    end
                    ITER = 0;
                    while abs(upper - lower) > 1e-3 && ITER < 100
                        ITER = ITER + 1;
                        mid = (lower + upper) / 2;
                        nuTempJudge(k) = mid;
                        [wCurrent, B] = wCurrentOMP(vCurrent, wCurrent, BCurrent, HatH1, C, nuTempJudge, ACurrent, K, Gt, Nt, UserRecJamPow, Pt, sigma2);
                        res = DualFunc(HatH1(:, :, k), B, C(:, :, k), vCurrent, wCurrent, k, UserRecJamPow, Pt, ACurrent, sigma2);
                        if res > bCurrent(k) - Gamma
                            lower = mid;
                        else
                            upper = mid;
                        end
                    end
                    nuTempJudge(k) = (lower + upper) / 2;
                end
            end
            nu = nuTempJudge;
            RES = DualFunc2(HatH1, B, C, vCurrent, wCurrent, k, UserRecJamPow, Pt, ACurrent, sigma2);
            RES = sum(RES .* (1 + nu));
            if abs(RES - RES2) < 1e-3
                break;
            end
        end
    end
end


function [wCurrent, B] = wCurrentOMP(vCurrent, wCurrent, BCurrent, HatH1, C, nu, ACurrent, K, Gt, Nt, UserRecJamPow, Pt, sigma2)

    A11 = zeros(K, K);
    A12 = zeros(K, K);
    Psi = zeros(K, Gt);
    gamma = 0;

    B = zeros(Gt, Nt);
    
    for i = 1 : K
        A11(i, i) = (1 + nu(i)) * real(ACurrent(1, 1, i));
        A12(i, i) = - (1 + nu(i)) * ACurrent(1, 2, i);
        Psi(i, :) = vCurrent(:, i)' * C(:, :, i)' * HatH1(:, :, i);
        gamma = gamma + (1 + nu(i)) * real(ACurrent(2, 2, i) * (UserRecJamPow(i) + sigma2 * norm(vCurrent(:, i), 2)^2) / Pt);
    end
    
    Phi = (A11)^(-0.5) * A12 * Psi;
    % [pos, wCurrent] = RLS_SOMP((A11)^(0.5), Phi, Gt, Nt, gamma);
    % [pos, WTemp] = Row_IHT(Phi, (A11)^(0.5), WInit, gamma, Nt, optsIHT);
    [pos, WTemp] = Nonconvex_ProxGrad(Phi, (A11)^(0.5), Nt, gamma);

    for n = 1 : Nt
        B(pos(n), n) = 1;
        wCurrent(n, :) = WTemp(pos(n), :);
    end

end

function rate = KUpdateSINR(vCurrent, wCurrent, HatH1, B, C, UserRecJamPow, Pt, Nt, K, sigma2)
    G = zeros(Nt, K);
    for k = 1 : K
        G(:, k) = B' * HatH1(:, :, k)' * C(:, :, k) * vCurrent(:, k);
    end

    wCurrent = sqrt(Pt) * wCurrent / norm(wCurrent, 'fro');

    He = abs(G' * wCurrent).^2;
    gr = zeros(K, 1);
    num = zeros(K, 1);
    for k = 1 : K
        tmp=He(k,:);
        num(k) = (sum(tmp)-tmp(k) + UserRecJamPow(k) + sigma2 * norm(vCurrent(:, k), 2)^2);
        gr(k)=tmp(k)/num(k);
    end
    % 采用ln
    rate = sum(log(1 + gr));

end

% ==== 针对第k个用户
function res = DualFunc(HatH1, B, C, vCurrent, wCurrent, k, UserRecJamPow, Pt, ACurrent, sigma2)
    
    res = real(ACurrent(1, 1, k)) + 2 * real(ACurrent(1, 2, k) * vCurrent(:,k)' * C' * HatH1 * B * wCurrent(:,k)) ... 
            + real(ACurrent(2, 2, k)) * (norm(vCurrent(:, k)' * C' * HatH1 * B * wCurrent, 2)^2 + ...
             norm(wCurrent(:), 2)^2 * (UserRecJamPow(k) + sigma2 * norm(vCurrent(:, k), 2)^2) / Pt);

end

function res = DualFunc2(HatH1, B, C, vCurrent, wCurrent, K, UserRecJamPow, Pt, ACurrent, sigma2)
    res = zeros(K, 1);
    for k = 1 : K
        res(k) = real(ACurrent(1, 1, k)) + 2 * real(ACurrent(1, 2, k) * vCurrent(:,k)' * C(:, :, k)' * HatH1(:, :, k) * B * wCurrent(:,k)) ... 
            + real(ACurrent(2, 2, k)) * (norm(vCurrent(:, k)' * C(:, :, k)' * HatH1(:, :, k) * B * wCurrent, 2)^2 + ...
             norm(wCurrent(:), 2)^2 * (UserRecJamPow(k) + sigma2 * norm(vCurrent(:, k), 2)^2) / Pt);
    end
end