function [wCurrent, B] = updatePrecoderWMMSEOMPPC(wCurrent, vCurrent, BCurrent, HatH1, C, Pt, UserRecJamPow, bCurrent, ACurrent, Gt, Gamma, Nt, K, sigma2)
    
    [wCurrent, B] = wCurrentWMMSEOMP(vCurrent, wCurrent, BCurrent, HatH1, C, K, Gt, Nt, UserRecJamPow, Pt, sigma2);
    

    Flag = 1;
    wCurrent = sqrt(Pt) * wCurrent / norm(wCurrent, 'fro');
    G = zeros(Nt, K);
    for k = 1 : K
        G(:, k) = B' * HatH1(:, :, k)' * C(:, :, k) * vCurrent(:, k);
    end

    He = abs(G' * wCurrent).^2;
    for k = 1 : K
        tmp=He(k,:);
        num = (sum(tmp)-tmp(k) + UserRecJamPow(k) + sigma2*norm(vCurrent(:, k), 2)^2);
        gr = log(1 +tmp(k)/num);
        if gr < Gamma
            Flag = 1;
        end
    end
    
    if Flag == 0
        wCurDirec = zeros(Nt, K);
        UserPow = zeros(K, 1);
        for k = 1 : K
            h = B' * HatH1(:, :, k)' * C(:, :, k) * vCurrent(:, k);
            UserPow(k) = norm(wCurrent(:, k), 2)^2;
            wCurDirec(:, k) = wCurrent(:, k) / norm(wCurrent(:, k), 2); %  * exp(-1j*angle(h' * wCurrent(:, k)))
        end
        iter = 0;
        PCIter = 0;
        % PCCurrent = 0.001;
        while(1)
            iter = iter + 1;
            
            [y, z] = KUpdateFrac(vCurrent, UserPow, wCurDirec, HatH1, B, C, UserRecJamPow, Pt, Nt, K, sigma2);
            cvx_begin quiet
                cvx_solver mosek
                variable UserPow(K, 1)
                G = zeros(Nt, K);
                for k = 1 : K
                    G(:, k) = B' * HatH1(:, :, k)' * C(:, :, k) * vCurrent(:, k);
                end
    
                He = abs(G' * wCurDirec).^2;
    
                MinMaxObj = 0;
                for k = 1 : K
                    MinMaxObj = MinMaxObj + z(k)^2 * sum(UserPow .* He(k,:).') - 2*z(k)*sqrt((1+y(k)) * UserPow(k) * He(k, k));
                end
    
                minimize MinMaxObj;
    
                subject to
                    for k = 1 : K
    
                        lhsA = - z(k)^2 * sum(UserPow .* He(k,:).') + 2*z(k)*sqrt((1+y(k)) * UserPow(k) * He(k, k));
    
                        rhsA = Gamma/log(2) - log2(1+y(k)) + y(k) + z(k)^2 * (UserRecJamPow(k) + sigma2);
    
                        lhsA >= rhsA;
    
                    end
                    sum(UserPow) <= Pt;
            cvx_end
            PCCurrent = MinMaxObj;
            if any(isnan(UserPow))
                break;
            else
                if norm(PCIter - PCCurrent) < 1e-3
                    for k = 1 : K
                        wCurrent(:, k) = sqrt(UserPow(k)) * wCurDirec(:, k);
                    end
                    break;
                end
                PCIter = PCCurrent;
            end
        end
    end

    
end


function [wCurrent, B] = wCurrentWMMSEOMP(vCurrent, wCurrent, BCurrent, HatH1, C, K, Gt, Nt, UserRecJamPow, Pt, sigma2)

    Omega = zeros(K, K);
    Psi = zeros(K, Gt);
    gamma = 0;

    B = zeros(Gt, Nt);
    
    for i = 1 : K
        Omega(i, i) = (1 - wCurrent(:, i)' * BCurrent' * HatH1(:, :, i)' * C(:, :, i) * vCurrent(:, i))^(-1);
        Psi(i, :) = vCurrent(:, i)' * C(:, :, i)' * HatH1(:, :, i);
        gamma = gamma + real(Omega(i, i) * (UserRecJamPow(i) + sigma2 * norm(vCurrent(:, i), 2)^2) / Pt);
    end
    
    Phi = (Omega)^(0.5) * Psi;
    [pos, wCurrent] = RLS_SOMP((Omega)^(0.5), Phi, Gt, Nt, gamma);
    % [pos, WTemp] = Row_IHT(Phi, (A11)^(0.5), WInit, gamma, Nt, optsIHT);
    % [pos, WTemp] = Nonconvex_ProxGrad(Phi, (A11)^(0.5), Nt, gamma);

    for n = 1 : Nt
        B(pos(n), n) = 1;
        % wCurrent(n, :) = WTemp(pos(n), :);
    end

end

function [y, z] = KUpdateFrac(vCurrent, UserPow, wCurDirec, HatH1, B, C, UserRecJamPow, Pt, Nt, K, sigma2)
    y = zeros(K, 1);
    z = zeros(K, 1);
    W = zeros(Nt, K);
    for k = 1 : K
        W(:, k) = sqrt(UserPow(k)) * wCurDirec(:, k);
    end

    for k = 1 : K
        h = B' * HatH1(:, :, k)' * C(:, :, k) * vCurrent(:, k);
        y(k) = abs(h' * W(:, k))^2 / (real(h' * W(:, setdiff(1 : K, k)) * W(:, setdiff(1 : K, k))' * h) + UserRecJamPow(k) + sigma2);
    end

    for k = 1 : K
        h = B' * HatH1(:, :, k)' * C(:, :, k) * vCurrent(:, k);
        z(k) = sqrt(1+y(k)) * real(h' * W(:, k)) / (real(h' * (W * W') * h) + UserRecJamPow(k) + sigma2);
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