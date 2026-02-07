function SINR = K_Update_SINR(vCurrent, wCurrent, H1, H2, Nt, K, P, sigma2)
    G = zeros(Nt, K);
    for k = 1 : K
        G(:, k) = H1(:, :, k)' * vCurrent(:, k);
    end

    He = abs(G' * wCurrent).^2;
    gr = zeros(K, 1);
    num = zeros(K, 1);
    for k = 1 : K
        tmp=He(k,:);
        num(k) = (sum(tmp)-tmp(k) + P * norm(vCurrent(:, k)' * squeeze(H2(:, k, :)), 2)^2 + sigma2*norm(vCurrent(:, k), 2)^2);
        gr(k)=tmp(k)/num(k);
    end
    
    SINR = sum(log2(1 + gr));

end