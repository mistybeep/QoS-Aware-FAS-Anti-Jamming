function [ACurrent, bCurrent] = GenMinMax(K, wCurrent, vCurrent, H1, UserRecJamPow, Pt, sigma2)
    ACurrent = zeros(2, 2, K);
    bCurrent = zeros(K, 1);
    for k = 1 : K
        D = [1; 0];
        BCurrent = [1, wCurrent(:, k)' * H1(:, :, k)' * vCurrent(:, k);
                vCurrent(:, k)' * H1(:, :, k) * wCurrent(:, k), norm(vCurrent(:, k)' * H1(:, :, k) * wCurrent, 2)^2 + norm(wCurrent, 'fro')^2 * (UserRecJamPow(k) + sigma2*norm(vCurrent(:, k), 2)^2) / Pt];

        ACurrent(:, :, k) = BCurrent^(-1) * D * (D' * BCurrent^(-1) * D)^(-1) * D' * BCurrent^(-1);

        bCurrent(k) = real(log(D' * BCurrent^(-1) * D) + trace(ACurrent(:, :, k) * BCurrent));

        % log(real(D'*BCurrent^(-1)*D))
    end

end