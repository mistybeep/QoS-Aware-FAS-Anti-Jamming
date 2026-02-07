clear all;
close all;

rng(42);

K = 3;
R = 2;
L = 8;
sigma = sqrt(1e-10);
data = load('dataset.mat');

rho = - 30;
alpha = 2.6;
Sigma1_TEST = zeros(L, L, K, 100);
Sigma2_TEST = zeros(L, L, K, R, 100);
for i = 1 : 100
    for k = 1:K
        Sigma1_TEST(:, :, k, i) = GenSigma(L, rho, alpha, data.distanceTransmit(k), sigma);
        for r = 1 : R
            Sigma2_TEST(:, :, k, r, i) = GenSigma(L, rho, alpha, data.distanceJammer(k, r), sigma);
        end
    end
end
save 'Sigma.mat' Sigma1_TEST Sigma2_TEST