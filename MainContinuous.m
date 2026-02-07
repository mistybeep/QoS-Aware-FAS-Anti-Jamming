clear all;
close all;

rng(42);

% === Parameters ===
Nt = 16;
Nr = 9;
R = 2;
L = 8;
K = 3;
Pt = db2pow(10 - 30);
P = Pt / (R * 10^(-2));
Gamma = 1 * log(2);
rho = - 30;
alpha = 2.2;
fc = 3e9;
lambda = 3e8/fc;
sigma2 = 1;
data = load('dataset.mat');
err = 4 / 180 * pi;
rhoSize = 3;
t1 = GenUPA(1, lambda);
% load("coord.mat", "t", "r");
t = GenUPA(Nt, lambda);
r = GenUPA(Nr, lambda);
K_r = cell(K, 1);
Q1 = 10;
Q2 = 10;

load Sigma.mat Sigma1_TEST Sigma2_TEST;

% ====== Init =======
wCurrent = randn(Nt, K) + 1j*randn(Nt, K);
wCurrent = sqrt(Pt) * wCurrent / norm(wCurrent, 'fro');

vCurrent = randn(Nr, K) + 1j*randn(Nr, K);
vCurrent = vCurrent / norm(vCurrent, 2);

% tic

rhoSizeSim = 4 : 8;
% --- 并行设置 ---
% ✅
% if isempty(gcp('nocreate'))
%     numWorkers = 12;
%     parpool('local', numWorkers);
%     addAttachedFiles(gcp, {'D:\Software\MATLAB\R2023b\bin\cvx'});
% else
%     fprintf('并行池已存在，直接使用\n');
% end

Sigma1_cell = cell(100, 1);
Sigma2_cell = cell(100, 1);
for sim = 1 : 100
    Sigma1_cell{sim} = Sigma1_TEST(:, :, :, sim);
    Sigma2_cell{sim} = Sigma2_TEST(:, :, :, :, sim);
end

% load("t.mat", "t");
% load("r.mat", "r");
rPos = r;

SimuRes = zeros(length(rhoSizeSim), 16);
for n = 1 : length(rhoSizeSim)    
    rhoSize = rhoSizeSim(n);
    % t = GenPosCurrent(rhoSize, lambda, Nt);
    % r = GenPosCurrent(rhoSize, lambda, Nr);
    
    temp_rates = zeros(100, 16);
    tic
    for sim = 1 : 100
        fprintf('  运行仿真 %d/%d...\n', sim, 100);
        result = FuncDiscrete(data, Nt, Nr, wCurrent, vCurrent, t, rPos, t1, Pt, P, Sigma1_cell{sim}, Sigma2_cell{sim}, err, rhoSize, lambda, Gamma, L, K, R, Q1, Q2, sigma2);
        temp_rates(sim, :) = result;
    end
    toc
    SimuRes(n, :) = mean(temp_rates);
end
% save("ContinuousSizeResult.mat", "SimuRes");
% disp(mean(temp_rates))
% plot(mean(temp_rates))
% toc


function res = FuncDiscrete(data, Nt, Nr, wCurrent, vCurrent, t, rPos, t1, Pt, P, Sigma1, Sigma2, err, rhoSize, lambda, Gamma, L, K, R, Q1, Q2, sigma2)
    
    K_r = cell(K, 1);
    H1 = zeros(Nr, Nt, K);
    F1 = zeros(L, Nr, K);
    G1 = zeros(L, Nt, K);
    
    % Sigma2 = zeros(L, L, K);
    H2 = zeros(Nr, K, R);
    F2 = zeros(L, Nr, K, R);

    Sigma1 = sqrt(sigma2) * Sigma1;
    Sigma2 = sqrt(sigma2) * Sigma2;
    
    thetaJammer = data.thetaJamReceive - err / 10;
    phiJammer = data.phiJamReceive + err / 10;
    
    for k = 1 : K
        K_r{k} = rPos;
        [F1(:, :, k), H1(:, :, k), G1(:, :, k)] = GenChannel(t, K_r{k}, lambda, data.thetaTransmit(:, k), data.phiTransmit(:, k), ...
                                                                    data.thetaReceive(:, k), data.phiReceive(:, k), Sigma1(:, :, k));
        % === Ture Channel ===
        for r = 1 : R
            [F2(:, :, k, r), H2(:, k, r), ~] = GenChannel(t1, K_r{k}, lambda, [], [], data.thetaJamReceive(:, k, r), data.phiJamReceive(:, k, r), Sigma2(:, :, k, r));
        end
    end
    
    % === 离散角度 ===
    thetaJammerDiscretization = zeros(L, K, R, Q1*Q2);
    phiJammerDiscretization = zeros(L, K, R, Q1*Q2);
    Index = 1;
    for p = 1 : Q1
        for q = 1 : Q2
            thetaJammerDiscretization(:, :, :, Index) = thetaJammer - err/2 + (p-1)*err/(Q1-1);
            phiJammerDiscretization(:, :, :, Index) = phiJammer - err/2 + (q-1)*err/(Q2-1);
            Index = Index + 1;
        end
    end
    
    SINR = K_Update_SINR(vCurrent, wCurrent, H1, H2, Nt, K, P, sigma2);
   
    % SINR = K_Update_SINR(vCurrent, wCurrent, H1, H2, P, Nt, K);
    iterD = 15;
    SINRIter = zeros(iterD + 1, 1);
    SINRIter(1) = SINR;
    i = 0;
    while i < iterD
        % tic
        i = i + 1;
        vCurrent = updateCombiner(vCurrent, wCurrent, H1, K_r, thetaJammerDiscretization, phiJammerDiscretization, Sigma2, Pt, P, Nr, K, L, R, Q1, Q2, lambda, sigma2);
        % [K_r, F1, H1] = updateRecPosition(data, K_r, wCurrent, vCurrent, Sigma1, G1, Sigma2, thetaJammerDiscretization, phiJammerDiscretization, Pt, P, Nt, Nr, K, L, R, Q1, Q2, lambda, rhoSize, sigma2);
        UserRecJamPow = zeros(K, 1);
        for k = 1 : K
            JTemp = zeros(Q1*Q2, 1);
            G2 = ones(L, 1);
            for r = 1 : R
                for q = 1 : Q1 * Q2
                    F = create_F(K_r{k}, thetaJammerDiscretization(:, k, r, q), phiJammerDiscretization(:, k, r, q), L, Nr, lambda);
                    JTemp(q) = JTemp(q) + P * norm(vCurrent(:, k)' * F' * Sigma2(:, :, k, r) * G2, 2)^2 / (Q1*Q2);
                end
            end
            UserRecJamPow(k) = sum(JTemp);
        end
        % SINR = K_Update_SINR_Iter(vCurrent, wCurrent, H1, UserRecJamPow, P, Nt, K)
        
        [ACurrent, bCurrent] = GenMinMax(K, wCurrent, vCurrent, H1, UserRecJamPow, Pt, sigma2);

        [wCurrent, Flag] = updatePrecoder(wCurrent, vCurrent, H1, UserRecJamPow, bCurrent, ACurrent, Nt, K, Pt, Gamma, sigma2);
        % [t, H1, G1] = updateTransPosition(data, t, wCurrent, vCurrent, Sigma1, F1, UserRecJamPow, ACurrent, bCurrent, Pt, Nt, Nr, K, L, lambda, rhoSize, Flag, Gamma);

        wCurrent = sqrt(Pt) * wCurrent / norm(wCurrent, 'fro');
        % SINR = K_Update_SINR(vCurrent, wCurrent, H1, H2, Nt, K, P, sigma2);
        % scatter(t(1, :), t(2, :), 80, 'filled', 'r');
        % scatter(K_r{1}(1, :), K_r{1}(2, :), 80, 'filled', 'r');
        % figure
        % scatter(K_r{2}(1, :), K_r{2}(2, :), 80, 'filled', 'r');
        % figure
        % scatter(K_r{3}(1, :), K_r{3}(2, :), 80, 'filled', 'r');
        for k = 1 : K
            % === Ture Channel ===
            for r = 1 : R
                [~, H2(:, k, r), ~] = GenChannel(t1, K_r{k}, lambda, [], [], data.thetaJamReceive(:, k, r), data.phiJamReceive(:, k, r), Sigma2(:, :, k, r));
            end
        end
        
        % SINR = K_Update_SINR_Iter(vCurrent, wCurrent, H1, UserRecJamPow, P, Nt, K);
        SINR = K_Update_SINR(vCurrent, wCurrent, H1, H2, Nt, K, P, sigma2);
        SINRIter(i + 1, 1) = SINR;
        if abs(SINRIter(i + 1, 1) - SINRIter(i, 1)) < 1e-3
            SINRIter(i + 2: end, 1) = SINRIter(i + 1, 1);
            break;
        end
        % toc
    end
    % toc
    RecPatternContinus(data, 1, vCurrent, K_r{1}, R, Nr, lambda)
    TransPatternContinous(data, 1, K, wCurrent, t, Nt, lambda)
    % plot(SINRIter)
    res = SINRIter;

end