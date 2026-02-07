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
% if isempty(gcp('nocreate'))
%     numWorkers = 10;
%     parpool('local', numWorkers);
%     % addAttachedFiles(gcp, {'D:\Software\MATLAB\R2023b\bin\cvx'});
% else
%     fprintf('并行池已存在，直接使用\n');
% end

rhoSizeSim = 4 : 8;

Sigma1_cell = cell(100, 1);
Sigma2_cell = cell(100, 1);
for sim = 1 : 100
    Sigma1_cell{sim} = Sigma1_TEST(:, :, :, sim);
    Sigma2_cell{sim} = Sigma2_TEST(:, :, :, :, sim);
end

SimuRes = zeros(length(rhoSizeSim), 16);
for n = 1 : length(rhoSizeSim)    
    rhoSize = rhoSizeSim(n);
    % t = GenPosCurrent(rhoSize, lambda, Nt);
    % r = GenPosCurrent(rhoSize, lambda, Nr);
    
    temp_rates = zeros(100, 16);
    tic
    for sim = 1 : 100
        fprintf('  运行仿真 %d/%d...\n', sim, 100);
        result = FuncDiscrete(data, Nt, Nr, wCurrent, vCurrent, t1, Pt, P, Sigma1_cell{sim}, Sigma2_cell{sim}, err, rhoSize, lambda, Gamma, L, K, R, Q1, Q2, sigma2);
        temp_rates(sim, :) = result;
    end
    toc
    SimuRes(n, :) = mean(temp_rates);
end
save("DiscreteSizeResult.mat", "SimuRes");
% disp(mean(temp_rates))
% plot(mean(temp_rates))
% toc


function res = FuncDiscrete(data, Nt, Nr, wCurrent, vCurrent, t1, Pt, P, Sigma1, Sigma2, err, rhoSize, lambda, Gamma, L, K, R, Q1, Q2, sigma2)
    Sigma1 = sqrt(sigma2) * Sigma1;
    Sigma2 = sqrt(sigma2) * Sigma2;
    
    thetaJammer = data.thetaJamReceive - err / 10;
    phiJammer = data.phiJamReceive + err / 10;

    Gr = (2 * rhoSize + 1)^2;
    Gt = (2 * rhoSize + 1)^2;
    RPos = zeros(2, Gr);
    T = zeros(2, Gt);
    xLow =  - rhoSize * lambda / 2;
    yLow = - rhoSize * lambda / 2;
    fla = 0;
    for gx=1:sqrt(Gr)
        for gy=1:sqrt(Gr)
            fla=fla+1;
            RPos(:, fla) = [(gx-1)*lambda/2 + xLow; (gy-1)*lambda/2+yLow];
        end
    end
    xLow =  - rhoSize * lambda / 2;
    yLow = - rhoSize * lambda / 2;
    fla = 0;
    for gx=1:sqrt(Gt)
        for gy=1:sqrt(Gt)
            fla=fla+1;
            T(:, fla) = [(gx-1)*lambda/2+xLow; (gy-1)*lambda/2+yLow];
        end
    end

    HatH1 = zeros(Gr, Gt, K);
    B = zeros(Gt, Nt);
    Col = vec(extract7x7Square(rhoSize, 7));
    % Col = [1:7, sqrt(Gt) + 1: sqrt(Gt) + 7, 2*sqrt(Gt)+1: 2*sqrt(Gt) + Nt - 14];
    for n = 1 : Nt
        B(Col(n), n) = 1;
    end
    HatH2 = zeros(Gr, K, R, Q1 * Q2);
    HatH2True = zeros(Gr, K, R);
    C = zeros(Gr, Nr, K);
    for k = 1 : K
        % C(:, :, k) = [eye(Nr); zeros(Gr-Nr, Nr)];
        for n = 1 : Nr
            C(Col(n), n, k) = 1;
        end
        [~, HatH1(:, :, k), ~] = GenChannel(T, RPos, lambda, data.thetaTransmit(:, k), data.phiTransmit(:, k), ...
                                   data.thetaReceive(:, k), data.phiReceive(:, k), Sigma1(:, :, k));
        
        % Cur_HH = zeros(Gr, Gr);
        for r = 1 : R
            [~, HatH2True(:, k, r), ~] = GenChannel(t1, RPos, lambda, [], [], ...
                                   data.thetaJamReceive(:, k, r), data.phiJamReceive(:, k, r), Sigma2(:, :, k, r));
            Index = 1;
            for p = 1 : Q1
                theta = thetaJammer(:, k, r) - err/2 + (p-1)*err/(Q1-1);
                for q = 1 : Q2
                    phi = phiJammer(:, k, r) - err/2 + (q-1)*err/(Q2-1);
                    [~, C_H, ~] = GenChannel(t1, RPos, lambda, [], [], ...
                        theta, phi, Sigma2(:, :, k, r));
                    % [~, C_H, ~] = GenChannel(t1, RPos, lambda, [], [], ...
                    %     thetaJammer(:, k, r), phiJammer(:, k, r), Sigma2(:, :, k, r));
                    HatH2(:, k, r, Index) = C_H;
                    Index = Index + 1;
                    % Cur_HH = Cur_HH + C_H * C_H' / (Q1*Q2);
                end
            end
        end
    end
    
    H1 = zeros(Nr, Nt, K);
    
    H2 = zeros(Nr, K, R);
    
    for k = 1 : K
        H1(:, :, k) = C(:, :, k)' * HatH1(:, :, k) * B;
        % === Ture Channel ===
        for r = 1 : R
            H2(:, k, r) = C(:, :, k)' * HatH2True(:, k, r);
        end
    end
    
    SINR = K_Update_SINR(vCurrent, wCurrent, H1, H2, Nt, K, P, sigma2);
   
    iterD = 15;
    SINRIter = zeros(iterD + 1, 1);
    SINRIter(1) = SINR;
    i = 0;
    while i < iterD
        % tic
        i = i + 1;
        [vCurrent, C] = updateCombinerOMP(vCurrent, HatH1, HatH2, B, C, wCurrent, Pt, P, Gr, Nr, R, Q1, Q2, sigma2);
        UserRecJamPow = zeros(K, 1);
        for k = 1 : K
            UserRecJamPow(k) = P * norm(vCurrent(:, k)' * C(:, :, k)' * reshape(squeeze(HatH2(:, k, :, :)), Gr, R*Q1*Q2), 2)^2 / (Q1 * Q2);
        end
        for k = 1 : K
            H1(:, :, k) = C(:, :, k)' * HatH1(:, :, k) * B;
        end
        
        [ACurrent, bCurrent] = GenMinMax(K, wCurrent, vCurrent, H1, UserRecJamPow, Pt, sigma2);

        [wCurrent, B] = updatePrecoderOMPPC(wCurrent, vCurrent, B, HatH1, C, Pt, UserRecJamPow, bCurrent, ACurrent, Gt, Gamma, Nt, K, sigma2);
    
        for k = 1 : K
            H1(:, :, k) = C(:, :, k)' * HatH1(:, :, k) * B;
            % === Ture Channel ===
            for r = 1 : R
                H2(:, k, r) = C(:, :, k)' * HatH2True(:, k, r);
            end
        end
        
        SINR = K_Update_SINR(vCurrent, wCurrent, H1, H2, Nt, K, P, sigma2);
        SINRIter(i + 1, 1) = SINR;
        % if abs(SINRIter(i + 1, 1) - SINRIter(i, 1)) < 1e-3
        %     SINRIter(i + 2: end, 1) = SINRIter(i + 1, 1);
        %     break;
        % end
        % toc
    end
    % toc
    RecPattern(data, 1, vCurrent, C, RPos, R, Nr, lambda)
    TransPattern(data, 1, K, wCurrent, B, T, Nt, lambda)
    % RecPositionChannelGain(data, 1, B, T, C, RPos, R, Sigma1 * sqrt(1e-10), Nt, Nr, rhoSize, L, lambda)
    % for r = 1 : R
    %     RecPositionJammerChannelGain(data, 1, C, RPos, r, Sigma2 * sqrt(1e-10), Nr, rhoSize, L, lambda)
    % end

    % plot(SINRIter)
    res = SINRIter;

end

function smallSquare = extract7x7Square(n, m)
    % 检查n是否大于等于3
    if n < 3
        error('n必须大于等于3');
    end
    
    % 总网格大小
    totalSize = 2*n + 1;

    Size = (m - 1) / 2;
    
    % 创建列优先的坐标矩阵
    coordinates = reshape(1:totalSize^2, totalSize, totalSize);
    % coordinates = flipud(coordinates);
    % 计算7x7区域的中心位置
    center = n + 1;  % 网格中心坐标
    
    % 提取7x7区域（行和列的范围）
    startRow = center - Size;
    endRow = center + Size;
    startCol = center - Size;
    endCol = center + Size;
    
    % 提取小正方形
    smallSquare = coordinates(startRow:endRow, startCol:endCol);
end