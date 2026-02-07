function GenPos(K, L, TranmitPos, UserAngle, JamPos, FileName)
    rng(42)
    UserDis = 100;  % 距离 单位: m

    R = size(JamPos, 2);

    sigma_aod = 5 * pi / 180;

    UserPos = zeros(3, K);
    for i = 1 : K
        theta = deg2rad(UserAngle(i, 1));
        phi = deg2rad(UserAngle(i, 2));

        % 计算相对于Tx_pos的坐标
        [x, y, z] = sph2cart(phi, theta, UserDis);

        % 转换为绝对坐标
        UserPos(1, i) = TranmitPos(1) + x;
        UserPos(2, i) = TranmitPos(2) + y;
        UserPos(3, i) = TranmitPos(3) + z;
        
    end
    thetaTransmit = zeros(L, K);
    phiTransmit = zeros(L, K);
    distanceTransmit = zeros(K, 1);

    distanceJammer = zeros(K, R);

    thetaReceive = zeros(L, K); 
    phiReceive = zeros(L, K);
    thetaJamReceive = zeros(L, K, R);
    phiJamReceive = zeros(L, K, R);

    theta_aod  = sqrt(sigma_aod)*randn(K, L - 1);
    phi_aod = sqrt(sigma_aod)*randn(K, L - 1);

    thetajam_aod  = sqrt(sigma_aod + pi/180) * randn(K, R, L - 1);
    phijam_aod = sqrt(sigma_aod + pi/180) * randn(K, R, L - 1);

    for k = 1:K

        [phiTransmit(1, k), thetaTransmit(1, k), ~] = cart2sph(UserPos(1, k) - TranmitPos(1), UserPos(2, k) - TranmitPos(2), UserPos(3, k) - TranmitPos(3));

        [phiReceive(1, k), thetaReceive(1, k), ~] = cart2sph(- UserPos(1, k) + TranmitPos(1), - UserPos(2, k) + TranmitPos(2), - UserPos(3, k) + TranmitPos(3));

        distanceTransmit(k) = UserDis; % 直射路径距离

        for r = 1 : R
            [phiJamReceive(1, k, r), thetaJamReceive(1, k, r), distanceJammer(k, r)] = cart2sph(- UserPos(1, k) + JamPos(1, r), - UserPos(2, k) + JamPos(2, r), - UserPos(3, k) + JamPos(3, r));
        end
        
        for l = 2 : L

            phiTransmit(l, k) = phiTransmit(1, k) + phi_aod(k, l - 1);

            thetaTransmit(l, k) = thetaTransmit(1, k) + theta_aod(k, l - 1);
            
            [x, y, z] = sph2cart(phiTransmit(l, k), thetaTransmit(l, k), UserDis / 2 + 10);

            [phiReceive(l, k), thetaReceive(l, k), ~] = cart2sph(x - UserPos(1, k), y - UserPos(2, k), z - UserPos(3, k));
            
            for r = 1 : R
                thetaJamReceive(l, k, r) = thetaJamReceive(1, k, r) + thetajam_aod(k, r, l - 1);

                phiJamReceive(l, k, r) = phiJamReceive(1, k, r) + phijam_aod(k, r, l - 1);
            end
            
        end
           
    end
    save(FileName, 'distanceTransmit', 'distanceJammer', 'thetaTransmit', 'phiTransmit', 'thetaReceive', 'phiReceive', 'thetaJamReceive', 'phiJamReceive');
end