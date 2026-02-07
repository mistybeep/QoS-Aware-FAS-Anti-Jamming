% 参数设置
K = 3; % 用户数量
L = 8; % 总路径数（包括直射路径和散射路径）
TranmitPos = [0, 0, 30]';  % Tx-DPSM位置 (x,y,z) 单位: m
JamPos = [4, -13.5, 48.5;
          -0.5, -8.5, 66]';  % 干扰源位置 (x,y,z) 单位: m
% TranmitPos = [0, 0, 60]';  % Tx-DPSM位置 (x,y,z) 单位: m
% JamPos = [-16.48, -14.66, 48.68;
%           12.38, -18.76, 66.46]';  % 干扰源位置 (x,y,z) 单位: m
FileName = 'dataset.mat'; % 输出文件名
t = JamPos - TranmitPos;

UserAngle = [-40, 30; -30, 50; -20, 70]; % theta 和 phi

% 调用函数
GenPos(K, L, TranmitPos, UserAngle, JamPos, FileName)

% load('communication_param.mat')