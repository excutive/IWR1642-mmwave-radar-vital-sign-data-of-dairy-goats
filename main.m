% main 主脚本：读一个 DCA1000 bin 文件，跑完整信号处理流水线，输出呼吸率
% 对应论文的三阶段框架：体动检测与恢复 -> 自适应窄带滤波 -> AIMF-EMD 提取呼吸分量
clc; clear; close all;

%% 数据文件路径（改成你要处理的数据文件）
Filename = '../Data/adc_data_9.5.bin';

%% 参数设置
Q_hi = 3;        % 高迟滞系数
Q_lo = 1.5;      % 低迟滞系数
delta = 0.25;    % 自适应滤波半带宽（Hz）
eps_E = 0.01;    % 能量一致性除零小量
eps_f = 0.001;   % 谱对齐除零小量
eps_p = 1e-3;    % 周期性除零小量
% 以下参数正文未给出，为 Fig.6 伪代码中的可调参数
w = [1 1 1 1 -1 1 1];   % 7 个谱特征的权重（峰度随运动下降，取负权）
stftWin = 64;           % STFT 帧长
stftHop = 8;            % STFT 步长
movWin = 9;             % 稳健标准化的滑动中位数窗口

%% 流水线
% 1. 读取 bin，去交错取 RX1，慢时间抽稀
[process_adc, fs] = readDCA1000(Filename);
% 2. 距离维 FFT、杂波抑制、最大能量选距离单元、相位提取
angle_fft_last2 = rangePhaseExtract(process_adc);
% 3. 体动检测（7 谱特征 + 双阈值迟滞）
motion_mask = detectMotion(angle_fft_last2, fs, Q_hi, Q_lo, w, stftWin, stftHop, movWin);
% 4. 体动恢复（频-相轨迹重建）
angle_fft_clean = restoreMotionPhase(angle_fft_last2, motion_mask, fs);
% 5. 自适应窄带滤波
[breath_data, f0] = adaptiveNarrowband(angle_fft_clean, fs, delta);
% 6. AIMF-EMD 提取呼吸分量并求呼吸率
[f_resp, RR] = aimfEmd(breath_data, fs, f0, eps_E, eps_f, eps_p);

%% 输出结果
fprintf('主频: %.3f Hz, 通带: %.3f-%.3f Hz\n', f0, max(0.10, f0-delta), min(2.00, f0+delta));
fprintf('呼吸频率: %.3f Hz (%.2f 次/分钟)\n', f_resp, RR);
