function [breath_data, f0, f_L, f_H] = adaptiveNarrowband(angle_fft_clean, fs, delta, window_sizes, butter_order)
% adaptiveNarrowband 自适应窄带滤波
% 平滑窗 SNR 最大化 -> 主频 -> 固定半带宽 -> Butterworth 带通
% 输入：
%   angle_fft_clean - 运动恢复后的相位差分信号（列向量）
%   fs              - 慢时间采样率（Hz）
%   delta           - 半带宽（论文取 0.25 Hz）
%   window_sizes    - 候选平滑窗长度
%   butter_order    - Butterworth 阶数
% 输出：
%   breath_data     - 滤波后的呼吸信号
%   f0              - 主频（Hz）
%   f_L, f_H        - 通带截止频率（Hz）

if nargin < 3 || isempty(delta), delta = 0.25; end
if nargin < 4 || isempty(window_sizes), window_sizes = [5 10 15 20]; end
if nargin < 5 || isempty(butter_order), butter_order = 4; end

angle_fft_clean = angle_fft_clean(:);

%% 平滑窗长度选择
best_snr = 0;
best_phi = angle_fft_clean;
for ws = window_sizes
    phi_w = smoothdata(angle_fft_clean, 'movmean', ws);
    signal_var = var(phi_w);
    noise_var = var(angle_fft_clean - phi_w);
    if noise_var > 0
        snr_temp = signal_var / noise_var;
        if snr_temp > best_snr
            best_snr = snr_temp;
            best_phi = phi_w;
        end
    end
end
phi = best_phi;

%% 主频
N1 = numel(phi);
phi_fft = abs(fft(phi));
f_axis = (0:N1-1) * (fs/N1);
band = (f_axis >= 0.1) & (f_axis <= 2.0);
f_band = f_axis(band);
fft_band = phi_fft(band);
[~, idx] = max(fft_band);
f0 = f_band(idx);

%% 半带宽与生理边界
f_L = max(0.10, f0 - delta);
f_H = min(2.00, f0 + delta);

%% 4 阶 Butterworth 带通 + 零相位滤波
[b, a] = butter(butter_order, [f_L, f_H]/(fs/2), 'bandpass');
breath_data = filtfilt(b, a, phi);
end
